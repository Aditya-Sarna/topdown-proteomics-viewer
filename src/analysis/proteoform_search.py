"""
Targeted proteoform search engine.
Inspired by the FLASHTnT one-spectrum-one-protein strategy:
  1. Generate candidate proteoforms (full-length + truncations + modifications)
  2. Calculate theoretical fragment ions for each
  3. Match against the MS2 spectrum
  4. Score and rank results
"""

from typing import List, Dict, Optional
import numpy as np

from ..data.models import Spectrum, Proteoform, Modification, SearchResult, FragmentIon
from ..data.amino_acids import PTM_DATABASE, PTM_TARGET_RESIDUES
from .fragment_ions import calc_ions
from .peak_matching import match_peaks, sequence_coverage_pct
from .mass_utils import calc_sequence_mass, mz_to_mass, ppm_error




def _truncation_candidates(sequence: str,
                            max_trunc: int = 10) -> List[tuple]:
    """Return (seq, start, end) for N/C-terminal truncations."""
    cands = []
    n = len(sequence)
    # N-terminal truncations  (lose first 1..max_trunc residues)
    for k in range(1, min(max_trunc + 1, n // 3)):
        cands.append((sequence[k:], k + 1, n))
    # C-terminal truncations
    for k in range(1, min(max_trunc + 1, n // 3)):
        cands.append((sequence[:-k], 1, n - k))
    # Both ends
    for k in range(1, min(6, n // 4)):
        for j in range(1, min(6, n // 4)):
            cands.append((sequence[k:-j], k + 1, n - j))
    return cands


def _mod_candidates(base_cands: List[tuple],
                    variable_mods: List[str],
                    max_mods: int = 2) -> List[tuple]:
    """Enumerate single/double modification combinations on the top candidates."""
    results = []
    for seq, start, end in base_cands:
        if len(seq) < 3:
            continue
        # single mods
        for mod_name in variable_mods:
            if mod_name not in PTM_DATABASE:
                continue
            shift    = PTM_DATABASE[mod_name]
            targets  = PTM_TARGET_RESIDUES.get(mod_name, 'ACDEFGHIKLMNPQRSTVWY')
            for pos_0, aa in enumerate(seq):
                if targets and aa not in targets:
                    continue
                mod = Modification(position=pos_0 + 1, name=mod_name,
                                   mass_shift=shift, residue=aa)
                results.append((seq, start, end, [mod]))
        if max_mods >= 2:
            for i, mod1_name in enumerate(variable_mods):
                for mod2_name in variable_mods[i:]:
                    if mod1_name not in PTM_DATABASE or mod2_name not in PTM_DATABASE:
                        continue
                    s1, s2 = PTM_DATABASE[mod1_name], PTM_DATABASE[mod2_name]
                    t1 = PTM_TARGET_RESIDUES.get(mod1_name, 'ACDEFGHIKLMNPQRSTVWY')
                    t2 = PTM_TARGET_RESIDUES.get(mod2_name, 'ACDEFGHIKLMNPQRSTVWY')
                    pos1 = next((p + 1 for p, aa in enumerate(seq) if not t1 or aa in t1), None)
                    pos2 = next((p + 1 for p, aa in enumerate(seq) if (not t2 or aa in t2) and p + 1 != pos1), None)
                    if pos1 and pos2:
                        mods = [
                            Modification(pos1, mod1_name, s1, seq[pos1-1]),
                            Modification(pos2, mod2_name, s2, seq[pos2-1]),
                        ]
                        results.append((seq, start, end, mods))
    return results


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def _score(n_matched: int, n_b: int, n_y: int, n_c: int, n_z: int,
           frac: float, obs_mass: float, th_mass: float) -> float:
    """Simple hyperscore-inspired scoring function."""
    import math
    ion_score = (math.log(n_b + 1) + math.log(n_y + 1) +
                 math.log(n_c + 1) + math.log(n_z + 1))
    mass_penalty = 1.0 / (1.0 + abs(ppm_error(obs_mass, th_mass)) / 5.0) if th_mass else 1.0
    return round(ion_score * frac * 100.0 * mass_penalty, 3)


# ---------------------------------------------------------------------------
# Main search function
# ---------------------------------------------------------------------------

def run_targeted_search(spectrum: Spectrum,
                         protein_sequence: str,
                         protein_name: str = "Target Protein",
                         ion_types: Optional[List[str]] = None,
                         max_charge: int = 4,
                         tolerance_ppm: float = 10.0,
                         search_truncations: bool = True,
                         search_modifications: bool = True,
                         variable_mods: Optional[List[str]] = None,
                         max_mods: int = 1,
                         progress_cb=None) -> List[SearchResult]:
    """
    One-spectrum-one-protein targeted search.

    Returns a ranked list of SearchResult (best first).
    """
    if ion_types is None:
        ion_types = ['b', 'y', 'c', 'z']
    if variable_mods is None:
        variable_mods = ['Phosphorylation', 'Acetylation', 'Oxidation']

    # Derive observed mass from precursor
    obs_mass = 0.0
    if spectrum.precursor_mass > 0:
        obs_mass = spectrum.precursor_mass
    elif spectrum.precursor_mz > 0 and spectrum.precursor_charge > 0:
        obs_mass = mz_to_mass(spectrum.precursor_mz, spectrum.precursor_charge)

    # --- Build candidates ---
    base_cands = [(protein_sequence, 1, len(protein_sequence))]
    if search_truncations:
        base_cands += _truncation_candidates(protein_sequence)

    full_cands = [(s, st, en, []) for s, st, en in base_cands]

    if search_modifications:
        mod_cands = _mod_candidates(base_cands[:8], variable_mods, max_mods)
        full_cands += [(s, st, en, mods) for s, st, en, mods in mod_cands[:60]]

    # Limit total candidates
    full_cands = full_cands[:100]

    # --- Score each candidate ---
    results: List[SearchResult] = []
    n_cands = len(full_cands)

    for cand_idx, (seq, start, end, mods) in enumerate(full_cands):
        if progress_cb and n_cands > 0:
            pct = int(15 + (cand_idx / n_cands) * 70)   # 15 → 85 %
            progress_cb(pct, f'{pct}%')
        if len(seq) < 3:
            continue
        mod_map: Dict[int, float] = {m.position: m.mass_shift for m in mods}
        th_mass = calc_sequence_mass(seq, mods)

        ions = calc_ions(seq, ion_types, mod_map, max_charge)
        if not ions:
            continue

        matched = match_peaks(ions, spectrum, tolerance_ppm)
        n_total   = len(matched)
        n_matched = n_b = n_y = n_c = n_z = 0
        for _ion in matched:
            if _ion.matched:
                n_matched += 1
                t = _ion.ion_type
                if t == 'b':   n_b += 1
                elif t == 'y': n_y += 1
                elif t == 'c': n_c += 1
                elif t == 'z': n_z += 1
        if n_total == 0 or n_matched == 0:
            continue

        frac = n_matched / n_total

        score = _score(n_matched, n_b, n_y, n_c, n_z, frac, obs_mass, th_mass)
        cov   = sequence_coverage_pct(seq, matched)

        err_da  = round(obs_mass - th_mass, 4) if obs_mass else 0.0
        err_ppm = round(ppm_error(obs_mass, th_mass), 2) if obs_mass and th_mass else 0.0

        pf = Proteoform(
            sequence=seq,
            protein_name=protein_name,
            start_pos=start,
            end_pos=end,
            modifications=mods,
            theoretical_mass=round(th_mass, 4),
            observed_mass=round(obs_mass, 4),
            mass_error_da=err_da,
            mass_error_ppm=err_ppm,
            score=score,
            matched_ions=n_matched,
            total_ions=n_total,
        )
        results.append(SearchResult(
            proteoform=pf,
            fragment_ions=matched,
            matched_b=n_b, matched_y=n_y,
            matched_c=n_c, matched_z=n_z,
            sequence_coverage=round(cov, 1),
        ))

    results.sort(key=lambda r: r.proteoform.score, reverse=True)
    return results[:20]
