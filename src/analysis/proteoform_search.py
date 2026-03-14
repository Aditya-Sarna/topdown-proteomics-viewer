"""
Targeted proteoform search engine.
Inspired by the FLASHTnT one-spectrum-one-protein strategy:
  1. Generate candidate proteoforms (full-length + truncations + modifications)
  2. Generate reversed-sequence decoys for FDR estimation
  3. Calculate theoretical fragment ions for each
  4. Match against the MS2 spectrum
  5. Score, compute e-values (binomial model), assign q-values (target-decoy)
  6. Return ranked target results with statistical annotations
"""

from typing import List, Dict, Optional, Tuple
import math
import numpy as np
from scipy.stats import binom as _binom

from ..data.models import Spectrum, Proteoform, Modification, SearchResult, FragmentIon
from ..data.amino_acids import PTM_DATABASE, PTM_TARGET_RESIDUES, AA_MASSES as _AA_MASSES
_VALID_AA = frozenset(_AA_MASSES)
from .fragment_ions import calc_ions
from .peak_matching import match_peaks, sequence_coverage_pct
from .mass_utils import calc_sequence_mass, mz_to_mass, ppm_error


# ---------------------------------------------------------------------------
# Candidate generation
# ---------------------------------------------------------------------------

def _truncation_candidates(sequence: str, max_trunc: int = 10) -> List[tuple]:
    """Return (seq, start, end) for N/C-terminal truncations."""
    cands = []
    n = len(sequence)
    for k in range(1, min(max_trunc + 1, n // 3)):
        cands.append((sequence[k:], k + 1, n))
    for k in range(1, min(max_trunc + 1, n // 3)):
        cands.append((sequence[:-k], 1, n - k))
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
        for mod_name in variable_mods:
            if mod_name not in PTM_DATABASE:
                continue
            shift   = PTM_DATABASE[mod_name]
            targets = PTM_TARGET_RESIDUES.get(mod_name, 'ACDEFGHIKLMNPQRSTVWY')
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
                    pos2 = next((p + 1 for p, aa in enumerate(seq)
                                 if (not t2 or aa in t2) and p + 1 != pos1), None)
                    if pos1 and pos2:
                        mods = [
                            Modification(pos1, mod1_name, s1, seq[pos1 - 1]),
                            Modification(pos2, mod2_name, s2, seq[pos2 - 1]),
                        ]
                        results.append((seq, start, end, mods))
    return results


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def _score(n_matched: int, n_b: int, n_y: int, n_c: int, n_z: int,
           frac: float, obs_mass: float, th_mass: float) -> float:
    """Hyperscore-inspired scoring function."""
    ion_score = (math.log(n_b + 1) + math.log(n_y + 1) +
                 math.log(n_c + 1) + math.log(n_z + 1))
    mass_penalty = 1.0 / (1.0 + abs(ppm_error(obs_mass, th_mass)) / 5.0) if th_mass else 1.0
    return round(ion_score * frac * 100.0 * mass_penalty, 3)


# ---------------------------------------------------------------------------
# E-value  (binomial random-match model, X!Tandem-style)
# ---------------------------------------------------------------------------

def _calc_evalue(n_matched: int, n_total_ions: int,
                 spectrum: Spectrum, tolerance_ppm: float,
                 n_candidates: int) -> float:
    """
    E-value = (number of candidates) × P(X ≥ n_matched)
    where X ~ Binomial(n_total_ions, p_rand) and p_rand is the probability
    that a single theoretical ion matches a random peak by chance.
    """
    if n_total_ions == 0 or n_matched == 0 or n_candidates == 0:
        return float(n_candidates)

    mz_arr = spectrum.mz_array
    n_peaks = len(mz_arr)
    if n_peaks == 0:
        return float(n_candidates)

    mz_range = float(mz_arr.max() - mz_arr.min())
    if mz_range <= 0:
        return float(n_candidates)

    median_mz = float(np.median(mz_arr))
    tol_da    = median_mz * tolerance_ppm / 1e6
    p_rand    = min(n_peaks * 2.0 * tol_da / mz_range, 0.9999)

    # P(X >= n_matched) = survival function of Binomial
    p_value = float(_binom.sf(n_matched - 1, n_total_ions, p_rand))
    return max(n_candidates * p_value, 1e-300)


# ---------------------------------------------------------------------------
# FDR / q-value  (target-decoy competition, Elias & Gygi 2007)
# ---------------------------------------------------------------------------

def _assign_qvalues(target_results: List[SearchResult],
                    decoy_results:  List[SearchResult]) -> None:
    """
    Mutates target_results in-place to set .q_value on each entry.

    Algorithm:
      1. Merge targets + decoys, sort by score descending.
      2. Walk the sorted list; at each position:
           FDR = #decoys_seen / max(#targets_seen, 1)
      3. q-value = running minimum FDR from that position onward
         (ensures monotonicity).
    """
    if not target_results:
        return

    all_entries: List[Tuple[SearchResult, bool]] = (
        [(r, False) for r in target_results] +
        [(r, True)  for r in decoy_results]
    )
    all_entries.sort(key=lambda x: x[0].proteoform.score, reverse=True)

    n_t = n_d = 0
    raw_fdrs: List[float] = []
    for r, is_decoy in all_entries:
        if is_decoy:
            n_d += 1
        else:
            n_t += 1
        raw_fdrs.append(n_d / n_t if n_t > 0 else 1.0)

    # Backward running minimum (monotone q-values)
    min_fdr = 1.0
    qvals: List[float] = []
    for f in reversed(raw_fdrs):
        min_fdr = min(min_fdr, f)
        qvals.append(min_fdr)
    qvals.reverse()

    # Assign q-value to each target result using object identity
    target_qval: Dict[int, float] = {
        id(r): qvals[i]
        for i, (r, is_decoy) in enumerate(all_entries)
        if not is_decoy
    }
    for r in target_results:
        r.q_value = round(target_qval.get(id(r), 1.0), 4)


# ---------------------------------------------------------------------------
# Core candidate-scoring helper (shared by target and decoy passes)
# ---------------------------------------------------------------------------

def _score_candidates(candidates: List[Tuple],
                      spectrum: Spectrum,
                      ion_types: List[str],
                      max_charge: int,
                      tolerance_ppm: float,
                      obs_mass: float,
                      protein_name: str,
                      is_decoy: bool) -> List[SearchResult]:
    results: List[SearchResult] = []
    for seq, start, end, mods in candidates:
        if len(seq) < 3:
            continue
        mod_map: Dict[int, float] = {m.position: m.mass_shift for m in mods}
        th_mass = calc_sequence_mass(seq, mods)
        ions    = calc_ions(seq, ion_types, mod_map, max_charge)
        if not ions:
            continue

        matched   = match_peaks(ions, spectrum, tolerance_ppm)
        n_total   = len(matched)
        n_matched = n_b = n_y = n_c = n_z = 0
        for _ion in matched:
            if _ion.matched:
                n_matched += 1
                t = _ion.ion_type
                if   t == 'b': n_b += 1
                elif t == 'y': n_y += 1
                elif t == 'c': n_c += 1
                elif t == 'z': n_z += 1
        if n_total == 0 or n_matched == 0:
            continue

        frac  = n_matched / n_total
        score = _score(n_matched, n_b, n_y, n_c, n_z, frac, obs_mass, th_mass)
        cov   = sequence_coverage_pct(seq, matched)

        err_da  = round(obs_mass - th_mass, 4) if obs_mass else 0.0
        err_ppm = round(ppm_error(obs_mass, th_mass), 2) if obs_mass and th_mass else 0.0

        pf = Proteoform(
            sequence=seq,
            protein_name=protein_name + (' [DECOY]' if is_decoy else ''),
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
            is_decoy=is_decoy,
        ))
    return results


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
                        max_mods: int = 1) -> List[SearchResult]:
    """
    One-spectrum-one-protein targeted search with statistical scoring.

    Returns target results ranked by score, each annotated with:
      - e_value  : expected number of false positives with this score or better
      - q_value  : minimum FDR at this score threshold (target-decoy)
    """
    if ion_types is None:
        ion_types = ['b', 'y', 'c', 'z']
    if variable_mods is None:
        variable_mods = ['Phosphorylation', 'Acetylation', 'Oxidation']

    obs_mass = 0.0
    if spectrum.precursor_mass > 0:
        obs_mass = spectrum.precursor_mass
    elif spectrum.precursor_mz > 0 and spectrum.precursor_charge > 0:
        obs_mass = mz_to_mass(spectrum.precursor_mz, spectrum.precursor_charge)

    # ── Build target candidates ────────────────────────────────────────────
    base_cands = [(protein_sequence, 1, len(protein_sequence))]
    if search_truncations:
        base_cands += _truncation_candidates(protein_sequence)

    target_cands = [(s, st, en, []) for s, st, en in base_cands]
    if search_modifications:
        mod_cands = _mod_candidates(base_cands[:8], variable_mods, max_mods)
        target_cands += [(s, st, en, mods) for s, st, en, mods in mod_cands[:60]]
    target_cands = target_cands[:100]

    # ── Build decoy candidates (reversed sequence) ─────────────────────────
    decoy_seq   = protein_sequence[::-1]
    decoy_base  = [(decoy_seq, 1, len(decoy_seq))]
    if search_truncations:
        decoy_base += _truncation_candidates(decoy_seq)
    decoy_cands = [(s, st, en, []) for s, st, en in decoy_base]
    if search_modifications:
        d_mod = _mod_candidates(decoy_base[:8], variable_mods, max_mods)
        decoy_cands += [(s, st, en, mds) for s, st, en, mds in d_mod[:60]]
    decoy_cands = decoy_cands[:100]

    n_total_candidates = len(target_cands) + len(decoy_cands)

    # ── Score both passes ──────────────────────────────────────────────────
    target_results = _score_candidates(target_cands, spectrum, ion_types,
                                       max_charge, tolerance_ppm, obs_mass,
                                       protein_name, is_decoy=False)
    decoy_results  = _score_candidates(decoy_cands, spectrum, ion_types,
                                       max_charge, tolerance_ppm, obs_mass,
                                       protein_name, is_decoy=True)

    # ── E-values ───────────────────────────────────────────────────────────
    for r in target_results + decoy_results:
        r.e_value = round(_calc_evalue(
            r.proteoform.matched_ions,
            r.proteoform.total_ions,
            spectrum, tolerance_ppm,
            n_total_candidates,
        ), 6)

    # ── q-values (target-decoy FDR) ────────────────────────────────────────
    target_results.sort(key=lambda r: r.proteoform.score, reverse=True)
    decoy_results.sort( key=lambda r: r.proteoform.score, reverse=True)
    _assign_qvalues(target_results, decoy_results)

    return target_results[:20]


# ---------------------------------------------------------------------------
# Database search  (multi-protein FASTA mode)
# ---------------------------------------------------------------------------

_DB_SEARCH_MAX_PROTEINS  = 2000  # hard cap per run (mass-filtered first)
_DB_SEARCH_MAX_CANDS_PER = 30    # candidates per protein (speed vs recall)
_DB_MASS_TOL_DA          = 3.0   # pre-filter window (Da) around precursor mass


def _passes_mass_filter(seq: str, obs_mass: float, tol_da: float) -> bool:
    """Quick mass check: does the full-length sequence mass lie within tol_da of obs_mass?

    Uses a fast residue sum (no Modification objects).  The window is wide
    enough to still include truncated forms — the per-aa mass range is used
    to widen the window by (max_trunc residues × heaviest_aa mass).
    """
    if obs_mass <= 0:
        return True   # no precursor info — don't filter
    full_mass = sum(_AA_MASSES.get(c, 0.0) for c in seq) + 18.01056  # +water
    # allow truncations: up to 4 residues removed can shift mass by ~600 Da
    # allow mods: common mods shift by < 160 Da
    margin = tol_da + 4 * 186.08 + 160.0   # ~1$10 Da generous window
    return abs(full_mass - obs_mass) <= margin


def run_database_search(
        spectrum: Spectrum,
        proteins: List[tuple],          # [(name, sequence), ...]
        ion_types: Optional[List[str]] = None,
        max_charge: int = 4,
        tolerance_ppm: float = 10.0,
        search_truncations: bool = True,
        top_n: int = 25) -> List[SearchResult]:
    """
    Search a spectrum against every protein in a FASTA-derived list.

    Strategy
    --------
    1. Pre-filter proteins by precursor mass (±_DB_MASS_TOL_DA + truncation
       margin) so only plausible proteins are scored.
    2. For each surviving protein: generate full-length + truncation
       candidates (no variable mods — speed) plus reversed-sequence decoys.
    3. E-values and q-values are computed globally across *all* proteins so
       that the FDR estimate reflects the full competition space.

    Returns
    -------
    At most *top_n* target SearchResult objects, ranked by score, each
    annotated with cross-database e_value and q_value.
    """
    if ion_types is None:
        ion_types = ['b', 'y', 'c', 'z']

    obs_mass = 0.0
    if spectrum.precursor_mass > 0:
        obs_mass = spectrum.precursor_mass
    elif spectrum.precursor_mz > 0 and spectrum.precursor_charge > 0:
        obs_mass = mz_to_mass(spectrum.precursor_mz, spectrum.precursor_charge)

    all_target: List[SearchResult] = []
    all_decoy:  List[SearchResult] = []
    total_cands = 0

    for prot_name, prot_seq in proteins[:_DB_SEARCH_MAX_PROTEINS]:
        seq = ''.join(c for c in prot_seq.upper() if c in _VALID_AA)
        if len(seq) < 5:
            continue

        # ── Precursor mass pre-filter ──────────────────────────────────────
        if not _passes_mass_filter(seq, obs_mass, _DB_MASS_TOL_DA):
            continue

        base = [(seq, 1, len(seq))]
        if search_truncations:
            base += _truncation_candidates(seq, max_trunc=4)
        target_cands = [(s, st, en, []) for s, st, en in base[:_DB_SEARCH_MAX_CANDS_PER]]

        decoy_seq = seq[::-1]
        decoy_base = [(decoy_seq, 1, len(decoy_seq))]
        if search_truncations:
            decoy_base += _truncation_candidates(decoy_seq, max_trunc=4)
        decoy_cands = [(s, st, en, []) for s, st, en in decoy_base[:_DB_SEARCH_MAX_CANDS_PER]]

        total_cands += len(target_cands) + len(decoy_cands)

        all_target.extend(_score_candidates(
            target_cands, spectrum, ion_types, max_charge, tolerance_ppm,
            obs_mass, prot_name, is_decoy=False))
        all_decoy.extend(_score_candidates(
            decoy_cands, spectrum, ion_types, max_charge, tolerance_ppm,
            obs_mass, prot_name, is_decoy=True))

    if not all_target:
        return []

    # E-values: use total candidate count so the null model reflects full DB size
    for r in all_target + all_decoy:
        r.e_value = round(_calc_evalue(
            r.proteoform.matched_ions,
            r.proteoform.total_ions,
            spectrum, tolerance_ppm,
            max(total_cands, 1),
        ), 6)

    # Global target-decoy q-values
    all_target.sort(key=lambda r: r.proteoform.score, reverse=True)
    all_decoy.sort( key=lambda r: r.proteoform.score, reverse=True)
    _assign_qvalues(all_target, all_decoy)

    return all_target[:top_n]

