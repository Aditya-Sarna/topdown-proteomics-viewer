from typing import List, Dict, Optional
from ..data.amino_acids import AA_MASSES, WATER, PROTON, NH3, CO
from ..data.models import FragmentIon


# ---------------------------------------------------------------------------
# OpenMS-backed theoretical spectrum (primary when pyopenms is available)
# ---------------------------------------------------------------------------

def _calc_ions_openms(sequence: str,
                      ion_types: List[str],
                      mod_map: Optional[Dict[int, float]] = None,
                      max_charge: int = 3) -> Optional[List[FragmentIon]]:
    """
    Generate fragment ions via pyopenms.TheoreticalSpectrumGenerator.
    Returns None if pyopenms is unavailable so the caller can fall back.
    """
    try:
        # import pyopenms  # temporarily disabled
        raise ImportError("pyopenms disabled")
    except ImportError:
        return None

    mod_map = mod_map or {}

    # Build an AASequence, injecting residue mass shifts as unnamed mods
    aa_seq = pyopenms.AASequence.fromString(sequence)
    for pos_1based, shift in mod_map.items():
        idx = pos_1based - 1
        if 0 <= idx < len(sequence):
            res = aa_seq.getResidue(idx)
            mod = pyopenms.ResidueModification()
            mod.setDiffMonoMass(shift)
            aa_seq.setModification(idx, mod)

    tsg       = pyopenms.TheoreticalSpectrumGenerator()
    params    = tsg.getParameters()
    ion_flags = {
        'b': b'add_b_ions',
        'y': b'add_y_ions',
        'a': b'add_a_ions',
        'c': b'add_c_ions',
        'z': b'add_z_ions',
    }
    for key in ion_flags.values():
        params.setValue(key, b'false')
    for itype in ion_types:
        if itype in ion_flags:
            params.setValue(ion_flags[itype], b'true')
    params.setValue(b'add_metainfo', b'true')
    params.setValue(b'max_charge', max_charge)
    tsg.setParameters(params)

    spec = pyopenms.MSSpectrum()
    tsg.getSpectrum(spec, aa_seq, 1, max_charge)

    ions: List[FragmentIon] = []
    for peak in spec:
        ann_list = []
        peak.getMetaValue(b'IonName', ann_list)
        mz = peak.getMZ()
        if not ann_list:
            continue
        ann   = ann_list[0].decode('utf-8') if isinstance(ann_list[0], bytes) else str(ann_list[0])
        # annotation format: "b3+2" / "y12+" / "c5++"
        import re
        m = re.match(r'([abcxyz])(\d+)(\+{1,4}|\^?\d+)', ann)
        if not m:
            continue
        itype    = m.group(1)
        position = int(m.group(2))
        charge_s = m.group(3)
        charge   = charge_s.count('+') if '+' in charge_s else int(charge_s.lstrip('^'))
        charge   = max(charge, 1)
        mass     = mz * charge - charge * PROTON
        if itype in ('b', 'a', 'c'):
            frag_seq = sequence[:position]
        else:
            frag_seq = sequence[len(sequence) - position:]
        ions.append(FragmentIon(
            ion_type=itype, position=position, charge=charge,
            mz=round(mz, 6), mass=round(mass, 6), sequence=frag_seq,
        ))
    return ions if ions else None


def _b_series_masses(sequence: str, mod_map: Dict[int, float]) -> List[float]:
    """Cumulative N-terminal residue masses (b-ion neutral masses)."""
    masses, cum = [], 0.0
    for i, aa in enumerate(sequence[:-1], start=1):
        cum += AA_MASSES.get(aa, 0.0) + mod_map.get(i, 0.0)
        masses.append(cum)
    return masses


def _y_series_masses(sequence: str, mod_map: Dict[int, float]) -> List[float]:
    """Cumulative C-terminal residue masses (y-ion neutral masses)."""
    n = len(sequence)
    masses, cum = [], WATER
    for i in range(n - 1, 0, -1):
        aa  = sequence[i]
        pos = i + 1   # 1-based
        cum += AA_MASSES.get(aa, 0.0) + mod_map.get(pos, 0.0)
        masses.append(cum)
    return masses


def calc_ions(sequence: str,
              ion_types: List[str],
              mod_map: Optional[Dict[int, float]] = None,
              max_charge: int = 3) -> List[FragmentIon]:
    """
    Calculate fragment ions for a given sequence.
    Uses pyopenms.TheoreticalSpectrumGenerator when available; falls back to
    the built-in recursive implementation.
    """
    result = _calc_ions_openms(sequence, ion_types, mod_map, max_charge)
    if result is not None:
        return result
    # ── pure-Python fallback ────────────────────────────────────────────────
    if mod_map is None:
        mod_map = {}

    b_masses = _b_series_masses(sequence, mod_map)
    y_masses = _y_series_masses(sequence, mod_map)

    offsets = {
        'b': (b_masses, 0.0,   sequence, True),   # N-terminal
        'a': (b_masses, -CO,   sequence, True),
        'c': (b_masses, +NH3,  sequence, True),
        'y': (y_masses, 0.0,   sequence[::-1], False),
        'z': (y_masses, -NH3,  sequence[::-1], False),
    }

    ions: List[FragmentIon] = []
    for itype in ion_types:
        if itype not in offsets:
            continue
        mass_list, offset, seq_ref, n_term = offsets[itype]
        n = len(sequence)
        for pos_idx, neutral_mass in enumerate(mass_list):
            position = pos_idx + 1   # ordinal (fragment size)
            mass = neutral_mass + offset
            if mass <= 0:
                continue
            # Build fragment sequence label
            if n_term:
                frag_seq = sequence[:position]
            else:
                frag_seq = sequence[n - position:]

            for z in range(1, max_charge + 1):
                mz = (mass + z * PROTON) / z
                if mz < 50:
                    continue
                ions.append(FragmentIon(
                    ion_type=itype,
                    position=position,
                    charge=z,
                    mz=round(mz, 6),
                    mass=round(mass, 6),
                    sequence=frag_seq,
                ))
    return ions


# ---------------------------------------------------------------------------
# Internal fragment ions (b-type, singly charged, two backbone cleavages)
# ---------------------------------------------------------------------------

def calc_internal_ions(
        sequence: str,
        obs_mz,
        tolerance_ppm: float = 20.0,
        max_frag_length: int = 15,
):
    """
    Compute b-type internal fragment ions and match them against an observed
    spectrum (obs_mz, already sorted).

    Returns a list of tuples:
        (start_0idx, end_0idx_exclusive, mz_theoretical, matched,
         observed_mz, ppm_error)

    Only internal fragments (excluding N- and C-terminal residues) of length
    2..max_frag_length are considered.
    """
    import numpy as np

    n = len(sequence)
    if n < 3:
        return []

    # ── Pass 1: enumerate all (start, end) pairs and accumulate masses ──────
    # This avoids creating one numpy array per fragment (was ~1000+ allocs).
    aa_vals = [AA_MASSES.get(c, 0.0) for c in sequence]
    starts, ends, th_mzs = [], [], []
    for i in range(1, n - 1):
        cum = 0.0
        limit = min(i + max_frag_length, n - 1)
        for j in range(i, limit):
            cum += aa_vals[j]
            if j - i + 1 >= 2:          # min fragment length = 2
                starts.append(i)
                ends.append(j + 1)
                th_mzs.append(cum + PROTON)   # b-type: no water

    if not th_mzs:
        return []

    th_arr = np.array(th_mzs, dtype=np.float64)
    obs_arr = np.asarray(obs_mz, dtype=np.float64)

    if len(obs_arr) == 0:
        return [(s, e, float(m), False, 0.0, 0.0)
                for s, e, m in zip(starts, ends, th_mzs)]

    # ── Pass 2: single batch binary-search match (replaces per-fragment argmin) ─
    idxs  = np.searchsorted(obs_arr, th_arr)
    i_l   = np.clip(idxs - 1, 0, len(obs_arr) - 1)
    i_r   = np.clip(idxs,     0, len(obs_arr) - 1)
    dl    = np.abs(obs_arr[i_l] - th_arr)
    dr    = np.abs(obs_arr[i_r] - th_arr)
    left  = dl <= dr
    best_diff = np.where(left, dl, dr)
    best_obs  = np.where(left, obs_arr[i_l], obs_arr[i_r])
    tols      = th_arr * tolerance_ppm / 1e6
    matched   = best_diff <= tols
    ppm_errs  = best_diff / th_arr * 1e6

    return [
        (starts[k], ends[k], float(th_arr[k]), bool(matched[k]),
         float(best_obs[k]) if matched[k] else 0.0,
         float(ppm_errs[k]) if matched[k] else 0.0)
        for k in range(len(starts))
    ]
