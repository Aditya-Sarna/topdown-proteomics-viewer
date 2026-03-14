from typing import List, Dict, Tuple
import numpy as np

from ..data.models import FragmentIon, Spectrum


def match_peaks(theoretical_ions: List[FragmentIon],
                spectrum: Spectrum,
                tolerance_ppm: float = 10.0) -> List[FragmentIon]:
    """
    Match theoretical fragment ions to observed peaks within tolerance.
    Uses binary search (O(n log m)) instead of a per-ion linear scan (O(n*m)).
    Returns a copy of each ion with .matched / .observed_mz / .mass_error_* filled in.
    """
    if not theoretical_ions:
        return []

    obs_sorted = np.sort(np.asarray(spectrum.mz_array, dtype=np.float64))
    n_obs = len(obs_sorted)

    # Empty spectrum — return all ions as unmatched
    if n_obs == 0:
        return [
            FragmentIon(
                ion_type=ion.ion_type, position=ion.position,
                charge=ion.charge, mz=ion.mz, mass=ion.mass,
                sequence=ion.sequence, matched=False,
            )
            for ion in theoretical_ions
        ]

    th_mzs = np.array([ion.mz for ion in theoretical_ions], dtype=np.float64)
    tols   = th_mzs * tolerance_ppm / 1e6

    # Binary search: find insertion points for all theoretical m/z values at once
    idxs    = np.searchsorted(obs_sorted, th_mzs)
    i_left  = np.clip(idxs - 1, 0, n_obs - 1)
    i_right = np.clip(idxs,     0, n_obs - 1)

    diff_left  = np.abs(obs_sorted[i_left]  - th_mzs)
    diff_right = np.abs(obs_sorted[i_right] - th_mzs)

    use_left   = diff_left <= diff_right
    best_diff  = np.where(use_left, diff_left,          diff_right)
    best_obs   = np.where(use_left, obs_sorted[i_left], obs_sorted[i_right])
    is_matched = best_diff <= tols

    matched: List[FragmentIon] = []
    for i, ion in enumerate(theoretical_ions):
        if is_matched[i]:
            err_da  = float(best_obs[i] - th_mzs[i])
            err_ppm = err_da / th_mzs[i] * 1e6
            matched.append(FragmentIon(
                ion_type=ion.ion_type, position=ion.position,
                charge=ion.charge, mz=ion.mz, mass=ion.mass,
                sequence=ion.sequence, matched=True,
                observed_mz=float(best_obs[i]),
                mass_error_da=round(err_da, 6),
                mass_error_ppm=round(err_ppm, 3),
            ))
        else:
            matched.append(FragmentIon(
                ion_type=ion.ion_type, position=ion.position,
                charge=ion.charge, mz=ion.mz, mass=ion.mass,
                sequence=ion.sequence, matched=False,
            ))
    return matched


def coverage_map(sequence: str,
                 matched_ions: List[FragmentIon]) -> Dict[int, List[str]]:
    """
    Return a dict mapping 0-based residue index -> list of ion types covering it.
    N-terminal ions (b/c/a) cover positions 0 .. position-1.
    C-terminal ions (y/z)   cover positions n-position .. n-1.
    """
    n    = len(sequence)
    cmap = {i: [] for i in range(n)}
    for ion in matched_ions:
        if not ion.matched:
            continue
        if ion.ion_type in ('b', 'c', 'a'):
            for p in range(ion.position):
                cmap[p].append(ion.ion_type)
        elif ion.ion_type in ('y', 'z'):
            for p in range(n - ion.position, n):
                cmap[p].append(ion.ion_type)
    return cmap


def sequence_coverage_pct(sequence: str,
                           matched_ions: List[FragmentIon]) -> float:
    n = len(sequence)
    if n == 0:
        return 0.0
    covered = np.zeros(n, dtype=bool)
    for ion in matched_ions:
        if not ion.matched:
            continue
        if ion.ion_type in ('b', 'c', 'a'):
            covered[:ion.position] = True
        elif ion.ion_type in ('y', 'z'):
            covered[n - ion.position:] = True
    return float(covered.sum()) / n * 100.0
