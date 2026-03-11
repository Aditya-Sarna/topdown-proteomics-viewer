from typing import List, Dict, Tuple
import numpy as np

from ..data.models import FragmentIon, Spectrum


def match_peaks(theoretical_ions: List[FragmentIon],
                spectrum: Spectrum,
                tolerance_ppm: float = 10.0) -> List[FragmentIon]:
    """
    Match theoretical fragment ions to observed peaks within tolerance.
    Returns a copy of each ion with .matched / .observed_mz / .mass_error_* filled in.
    """
    obs_mz  = spectrum.mz_array
    matched: List[FragmentIon] = []

    for ion in theoretical_ions:
        th_mz = ion.mz
        tol   = th_mz * tolerance_ppm / 1e6
        diffs  = np.abs(obs_mz - th_mz)
        idx    = int(np.argmin(diffs))

        if diffs[idx] <= tol:
            err_da  = float(obs_mz[idx] - th_mz)
            err_ppm = err_da / th_mz * 1e6
            matched.append(FragmentIon(
                ion_type=ion.ion_type, position=ion.position,
                charge=ion.charge, mz=ion.mz, mass=ion.mass,
                sequence=ion.sequence, matched=True,
                observed_mz=float(obs_mz[idx]),
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
    cmap = coverage_map(sequence, matched_ions)
    covered = sum(1 for v in cmap.values() if v)
    return covered / len(sequence) * 100.0 if sequence else 0.0
