from typing import List, Dict, Optional
from ..data.amino_acids import AA_MASSES, WATER, PROTON, NH3, CO
from ..data.models import FragmentIon


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

    Parameters
    ----------
    sequence   : protein/peptide sequence string
    ion_types  : list of ion types, e.g. ['b', 'y', 'c', 'z']
    mod_map    : dict mapping 1-based position -> mass shift
    max_charge : maximum charge state to calculate

    Returns
    -------
    List of FragmentIon objects
    """
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
