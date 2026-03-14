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
        import pyopenms
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
