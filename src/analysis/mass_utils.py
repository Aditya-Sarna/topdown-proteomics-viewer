from typing import List, Dict, Optional, Tuple
import math
import numpy as np

from ..data.amino_acids import AA_MASSES, WATER, PROTON, PTM_DATABASE
from ..data.models import Modification


def calc_sequence_mass(sequence: str,
                        modifications: Optional[List[Modification]] = None,
                        nterm_mod: float = 0.0,
                        cterm_mod: float = 0.0) -> float:
    """
    Monoisotopic neutral mass of a proteoform.
    Uses pyopenms.AASequence when available for higher accuracy.
    """
    # pyopenms temporarily disabled
    # try:
    #     import pyopenms
    #     aa_seq = pyopenms.AASequence.fromString(sequence)
    #     mass = aa_seq.getMonoWeight()   # already includes H2O terminal
    #     mass += nterm_mod + cterm_mod
    #     if modifications:
    #         for mod in modifications:
    #             mass += mod.mass_shift
    #     return mass
    # except Exception:
    #     pass  # fall back to manual calculation

    mass = WATER + nterm_mod + cterm_mod
    for aa in sequence:
        mass += AA_MASSES.get(aa, 0.0)
    if modifications:
        for mod in modifications:
            mass += mod.mass_shift
    return mass


def mass_to_mz(mass: float, charge: int) -> float:
    return (mass + charge * PROTON) / charge


def mz_to_mass(mz: float, charge: int) -> float:
    return mz * charge - charge * PROTON


def ppm_error(observed: float, theoretical: float) -> float:
    if theoretical == 0:
        return 0.0
    return (observed - theoretical) / theoretical * 1e6


def da_error(observed: float, theoretical: float) -> float:
    return observed - theoretical


def possible_charges(mass: float,
                     min_mz: float = 200.0,
                     max_mz: float = 2000.0) -> List[int]:
    charges = []
    for z in range(1, 200):
        mz = mass_to_mz(mass, z)
        if mz < min_mz:
            break
        if mz <= max_mz:
            charges.append(z)
    return charges


def suggest_modifications(mass_diff: float,
                           tolerance_da: float = 0.05) -> List[Dict]:
    """Return PTMs whose mass shift is within tolerance_da of mass_diff."""
    suggestions = []
    for name, shift in PTM_DATABASE.items():
        if abs(abs(mass_diff) - abs(shift)) <= tolerance_da:
            suggestions.append({
                'name': name,
                'mass_shift': round(shift, 5),
                'residual_da': round(mass_diff - shift, 5),
            })
    # Two-modification combinations
    items = list(PTM_DATABASE.items())
    for i in range(len(items)):
        for j in range(i, len(items)):
            combined = items[i][1] + items[j][1]
            if abs(mass_diff - combined) <= tolerance_da:
                suggestions.append({
                    'name': f"{items[i][0]} + {items[j][0]}",
                    'mass_shift': round(combined, 5),
                    'residual_da': round(mass_diff - combined, 5),
                })
    suggestions.sort(key=lambda x: abs(x['residual_da']))
    return suggestions[:12]


def generate_isotope_envelope(mass: float,
                               charge: int,
                               n_isotopes: int = 5) -> Tuple[np.ndarray, np.ndarray]:
    """Approximate isotope envelope using the averagine model."""
    mz_base    = mass_to_mz(mass, charge)
    mz_spacing = 1.003355 / charge
    mzs        = np.array([mz_base + i * mz_spacing for i in range(n_isotopes)])

    # Simplified Poisson distribution (averagine ~ 1 C per 112 Da)
    lam = (mass / 112.0) * 0.0107
    intensities = np.array([
        np.exp(-lam) * (lam ** k) / max(math.factorial(k), 1)
        for k in range(n_isotopes)
    ], dtype=float)
    if intensities.max() > 0:
        intensities = intensities / intensities.max() * 100.0
    return mzs, intensities
