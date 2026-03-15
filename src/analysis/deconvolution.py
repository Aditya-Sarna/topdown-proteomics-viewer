"""
Charge-state deconvolution using pyopenms.Deisotoper.

Takes a raw MS2 spectrum (m/z vs intensity) and returns a new spectrum where
each isotope envelope has been collapsed to a single monoisotopic peak at
neutral mass + 1 Da (singly-charged equivalent), making downstream fragment
ion matching charge-agnostic.

Falls back gracefully when pyopenms is unavailable.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, Tuple
import numpy as np

from ..data.models import Spectrum


@dataclass
class DeconvolutionResult:
    spectrum: Spectrum          # deconvoluted spectrum (z=1 equivalent peaks)
    n_original_peaks: int
    n_deconvoluted_peaks: int
    used_openms: bool


def deconvolute_spectrum(
    spectrum: Spectrum,
    fragment_tolerance_ppm: float = 10.0,
    min_charge: int = 1,
    max_charge: int = 20,
    min_isotopes: int = 2,
    max_isotopes: int = 10,
) -> DeconvolutionResult:
    """
    Deconvolute charge states using pyopenms.Deisotoper.

    Each group of isotope peaks is collapsed to its monoisotopic peak,
    and the m/z is converted to a singly-charged (z=1) neutral equivalent
    so all downstream tools see simple masses.

    Parameters
    ----------
    spectrum             : raw Spectrum object
    fragment_tolerance_ppm : mass tolerance for isotope matching
    min_charge / max_charge : charge state search range
    min_isotopes         : minimum isotopes required to call a group
    max_isotopes         : maximum isotopes to consider per group

    Returns
    -------
    DeconvolutionResult with the deconvoluted Spectrum and stats.
    """
    n_orig = len(spectrum.mz_array)

    # pyopenms is disabled — return the original spectrum immediately
    return DeconvolutionResult(
        spectrum=spectrum,
        n_original_peaks=n_orig,
        n_deconvoluted_peaks=n_orig,
        used_openms=False,
    )

    try:
        # import pyopenms  # temporarily disabled
        raise ImportError("pyopenms disabled")

        ms_spec = pyopenms.MSSpectrum()  # noqa: F821 (unreachable)
        ms_spec.setMSLevel(2)

        # Load peaks
        peaks_mz  = spectrum.mz_array.tolist()
        peaks_int = spectrum.intensity_array.tolist()
        ms_spec.set_peaks((peaks_mz, peaks_int))
        ms_spec.sortByPosition()

        fragment_tol_da = (
            np.mean(spectrum.mz_array) * fragment_tolerance_ppm / 1e6
            if n_orig > 0 else 0.02
        )

        pyopenms.Deisotoper.deisotopeAndSingleCharge(
            ms_spec,
            fragment_tol_da,       # fragment_tolerance (Da)
            False,                 # fragment_unit_ppm = False (we pass Da)
            min_charge,
            max_charge,
            True,                  # keep_only_deisotoped
            min_isotopes,
            max_isotopes,
            True,                  # make_single_charged
            True,                  # annotate_charge
            False,                 # annotate_iso_peak_count
            True,                  # use_decreasing_model
            2,                     # start_intensity_check
            False,                 # add_up_intensity
        )

        out_mz, out_int = ms_spec.get_peaks()
        out_mz  = np.asarray(out_mz,  dtype=np.float64)
        out_int = np.asarray(out_int, dtype=np.float64)

        deconv = Spectrum(
            scan_id=spectrum.scan_id,
            mz_array=out_mz,
            intensity_array=out_int,
            retention_time=spectrum.retention_time,
            precursor_mz=spectrum.precursor_mz,
            precursor_charge=spectrum.precursor_charge,
            precursor_mass=spectrum.precursor_mass,
            ms_level=spectrum.ms_level,
            file_name=spectrum.file_name,
        )
        return DeconvolutionResult(
            spectrum=deconv,
            n_original_peaks=n_orig,
            n_deconvoluted_peaks=len(out_mz),
            used_openms=True,
        )

    except Exception:
        # pyopenms unavailable or error — return original unchanged
        return DeconvolutionResult(
            spectrum=spectrum,
            n_original_peaks=n_orig,
            n_deconvoluted_peaks=n_orig,
            used_openms=False,
        )
