"""
Unit and integration tests for the core analysis layer.
Run with:  python -m pytest tests/ -v
"""
import math
import numpy as np
import pytest

from src.data.amino_acids import AA_MASSES, PTM_DATABASE, PTM_TARGET_RESIDUES, WATER
from src.data.models import Spectrum, FragmentIon, Modification
from src.data.parsers import parse_fasta, generate_demo_spectrum, UBIQUITIN
from src.analysis.fragment_ions import calc_ions
from src.analysis.mass_utils import (
    calc_sequence_mass, mz_to_mass, mass_to_mz, ppm_error,
)
from src.analysis.peak_matching import match_peaks, sequence_coverage_pct
from src.analysis.proteoform_search import (
    run_targeted_search,
    run_database_search,
    _passes_mass_filter,
)


# ─────────────────────────────────────────────────────────────────────────────
# Fixtures
# ─────────────────────────────────────────────────────────────────────────────

UBIQ_10 = "MQIFVKTLTG"   # first 10 residues of ubiquitin

def _make_spectrum(mz_values, intensities=None, precursor_mz=0.0,
                   precursor_charge=0, precursor_mass=0.0):
    if intensities is None:
        intensities = [1000.0] * len(mz_values)
    return Spectrum(
        scan_id="test",
        mz_array=np.array(sorted(mz_values), dtype=float),
        intensity_array=np.array(intensities, dtype=float),
        precursor_mz=precursor_mz,
        precursor_charge=precursor_charge,
        precursor_mass=precursor_mass,
    )


# ─────────────────────────────────────────────────────────────────────────────
# amino_acids
# ─────────────────────────────────────────────────────────────────────────────

class TestAminoAcids:
    def test_all_standard_aas_present(self):
        for aa in "ACDEFGHIKLMNPQRSTVWY":
            assert aa in AA_MASSES, f"Missing standard AA: {aa}"

    def test_ptm_database_not_empty(self):
        assert len(PTM_DATABASE) >= 30, "PTM database too small"

    def test_ptm_has_target_residues_entry(self):
        for name in PTM_DATABASE:
            assert name in PTM_TARGET_RESIDUES, (
                f"PTM '{name}' missing from PTM_TARGET_RESIDUES"
            )

    def test_water_mass(self):
        assert abs(WATER - 18.01056) < 1e-4


# ─────────────────────────────────────────────────────────────────────────────
# mass_utils
# ─────────────────────────────────────────────────────────────────────────────

class TestMassUtils:
    def test_calc_sequence_mass_glycine(self):
        # Gly: 57.02146 + H2O = 75.03203
        mass = calc_sequence_mass("G")
        assert abs(mass - 75.03203) < 1e-4

    def test_calc_sequence_mass_dipeptide(self):
        # GG = 57.02146*2 + 18.01056 = 132.05348
        mass = calc_sequence_mass("GG")
        assert abs(mass - 132.05348) < 1e-4

    def test_calc_sequence_mass_with_mod(self):
        mod = Modification(position=1, name="Phosphorylation",
                           mass_shift=79.96633, residue="S")
        seq = "SG"
        base = calc_sequence_mass(seq)
        modded = calc_sequence_mass(seq, modifications=[mod])
        assert abs(modded - base - 79.96633) < 1e-4

    def test_mz_roundtrip(self):
        mass = 8565.0
        for z in [1, 4, 10, 25]:
            mz = mass_to_mz(mass, z)
            recovered = mz_to_mass(mz, z)
            assert abs(recovered - mass) < 1e-6, f"Roundtrip failed at z={z}"

    def test_ppm_error_zero_theoretical(self):
        # Should not raise
        assert ppm_error(100.0, 0.0) == 0.0

    def test_ppm_error_exact(self):
        assert ppm_error(100.0, 100.0) == 0.0

    def test_ppm_error_sign(self):
        assert ppm_error(100.001, 100.0) > 0
        assert ppm_error(99.999, 100.0) < 0


# ─────────────────────────────────────────────────────────────────────────────
# fragment_ions
# ─────────────────────────────────────────────────────────────────────────────

class TestFragmentIons:
    def test_b_ions_count(self):
        seq = UBIQ_10
        ions = calc_ions(seq, ['b'], max_charge=1)
        # b-series: positions 1 to n-1  (one per residue except last)
        b_ions = [i for i in ions if i.ion_type == 'b' and i.charge == 1]
        assert len(b_ions) == len(seq) - 1

    def test_y_ions_count(self):
        seq = UBIQ_10
        ions = calc_ions(seq, ['y'], max_charge=1)
        y_ions = [i for i in ions if i.ion_type == 'y' and i.charge == 1]
        assert len(y_ions) == len(seq) - 1

    def test_ion_types_respected(self):
        ions = calc_ions(UBIQ_10, ['b', 'y'], max_charge=1)
        types = {i.ion_type for i in ions}
        assert types == {'b', 'y'}

    def test_multiply_charged_ions(self):
        ions = calc_ions(UBIQ_10, ['b'], max_charge=3)
        charges = {i.charge for i in ions}
        assert charges == {1, 2, 3}

    def test_mod_shifts_ion_mass(self):
        seq = "ACDE"
        base_ions  = calc_ions(seq, ['b'], max_charge=1, mod_map={})
        mod_map    = {2: 79.96633}   # phospho on C (pos 2)
        mod_ions   = calc_ions(seq, ['b'], max_charge=1, mod_map=mod_map)
        # b1 (A only) should be unaffected; b2 (AC) should shift by ~80
        base_b1 = next(i for i in base_ions  if i.position == 1 and i.charge == 1)
        mod_b1  = next(i for i in mod_ions   if i.position == 1 and i.charge == 1)
        base_b2 = next(i for i in base_ions  if i.position == 2 and i.charge == 1)
        mod_b2  = next(i for i in mod_ions   if i.position == 2 and i.charge == 1)
        assert abs(base_b1.mz - mod_b1.mz) < 1e-5
        assert abs((mod_b2.mz - base_b2.mz) - 79.96633) < 1e-3

    def test_empty_sequence(self):
        # Should not raise, just return no ions
        ions = calc_ions("", ['b', 'y'], max_charge=2)
        assert ions == []


# ─────────────────────────────────────────────────────────────────────────────
# peak_matching
# ─────────────────────────────────────────────────────────────────────────────

class TestPeakMatching:
    def _perfect_spectrum(self, seq, ion_types=('b', 'y'), max_charge=1):
        """Return a spectrum containing every theoretical ion exactly."""
        ions = calc_ions(seq, list(ion_types), max_charge=max_charge)
        mzs  = [ion.mz for ion in ions]
        return _make_spectrum(mzs), ions

    def test_all_match_perfect_spectrum(self):
        seq = UBIQ_10
        spec, ions = self._perfect_spectrum(seq)
        matched = match_peaks(ions, spec, tolerance_ppm=5.0)
        n_matched = sum(1 for m in matched if m.matched)
        assert n_matched == len(ions)

    def test_no_match_empty_spectrum(self):
        ions = calc_ions(UBIQ_10, ['b', 'y'], max_charge=1)
        spec = _make_spectrum([])
        matched = match_peaks(ions, spec, tolerance_ppm=5.0)
        assert all(not m.matched for m in matched)

    def test_tolerance_respected(self):
        # Place a peak 100 ppm away from the first b-ion — should NOT match at 5 ppm
        ions = calc_ions(UBIQ_10, ['b'], max_charge=1)
        b1 = next(i for i in ions if i.position == 1 and i.charge == 1)
        far_mz = b1.mz * (1 + 100e-6)   # 100 ppm offset
        spec = _make_spectrum([far_mz])
        matched = match_peaks(ions, spec, tolerance_ppm=5.0)
        b1_matched = next(m for m in matched if m.position == 1 and m.charge == 1
                          and m.ion_type == 'b')
        assert not b1_matched.matched

    def test_ppm_error_recorded(self):
        ions = calc_ions(UBIQ_10, ['b'], max_charge=1)
        b1 = next(i for i in ions if i.position == 1 and i.charge == 1)
        offset_ppm = 3.0
        obs_mz = b1.mz * (1 + offset_ppm * 1e-6)
        spec = _make_spectrum([obs_mz])
        matched = match_peaks(ions, spec, tolerance_ppm=5.0)
        b1_m = next(m for m in matched if m.ion_type == 'b' and m.position == 1
                    and m.charge == 1)
        assert b1_m.matched
        assert abs(b1_m.mass_error_ppm - offset_ppm) < 0.1

    def test_sequence_coverage_full(self):
        seq = UBIQ_10
        spec, ions = self._perfect_spectrum(seq)
        matched = match_peaks(ions, spec, tolerance_ppm=5.0)
        cov = sequence_coverage_pct(seq, matched)
        assert cov == 100.0

    def test_sequence_coverage_zero(self):
        ions = calc_ions(UBIQ_10, ['b', 'y'], max_charge=1)
        spec = _make_spectrum([])
        matched = match_peaks(ions, spec, tolerance_ppm=5.0)
        cov = sequence_coverage_pct(UBIQ_10, matched)
        assert cov == 0.0

    def test_returns_same_count_as_input(self):
        ions = calc_ions(UBIQ_10, ['b', 'y'], max_charge=2)
        spec = _make_spectrum([999.0])
        matched = match_peaks(ions, spec, tolerance_ppm=5.0)
        assert len(matched) == len(ions)


# ─────────────────────────────────────────────────────────────────────────────
# parse_fasta
# ─────────────────────────────────────────────────────────────────────────────

class TestParseFasta:
    def test_single_entry(self):
        txt = ">sp|P1|gene Protein\nMQIFVKTLTGK\n"
        result = parse_fasta(txt)
        assert len(result) == 1
        assert result[0][0] == "sp|P1|gene"
        assert result[0][1] == "MQIFVKTLTGK"

    def test_multi_entry(self):
        txt = ">P1\nACDEF\n>P2\nMNPQR\n"
        result = parse_fasta(txt)
        assert len(result) == 2
        assert result[0][1] == "ACDEF"
        assert result[1][1] == "MNPQR"

    def test_multiline_sequence(self):
        txt = ">P1\nACDE\nFGHI\nKLMN\n"
        result = parse_fasta(txt)
        assert result[0][1] == "ACDEFGHIKLMN"

    def test_strips_numbers_and_whitespace(self):
        txt = ">P1\n10 MQIFVKTLTG 20\n"
        result = parse_fasta(txt)
        assert result[0][1] == "MQIFVKTLTG"

    def test_empty_string(self):
        assert parse_fasta("") == []

    def test_no_header(self):
        # Lines without '>' are ignored on their own
        assert parse_fasta("ACDEFGHIK\n") == []

    def test_uppercase(self):
        txt = ">P1\nacdef\n"
        result = parse_fasta(txt)
        assert result[0][1] == "ACDEF"


# ─────────────────────────────────────────────────────────────────────────────
# proteoform_search — targeted
# ─────────────────────────────────────────────────────────────────────────────

class TestTargetedSearch:
    @pytest.fixture(scope="class")
    def ubiq_results(self):
        spec = generate_demo_spectrum()
        return run_targeted_search(
            spec, UBIQUITIN, "Ubiquitin",
            ion_types=['b', 'y'], max_charge=3,
            tolerance_ppm=10, search_truncations=True,
            search_modifications=False,
        )

    def test_returns_results(self, ubiq_results):
        assert len(ubiq_results) > 0

    def test_top_hit_is_full_length(self, ubiq_results):
        top = ubiq_results[0]
        assert top.proteoform.start_pos == 1
        assert top.proteoform.end_pos == len(UBIQUITIN)

    def test_e_value_is_positive(self, ubiq_results):
        # e_value may round to 0.0 for extremely strong matches but must never be negative
        for r in ubiq_results:
            assert r.e_value >= 0

    def test_q_value_in_range(self, ubiq_results):
        for r in ubiq_results:
            assert 0.0 <= r.q_value <= 1.0

    def test_results_sorted_by_score_desc(self, ubiq_results):
        scores = [r.proteoform.score for r in ubiq_results]
        assert scores == sorted(scores, reverse=True)

    def test_no_decoys_in_output(self, ubiq_results):
        assert all(not r.is_decoy for r in ubiq_results)

    def test_coverage_in_range(self, ubiq_results):
        for r in ubiq_results:
            assert 0.0 <= r.sequence_coverage <= 100.0

    def test_empty_spectrum_returns_empty(self):
        spec = _make_spectrum([])
        results = run_targeted_search(spec, UBIQUITIN, "Ubiquitin",
                                      search_truncations=False,
                                      search_modifications=False)
        assert results == []


# ─────────────────────────────────────────────────────────────────────────────
# proteoform_search — database
# ─────────────────────────────────────────────────────────────────────────────

class TestDatabaseSearch:
    _PROTEINS = [
        ("Ubiquitin",  UBIQUITIN),
        ("Decoy_Prot", "ACDEFGHIKLMNPQRSTVWY" * 4),  # unlikely to match well
    ]

    @pytest.fixture(scope="class")
    def db_results(self):
        spec = generate_demo_spectrum()
        return run_database_search(
            spec, self._PROTEINS,
            ion_types=['b', 'y'], max_charge=3,
            tolerance_ppm=10, top_n=10,
        )

    def test_returns_results(self, db_results):
        assert len(db_results) > 0

    def test_top_hit_is_ubiquitin(self, db_results):
        assert "Ubiquitin" in db_results[0].proteoform.protein_name

    def test_no_decoys_in_output(self, db_results):
        assert all(not r.is_decoy for r in db_results)

    def test_e_values_positive(self, db_results):
        # e_value may round to 0.0 for extremely strong matches but must never be negative
        for r in db_results:
            assert r.e_value >= 0

    def test_q_values_in_range(self, db_results):
        for r in db_results:
            assert 0.0 <= r.q_value <= 1.0


# ─────────────────────────────────────────────────────────────────────────────
# mass pre-filter
# ─────────────────────────────────────────────────────────────────────────────

class TestMassFilter:
    def test_passes_exact_match(self):
        seq = "MQIFVKTLTG"
        mass = calc_sequence_mass(seq)
        assert _passes_mass_filter(seq, mass, tol_da=3.0)

    def test_passes_no_precursor(self):
        # obs_mass == 0 means unknown → always pass
        assert _passes_mass_filter("ACDEFG", 0.0, tol_da=3.0)

    def test_rejects_wrong_protein(self):
        # Ubiquitin ~8565 Da; obs_mass far away should fail
        assert not _passes_mass_filter(UBIQUITIN, 100.0, tol_da=3.0)

    def test_passes_within_truncation_window(self):
        # The full sequence passes even if obs_mass is a few residues lighter
        seq = "MQIFVKTLTGKTITLEVEPSDTIENVKAK"   # 29 aa ~3200 Da
        mass = calc_sequence_mass(seq)
        # Obs mass matches a 4-residue N-terminal truncation → should still pass
        trunc_mass = calc_sequence_mass(seq[4:])
        assert _passes_mass_filter(seq, trunc_mass, tol_da=3.0)
