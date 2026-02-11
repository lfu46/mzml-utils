"""Tests for the fragments module."""

import numpy as np
import pytest

from mzml_utils.fragments import (
    FragmentCalculator,
    TheoreticalIon,
    MatchedIon,
    match_peaks,
    calculate_false_match_rate,
    calculate_annotation_statistics,
)
from mzml_utils.constants import PROTON, H2O, AA_MASSES


@pytest.fixture
def simple_calc():
    """FragmentCalculator for a simple unmodified peptide: ACDK (charge 2)."""
    return FragmentCalculator("ACDK", [], precursor_charge=2)


@pytest.fixture
def glyco_calc():
    """FragmentCalculator for TMT-labelled O-GlcNAc glycopeptide."""
    mods = [
        {'position': 0, 'residue': 'N-term', 'mass': 229.1629},
        {'position': 4, 'residue': 'S', 'mass': 528.2859},
    ]
    return FragmentCalculator("AGYSQGATQYTQAQQTR", mods, precursor_charge=3)


class TestFragmentCalculatorInit:
    def test_length(self, simple_calc):
        assert simple_calc.length == 4

    def test_precursor_mass(self, simple_calc):
        expected = AA_MASSES['A'] + AA_MASSES['C'] + AA_MASSES['D'] + AA_MASSES['K'] + H2O
        assert simple_calc.precursor_mass == pytest.approx(expected, abs=0.01)

    def test_precursor_mz(self, simple_calc):
        expected = (simple_calc.precursor_mass + 2 * PROTON) / 2
        assert simple_calc.precursor_mz == pytest.approx(expected, abs=0.001)


class TestBYIons:
    def test_b_ion_count(self, simple_calc):
        b_ions = simple_calc.calculate_b_ions(charges=[1])
        assert len(b_ions) == 3  # b1, b2, b3

    def test_y_ion_count(self, simple_calc):
        y_ions = simple_calc.calculate_y_ions(charges=[1])
        assert len(y_ions) == 3  # y1, y2, y3

    def test_b1_mass(self, simple_calc):
        b_ions = simple_calc.calculate_b_ions(charges=[1])
        b1 = [i for i in b_ions if i.ion_number == 1][0]
        expected = AA_MASSES['A'] + PROTON
        assert b1.mz == pytest.approx(expected, abs=0.01)

    def test_y1_mass(self, simple_calc):
        y_ions = simple_calc.calculate_y_ions(charges=[1])
        y1 = [i for i in y_ions if i.ion_number == 1][0]
        expected = AA_MASSES['K'] + H2O + PROTON
        assert y1.mz == pytest.approx(expected, abs=0.01)

    def test_charge_states(self, simple_calc):
        b_ions = simple_calc.calculate_b_ions(charges=[1, 2])
        # 3 positions * 2 charges = 6
        assert len(b_ions) == 6


class TestCZIons:
    def test_c_ion_count(self, simple_calc):
        c_ions = simple_calc.calculate_c_ions(charges=[1])
        assert len(c_ions) == 3

    def test_z_ion_count(self, simple_calc):
        z_ions = simple_calc.calculate_z_ions(charges=[1])
        assert len(z_ions) == 3


class TestGlycoIons:
    def test_Y_ions(self, glyco_calc):
        y_ions = glyco_calc.calculate_Y_ions(charges=[1])
        assert len(y_ions) > 0
        annotations = [i.annotation for i in y_ions]
        assert any('Y0' in a for a in annotations)
        assert any('Y1' in a for a in annotations)

    def test_oxonium_ions(self, glyco_calc):
        ox = glyco_calc.calculate_oxonium_ions()
        assert len(ox) == 7
        assert all(i.charge == 1 for i in ox)


class TestCalculateAllIons:
    def test_keys(self, glyco_calc):
        result = glyco_calc.calculate_all_ions()
        assert 'b' in result
        assert 'y' in result
        assert 'c' in result
        assert 'z' in result
        assert 'Y' in result
        assert 'oxonium' in result
        assert 'b_NL' in result

    def test_flat(self, glyco_calc):
        flat = glyco_calc.get_all_ions_flat()
        assert len(flat) > 50  # should be many ions


class TestMatchPeaks:
    def test_match(self, simple_calc):
        b_ions = simple_calc.calculate_b_ions(charges=[1])
        # Create fake spectrum with exact theoretical m/z
        exp_mz = np.array([b.mz for b in b_ions])
        exp_int = np.array([1000.0] * len(b_ions))

        matched = match_peaks(b_ions, exp_mz, exp_int, tolerance_ppm=20.0,
                              match_isotopes=False)
        assert len(matched) == len(b_ions)
        for m in matched:
            assert abs(m.mass_error_ppm) < 1.0

    def test_no_match(self, simple_calc):
        b_ions = simple_calc.calculate_b_ions(charges=[1])
        exp_mz = np.array([1.0, 2.0, 3.0])
        exp_int = np.array([100.0, 200.0, 300.0])
        matched = match_peaks(b_ions, exp_mz, exp_int, tolerance_ppm=20.0,
                              match_isotopes=False)
        assert len(matched) == 0


class TestFalseMatchRate:
    def test_perfect_spectrum(self, simple_calc):
        ions = simple_calc.get_all_ions_flat()
        exp_mz = np.array([i.mz for i in ions])
        exp_int = np.array([1000.0] * len(ions))
        fmr = calculate_false_match_rate(ions, exp_mz, exp_int, tolerance_ppm=20.0)
        # Perfect match should have low FMR
        assert fmr.fmr_peaks < 0.5
        assert fmr.matched_peaks > 0

    def test_empty(self, simple_calc):
        ions = simple_calc.get_all_ions_flat()
        fmr = calculate_false_match_rate(ions, np.array([]), np.array([]))
        assert fmr.matched_peaks == 0
        assert fmr.fmr_peaks == 0.0


class TestAnnotationStatistics:
    def test_full_coverage(self, simple_calc):
        b_ions = simple_calc.calculate_b_ions(charges=[1])
        y_ions = simple_calc.calculate_y_ions(charges=[1])
        all_theo = b_ions + y_ions

        exp_mz = np.array([i.mz for i in all_theo])
        exp_int = np.array([1000.0] * len(all_theo))

        matched = match_peaks(all_theo, exp_mz, exp_int, tolerance_ppm=20.0,
                              match_isotopes=False)
        stats = calculate_annotation_statistics(
            matched, all_theo, exp_mz, exp_int, simple_calc.length)

        assert stats['sequence_coverage'] == pytest.approx(1.0)
        assert stats['fragments_found'] == pytest.approx(1.0)
