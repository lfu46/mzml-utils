"""Tests for the ions module."""

import numpy as np
import pytest

from mzml_utils.ions import (
    ppm_error,
    da_error,
    within_tolerance,
    find_ion,
    search_ions,
    ion_rank,
)


class TestPpmError:
    def test_zero_error(self):
        assert ppm_error(500.0, 500.0) == pytest.approx(0.0)

    def test_positive_error(self):
        # 500.001 vs 500.0 => 2 ppm
        assert ppm_error(500.001, 500.0) == pytest.approx(2.0)

    def test_negative_error(self):
        assert ppm_error(499.999, 500.0) == pytest.approx(-2.0)


class TestDaError:
    def test_zero(self):
        assert da_error(204.0864, 204.0864) == pytest.approx(0.0)

    def test_positive(self):
        assert da_error(204.09, 204.0864) == pytest.approx(0.0036, abs=1e-4)


class TestWithinTolerance:
    def test_ppm_within(self):
        assert within_tolerance(500.005, 500.0, 20, 'ppm') is True

    def test_ppm_outside(self):
        assert within_tolerance(500.1, 500.0, 20, 'ppm') is False

    def test_da_within(self):
        assert within_tolerance(204.1, 204.0864, 0.1, 'Da') is True

    def test_da_outside(self):
        assert within_tolerance(204.3, 204.0864, 0.1, 'Da') is False

    def test_invalid_unit(self):
        with pytest.raises(ValueError):
            within_tolerance(500.0, 500.0, 10, 'amu')


class TestFindIon:
    @pytest.fixture
    def spectrum(self):
        mz = np.array([100.0, 204.087, 300.0, 500.0])
        intensity = np.array([1000.0, 5000.0, 2000.0, 3000.0])
        return mz, intensity

    def test_found(self, spectrum):
        mz, ints = spectrum
        result = find_ion(mz, ints, 204.0864, tolerance=0.1, unit='Da')
        assert result is not None
        assert result['mz'] == pytest.approx(204.087)
        assert result['intensity'] == pytest.approx(5000.0)
        assert abs(result['ppm_error']) < 50

    def test_not_found(self, spectrum):
        mz, ints = spectrum
        result = find_ion(mz, ints, 600.0, tolerance=0.1, unit='Da')
        assert result is None

    def test_empty_spectrum(self):
        result = find_ion(np.array([]), np.array([]), 204.0864, 0.1, 'Da')
        assert result is None

    def test_ppm_tolerance(self, spectrum):
        mz, ints = spectrum
        result = find_ion(mz, ints, 204.0864, tolerance=20, unit='ppm')
        assert result is not None

    def test_ppm_tolerance_too_tight(self, spectrum):
        mz, ints = spectrum
        # 204.087 vs 204.0864: ~2.9 ppm, so 1 ppm should miss
        result = find_ion(mz, ints, 204.0864, tolerance=1, unit='ppm')
        assert result is None


class TestSearchIons:
    def test_multiple_ions(self):
        mz = np.array([138.055, 144.065, 168.065, 186.076, 204.086, 500.0])
        ints = np.array([100, 200, 300, 400, 500, 1000])
        ion_list = {
            'HexNAc': 204.0864,
            'HexNAc-H2O': 186.0760,
            'missing': 700.0,
        }
        results = search_ions(mz, ints, ion_list, tolerance=0.1, unit='Da')
        assert results['HexNAc'] is not None
        assert results['HexNAc-H2O'] is not None
        assert results['missing'] is None

    def test_rel_threshold(self):
        mz = np.array([204.086, 500.0])
        ints = np.array([10.0, 1000.0])
        results = search_ions(mz, ints, {'HexNAc': 204.0864},
                              tolerance=0.1, unit='Da', rel_threshold=5.0)
        assert results['HexNAc'] is not None
        assert results['HexNAc']['above_threshold'] is False
        assert results['HexNAc']['rel_intensity'] == pytest.approx(1.0)


class TestIonRank:
    def test_rank(self):
        mz = np.array([100.0, 204.086, 300.0, 500.0])
        ints = np.array([1000.0, 5000.0, 2000.0, 3000.0])
        rank = ion_rank(mz, ints, 204.0864, tolerance=0.1, unit='Da')
        assert rank == 1  # most intense

    def test_not_found(self):
        mz = np.array([100.0, 200.0])
        ints = np.array([1000.0, 2000.0])
        assert ion_rank(mz, ints, 500.0, tolerance=0.1, unit='Da') is None
