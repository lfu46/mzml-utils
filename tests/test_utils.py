"""Tests for the utils module."""

import os
import tempfile
import pytest

from mzml_utils.utils import parse_spectrum_id, find_mzml_file, parse_modifications


class TestParseSpectrumId:
    def test_normal(self):
        fname, scan, charge = parse_spectrum_id("my_experiment.12345.12345.3")
        assert fname == "my_experiment"
        assert scan == 12345
        assert charge == 3

    def test_dots_in_name(self):
        fname, scan, charge = parse_spectrum_id("sample.01.raw.999.999.2")
        assert fname == "sample.01.raw"
        assert scan == 999
        assert charge == 2

    def test_invalid(self):
        assert parse_spectrum_id("bad_string") == (None, None, None)

    def test_empty(self):
        assert parse_spectrum_id("") == (None, None, None)


class TestFindMzmlFile:
    def test_finds_calibrated(self, tmp_path):
        # Create a fake calibrated mzML
        (tmp_path / "sample_calibrated.mzML").write_text("")
        result = find_mzml_file("sample", str(tmp_path))
        assert result is not None
        assert result.endswith("sample_calibrated.mzML")

    def test_finds_raw(self, tmp_path):
        (tmp_path / "sample.mzML").write_text("")
        result = find_mzml_file("sample", str(tmp_path))
        assert result is not None

    def test_not_found(self, tmp_path):
        result = find_mzml_file("nonexistent", str(tmp_path))
        assert result is None

    def test_custom_patterns(self, tmp_path):
        (tmp_path / "sample_custom.mzML").write_text("")
        result = find_mzml_file("sample", str(tmp_path),
                                patterns=["{name}_custom.mzML"])
        assert result is not None


class TestParseModifications:
    def test_nterm_and_residue(self):
        mods = parse_modifications("N-term(229.1629),4S(528.2859),19K(229.1629)")
        assert len(mods) == 3
        assert mods[0] == {'position': 0, 'residue': 'N-term', 'mass': 229.1629}
        assert mods[1] == {'position': 4, 'residue': 'S', 'mass': 528.2859}
        assert mods[2] == {'position': 19, 'residue': 'K', 'mass': 229.1629}

    def test_empty(self):
        assert parse_modifications("") == []
        assert parse_modifications("nan") == []
        assert parse_modifications(None) == []

    def test_cterm(self):
        mods = parse_modifications("C-term(100.0)")
        assert len(mods) == 1
        assert mods[0]['position'] == -1

    def test_single(self):
        mods = parse_modifications("5T(203.0794)")
        assert len(mods) == 1
        assert mods[0] == {'position': 5, 'residue': 'T', 'mass': 203.0794}
