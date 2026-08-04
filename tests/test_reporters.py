"""Tests for isobaric reporter-ion extraction."""
import numpy as np
import pytest

from mzml_utils import (
    extract_reporters,
    extract_reporters_with_fallback,
    REPORTER_SETS,
    TMT_6PLEX,
    TMT_10PLEX,
    TMT_16PLEX,
)


class FakeSpectrum:
    def __init__(self, mz, intensity, ms_level=2):
        self.mz = np.asarray(mz, dtype=float)
        self.intensity = np.asarray(intensity, dtype=float)
        self.ms_level = ms_level


class FakeReader:
    """Minimal stand-in for the object open_spectra() returns."""

    def __init__(self, spectra):
        self._spectra = spectra

    def get_spectrum(self, scan):
        if scan not in self._spectra:
            raise KeyError(scan)
        return self._spectra[scan]


def _reporter_spectrum(values, channels=TMT_6PLEX):
    mz, inten = [], []
    for (_, theo), v in zip(channels.items(), values):
        mz.append(theo)
        inten.append(v)
    return FakeSpectrum(mz, inten)


def test_extracts_all_six_channels():
    vals = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0]
    reader = FakeReader({10: _reporter_spectrum(vals)})
    out = extract_reporters(reader, 10, "tmt6")
    assert list(out.values()) == vals
    assert list(out) == list(TMT_6PLEX)


def test_missing_channel_records_zero():
    # only the first channel present
    reader = FakeReader({1: FakeSpectrum([TMT_6PLEX["126"]], [42.0])})
    out = extract_reporters(reader, 1, "tmt6")
    assert out["126"] == 42.0
    assert all(out[c] == 0.0 for c in list(TMT_6PLEX)[1:])


def test_missing_value_is_configurable():
    reader = FakeReader({1: FakeSpectrum([TMT_6PLEX["126"]], [42.0])})
    out = extract_reporters(reader, 1, "tmt6", missing=float("nan"))
    assert out["126"] == 42.0
    assert np.isnan(out["127"])


def test_takes_max_within_window_not_sum():
    theo = TMT_6PLEX["126"]
    reader = FakeReader({1: FakeSpectrum([theo - 0.003, theo + 0.004], [10.0, 70.0])})
    out = extract_reporters(reader, 1, "tmt6")
    assert out["126"] == 70.0


def test_tolerance_excludes_far_peaks():
    theo = TMT_6PLEX["126"]
    reader = FakeReader({1: FakeSpectrum([theo + 0.05], [999.0])})
    assert extract_reporters(reader, 1, "tmt6")["126"] == 0.0


def test_ppm_unit():
    theo = TMT_6PLEX["126"]
    off = theo * 5e-6                       # 5 ppm away
    reader = FakeReader({1: FakeSpectrum([theo + off], [5.0])})
    assert extract_reporters(reader, 1, "tmt6", tolerance=20, unit="ppm")["126"] == 5.0
    assert extract_reporters(reader, 1, "tmt6", tolerance=1, unit="ppm")["126"] == 0.0


def test_empty_and_unreadable_scans_return_none():
    reader = FakeReader({1: FakeSpectrum([], [])})
    assert extract_reporters(reader, 1, "tmt6") is None
    assert extract_reporters(reader, 999, "tmt6") is None


def test_unknown_reporter_set_raises():
    reader = FakeReader({1: _reporter_spectrum([1, 2, 3, 4, 5, 6])})
    with pytest.raises(ValueError, match="unknown reporter set"):
        extract_reporters(reader, 1, "tmt7")


def test_custom_channel_map():
    custom = {"a": 126.127726, "b": 127.131081}
    reader = FakeReader({1: FakeSpectrum([126.127726, 127.131081], [3.0, 4.0])})
    assert extract_reporters(reader, 1, custom) == {"a": 3.0, "b": 4.0}


@pytest.mark.parametrize("name,n", [("tmt6", 6), ("tmt10", 10), ("tmt11", 11),
                                    ("tmt16", 16), ("tmt18", 18), ("itraq4", 4)])
def test_reporter_set_sizes(name, n):
    assert len(REPORTER_SETS[name]) == n


def test_plex_sets_are_nested_and_monotonic():
    assert set(TMT_10PLEX) < set(TMT_16PLEX)
    for s in (TMT_6PLEX, TMT_10PLEX, TMT_16PLEX):
        vals = list(s.values())
        assert vals == sorted(vals), "channel m/z should be ascending"


# --- the EThcD -> paired-HCD rescue -------------------------------------------------------


def test_fallback_used_when_primary_is_all_zero():
    reader = FakeReader({
        5: _reporter_spectrum([0.0] * 6),            # EThcD: no reporter signal
        4: _reporter_spectrum([1, 2, 3, 4, 5, 6]),   # paired HCD: good signal
    })
    out, src = extract_reporters_with_fallback(reader, 5, 4, "tmt6")
    assert src == "fallback"
    assert list(out.values()) == [1, 2, 3, 4, 5, 6]


def test_primary_kept_when_it_has_signal():
    reader = FakeReader({
        5: _reporter_spectrum([9, 9, 9, 9, 9, 9]),
        4: _reporter_spectrum([1, 2, 3, 4, 5, 6]),
    })
    out, src = extract_reporters_with_fallback(reader, 5, 4, "tmt6")
    assert src == "primary"
    assert list(out.values()) == [9] * 6


def test_no_fallback_scan_available():
    reader = FakeReader({5: _reporter_spectrum([0.0] * 6)})
    out, src = extract_reporters_with_fallback(reader, 5, None, "tmt6")
    assert src == "primary"
    assert all(v == 0.0 for v in out.values())


def test_both_empty_returns_none():
    reader = FakeReader({5: FakeSpectrum([], []), 4: FakeSpectrum([], [])})
    out, src = extract_reporters_with_fallback(reader, 5, 4, "tmt6")
    assert out is None and src == "none"
