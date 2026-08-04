"""Isobaric (TMT / iTRAQ) reporter-ion extraction.

Promoted from OGlycoTM/data_analysis/fix_ogalnac_site_quant.py, which needed it to rescue
quantification for EThcD PSMs: EThcD reporter yield is poor, so the reporters are taken from
the paired HCD scan of the same precursor instead. That rescue is real methodology and had no
equivalent anywhere in the shared libraries -- only a private copy in one project's script.

Pairing itself lives next door in `pairing.pair_hcd_ethcd`; this module only reads reporters
out of whatever scan you hand it.
"""
from __future__ import annotations

from typing import Dict, Iterable, Optional, Sequence

import numpy as np

# Monoisotopic reporter m/z. TMT values are the standard Thermo channel masses; the 6-plex set
# is the subset used by OGlycoTM (126, 127N, 128C, 129N, 130C, 131N).
TMT_6PLEX = {
    "126":  126.127726,
    "127":  127.131081,
    "128":  128.134436,
    "129":  129.137790,
    "130":  130.141145,
    "131":  131.138180,
}

TMT_10PLEX = {
    "126":  126.127726, "127N": 127.124761, "127C": 127.131081, "128N": 128.128116,
    "128C": 128.134436, "129N": 129.131471, "129C": 129.137790, "130N": 130.134825,
    "130C": 130.141145, "131N": 131.138180,
}

TMT_11PLEX = {**TMT_10PLEX, "131C": 131.144499}

TMT_16PLEX = {
    **TMT_11PLEX,
    "132N": 132.141535, "132C": 132.147855, "133N": 133.144890,
    "133C": 133.151210, "134N": 134.148245,
}

TMT_18PLEX = {**TMT_16PLEX, "134C": 134.154565, "135N": 135.151600}

ITRAQ_4PLEX = {"114": 114.111228, "115": 115.108263, "116": 116.111618, "117": 117.114973}

REPORTER_SETS = {
    "tmt6": TMT_6PLEX, "tmt10": TMT_10PLEX, "tmt11": TMT_11PLEX,
    "tmt16": TMT_16PLEX, "tmt18": TMT_18PLEX, "itraq4": ITRAQ_4PLEX,
}

DEFAULT_TOL_DA = 0.01  # 10 mDa, the tolerance the OGlycoTM rescue used


def extract_reporters(
    source,
    scan_num: int,
    reporters: Dict[str, float] | str = "tmt6",
    *,
    tolerance: float = DEFAULT_TOL_DA,
    unit: str = "Da",
    missing: float = 0.0,
) -> Optional[Dict[str, float]]:
    """Reporter-ion intensities from one scan.

    Takes the most intense peak inside the window per channel, matching the original
    implementation -- with a 10 mDa window the reporter region is not usually crowded, but a
    max is safer than a sum if a neighbouring channel's isotope leaks in.

    Args:
        source:     an open reader (anything with ``.get_spectrum(scan)``) -- e.g. the object
                    returned by ``mzml_utils.open_spectra(path)``.
        scan_num:   scan number to read.
        reporters:  channel name -> theoretical m/z, or a key of ``REPORTER_SETS``.
        tolerance:  match window, in ``unit``.
        unit:       ``"Da"`` or ``"ppm"``.
        missing:    value recorded for a channel with no peak in the window.

    Returns:
        ``{channel: intensity}``, or ``None`` if the scan is unreadable or empty.
    """
    if isinstance(reporters, str):
        try:
            reporters = REPORTER_SETS[reporters]
        except KeyError:
            raise ValueError(
                f"unknown reporter set {reporters!r}; choose from {sorted(REPORTER_SETS)}"
            ) from None

    try:
        spec = source.get_spectrum(scan_num)
    except Exception:
        return None
    if spec is None:
        return None

    mz = np.asarray(spec.mz)
    intensity = np.asarray(spec.intensity)
    if mz.size == 0:
        return None

    out: Dict[str, float] = {}
    for channel, theo_mz in reporters.items():
        win = tolerance * theo_mz / 1e6 if unit.lower() == "ppm" else tolerance
        mask = np.abs(mz - theo_mz) < win
        out[channel] = float(intensity[mask].max()) if mask.any() else missing
    return out


def extract_reporters_with_fallback(
    source,
    scan_num: int,
    fallback_scan: Optional[int],
    reporters: Dict[str, float] | str = "tmt6",
    **kwargs,
) -> tuple[Optional[Dict[str, float]], str]:
    """Reporters from ``scan_num``, falling back to ``fallback_scan`` when they are unusable.

    The EThcD rescue: ETD-type scans often carry little or no reporter signal, so quantification
    is taken from the paired HCD scan of the same precursor (find it with
    ``pairing.pair_hcd_ethcd``). Returns the values and which scan they came from, so the
    provenance can be recorded rather than silently lost.

    Returns:
        ``(reporters_or_None, source_label)`` where ``source_label`` is ``"primary"``,
        ``"fallback"`` or ``"none"``.
    """
    primary = extract_reporters(source, scan_num, reporters, **kwargs)
    if primary and any(v > 0 for v in primary.values()):
        return primary, "primary"

    if fallback_scan is not None:
        secondary = extract_reporters(source, fallback_scan, reporters, **kwargs)
        if secondary and any(v > 0 for v in secondary.values()):
            return secondary, "fallback"

    return (primary if primary else None), ("primary" if primary else "none")
