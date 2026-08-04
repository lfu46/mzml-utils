"""Canonical match layer.

Unifies the theoretical-to-experimental ion-matching policy used by validation
and annotation callers. Each caller historically reproduced the same envelope
around :func:`mzml_utils.match_peaks`: (deisotope decision → primary match →
optional raw fallback → optional charge-reduced exclusion). Centralising the
envelope keeps validator/annotator paths bit-identical and makes divergences
traceable to a single helper.

Invariant: with default flags (``raw_fallback=False``, ``exclude_mzs=None``),
``canonical_match(ions, ctx)`` is byte-identical to
``match_peaks(ions, ctx.exp_mz, ctx.exp_intensity, ctx.tolerance_ppm,
match_isotopes=not ctx.deisotope_applied)``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Optional

import numpy as np

from .deisotope import deisotope
from .fragments import match_peaks as _default_match_peaks
from .fragments import TheoreticalIon, MatchedIon


@dataclass
class MatchContext:
    """Carries both the processed (possibly deisotoped) spectrum and the raw
    spectrum so callers can opt into raw-fallback matching per ion set.

    ``exp_mz/exp_intensity`` — what ``match_peaks`` sees on the primary pass.
    ``raw_mz/raw_intensity`` — always the original arrays (pre-deisotope). When
        ``deisotope_applied`` is False these equal the exp arrays by reference.
    """
    exp_mz: np.ndarray
    exp_intensity: np.ndarray
    raw_mz: np.ndarray
    raw_intensity: np.ndarray
    tolerance_ppm: float
    deisotope_applied: bool


def from_spectrum(
    spectrum,
    tolerance_ppm: float = 20.0,
    deisotope_spectrum: bool = True,
    max_charge: Optional[int] = None,
) -> MatchContext:
    """Build a MatchContext from a ``mzml_utils.Spectrum``.

    If ``deisotope_spectrum`` is True, ``exp_mz/exp_intensity`` are the
    deisotoped arrays; ``raw_mz/raw_intensity`` retain the original.
    """
    raw_mz = spectrum.mz
    raw_int = spectrum.intensity
    if deisotope_spectrum and len(raw_mz) > 0:
        exp_mz, exp_int = deisotope(raw_mz, raw_int, max_charge=max_charge or 4)
    else:
        exp_mz, exp_int = raw_mz, raw_int
    return MatchContext(
        exp_mz=exp_mz, exp_intensity=exp_int,
        raw_mz=raw_mz, raw_intensity=raw_int,
        tolerance_ppm=tolerance_ppm,
        deisotope_applied=bool(deisotope_spectrum and len(raw_mz) > 0),
    )


def from_arrays(
    exp_mz: np.ndarray,
    exp_intensity: np.ndarray,
    tolerance_ppm: float = 20.0,
    raw_mz: Optional[np.ndarray] = None,
    raw_intensity: Optional[np.ndarray] = None,
    deisotope_applied: bool = False,
) -> MatchContext:
    """Build a MatchContext directly from arrays (for callers that have already
    deisotoped or are providing processed+raw pairs themselves, e.g., the
    spectrum annotator).
    """
    if raw_mz is None:
        raw_mz = exp_mz
    if raw_intensity is None:
        raw_intensity = exp_intensity
    return MatchContext(
        exp_mz=exp_mz, exp_intensity=exp_intensity,
        raw_mz=raw_mz, raw_intensity=raw_intensity,
        tolerance_ppm=tolerance_ppm,
        deisotope_applied=deisotope_applied,
    )


def canonical_match(
    ions: List[TheoreticalIon],
    ctx: MatchContext,
    *,
    raw_fallback: bool = False,
    exclude_mzs: Optional[List[float]] = None,
    exclude_tol_ppm: float = 15.0,
    match_fn: Optional[Callable] = None,
) -> List[MatchedIon]:
    """Canonical theoretical→experimental match.

    Steps (in order):
      1. Primary match against ctx.exp_mz/exp_intensity via ``match_fn``
         (defaults to ``mzml_utils.fragments.match_peaks``). When
         ``ctx.deisotope_applied`` is True, ``match_isotopes=False`` is passed
         (deisotoped peaks already represent monoisotopic only); else
         ``match_isotopes=True`` (raw spectrum; let the matcher walk isotopes).
      2. If ``raw_fallback`` is True, a second pass against the raw spectrum
         catches ions whose monoisotopic peak was removed by the deisotoper.
         De-duped by (ion_type, ion_number, charge) against primary matches.
      3. If ``exclude_mzs`` is provided, any match whose observed m/z is
         within ``exclude_tol_ppm`` of an exclusion target is dropped.

    ``match_fn`` lets callers plug in a project-local ``match_peaks``
    implementation (e.g., the annotator's extended ``MatchedIon`` subclass
    with ``has_modification``). The function must accept the same signature
    as :func:`mzml_utils.fragments.match_peaks`.

    Returns a list of ``MatchedIon`` (or caller-specific subclass) in the
    same order the underlying matcher produces.
    """
    if len(ions) == 0 or len(ctx.exp_mz) == 0:
        return []

    mp = match_fn or _default_match_peaks

    primary = mp(
        ions, ctx.exp_mz, ctx.exp_intensity,
        tolerance_ppm=ctx.tolerance_ppm,
        match_isotopes=not ctx.deisotope_applied,
    )

    if raw_fallback and ctx.deisotope_applied:
        matched_keys = {(m.ion_type, m.ion_number, m.charge) for m in primary}
        unmatched_theo = [
            t for t in ions
            if (t.ion_type, t.ion_number, t.charge) not in matched_keys
        ]
        if unmatched_theo:
            raw_hits = mp(
                unmatched_theo, ctx.raw_mz, ctx.raw_intensity,
                tolerance_ppm=ctx.tolerance_ppm,
                match_isotopes=False,
            )
            primary.extend(raw_hits)

    if exclude_mzs:
        kept: List[MatchedIon] = []
        for m in primary:
            drop = False
            for cr_mz in exclude_mzs:
                if cr_mz <= 0:
                    continue
                if abs(m.mz - cr_mz) / cr_mz * 1e6 < exclude_tol_ppm:
                    drop = True
                    break
            if not drop:
                kept.append(m)
        return kept

    return primary
