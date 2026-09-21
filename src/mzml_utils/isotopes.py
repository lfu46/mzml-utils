"""
Exact-composition isotope distributions.

``deisotope`` models an isotope envelope with a one-parameter Poisson averagine,
which is right for an unknown CHNOS peptide and cannot express anything else: a
halogen tag (the 79Br/81Br doublet puts M+2 above M), selenium, a heavy-isotope
label, or a glycan-heavy composition. This module computes the envelope from
the elemental composition itself.

Peaks are aggregated by nominal mass, which is what a peptide-scale envelope
looks like at ordinary resolving power: the 81Br and the 13C2 species of M+2 are
one peak. Each aggregated peak carries its abundance-weighted centroid mass, so
the spacing between peaks is the real one (M+1 of a brominated peptide is the
13C peak at +1.00335, not a fictitious 80Br at +1.00019).

Element isotope masses and abundances come from ``pyteomics.mass.nist_mass``;
nothing is typed in here.

Typical use::

    from pyteomics import mass
    from mzml_utils import isotope_distribution

    comp = mass.Composition(sequence="TPENFPSK") + mass.Composition({"Br": 1, "H": -1})
    dist = isotope_distribution(comp, n_peaks=3, charge=2)
    dist.abundance      # array([0.383, 0.188, 0.429])
    dist.mz             # centroid m/z of M, M+1, M+2 at 2+
"""

from __future__ import annotations

import re
import numpy as np
from dataclasses import dataclass
from typing import Mapping, Tuple

from pyteomics.mass import nist_mass

from .constants import PROTON

_ELEMENT_KEY = re.compile(r"^([A-Za-z+*]+?)(?:\[(\d+)\])?$")


@dataclass
class IsotopeDistribution:
    """Aggregated isotope envelope of one elemental composition.

    Peak 0 is the species made of the lightest stable isotope of every element:
    the monoisotopic peak for C, H, N, O, S, P and the halogens. (For B, Se or Fe
    the lightest isotope is not the most abundant one, so peak 0 is not the
    tallest-isotope peak there.)
    """

    mass: np.ndarray
    """Neutral abundance-weighted centroid mass of each nominal peak, in Da."""
    abundance: np.ndarray
    """Relative abundance of each peak, normalised to sum to 1 over the peaks returned."""
    charge: int = 0

    @property
    def n_peaks(self) -> int:
        return int(len(self.mass))

    @property
    def mz(self) -> np.ndarray:
        """Centroid m/z of each peak as [M+zH]z+; the neutral mass when charge is 0."""
        if not self.charge:
            return self.mass.copy()
        return (self.mass + self.charge * PROTON) / self.charge

    @property
    def spacing(self) -> np.ndarray:
        """Neutral mass of each peak relative to peak 0, in Da."""
        return self.mass - self.mass[0]


def _element_polynomial(key: str, n_peaks: int) -> Tuple[np.ndarray, np.ndarray]:
    """(probability, probability-weighted mass) per nominal offset for ONE atom."""
    m = _ELEMENT_KEY.match(key)
    if m is None or m.group(1) not in nist_mass:
        raise ValueError(f"Unknown element '{key}' in composition")
    element, fixed = m.group(1), m.group(2)
    table = nist_mass[element]
    prob = np.zeros(n_peaks)
    wmass = np.zeros(n_peaks)
    if fixed is not None:
        # 'C[13]' -- one named isotope, no distribution.
        if int(fixed) not in table:
            raise ValueError(f"Unknown isotope '{key}' in composition")
        prob[0], wmass[0] = 1.0, table[int(fixed)][0]
        return prob, wmass
    natural = {a: mp for a, mp in table.items() if a != 0 and mp[1] > 0}
    if not natural:
        # 'e*': a fixed mass, not an element with isotopes.
        prob[0], wmass[0] = 1.0, table[0][0]
        return prob, wmass
    lightest = min(natural)
    for a, (iso_mass, p) in natural.items():
        k = a - lightest
        if k < n_peaks:
            prob[k] += p
            wmass[k] += p * iso_mass
    return prob, wmass


def _convolve(a: Tuple[np.ndarray, np.ndarray], b: Tuple[np.ndarray, np.ndarray],
              n_peaks: int) -> Tuple[np.ndarray, np.ndarray]:
    """Combine two independent atom groups, truncated to n_peaks offsets.

    The masses add, so prob_ab * (m_a + m_b) = (p_a m_a) p_b + p_a (p_b m_b).
    Truncation is exact: offset k only ever depends on offsets <= k.
    """
    pa, wa = a
    pb, wb = b
    prob = np.convolve(pa, pb)[:n_peaks]
    wmass = (np.convolve(wa, pb) + np.convolve(pa, wb))[:n_peaks]
    return prob, wmass


def _power(base: Tuple[np.ndarray, np.ndarray], count: int,
           n_peaks: int) -> Tuple[np.ndarray, np.ndarray]:
    """base ** count by repeated squaring."""
    result = (np.eye(1, n_peaks)[0], np.zeros(n_peaks))
    while count:
        if count & 1:
            result = _convolve(result, base, n_peaks)
        count >>= 1
        if count:
            base = _convolve(base, base, n_peaks)
    return result


def isotope_distribution(composition: Mapping[str, int], n_peaks: int = 6,
                         charge: int = 0) -> IsotopeDistribution:
    """Isotope envelope of an elemental composition, aggregated by nominal mass.

    Parameters
    ----------
    composition : mapping of element -> count
        A dict or a ``pyteomics.mass.Composition`` for the NEUTRAL species, e.g.
        ``{'C': 41, 'H': 61, 'N': 10, 'O': 14, 'Br': 1}``. A key such as
        ``'C[13]'`` fixes that many atoms to one isotope. Counts must be >= 0.
    n_peaks : int
        Number of peaks to return (M+0 through M+(n-1)).
    charge : int
        Charge used by ``IsotopeDistribution.mz``; 0 leaves neutral masses.

    Returns
    -------
    IsotopeDistribution
        Centroid mass and abundance of each peak. Abundances are normalised over
        the peaks returned, the same convention as the averagine in ``deisotope``.
    """
    if n_peaks < 1:
        raise ValueError("n_peaks must be >= 1")
    total = (np.eye(1, n_peaks)[0], np.zeros(n_peaks))
    for key, count in composition.items():
        count = int(count)
        if count < 0:
            raise ValueError(f"Negative count for '{key}': pass the net composition")
        if count:
            total = _convolve(total, _power(_element_polynomial(key, n_peaks), count, n_peaks),
                              n_peaks)
    prob, wmass = total
    with np.errstate(invalid="ignore", divide="ignore"):
        centroid = np.where(prob > 0, wmass / prob, np.nan)
    norm = prob.sum()
    return IsotopeDistribution(mass=centroid, abundance=prob / norm if norm > 0 else prob,
                               charge=int(charge))
