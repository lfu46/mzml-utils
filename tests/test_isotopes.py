"""Tests for the isotopes module."""

from collections import defaultdict

import numpy as np
import pytest
from pyteomics import mass

from mzml_utils.constants import PROTON
from mzml_utils.isotopes import isotope_distribution, IsotopeDistribution


def _brute_force(composition, n_peaks):
    """Aggregate pyteomics' explicit isotopologue enumeration by nominal offset."""
    mono = mass.calculate_mass(composition=composition)
    prob = defaultdict(float)
    wmass = defaultdict(float)
    # isotope_threshold keeps 2H (1.15e-4) and 17O (3.8e-4) but drops the
    # zero-abundance nuclides, which would otherwise explode the enumeration.
    for iso, ab in mass.isotopologues(composition=composition, report_abundance=True,
                                      isotope_threshold=1e-6, overall_threshold=0.0):
        m = mass.calculate_mass(composition=iso)
        k = round(m - mono)
        prob[k] += ab
        wmass[k] += ab * m
    total = sum(prob[k] for k in range(n_peaks))
    return (np.array([prob[k] / total for k in range(n_peaks)]),
            np.array([wmass[k] / prob[k] for k in range(n_peaks)]))


BRIS_1 = mass.Composition(sequence="TPENFPSK") + mass.Composition({"Br": 1, "H": -1})


class TestIsotopeDistribution:
    def test_matches_brute_force_enumeration(self):
        comp = mass.Composition({"C": 3, "H": 5, "O": 2, "Br": 1})
        abundance, centroid = _brute_force(comp, 4)
        dist = isotope_distribution(comp, n_peaks=4)
        np.testing.assert_allclose(dist.abundance, abundance, rtol=1e-9)
        np.testing.assert_allclose(dist.mass, centroid, rtol=1e-12)

    def test_monoisotopic_mass(self):
        dist = isotope_distribution(BRIS_1, n_peaks=3)
        assert dist.mass[0] == pytest.approx(mass.calculate_mass(composition=BRIS_1), abs=1e-9)

    def test_bromine_doublet(self):
        # Published BrIS-1 abundances (Li et al. 2026, Supporting Data 2, brainpy).
        dist = isotope_distribution(BRIS_1, n_peaks=3)
        np.testing.assert_allclose(dist.abundance, [0.38276, 0.18844, 0.42880], atol=1e-5)
        assert dist.abundance[2] > dist.abundance[0] > dist.abundance[1]

    def test_m_plus_1_of_a_brominated_peptide_is_the_13c_peak(self):
        # The same table lists M+1 at +1.0001922 (a zero-abundance 80Br); it is 13C.
        dist = isotope_distribution(BRIS_1, n_peaks=3)
        assert dist.spacing[1] == pytest.approx(1.0029, abs=5e-4)
        assert dist.spacing[2] == pytest.approx(1.9990, abs=5e-4)

    def test_abundances_sum_to_one(self):
        assert isotope_distribution(BRIS_1, n_peaks=6).abundance.sum() == pytest.approx(1.0)

    def test_truncation_does_not_change_the_leading_peaks(self):
        short = isotope_distribution(BRIS_1, n_peaks=3)
        long = isotope_distribution(BRIS_1, n_peaks=8)
        np.testing.assert_allclose(short.mass, long.mass[:3])
        np.testing.assert_allclose(short.abundance / short.abundance[0],
                                   long.abundance[:3] / long.abundance[0])

    def test_mz(self):
        dist = isotope_distribution(BRIS_1, n_peaks=3, charge=2)
        np.testing.assert_allclose(dist.mz, (dist.mass + 2 * PROTON) / 2)
        np.testing.assert_allclose(isotope_distribution(BRIS_1, n_peaks=3).mz, dist.mass)

    def test_accepts_a_plain_dict(self):
        dist = isotope_distribution({"C": 6, "H": 12, "O": 6}, n_peaks=3)
        assert isinstance(dist, IsotopeDistribution)
        assert dist.n_peaks == 3
        assert dist.mass[0] == pytest.approx(180.06339, abs=1e-4)

    def test_fixed_isotope_key(self):
        light = isotope_distribution({"C": 6, "H": 12, "O": 6}, n_peaks=2)
        heavy = isotope_distribution({"C": 5, "C[13]": 1, "H": 12, "O": 6}, n_peaks=2)
        assert heavy.mass[0] - light.mass[0] == pytest.approx(1.0033548, abs=1e-6)
        assert heavy.abundance[1] < light.abundance[1]

    def test_zero_count_is_ignored(self):
        np.testing.assert_allclose(
            isotope_distribution({"C": 6, "H": 12, "O": 6, "Br": 0}, n_peaks=3).abundance,
            isotope_distribution({"C": 6, "H": 12, "O": 6}, n_peaks=3).abundance)

    def test_rejects_bad_input(self):
        with pytest.raises(ValueError):
            isotope_distribution({"C": 6, "H": -1})
        with pytest.raises(ValueError):
            isotope_distribution({"Xx": 1})
        with pytest.raises(ValueError):
            isotope_distribution({"C": 1}, n_peaks=0)
