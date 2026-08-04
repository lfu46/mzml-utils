"""Tests for the deisotope module."""

import math

import numpy as np
import pytest

from mzml_utils.constants import NEUTRON_MASS, PROTON
from mzml_utils.deisotope import (
    IsotopeCluster,
    _kl_score,
    _poisson_distribution,
    deisotope,
    score_ms1_envelope,
)


class TestPoissonDistribution:
    def test_small_mass(self):
        """Poisson distribution for small peptide (~500 Da)."""
        dist = _poisson_distribution(500.0, 4)
        assert len(dist) == 4
        assert abs(dist.sum() - 1.0) < 1e-6
        # For 500 Da, lambda=0.297, M+0 should dominate
        assert dist[0] > dist[1] > dist[2]

    def test_large_mass(self):
        """Poisson distribution for large glycopeptide (~3000 Da)."""
        dist = _poisson_distribution(3000.0, 6)
        assert len(dist) == 6
        assert abs(dist.sum() - 1.0) < 1e-6
        # For 3000 Da, lambda=1.782, M+1 should be close to M+0
        assert dist[1] > dist[0] * 0.5

    def test_zero_mass(self):
        """Edge case: zero mass."""
        dist = _poisson_distribution(0.0, 3)
        assert len(dist) == 3
        # lambda ~ 0, M+0 should be ~1.0
        assert dist[0] > 0.99


class TestKLScore:
    def test_perfect_match(self):
        """Identical distributions should give score = 1.0."""
        dist = np.array([0.5, 0.3, 0.2])
        score = _kl_score(dist, dist)
        assert score > 0.99

    def test_poor_match(self):
        """Very different distributions should give low score."""
        obs = np.array([0.1, 0.9])
        theo = np.array([0.9, 0.1])
        score = _kl_score(obs, theo)
        assert score < 0.5

    def test_zero_intensity(self):
        """All-zero observed should return 0."""
        obs = np.array([0.0, 0.0])
        theo = np.array([0.5, 0.5])
        score = _kl_score(obs, theo)
        assert score == 0.0


class TestDeisotope:
    def test_empty_spectrum(self):
        """Empty arrays returned unchanged."""
        mz, inten = deisotope(np.array([]), np.array([]))
        assert len(mz) == 0
        assert len(inten) == 0

    def test_singleton_peaks_retained(self):
        """Peaks with no isotope neighbors are kept."""
        mz = np.array([200.0, 400.0, 600.0])
        inten = np.array([1000.0, 500.0, 300.0])
        out_mz, out_inten = deisotope(mz, inten)
        assert len(out_mz) == 3

    def test_charge1_doublet(self):
        """M+1 peak at charge 1 should be removed."""
        mono_mz = 500.0
        m1_mz = mono_mz + NEUTRON_MASS  # M+1 at z=1
        mz = np.array([mono_mz, m1_mz])
        # Reasonable isotope ratio for ~500 Da peptide
        inten = np.array([1000.0, 300.0])
        out_mz, out_inten = deisotope(mz, inten, max_charge=2)
        assert len(out_mz) == 1
        assert abs(out_mz[0] - mono_mz) < 0.001

    def test_charge2_triplet(self):
        """M+1 and M+2 peaks at charge 2 should be removed."""
        mono_mz = 500.0
        spacing = NEUTRON_MASS / 2
        mz = np.array([mono_mz, mono_mz + spacing, mono_mz + 2 * spacing])
        # Averagine-like ratios for ~1000 Da neutral mass
        inten = np.array([1000.0, 600.0, 200.0])
        out_mz, out_inten = deisotope(mz, inten, max_charge=3)
        assert len(out_mz) == 1
        assert abs(out_mz[0] - mono_mz) < 0.001

    def test_mixed_peaks(self):
        """Isotope peaks removed while unrelated peaks kept."""
        mono_mz = 500.0
        m1_mz = mono_mz + NEUTRON_MASS
        unrelated = 700.0
        mz = np.array([mono_mz, m1_mz, unrelated])
        inten = np.array([1000.0, 300.0, 800.0])
        out_mz, out_inten = deisotope(mz, inten)
        assert len(out_mz) == 2
        assert abs(out_mz[0] - mono_mz) < 0.001
        assert abs(out_mz[1] - unrelated) < 0.001

    def test_return_clusters(self):
        """return_clusters=True returns IsotopeCluster objects."""
        mono_mz = 500.0
        m1_mz = mono_mz + NEUTRON_MASS
        mz = np.array([mono_mz, m1_mz])
        inten = np.array([1000.0, 300.0])
        out_mz, out_inten, clusters = deisotope(mz, inten, return_clusters=True)
        assert len(clusters) >= 1
        c = clusters[0]
        assert isinstance(c, IsotopeCluster)
        assert c.charge >= 1
        assert c.n_peaks == 2
        assert abs(c.mono_mz - mono_mz) < 0.001
        assert c.score > 0

    def test_noise_peak_not_removed(self):
        """Peak at isotope spacing but wrong intensity ratio is not removed."""
        mono_mz = 500.0
        m1_mz = mono_mz + NEUTRON_MASS
        mz = np.array([mono_mz, m1_mz])
        # M+1 is 50x more intense than M+0 — not a valid isotope pattern
        inten = np.array([10.0, 5000.0])
        out_mz, out_inten = deisotope(mz, inten)
        # Both peaks should be retained (bad averagine score)
        assert len(out_mz) == 2

    def test_two_independent_clusters(self):
        """Two separate isotope clusters both processed correctly."""
        # Cluster A: z=1 at 400 Da
        a_mono = 400.0
        a_m1 = a_mono + NEUTRON_MASS
        # Cluster B: z=1 at 800 Da
        b_mono = 800.0
        b_m1 = b_mono + NEUTRON_MASS
        mz = np.array([a_mono, a_m1, b_mono, b_m1])
        inten = np.array([1000.0, 250.0, 800.0, 350.0])
        out_mz, out_inten = deisotope(mz, inten)
        # Should keep both monoisotopic, remove both M+1
        assert len(out_mz) == 2
        assert abs(out_mz[0] - a_mono) < 0.001
        assert abs(out_mz[1] - b_mono) < 0.001

    def test_performance(self):
        """2000 peaks should deisotope in under 2 seconds."""
        import time
        rng = np.random.default_rng(42)
        mz = np.sort(rng.uniform(100, 2000, 2000))
        inten = rng.uniform(100, 10000, 2000)
        start = time.perf_counter()
        deisotope(mz, inten, max_charge=4)
        elapsed = time.perf_counter() - start
        assert elapsed < 2.0, f"Deisotoping took {elapsed:.2f}s (>2s)"


class TestScoreMS1Envelope:
    def _make_envelope(self, mono_mz, charge, neutral_mass, n_peaks=5):
        """Create a synthetic isotope envelope matching averagine."""
        spacing = NEUTRON_MASS / charge
        dist = _poisson_distribution(neutral_mass, n_peaks)
        mz = np.array([mono_mz + i * spacing for i in range(n_peaks)])
        inten = dist * 10000
        return mz, inten

    def test_clean_envelope(self):
        """Well-behaved isotope envelope scores high."""
        mono_mz = 800.0
        charge = 2
        neutral_mass = mono_mz * charge - charge * PROTON
        mz, inten = self._make_envelope(mono_mz, charge, neutral_mass)
        cluster = score_ms1_envelope(mz, inten, mono_mz, charge)
        assert cluster.score > 0.8
        assert abs(cluster.mono_mz - mono_mz) < 0.01
        assert cluster.n_peaks == 5

    def test_wrong_monoisotopic(self):
        """Reported precursor at M+1 — should find true M+0."""
        mono_mz = 800.0
        charge = 2
        neutral_mass = mono_mz * charge - charge * PROTON
        spacing = NEUTRON_MASS / charge
        mz, inten = self._make_envelope(mono_mz, charge, neutral_mass)
        # Report M+1 as precursor
        reported_mz = mono_mz + spacing
        cluster = score_ms1_envelope(mz, inten, reported_mz, charge)
        assert abs(cluster.mono_mz - mono_mz) < 0.01

    def test_empty_spectrum(self):
        """Empty MS1 spectrum returns default cluster."""
        cluster = score_ms1_envelope(np.array([]), np.array([]), 800.0, 2)
        assert cluster.score == 0.0
        assert cluster.n_peaks == 0

    def test_single_peak(self):
        """Only one peak found — score should be 0 (can't evaluate envelope)."""
        mz = np.array([800.0])
        inten = np.array([10000.0])
        cluster = score_ms1_envelope(mz, inten, 800.0, 2)
        assert cluster.n_peaks == 1
        assert cluster.score == 0.0  # can't score with 1 peak
