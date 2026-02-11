"""Tests for the pairing module."""

import pytest

from mzml_utils.pairing import match_precursors, group_ms1_cycles, pair_hcd_ethcd


class TestMatchPrecursors:
    def test_exact_match(self):
        matched, ppm, offset = match_precursors(500.0, 500.0, 2, tolerance_ppm=10.0)
        assert matched is True
        assert ppm == pytest.approx(0.0)
        assert offset == 0

    def test_within_tolerance(self):
        matched, ppm, offset = match_precursors(500.001, 500.0, 2, tolerance_ppm=10.0)
        assert matched is True
        assert ppm < 10.0

    def test_outside_tolerance(self):
        # 500.0 vs 510.0 at charge 2 is too far even with isotope offsets
        matched, ppm, offset = match_precursors(500.0, 510.0, 2, tolerance_ppm=10.0)
        assert matched is False

    def test_isotope_offset(self):
        # Shift by 1 isotope for charge 2: 1.003355/2 = 0.5017
        from mzml_utils.constants import NEUTRON_MASS
        shifted_mz = 500.0 - NEUTRON_MASS / 2
        matched, ppm, offset = match_precursors(500.0, shifted_mz, 2, tolerance_ppm=10.0)
        assert matched is True
        assert offset == 1


class TestGroupMs1Cycles:
    def test_basic_grouping(self):
        scans = [
            {'scan_num': 1, 'activation_type': 'MS1'},
            {'scan_num': 2, 'activation_type': 'HCD'},
            {'scan_num': 3, 'activation_type': 'HCD'},
            {'scan_num': 4, 'activation_type': 'EThcD'},
            {'scan_num': 5, 'activation_type': 'MS1'},
            {'scan_num': 6, 'activation_type': 'HCD'},
        ]
        cycles = group_ms1_cycles(scans)
        assert len(cycles) == 2
        assert cycles[0][0]['scan_num'] == 1
        assert len(cycles[0][1]) == 3  # 2 HCD + 1 EThcD
        assert cycles[1][0]['scan_num'] == 5
        assert len(cycles[1][1]) == 1

    def test_no_ms1(self):
        scans = [
            {'scan_num': 1, 'activation_type': 'HCD'},
            {'scan_num': 2, 'activation_type': 'EThcD'},
        ]
        cycles = group_ms1_cycles(scans)
        assert len(cycles) == 1
        assert cycles[0][0] is None
        assert len(cycles[0][1]) == 2

    def test_act_type_alias(self):
        """Should also accept 'act_type' key."""
        scans = [
            {'scan_num': 1, 'act_type': 'MS1'},
            {'scan_num': 2, 'act_type': 'HCD'},
        ]
        cycles = group_ms1_cycles(scans)
        assert len(cycles) == 1


class TestPairHcdEthcd:
    def test_basic_pairing(self):
        cycles = [
            (
                {'scan_num': 1, 'activation_type': 'MS1'},
                [
                    {'scan_num': 2, 'activation_type': 'HCD',
                     'precursor_mz': 500.25, 'precursor_charge': 2},
                    {'scan_num': 3, 'activation_type': 'HCD',
                     'precursor_mz': 600.30, 'precursor_charge': 3},
                    {'scan_num': 4, 'activation_type': 'EThcD',
                     'precursor_mz': 500.25, 'precursor_charge': 2},
                ],
            )
        ]
        paired = pair_hcd_ethcd(cycles)
        assert 2 in paired
        assert paired[2] == 4
        assert 3 not in paired

    def test_no_ethcd(self):
        cycles = [
            (
                {'scan_num': 1, 'activation_type': 'MS1'},
                [
                    {'scan_num': 2, 'activation_type': 'HCD',
                     'precursor_mz': 500.0, 'precursor_charge': 2},
                ],
            )
        ]
        paired = pair_hcd_ethcd(cycles)
        assert len(paired) == 0
