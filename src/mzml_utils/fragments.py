"""
Fragment ion calculator for peptide and glycopeptide MS/MS spectra.

Calculates theoretical b/y/c/z ions with modifications, glycan Y ions,
oxonium ions, neutral losses, and charge-reduced precursors. Includes
peak matching and false match rate estimation.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

from .constants import (
    AA_MASSES, MOD_MASSES, PROTON, H2O, NH3,
    ION_ADJUSTMENTS, NEUTRAL_LOSSES, OXONIUM_IONS,
    GLYCAN_Y_LOSSES, ISOTOPE_SPACING,
)


# =============================================================================
# Data classes
# =============================================================================

@dataclass
class Modification:
    """A modification on a peptide residue or terminus."""
    position: int       # 0 = N-term, -1 = C-term, 1-based for residues
    residue: str
    mass: float
    name: str = ""


@dataclass
class TheoreticalIon:
    """A theoretical fragment ion."""
    ion_type: str       # 'b', 'y', 'c', 'z', 'Y', 'oxonium', 'precursor'
    ion_number: int
    charge: int
    mz: float
    sequence: str
    neutral_loss: str = ""
    annotation: str = ""


@dataclass
class MatchedIon:
    """A theoretical ion matched to an experimental peak."""
    ion_type: str
    ion_number: int
    charge: int
    mz: float
    sequence: str
    neutral_loss: str = ""
    annotation: str = ""
    exp_mz: float = 0.0
    exp_intensity: float = 0.0
    mass_error_ppm: float = 0.0


@dataclass
class FalseMatchRate:
    """Results of false match rate estimation."""
    fmr_peaks: float
    fmr_intensity: float
    matched_peaks: int
    matched_intensity: float
    avg_random_peaks: float
    avg_random_intensity: float
    n_shifts: int


# =============================================================================
# Fragment calculator
# =============================================================================

class FragmentCalculator:
    """
    Calculate theoretical fragment ions for peptides and glycopeptides.

    Supports b, y, c, z, Y (glycan), oxonium, and charge-reduced precursor
    ions with neutral losses.

    Args:
        peptide: Amino acid sequence (single-letter codes).
        modifications: List of dicts with ``position``, ``residue``, ``mass``
            (and optionally ``name``). Use position 0 for N-term, -1 for C-term.
        precursor_charge: Precursor charge state.
        max_fragment_charge: Maximum fragment ion charge (default 2).

    Example::

        calc = FragmentCalculator(
            "AGYSQGATQYTQAQQTR",
            [{'position': 0, 'residue': 'N-term', 'mass': 229.1629},
             {'position': 4, 'residue': 'S', 'mass': 528.2859}],
            precursor_charge=3,
        )
        ions = calc.get_all_ions_flat()
    """

    def __init__(self,
                 peptide: str,
                 modifications: List[Dict],
                 precursor_charge: int,
                 max_fragment_charge: int = 2):
        self.peptide = peptide.upper()
        self.length = len(peptide)
        self.precursor_charge = precursor_charge
        self.max_fragment_charge = min(max_fragment_charge, precursor_charge - 1)

        self.modifications: List[Modification] = []
        self.mod_by_position: Dict[int, Modification] = {}

        for mod in modifications:
            m = Modification(
                position=mod['position'],
                residue=mod.get('residue', ''),
                mass=mod['mass'],
                name=mod.get('name', ''),
            )
            self.modifications.append(m)
            self.mod_by_position[m.position] = m

        self._calculate_residue_masses()
        self.precursor_mass = self._calculate_precursor_mass()
        self.precursor_mz = (self.precursor_mass + precursor_charge * PROTON) / precursor_charge

    # ----- internal helpers -----

    def _calculate_residue_masses(self):
        self.residue_masses: List[float] = []
        for i, aa in enumerate(self.peptide):
            mass = AA_MASSES.get(aa, 0.0)
            if (i + 1) in self.mod_by_position:
                mass += self.mod_by_position[i + 1].mass
            self.residue_masses.append(mass)

        self.nterm_mod_mass = self.mod_by_position.get(0, Modification(0, '', 0)).mass
        self.cterm_mod_mass = self.mod_by_position.get(-1, Modification(-1, '', 0)).mass

    def _calculate_precursor_mass(self) -> float:
        return sum(self.residue_masses) + H2O + self.nterm_mod_mass + self.cterm_mod_mass

    def _get_glycan_position(self) -> Optional[int]:
        for pos, mod in self.mod_by_position.items():
            if abs(mod.mass - MOD_MASSES['HexNAc_TMT']) < 0.01:
                return pos
        return None

    def _get_glycan_type(self) -> Optional[str]:
        for pos, mod in self.mod_by_position.items():
            if abs(mod.mass - MOD_MASSES['HexNAc_TMT']) < 0.1:
                return 'HexNAc_TMT'
            elif abs(mod.mass - MOD_MASSES['HexNAc']) < 0.1:
                return 'HexNAc'
        return None

    # ----- backbone ions -----

    def calculate_b_ions(self, charges: Optional[List[int]] = None) -> List[TheoreticalIon]:
        """Calculate b ions (N-terminal, HCD)."""
        if charges is None:
            charges = list(range(1, self.max_fragment_charge + 1))
        ions: List[TheoreticalIon] = []
        cumulative = self.nterm_mod_mass
        for i in range(self.length - 1):
            cumulative += self.residue_masses[i]
            neutral = cumulative + ION_ADJUSTMENTS['b'] - PROTON
            for z in charges:
                mz = (neutral + z * PROTON) / z
                ions.append(TheoreticalIon('b', i + 1, z, mz, self.peptide[:i + 1],
                                           annotation=f"b{i + 1}" + "+" * z))
        return ions

    def calculate_y_ions(self, charges: Optional[List[int]] = None) -> List[TheoreticalIon]:
        """Calculate y ions (C-terminal, HCD)."""
        if charges is None:
            charges = list(range(1, self.max_fragment_charge + 1))
        ions: List[TheoreticalIon] = []
        cumulative = self.cterm_mod_mass
        for i in range(self.length - 1):
            cumulative += self.residue_masses[self.length - 1 - i]
            neutral = cumulative + ION_ADJUSTMENTS['y'] - PROTON
            for z in charges:
                mz = (neutral + z * PROTON) / z
                ions.append(TheoreticalIon('y', i + 1, z, mz, self.peptide[-(i + 1):],
                                           annotation=f"y{i + 1}" + "+" * z))
        return ions

    def calculate_c_ions(self, charges: Optional[List[int]] = None) -> List[TheoreticalIon]:
        """Calculate c ions (N-terminal, ETD)."""
        if charges is None:
            charges = list(range(1, self.max_fragment_charge + 1))
        ions: List[TheoreticalIon] = []
        cumulative = self.nterm_mod_mass
        for i in range(self.length - 1):
            cumulative += self.residue_masses[i]
            neutral = cumulative + ION_ADJUSTMENTS['c'] - PROTON
            for z in charges:
                mz = (neutral + z * PROTON) / z
                ions.append(TheoreticalIon('c', i + 1, z, mz, self.peptide[:i + 1],
                                           annotation=f"c{i + 1}" + "+" * z))
        return ions

    def calculate_z_ions(self, charges: Optional[List[int]] = None) -> List[TheoreticalIon]:
        """Calculate z ions (C-terminal, ETD)."""
        if charges is None:
            charges = list(range(1, self.max_fragment_charge + 1))
        ions: List[TheoreticalIon] = []
        cumulative = self.cterm_mod_mass
        for i in range(self.length - 1):
            cumulative += self.residue_masses[self.length - 1 - i]
            neutral = cumulative + ION_ADJUSTMENTS['z'] - PROTON
            for z in charges:
                mz = (neutral + z * PROTON) / z
                ions.append(TheoreticalIon('z', i + 1, z, mz, self.peptide[-(i + 1):],
                                           annotation=f"z{i + 1}" + "+" * z))
        return ions

    # ----- glycan and special ions -----

    def calculate_Y_ions(self, charges: Optional[List[int]] = None) -> List[TheoreticalIon]:
        """Calculate glycan Y ions (peptide + partial/no glycan)."""
        if charges is None:
            charges = list(range(1, self.precursor_charge + 1))
        ions: List[TheoreticalIon] = []

        glycan_type = self._get_glycan_type()
        if glycan_type and glycan_type in GLYCAN_Y_LOSSES:
            for y_name, loss_mass in GLYCAN_Y_LOSSES[glycan_type].items():
                y_mass = self.precursor_mass - loss_mass
                for z in charges:
                    mz = (y_mass + z * PROTON) / z
                    ion_number = 0 if 'Y0' in y_name else 1
                    ann = f"{y_name} {z}+" if 'Y0' in y_name else f"Y* {z}+"
                    ions.append(TheoreticalIon('Y', ion_number, z, mz, self.peptide,
                                               annotation=ann))

        for z in charges:
            mz = (self.precursor_mass + z * PROTON) / z
            ions.append(TheoreticalIon(
                'Y', 1 if glycan_type else 0, z, mz, self.peptide,
                annotation=f"Y1 {z}+" if glycan_type else f"[M+{z}H]{z}+"))

        return ions

    def calculate_charge_reduced_precursor(self) -> List[TheoreticalIon]:
        """Calculate charge-reduced precursor species from ETD."""
        ions: List[TheoreticalIon] = []
        for z in range(1, self.precursor_charge):
            mz = (self.precursor_mass + z * PROTON) / z
            ions.append(TheoreticalIon(
                'precursor', 0, z, mz, self.peptide,
                annotation=f"[M+{self.precursor_charge}H]{z}+\u2022 (CR)"))
        return ions

    def calculate_precursor_isotopes(self, n_isotopes: int = 4) -> List[TheoreticalIon]:
        """Calculate precursor isotope peaks."""
        ions: List[TheoreticalIon] = []
        spacing = ISOTOPE_SPACING / self.precursor_charge
        for i in range(n_isotopes):
            mz = self.precursor_mz + i * spacing
            ions.append(TheoreticalIon(
                'precursor', i, self.precursor_charge, mz, self.peptide,
                annotation=f"[M+{self.precursor_charge}H]{self.precursor_charge}+ iso{i}"))
        return ions

    def calculate_oxonium_ions(self) -> List[TheoreticalIon]:
        """Calculate glycan oxonium (diagnostic) ions."""
        return [
            TheoreticalIon('oxonium', 0, 1, mz, '', annotation=name)
            for name, mz in OXONIUM_IONS.items()
        ]

    # ----- neutral losses -----

    def calculate_neutral_loss_ions(self,
                                    base_ions: List[TheoreticalIon],
                                    loss_types: Optional[List[str]] = None
                                    ) -> List[TheoreticalIon]:
        """Calculate neutral loss variants of base fragment ions."""
        if loss_types is None:
            loss_types = ['H2O', 'NH3']
        glycan_pos = self._get_glycan_position()
        ions: List[TheoreticalIon] = []

        for base in base_ions:
            for loss_type in loss_types:
                loss_mass = NEUTRAL_LOSSES.get(loss_type, 0)
                if loss_mass == 0:
                    continue
                if loss_type in ('HexNAc_TMT', 'HexNAc') and glycan_pos:
                    if base.ion_type in ('b', 'c') and base.ion_number < glycan_pos:
                        continue
                    if base.ion_type in ('y', 'z') and base.ion_number < (self.length - glycan_pos + 1):
                        continue
                new_mz = base.mz - loss_mass / base.charge
                if new_mz > 0:
                    ions.append(TheoreticalIon(
                        base.ion_type, base.ion_number, base.charge, new_mz,
                        base.sequence, neutral_loss=loss_type,
                        annotation=f"{base.annotation}-{loss_type}"))
        return ions

    # ----- convenience methods -----

    def calculate_all_ions(self,
                           include_neutral_losses: bool = True,
                           neutral_loss_types: Optional[List[str]] = None
                           ) -> Dict[str, List[TheoreticalIon]]:
        """Calculate all theoretical ions, organised by type."""
        if neutral_loss_types is None:
            neutral_loss_types = ['H2O', 'NH3', 'HexNAc_TMT', 'HexNAc']

        result: Dict[str, List[TheoreticalIon]] = {
            'b': self.calculate_b_ions(),
            'y': self.calculate_y_ions(),
            'c': self.calculate_c_ions(),
            'z': self.calculate_z_ions(),
            'Y': self.calculate_Y_ions(),
            'precursor': self.calculate_precursor_isotopes() + self.calculate_charge_reduced_precursor(),
            'oxonium': self.calculate_oxonium_ions(),
        }
        if include_neutral_losses:
            for ion_type in ('b', 'y', 'c', 'z'):
                result[f'{ion_type}_NL'] = self.calculate_neutral_loss_ions(
                    result[ion_type], neutral_loss_types)
        return result

    def get_all_ions_flat(self, **kwargs) -> List[TheoreticalIon]:
        """Get all ions as a flat list."""
        return [ion for ions in self.calculate_all_ions(**kwargs).values() for ion in ions]


# =============================================================================
# Peak matching
# =============================================================================

def match_peaks(theoretical_ions: List[TheoreticalIon],
                exp_mz: np.ndarray,
                exp_intensity: np.ndarray,
                tolerance_ppm: float = 20.0,
                match_isotopes: bool = True,
                max_isotope: int = 2) -> List[MatchedIon]:
    """
    Match experimental peaks to theoretical ions.

    Args:
        theoretical_ions: Theoretical ion list.
        exp_mz: Experimental m/z array.
        exp_intensity: Experimental intensity array.
        tolerance_ppm: Mass tolerance in ppm.
        match_isotopes: Also try M+1, M+2 isotope peaks.
        max_isotope: Maximum isotope offset.

    Returns:
        List of :class:`MatchedIon` objects.
    """
    matched: List[MatchedIon] = []
    used_peaks: Set[int] = set()

    def _sort_key(ion: TheoreticalIon):
        if ion.ion_type == 'Y':
            p = 0
        elif not ion.neutral_loss:
            p = 1
        else:
            p = 2
        return (p, ion.mz)

    for ion in sorted(theoretical_ions, key=_sort_key):
        offsets = list(range(max_isotope + 1)) if match_isotopes else [0]
        found = False
        for iso in offsets:
            if found:
                break
            target = ion.mz + iso * ISOTOPE_SPACING / ion.charge
            tol = target * tolerance_ppm / 1e6
            indices = np.where(np.abs(exp_mz - target) <= tol)[0]
            if len(indices) == 0:
                continue
            for idx in indices[np.argsort(np.abs(exp_mz[indices] - target))]:
                if idx in used_peaks:
                    continue
                err = (exp_mz[idx] - target) / target * 1e6
                ann = f"{ion.annotation}+{iso}" if iso > 0 else ion.annotation
                matched.append(MatchedIon(
                    ion_type=ion.ion_type,
                    ion_number=ion.ion_number,
                    charge=ion.charge,
                    mz=ion.mz,
                    sequence=ion.sequence,
                    neutral_loss=ion.neutral_loss,
                    annotation=ann,
                    exp_mz=exp_mz[idx],
                    exp_intensity=exp_intensity[idx],
                    mass_error_ppm=err,
                ))
                used_peaks.add(idx)
                found = True
                break
    return matched


# =============================================================================
# False match rate
# =============================================================================

def calculate_false_match_rate(
    theoretical_ions: List[TheoreticalIon],
    exp_mz: np.ndarray,
    exp_intensity: np.ndarray,
    tolerance_ppm: float = 20.0,
    shift_range: float = 25.0,
    shift_step: float = 1.0,
) -> FalseMatchRate:
    """
    Estimate false match rate via spectrum shifting (Schulte et al., Anal. Chem. 2025).

    The spectrum is shifted by pi +/- *shift_range* Th in *shift_step* increments.
    The pi offset prevents false matches from isotope patterns.

    Args:
        theoretical_ions: Theoretical ion list.
        exp_mz: Experimental m/z array.
        exp_intensity: Experimental intensity array.
        tolerance_ppm: Matching tolerance in ppm.
        shift_range: Range of shifts in Th.
        shift_step: Step size for shifts in Th.

    Returns:
        :class:`FalseMatchRate` with peak- and intensity-based FMR.
    """
    true_matched = match_peaks(theoretical_ions, exp_mz, exp_intensity,
                               tolerance_ppm, match_isotopes=False)
    true_count = len(true_matched)
    true_intensity = sum(m.exp_intensity for m in true_matched)

    if true_count == 0:
        return FalseMatchRate(0.0, 0.0, 0, 0.0, 0.0, 0.0, 0)

    shifts = np.arange(-shift_range, shift_range + shift_step, shift_step) + np.pi
    peak_counts: List[int] = []
    int_sums: List[float] = []

    for shift in shifts:
        shifted = match_peaks(theoretical_ions, exp_mz + shift, exp_intensity,
                              tolerance_ppm, match_isotopes=False)
        peak_counts.append(len(shifted))
        int_sums.append(sum(m.exp_intensity for m in shifted))

    avg_peaks = float(np.mean(peak_counts))
    avg_int = float(np.mean(int_sums))

    return FalseMatchRate(
        fmr_peaks=avg_peaks / true_count,
        fmr_intensity=avg_int / true_intensity if true_intensity > 0 else 0.0,
        matched_peaks=true_count,
        matched_intensity=true_intensity,
        avg_random_peaks=avg_peaks,
        avg_random_intensity=avg_int,
        n_shifts=len(shifts),
    )


# =============================================================================
# Annotation statistics
# =============================================================================

def calculate_annotation_statistics(
    matched_ions: List[MatchedIon],
    theoretical_ions: List[TheoreticalIon],
    exp_mz: np.ndarray,
    exp_intensity: np.ndarray,
    peptide_length: int,
) -> Dict:
    """
    Calculate comprehensive annotation statistics.

    Returns a dict with sequence_coverage, peaks_annotated,
    intensity_annotated, fragments_found, and detail strings.
    """
    total_bonds = peptide_length - 1
    n_term: Set[int] = set()
    c_term: Set[int] = set()

    for ion in matched_ions:
        if ion.ion_type in ('b', 'c') and ion.ion_number > 0:
            n_term.add(ion.ion_number)
        elif ion.ion_type in ('y', 'z') and ion.ion_number > 0:
            c_term.add(ion.ion_number)

    covered = set()
    for pos in n_term:
        covered.add(pos)
    for pos in c_term:
        covered.add(peptide_length - pos)

    seq_cov = len(covered) / total_bonds if total_bonds > 0 else 0.0

    matched_mzs = {m.exp_mz for m in matched_ions}
    peaks_ann = len(matched_mzs) / len(exp_mz) if len(exp_mz) > 0 else 0.0

    total_int = float(np.sum(exp_intensity))
    matched_int = sum(m.exp_intensity for m in matched_ions)
    int_ann = matched_int / total_int if total_int > 0 else 0.0

    theo_backbone = {(i.ion_type, i.ion_number, i.charge)
                     for i in theoretical_ions
                     if i.ion_type in ('b', 'y', 'c', 'z') and not i.neutral_loss}
    matched_backbone = {(m.ion_type, m.ion_number, m.charge)
                        for m in matched_ions
                        if m.ion_type in ('b', 'y', 'c', 'z') and not m.neutral_loss}
    frag_found = len(matched_backbone) / len(theo_backbone) if theo_backbone else 0.0

    return {
        'sequence_coverage': seq_cov,
        'sequence_coverage_bonds': f"{len(covered)}/{total_bonds}",
        'peaks_annotated': peaks_ann,
        'peaks_annotated_count': f"{len(matched_mzs)}/{len(exp_mz)}",
        'intensity_annotated': int_ann,
        'fragments_found': frag_found,
        'fragments_found_count': f"{len(matched_backbone)}/{len(theo_backbone)}",
        'n_term_coverage': n_term,
        'c_term_coverage': c_term,
    }
