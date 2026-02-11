"""
Physical constants, amino acid masses, modification masses, and common ion lists
for mass spectrometry calculations.

All masses are monoisotopic unless otherwise noted.
"""

# =============================================================================
# Physical constants
# =============================================================================

PROTON = 1.007276
"""Proton mass in Da."""

H2O = 18.010565
"""Water mass in Da."""

NH3 = 17.026549
"""Ammonia mass in Da."""

CO = 27.994915
"""Carbon monoxide mass in Da."""

ELECTRON = 0.000549
"""Electron mass in Da."""

NEUTRON_MASS = 1.003355
"""Neutron mass difference (C13 - C12) in Da."""

ISOTOPE_SPACING = NEUTRON_MASS
"""Mass spacing between isotope peaks (alias for NEUTRON_MASS)."""

# =============================================================================
# Amino acid residue masses (monoisotopic)
# =============================================================================

AA_MASSES = {
    'A': 71.03711,   # Alanine
    'R': 156.10111,  # Arginine
    'N': 114.04293,  # Asparagine
    'D': 115.02694,  # Aspartic acid
    'C': 103.00919,  # Cysteine
    'E': 129.04259,  # Glutamic acid
    'Q': 128.05858,  # Glutamine
    'G': 57.02146,   # Glycine
    'H': 137.05891,  # Histidine
    'I': 113.08406,  # Isoleucine
    'L': 113.08406,  # Leucine
    'K': 128.09496,  # Lysine
    'M': 131.04049,  # Methionine
    'F': 147.06841,  # Phenylalanine
    'P': 97.05276,   # Proline
    'S': 87.03203,   # Serine
    'T': 101.04768,  # Threonine
    'W': 186.07931,  # Tryptophan
    'Y': 163.06333,  # Tyrosine
    'V': 99.06841,   # Valine
}
"""Monoisotopic residue masses for the 20 standard amino acids."""

# =============================================================================
# Modification masses
# =============================================================================

MOD_MASSES = {
    'TMT6plex': 229.1629,
    'HexNAc': 203.0794,
    'HexNAc_TMT': 528.2859,
    'Carbamidomethyl': 57.02146,
    'Oxidation': 15.9949,
}
"""Common modification masses in Da."""

# =============================================================================
# Neutral loss masses
# =============================================================================

NEUTRAL_LOSSES = {
    'H2O': 18.010565,
    'NH3': 17.026549,
    'HexNAc_TMT': 528.2859,
    'HexNAc': 300.1308,
}
"""Common neutral loss masses in Da."""

# =============================================================================
# Oxonium (diagnostic) ions
# =============================================================================

OXONIUM_IONS = {
    'HexNAc_TMT': 529.2937,
    'HexNAc+': 300.1308,
    'HexNAc': 204.0864,
    'HexNAc-H2O': 186.0760,
    'HexNAc-2H2O': 168.0652,
    'HexNAc_frag1': 144.0652,
    'HexNAc_frag2': 138.0546,
}
"""Glycan oxonium (B-type diagnostic) ion m/z values."""

# =============================================================================
# Glycan Y ion loss definitions
# =============================================================================

GLYCAN_Y_LOSSES = {
    'HexNAc_TMT': {
        'Y0': 528.2859,
        'Y0-H2O': 528.2859 - 18.0106,
        'Y*': 203.0794,
    },
    'HexNAc': {
        'Y0': 203.0794,
        'Y0-H2O': 203.0794 - 18.0106,
    },
}
"""Mass losses from intact glycopeptide to produce Y ions."""

# =============================================================================
# Ion type mass adjustments
# =============================================================================

ION_ADJUSTMENTS = {
    'b': PROTON,
    'y': H2O + PROTON,
    'c': NH3 + PROTON,
    'z': -NH3 + H2O + PROTON,
}
"""Mass adjustments added to residue sum for each ion type."""
