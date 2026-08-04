"""
Protease definitions and in silico protein digestion.

Protease rules based on mzLib (Smith lab, UW-Madison) proteases.tsv.
Supports standard, semi-specific, and dual-enzyme digestion.
"""

import re
from dataclasses import dataclass, field
from typing import List, Dict, Set, Tuple, Optional


@dataclass
class CleavageRule:
    """A single cleavage rule: residue + position + exclusions."""
    residues: str          # amino acids that trigger cleavage (e.g., "KR")
    position: str          # 'C' = cleave C-terminal, 'N' = cleave N-terminal
    exclude_after: str = ""  # residues that block cleavage when following (e.g., "P")

    def cleaves_at(self, sequence: str, pos: int) -> bool:
        """Check if this rule cleaves at position pos (0-indexed).

        For C-terminal rules: cleaves after sequence[pos].
        For N-terminal rules: cleaves before sequence[pos].
        """
        if sequence[pos] not in self.residues:
            return False
        if self.position == 'C' and self.exclude_after:
            if pos + 1 < len(sequence) and sequence[pos + 1] in self.exclude_after:
                return False
        if self.position == 'N' and self.exclude_after:
            if pos > 0 and sequence[pos - 1] in self.exclude_after:
                return False
        return True


@dataclass
class Protease:
    """A protease with cleavage rules."""
    name: str
    rules: List[CleavageRule] = field(default_factory=list)
    specificity: str = "full"  # "full", "semi", "none"

    def find_cleavage_sites(self, sequence: str) -> List[int]:
        """Find all cleavage positions in a sequence.

        Returns list of positions where cleavage occurs.
        For C-terminal rules: position is after the residue (peptide ends here).
        For N-terminal rules: position is before the residue (peptide starts here).

        Returns positions as indices into sequence where a new peptide starts.
        """
        sites = set()
        for rule in self.rules:
            for i in range(len(sequence)):
                if rule.cleaves_at(sequence, i):
                    if rule.position == 'C':
                        sites.add(i + 1)  # new peptide starts after this residue
                    else:  # N-terminal
                        sites.add(i)  # new peptide starts at this residue
        return sorted(sites)

    def digest(self, sequence: str, missed_cleavages: int = 2,
               min_length: int = 6, max_length: int = 50,
               include_protein_termini: bool = True) -> List[Tuple[int, int, str]]:
        """In silico digest a protein sequence.

        Args:
            sequence: protein amino acid sequence
            missed_cleavages: max allowed missed cleavages
            min_length: minimum peptide length
            max_length: maximum peptide length
            include_protein_termini: include protein N/C-terminus as valid boundaries

        Returns:
            List of (start, end, peptide_sequence) tuples (0-indexed, end exclusive)
        """
        if self.specificity == 'none':
            return [(0, len(sequence), sequence)]

        sites = self.find_cleavage_sites(sequence)

        # Add protein termini
        boundaries = [0] + sites + [len(sequence)]
        boundaries = sorted(set(boundaries))

        peptides = []
        for i in range(len(boundaries) - 1):
            for j in range(i + 1, min(i + 2 + missed_cleavages, len(boundaries))):
                start = boundaries[i]
                end = boundaries[j]
                pep_len = end - start
                if pep_len < min_length or pep_len > max_length:
                    continue
                peptides.append((start, end, sequence[start:end]))

        return peptides


def _parse_motif(motif_str: str) -> List[CleavageRule]:
    """Parse mzLib-style motif string into CleavageRule objects.

    Syntax:
        K| = cleave C-terminal to K
        |D = cleave N-terminal to D
        K[P]| = cleave C-terminal to K, unless followed by P
        K|,R| = multiple rules (comma-separated)
        X = any amino acid
    """
    if not motif_str or motif_str.strip() == '':
        return []

    rules = []
    for part in motif_str.split(','):
        part = part.strip()
        if not part:
            continue

        # Extract exclusion (e.g., [P])
        exclude = ""
        exclude_match = re.search(r'\[([A-Z]+)\]', part)
        if exclude_match:
            exclude = exclude_match.group(1)
            part = part[:exclude_match.start()] + part[exclude_match.end():]

        # Determine position
        if '|' in part:
            before, after = part.split('|', 1)
            if before and not after:
                # K| = C-terminal cleavage
                residues = before.replace('X', 'ACDEFGHIKLMNPQRSTVWY')
                rules.append(CleavageRule(residues, 'C', exclude))
            elif after and not before:
                # |D = N-terminal cleavage
                residues = after.replace('X', 'ACDEFGHIKLMNPQRSTVWY')
                rules.append(CleavageRule(residues, 'N', exclude))
            else:
                # Complex motif like GPX|GPX — simplified handling
                residues = before[-1] if before else ''
                if residues:
                    rules.append(CleavageRule(residues, 'C', exclude))

    return rules


def dual_enzyme_digest(sequence: str,
                       enzyme1: Protease, enzyme2: Protease,
                       missed_cleavages: int = 6,
                       min_length: int = 6, max_length: int = 35,
                       mature_start: int = 0) -> List[Tuple[int, int, str, str]]:
    """Digest with two enzymes where enzyme1 defines N-terminus and enzyme2 defines C-terminus.

    This models IMPa (N-term S/T) + Trypsin (C-term K/R) dual digestion.

    Args:
        sequence: protein amino acid sequence
        enzyme1: enzyme defining N-terminal cleavage (e.g., IMPa)
        enzyme2: enzyme defining C-terminal cleavage (e.g., trypsin)
        missed_cleavages: max missed cleavages for enzyme2
        min_length: minimum peptide length
        max_length: maximum peptide length
        mature_start: 0-indexed mature protein start (after signal peptide)

    Returns:
        List of (start, end, peptide, source) tuples.
        source is 'enzymatic', 'nterm_mature', or 'tryptic'.
    """
    n = len(sequence)
    if n == 0:
        return []

    # Enzyme1 N-terminal sites
    nterm_sites = set(enzyme1.find_cleavage_sites(sequence))
    # Enzyme2 C-terminal sites
    cterm_sites = enzyme2.find_cleavage_sites(sequence)
    cterm_positions = sorted(set(cterm_sites + [n]))  # include protein C-term

    seen = set()
    peptides = []

    def _add(start, end, source):
        key = (start, end)
        if key in seen:
            return
        pep_len = end - start
        if pep_len < min_length or pep_len > max_length:
            return
        # Count missed cleavages (internal enzyme2 sites)
        internal_mc = sum(1 for c in cterm_sites if start < c < end)
        if internal_mc > missed_cleavages:
            return
        seen.add(key)
        peptides.append((start, end, sequence[start:end], source))

    # Type 1: Enzymatic (enzyme1 N-term + enzyme2 C-term)
    for nterm in nterm_sites:
        for cterm in cterm_positions:
            if cterm <= nterm:
                continue
            _add(nterm, cterm, 'enzymatic')

    # Type 2: Tryptic (after enzyme2 site + enzyme2 C-term)
    tryptic_starts = [0] + [c for c in cterm_sites if c < n]
    for start in tryptic_starts:
        if start in nterm_sites:
            continue  # already covered by enzymatic
        for cterm in cterm_positions:
            if cterm <= start:
                continue
            _add(start, cterm, 'tryptic')

    # Type 3: Mature N-terminal
    if mature_start > 0 and mature_start < n:
        starts = [mature_start]
        if sequence[mature_start] == 'M' and mature_start + 1 < n:
            starts.append(mature_start + 1)
        for start in starts:
            if start in nterm_sites:
                continue
            for cterm in cterm_positions:
                if cterm <= start:
                    continue
                _add(start, cterm, 'nterm_mature')

    return peptides


# --- Built-in protease definitions ---

PROTEASES: Dict[str, Protease] = {}


def _register(name, motif, specificity='full'):
    PROTEASES[name] = Protease(name, _parse_motif(motif), specificity)


# Standard proteases from mzLib proteases.tsv
_register('trypsin', 'K|,R|')
_register('trypsin/P', 'K[P]|,R[P]|')
_register('Lys-C', 'K|')
_register('Lys-C/P', 'K[P]|')
_register('Lys-N', '|K')
_register('Arg-C', 'R|')
_register('Asp-N', '|D')
_register('Glu-C', 'E|')
_register('Glu-C (D)', 'E|,D|')
_register('chymotrypsin', 'F[P]|,W[P]|,Y[P]|,L[P]|')
_register('CNBr', 'M|')
_register('elastase', 'A[P]|,V[P]|,S[P]|,G[P]|,L[P]|,I[P]|')
_register('ProAlanase', 'P|,A|')
_register('non-specific', 'X|', 'none')

# O-glycoproteases
_register('IMPa', '|S,|T')       # cleaves N-terminal to Ser/Thr
_register('OgpA', '|S,|T')       # cleaves N-terminal to glyco-Ser/Thr
_register('SmE', '|S,|T')        # cleaves N-terminal to glyco-Ser/Thr
_register('StcE', '|S,|T')       # cleaves at T/S*_X_T/S motif

# Dual enzyme combinations
_register('StcE-trypsin', '|S,|T,K|,R|')
