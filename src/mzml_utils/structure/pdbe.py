"""PDBe -- experimental structures and SIFTS residue mapping.

``best_structures`` is the single most useful call when you want a *real*
structure for an accession: PDBe returns them pre-ranked by coverage and
resolution, so no client-side ranking is needed.

Numbering traps, all of which produce silent off-by-N errors:

* ``unp_start``/``unp_end`` are **UniProt** numbering; ``start``/``end`` are
  **PDB** numbering.  They are not interchangeable.
* SIFTS mappings additionally distinguish ``residue_number`` (label seq id) from
  ``author_residue_number``.  Coordinate files and papers use the author
  numbering; most APIs return both.
* One row comes back **per chain**, so the same ``pdb_id`` repeats.
* ``resolution`` is ``None`` for NMR entries -- do not sort on it unguarded.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from . import http
from .records import EXPERIMENTAL, StructureFetchError, StructureRecord

BASE = "https://www.ebi.ac.uk/pdbe"
BEST_STRUCTURES_URL = BASE + "/api/mappings/best_structures/{accession}"
UNIPROT_MAPPING_URL = BASE + "/api/mappings/uniprot/{pdb_id}"
SUMMARY_URL = BASE + "/api/pdb/entry/summary/{pdb_id}"
#: The "updated" mmCIF carries the SIFTS annotations merged in.
UPDATED_CIF_URL = BASE + "/static/entry/{pdb_id}_updated.cif"


def best_structures(accession: str) -> List[Dict[str, Any]]:
    """Experimental structures for an accession, best first. ``[]`` if none."""
    status, payload = http.get_json(
        BEST_STRUCTURES_URL.format(accession=accession)
    )
    if status == 404:
        return []
    if status >= 400 or not isinstance(payload, dict):
        raise StructureFetchError(
            f"PDBe best_structures returned HTTP {status} for {accession}"
        )
    for key in (accession, accession.upper(), accession.lower()):
        if key in payload:
            return list(payload[key] or [])
    # Single-key payloads are sometimes keyed by the canonical accession.
    values = list(payload.values())
    return list(values[0] or []) if len(values) == 1 else []


def best_structure_records(accession: str) -> List[StructureRecord]:
    """``best_structures`` normalised to StructureRecords, deduplicated by entry.

    One record per PDB entry (the first/best chain), since a caller asking for a
    structure wants entries, not chains.  Use :func:`best_structures` for the
    full per-chain listing.
    """
    seen: set = set()
    out: List[StructureRecord] = []
    for hit in best_structures(accession):
        pdb_id = str(hit.get("pdb_id") or "").lower()
        if not pdb_id or pdb_id in seen:
            continue
        seen.add(pdb_id)
        out.append(
            StructureRecord(
                accession=accession,
                provider="PDBe",
                model_identifier=pdb_id,
                model_url=UPDATED_CIF_URL.format(pdb_id=pdb_id),
                model_category=EXPERIMENTAL,
                model_format="MMCIF",
                model_page_url=f"https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id}",
                coverage=hit.get("coverage"),
                uniprot_start=hit.get("unp_start"),
                uniprot_end=hit.get("unp_end"),
                experimental_method=hit.get("experimental_method") or "",
                resolution=hit.get("resolution"),
                raw=dict(hit),
            )
        )
    return out


def uniprot_mappings(pdb_id: str) -> Dict[str, Any]:
    """SIFTS UniProt mappings for a PDB entry, keyed by accession.

    Each mapping carries ``chain_id``, ``struct_asym_id``, ``entity_id``,
    ``identity``, ``unp_start``/``unp_end``, and ``start``/``end`` objects with
    both ``residue_number`` and ``author_residue_number``.
    """
    pdb_id = pdb_id.lower()
    status, payload = http.get_json(UNIPROT_MAPPING_URL.format(pdb_id=pdb_id))
    if status == 404:
        return {}
    if status >= 400 or not isinstance(payload, dict):
        raise StructureFetchError(
            f"PDBe mappings returned HTTP {status} for {pdb_id}"
        )
    return (payload.get(pdb_id) or {}).get("UniProt", {}) or {}


def map_uniprot_residue(
    pdb_id: str, accession: str, position: int
) -> List[Dict[str, Any]]:
    """Where one UniProt residue sits in a PDB entry, per chain.

    Returns rows of ``chain_id / entity_id / author_residue_number /
    residue_number``.  Only chains whose mapped range covers ``position`` are
    returned; an empty list means the residue is not resolved in that entry
    (disordered loops and cleaved propeptides are routinely absent).

    The offset is computed from the mapping's own start values rather than
    assumed, because PDB numbering is frequently neither 1-based nor contiguous
    with UniProt.
    """
    out: List[Dict[str, Any]] = []
    entry = uniprot_mappings(pdb_id).get(accession) or {}
    for m in entry.get("mappings") or []:
        unp_start, unp_end = m.get("unp_start"), m.get("unp_end")
        if unp_start is None or unp_end is None:
            continue
        if not (unp_start <= position <= unp_end):
            continue
        start = m.get("start") or {}
        offset = position - unp_start
        auth = start.get("author_residue_number")
        label = start.get("residue_number")
        out.append(
            {
                "pdb_id": pdb_id.lower(),
                "chain_id": m.get("chain_id", ""),
                "struct_asym_id": m.get("struct_asym_id", ""),
                "entity_id": m.get("entity_id"),
                "author_residue_number": None if auth is None else auth + offset,
                "residue_number": None if label is None else label + offset,
                "identity": m.get("identity"),
            }
        )
    return out


def summary(pdb_id: str) -> Optional[Dict[str, Any]]:
    """Entry summary (title, method, release date, assemblies)."""
    pdb_id = pdb_id.lower()
    status, payload = http.get_json(SUMMARY_URL.format(pdb_id=pdb_id))
    if status == 404:
        return None
    if status >= 400 or not isinstance(payload, dict):
        raise StructureFetchError(f"PDBe summary returned HTTP {status} for {pdb_id}")
    entries = payload.get(pdb_id) or []
    return entries[0] if entries else None
