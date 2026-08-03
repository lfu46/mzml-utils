"""AlphaFold DB.

Use :mod:`.beacons` to *resolve* a structure; use this module when you need
AlphaFold-specific metadata -- mean pLDDT, the PAE document, per-version history,
isoform entries.

Three facts this module exists to encode, each of which has already broken
hand-written code:

1. **Never hard-code the model version.**  v6 is current (released 2025-10-21);
   the v4 and v5 *files are deleted* -- ``AF-P08571-F1-model_v4.cif`` is a 404.
   Read ``latestVersion`` from the API, or take the URL 3D-Beacons hands you.
2. **v6 abolished the fragment scheme.**  There is no ``-F2-``/``-F3-`` any more:
   proteins over ~2700 aa are simply absent (titin ``Q8WZ42`` -> 404).  The route
   to those proteins is an isoform entry (``AF-{acc}-{n}-F1``) or an experimental
   structure, not a fragment suffix.
3. **``/api/prediction`` returns an array that can hold more than one model.**
   ``P08571`` also yields a checksum-keyed entry ``AF-0000000003510504``.  Filter
   on ``isUniProt``; never take ``[0]``.

The API is also mid-rename (announced sunset 2026-06-25, both field sets still
present as of 2026-08-03).  ``_pick`` reads the new name with the old as
fallback, so this keeps working whichever set disappears.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from . import http
from .records import (
    AB_INITIO,
    FetchStatus,
    ResolveResult,
    StructureFetchError,
    StructureRecord,
)

BASE = "https://alphafold.ebi.ac.uk"
PREDICTION_URL = BASE + "/api/prediction/{qualifier}"
UNIPROT_SUMMARY_URL = BASE + "/api/uniprot/summary/{qualifier}.json"
OPENAPI_URL = BASE + "/api/openapi.json"

#: Field renames announced by EBI. new name -> legacy name.
FIELD_ALIASES = {
    "modelEntityId": "entryId",
    "sequenceStart": "uniprotStart",
    "sequenceEnd": "uniprotEnd",
    "sequence": "uniprotSequence",
    "isUniProtReviewed": "isReviewed",
    "isUniProtReferenceProteome": "isReferenceProteome",
}


def _pick(entry: Dict[str, Any], key: str, default: Any = None) -> Any:
    """Read ``key``, falling back to its legacy name if the new one is absent."""
    if key in entry and entry[key] is not None:
        return entry[key]
    legacy = FIELD_ALIASES.get(key)
    if legacy and entry.get(legacy) is not None:
        return entry[legacy]
    # also handle being handed a legacy name while the API serves the new one
    for new, old in FIELD_ALIASES.items():
        if old == key and entry.get(new) is not None:
            return entry[new]
    return default


def _to_record(entry: Dict[str, Any]) -> StructureRecord:
    url = entry.get("cifUrl") or entry.get("pdbUrl") or ""
    version = entry.get("latestVersion")
    return StructureRecord(
        accession=entry.get("uniprotAccession") or "",
        provider="AlphaFold",
        model_identifier=_pick(entry, "modelEntityId", "") or "",
        model_url=url,
        model_category=AB_INITIO,
        model_format="MMCIF" if url.endswith(".cif") else "PDB",
        model_page_url=f"{BASE}/entry/{_pick(entry, 'modelEntityId', '')}",
        confidence_avg_local_score=entry.get("globalMetricValue"),
        confidence_type="pLDDT",
        uniprot_start=_pick(entry, "sequenceStart"),
        uniprot_end=_pick(entry, "sequenceEnd"),
        model_version=int(version) if isinstance(version, int) else None,
        raw=dict(entry),
    )


def prediction(qualifier: str, *, canonical_only: bool = True) -> List[Dict[str, Any]]:
    """Raw ``/api/prediction`` entries for a UniProt accession or AFDB id.

    ``canonical_only`` keeps the ``AF-{accession}-F1`` entry and drops the
    checksum-keyed twin.

    Note on how that filter is done: the obvious ``isUniProt`` flag **does not
    discriminate**.  Verified against P08571 -- both ``AF-P08571-F1`` and
    ``AF-0000000003510504`` carry ``isUniProt: True`` and
    ``uniprotAccession: "P08571"``; they differ only in the entry id (and in
    ``latestVersion``, 6 vs 1).  So the entry id is what is matched here.

    Returns ``[]`` when AlphaFold has no model (HTTP 404).  Raises
    ``StructureFetchError`` for a malformed identifier (HTTP 400) or a transport
    failure, so "no model" is never conflated with "bad request".
    """
    status, payload = http.get_json(PREDICTION_URL.format(qualifier=qualifier))
    if status == 404:
        return []
    if status == 400:
        raise StructureFetchError(
            f"AlphaFold rejected {qualifier!r} as a malformed identifier (HTTP 400)"
        )
    if status >= 400 or payload is None:
        raise StructureFetchError(
            f"AlphaFold returned HTTP {status} for {qualifier}"
        )
    entries = payload if isinstance(payload, list) else [payload]
    entries = [e for e in entries if isinstance(e, dict)]
    if canonical_only and entries:
        base = qualifier.split("-")[0].upper()
        canon = [
            e
            for e in entries
            if base in str(_pick(e, "modelEntityId", "") or "").upper()
        ]
        # Fall back to everything rather than to nothing: an AFDB id passed as
        # the qualifier legitimately has no accession in its entry id.
        return canon or entries
    return entries


def latest_version(qualifier: str) -> Optional[int]:
    """The current model version for this accession, read from the API."""
    entries = prediction(qualifier)
    for e in entries:
        v = e.get("latestVersion")
        if isinstance(v, int):
            return v
    return None


def records(qualifier: str) -> List[StructureRecord]:
    """AlphaFold models for an accession, as StructureRecords."""
    return [_to_record(e) for e in prediction(qualifier)]


def resolve(accession: str) -> ResolveResult:
    """Resolve via AlphaFold alone. Prefer :func:`.beacons.resolve` in new code."""
    try:
        recs = records(accession)
    except StructureFetchError as exc:
        bad = "malformed identifier" in str(exc)
        return ResolveResult(
            requested=accession,
            accession=accession,
            status=(
                FetchStatus.INVALID_ACCESSION if bad else FetchStatus.FETCH_FAILED
            ),
            detail=str(exc),
        )
    return ResolveResult(
        requested=accession,
        accession=recs[0].accession if recs else accession,
        records=recs,
        status=FetchStatus.RESOLVED if recs else FetchStatus.NO_MODEL_AVAILABLE,
        detail="" if recs else "AlphaFold has no model for this accession",
    )


def pae_url(qualifier: str) -> Optional[str]:
    """URL of the predicted-aligned-error JSON, at the version the API reports."""
    for e in prediction(qualifier):
        url = e.get("paeDocUrl")
        if url:
            return url
    return None


def isoform_entries(accession: str, max_isoform: int = 20) -> List[Dict[str, Any]]:
    """Probe isoform entries for a protein whose canonical model is absent.

    This is the supported route for proteins over the ~2700 aa ceiling: the
    canonical sequence has no model, but individual isoforms often do (nesprin-1
    ``Q8NF91`` has no canonical entry and six isoform entries, all under 1725 aa).

    Costs one request per probed isoform, so call it only after the canonical
    lookup has come back empty.
    """
    base = accession.split("-")[0]
    found: List[Dict[str, Any]] = []
    for n in range(1, max_isoform + 1):
        try:
            entries = prediction(f"{base}-{n}")
        except StructureFetchError:
            continue
        found.extend(entries)
    return found
