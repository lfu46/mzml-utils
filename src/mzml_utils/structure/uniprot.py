"""UniProtKB -- entries, sequence features, and accession liveness.

Two things here are easy to get wrong and both have bitten real pipelines:

* **A merged accession answers 303, not 404.**  ``Q9UEV2`` redirects to
  ``Q13683``.  A client with redirects disabled loses those proteins silently.
  :func:`accession_state` reports MERGED / DELETED / DEMERGED explicitly rather
  than letting a dead accession look like a protein with no structure.
* **Glycosylation type lives in ``description``, not ``type``.**  Every site has
  ``type == "Glycosylation"``; N- vs O-linked is a string match on
  ``"N-linked (GlcNAc..."`` vs ``"O-linked (GalNAc..."``.  Positions can also
  carry a ``modifier`` of ``UNKNOWN`` or ``OUTSIDE``, so an unchecked
  ``location.start.value`` is not necessarily a real residue number.
"""
from __future__ import annotations

from typing import Any, Dict, Iterator, List, Optional

from . import http
from .records import StructureFetchError

BASE = "https://rest.uniprot.org"
ENTRY_URL = BASE + "/uniprotkb/{accession}.json"
SEARCH_URL = BASE + "/uniprotkb/search"
STREAM_URL = BASE + "/uniprotkb/stream"

#: Feature types used for glyco/topology work (exact, case-sensitive).
GLYCOSYLATION = "Glycosylation"
SIGNAL = "Signal"
TRANSMEMBRANE = "Transmembrane"
TOPOLOGICAL_DOMAIN = "Topological domain"
DOMAIN = "Domain"
CHAIN = "Chain"
DISULFIDE = "Disulfide bond"

#: Position modifiers that mean the coordinate is not a definite residue number.
UNCERTAIN_MODIFIERS = frozenset({"UNKNOWN", "OUTSIDE", "UNSURE"})


def entry(accession: str, *, fields: Optional[str] = None) -> Optional[Dict[str, Any]]:
    """Fetch one UniProtKB entry. ``None`` when the accession does not exist.

    Redirects are followed, so a merged accession returns its successor's entry
    -- check :func:`accession_state` first when that matters.
    """
    params = {"fields": fields} if fields else None
    status, payload = http.get_json(
        ENTRY_URL.format(accession=accession), params=params
    )
    if status == 404:
        return None
    if status >= 400 or payload is None:
        raise StructureFetchError(f"UniProt returned HTTP {status} for {accession}")
    return payload


def accession_state(accession: str) -> Dict[str, Any]:
    """Whether an accession is live, and what replaced it if not.

    Queried through ``search`` rather than the entry endpoint precisely because
    search does *not* redirect -- it reports the inactive record itself.

    Returns ``{"accession", "active", "reason", "replaced_by": [...]}`` where
    ``reason`` is ``""``, ``MERGED``, ``DELETED`` or ``DEMERGED``.  A demerged
    accession maps to several successors.
    """
    status, payload = http.get_json(
        SEARCH_URL, params={"query": f"accession:{accession}", "format": "json"}
    )
    if status >= 400 or not isinstance(payload, dict):
        raise StructureFetchError(f"UniProt search returned HTTP {status} for {accession}")
    results = payload.get("results") or []
    if not results:
        return {"accession": accession, "active": False, "reason": "DELETED", "replaced_by": []}
    rec = results[0]
    if rec.get("entryType") != "Inactive":
        return {
            "accession": rec.get("primaryAccession", accession),
            "active": True,
            "reason": "",
            "replaced_by": [],
        }
    reason = (rec.get("inactiveReason") or {})
    return {
        "accession": rec.get("primaryAccession", accession),
        "active": False,
        "reason": reason.get("inactiveReasonType", "") or "",
        "replaced_by": list(reason.get("mergeDemergeTo") or []),
    }


def resolve_accession(accession: str) -> str:
    """Follow a merge to the live accession. Returns the input if already live."""
    state = accession_state(accession)
    if state["active"]:
        return state["accession"]
    replaced = state["replaced_by"]
    return replaced[0] if len(replaced) == 1 else accession


def features(
    accession: str,
    *,
    types: Optional[List[str]] = None,
    entry_json: Optional[Dict[str, Any]] = None,
) -> List[Dict[str, Any]]:
    """Sequence features, flattened to ``type / description / start / end / certain``.

    ``certain`` is False when either endpoint carries an ``UNKNOWN``/``OUTSIDE``
    modifier -- treat those positions as unusable rather than as residue numbers.
    """
    data = entry_json if entry_json is not None else entry(accession)
    if not data:
        return []
    wanted = set(types) if types else None
    out: List[Dict[str, Any]] = []
    for feat in data.get("features") or []:
        ftype = feat.get("type", "")
        if wanted and ftype not in wanted:
            continue
        loc = feat.get("location") or {}
        start, end = loc.get("start") or {}, loc.get("end") or {}
        smod = (start.get("modifier") or "").upper()
        emod = (end.get("modifier") or "").upper()
        out.append(
            {
                "type": ftype,
                "description": feat.get("description", "") or "",
                "start": start.get("value"),
                "end": end.get("value"),
                "certain": smod not in UNCERTAIN_MODIFIERS
                and emod not in UNCERTAIN_MODIFIERS,
                "feature_id": feat.get("featureId", "") or "",
                "evidences": feat.get("evidences") or [],
            }
        )
    return out


def glycosylation_sites(
    accession: str,
    *,
    linkage: str = "N",
    entry_json: Optional[Dict[str, Any]] = None,
) -> List[Dict[str, Any]]:
    """Annotated glycosylation sites. ``linkage`` is 'N', 'O' or 'any'.

    The linkage is read from ``description`` because ``type`` is
    ``"Glycosylation"`` for every site regardless of linkage.
    """
    if linkage not in ("N", "O", "any"):
        raise ValueError(f"linkage must be 'N', 'O' or 'any', got {linkage!r}")
    feats = features(accession, types=[GLYCOSYLATION], entry_json=entry_json)
    if linkage == "any":
        return feats
    prefix = f"{linkage}-linked"
    return [f for f in feats if f["description"].startswith(prefix)]


def cross_references(
    accession: str,
    database: Optional[str] = None,
    *,
    entry_json: Optional[Dict[str, Any]] = None,
) -> List[Dict[str, Any]]:
    """Cross-references, with ``properties`` flattened into a dict.

    The PDB entries carry a ``Chains`` property like ``"A=26-335"``, which gives
    chain-to-UniProt-range mapping without a separate SIFTS call.
    """
    data = entry_json if entry_json is not None else entry(accession)
    if not data:
        return []
    out = []
    for xref in data.get("uniProtKBCrossReferences") or []:
        if database and xref.get("database") != database:
            continue
        props = {p.get("key"): p.get("value") for p in xref.get("properties") or []}
        out.append(
            {
                "database": xref.get("database", ""),
                "id": xref.get("id", ""),
                "properties": props,
            }
        )
    return out


def pdb_ids(accession: str, *, entry_json: Optional[Dict[str, Any]] = None) -> List[str]:
    """PDB identifiers cross-referenced from a UniProt entry."""
    return [x["id"] for x in cross_references(accession, "PDB", entry_json=entry_json)]


def sequence(accession: str, *, entry_json: Optional[Dict[str, Any]] = None) -> str:
    """The canonical sequence."""
    data = entry_json if entry_json is not None else entry(accession)
    if not data:
        return ""
    return (data.get("sequence") or {}).get("value", "") or ""


def search(
    query: str,
    *,
    fields: Optional[str] = None,
    size: int = 500,
    max_pages: Optional[int] = None,
) -> Iterator[Dict[str, Any]]:
    """Paginated search, following the ``Link`` header cursor.

    Cursors are opaque -- follow the URL UniProt gives back, never construct
    ``?page=N``.  ``size`` is capped at 500 by the API.
    """
    params: Optional[Dict[str, Any]] = {
        "query": query,
        "size": min(size, 500),
        "format": "json",
    }
    if fields:
        params["fields"] = fields
    url = SEARCH_URL
    pages = 0
    while url:
        resp = http.request("GET", url, params=params)
        if not resp.ok:
            raise StructureFetchError(
                f"UniProt search returned HTTP {resp.status_code} for {query!r}"
            )
        payload = resp.json()
        for rec in payload.get("results") or []:
            yield rec
        pages += 1
        if max_pages is not None and pages >= max_pages:
            return
        nxt = resp.links.get("next") or {}
        url = nxt.get("url")
        params = None  # the cursor URL already carries every parameter
