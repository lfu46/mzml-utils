"""InterPro -- domain and family annotation by UniProt accession.

Two shape details matter:

* **Pagination is cursor-based.**  Follow the ``next`` URL verbatim until it is
  ``null``; constructing ``?page=N`` does not work.
* **Coordinates arrive as ``fragments``.**  A discontinuous domain is several
  fragments carrying a ``dc-status``.  Collapsing them to ``min(start)``/
  ``max(end)`` silently swallows the gap and reports a domain spanning residues
  it does not occupy, so :func:`domains` keeps the fragment list intact and only
  reports the envelope alongside it.
"""
from __future__ import annotations

from typing import Any, Dict, Iterator, List, Optional

from . import http
from .records import StructureFetchError

BASE = "https://www.ebi.ac.uk/interpro/api"
ENTRY_BY_PROTEIN_URL = BASE + "/entry/{db}/protein/uniprot/{accession}/"

DEFAULT_PAGE_SIZE = 100


def _paginate(url: str, params: Optional[Dict[str, Any]] = None) -> Iterator[Dict[str, Any]]:
    """Yield ``results`` across pages, following the cursor ``next`` link."""
    while url:
        status, payload = http.get_json(url, params=params)
        if status == 204:  # documented empty response
            return
        if status == 404:
            return
        if status >= 400 or not isinstance(payload, dict):
            raise StructureFetchError(f"InterPro returned HTTP {status} for {url}")
        for item in payload.get("results") or []:
            yield item
        url = payload.get("next")
        params = None  # the cursor URL is already complete


def entries(
    accession: str, *, db: str = "interpro", page_size: int = DEFAULT_PAGE_SIZE
) -> List[Dict[str, Any]]:
    """Raw InterPro entries matched to a protein.

    ``db`` is ``interpro`` for integrated entries, ``all`` for every member
    database, or a member name such as ``pfam``.
    """
    url = ENTRY_BY_PROTEIN_URL.format(db=db, accession=accession)
    return list(_paginate(url, {"page_size": page_size}))


def domains(
    accession: str, *, db: str = "interpro", types: Optional[List[str]] = None
) -> List[Dict[str, Any]]:
    """Domain annotations with coordinates, one row per matched location.

    Each row carries ``fragments`` (the authoritative coordinates) plus
    ``start``/``end`` as the envelope across them.  For a discontinuous domain
    the envelope is wider than the residues actually covered -- use
    ``fragments`` whenever residue-level accuracy matters.

    ``types`` filters InterPro's ``metadata.type`` (``domain``, ``family``,
    ``repeat``, ``homologous_superfamily``, ...).
    """
    wanted = {t.lower() for t in types} if types else None
    out: List[Dict[str, Any]] = []
    for item in entries(accession, db=db):
        meta = item.get("metadata") or {}
        etype = (meta.get("type") or "").lower()
        if wanted and etype not in wanted:
            continue
        for protein in item.get("proteins") or []:
            for loc in protein.get("entry_protein_locations") or []:
                frags = [
                    {
                        "start": f.get("start"),
                        "end": f.get("end"),
                        "dc_status": f.get("dc-status", ""),
                    }
                    for f in loc.get("fragments") or []
                ]
                starts = [f["start"] for f in frags if f["start"] is not None]
                ends = [f["end"] for f in frags if f["end"] is not None]
                if not starts or not ends:
                    continue
                out.append(
                    {
                        "accession": accession,
                        "entry_accession": meta.get("accession", ""),
                        "name": meta.get("name", ""),
                        "type": etype,
                        "source_database": meta.get("source_database", ""),
                        "integrated": meta.get("integrated"),
                        "protein_length": protein.get("protein_length"),
                        "in_alphafold": protein.get("in_alphafold"),
                        "start": min(starts),
                        "end": max(ends),
                        "fragments": frags,
                        "discontinuous": len(frags) > 1,
                        "score": loc.get("score"),
                        "representative": loc.get("representative"),
                    }
                )
    return out


def domains_at(accession: str, position: int, **kwargs: Any) -> List[Dict[str, Any]]:
    """Domains actually covering a residue, fragment-aware.

    Tests the fragments rather than the envelope, so a residue sitting in the
    gap of a discontinuous domain is correctly reported as *not* in it.
    """
    hits = []
    for d in domains(accession, **kwargs):
        for frag in d["fragments"]:
            s, e = frag["start"], frag["end"]
            if s is not None and e is not None and s <= position <= e:
                hits.append(d)
                break
    return hits


def in_alphafold(accession: str) -> Optional[bool]:
    """InterPro's own flag for whether AlphaFold covers this protein.

    A cheap pre-check before spending a request on the AlphaFold API.  ``None``
    when InterPro has no record of the protein at all.
    """
    for item in entries(accession, page_size=1):
        for protein in item.get("proteins") or []:
            flag = protein.get("in_alphafold")
            if flag is not None:
                return bool(flag)
    return None
