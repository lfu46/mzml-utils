"""3D-Beacons -- the federated structure lookup, and the preferred resolver.

One request by UniProt accession returns every model any member provider holds:
PDBe (experimental), AlphaFold, SWISS-MODEL, AlphaFill, PED.  Critically for this
package, **3D-Beacons emits the AlphaFold file URL with the correct version
already in it**, which is what makes hand-building
``AF-{acc}-F1-model_v{N}.cif`` unnecessary -- and hand-building is exactly what
broke when v4 and v5 files were deleted from the AlphaFold DB.

Gotcha that costs a wrong answer if missed: an accession with no models returns
an **empty object with HTTP 200**, not a 404.  Test for emptiness.
"""
from __future__ import annotations

import re
from typing import Any, Dict, Iterable, List, Optional

from . import http
from .records import (
    AB_INITIO,
    EXPERIMENTAL,
    FetchStatus,
    ResolveResult,
    StructureFetchError,
    StructureRecord,
    TEMPLATE_BASED,
)

BASE = "https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api"
SUMMARY_URL = BASE + "/uniprot/summary/{qualifier}.json"
BATCH_SUMMARY_URL = BASE + "/uniprot/summary"
HEALTH_URL = BASE + "/health/"

#: AF-P08571-F1-model_v6.cif -> 6
_AF_VERSION_RE = re.compile(r"-model_v(\d+)\.(?:cif|pdb|bcif)", re.IGNORECASE)


def parse_model_version(url: str) -> Optional[int]:
    """Pull the AlphaFold model version out of a file URL, if present."""
    m = _AF_VERSION_RE.search(url or "")
    return int(m.group(1)) if m else None


def _to_record(accession: str, entry: Dict[str, Any]) -> Optional[StructureRecord]:
    """Normalise one ``structures[]`` element into a StructureRecord.

    Accepts both the documented ``{"summary": {...}}`` shape and a flat shape,
    because provider deployments have shipped both.
    """
    summary = entry.get("summary") if isinstance(entry, dict) else None
    if not isinstance(summary, dict):
        summary = entry if isinstance(entry, dict) else None
    if not summary:
        return None
    url = summary.get("model_url") or ""
    if not url:
        return None
    return StructureRecord(
        accession=accession,
        provider=summary.get("provider") or "",
        model_identifier=summary.get("model_identifier") or "",
        model_url=url,
        model_category=summary.get("model_category") or "",
        model_format=summary.get("model_format") or "",
        model_page_url=summary.get("model_page_url") or "",
        confidence_avg_local_score=summary.get("confidence_avg_local_score"),
        confidence_type=summary.get("confidence_type") or "",
        coverage=summary.get("coverage"),
        uniprot_start=summary.get("uniprot_start"),
        uniprot_end=summary.get("uniprot_end"),
        experimental_method=summary.get("experimental_method") or "",
        resolution=summary.get("resolution"),
        sequence_identity=summary.get("sequence_identity"),
        model_version=parse_model_version(url),
        raw=dict(summary),
    )


def parse_summary(accession: str, payload: Any) -> List[StructureRecord]:
    """Turn a summary payload into records. Empty payload -> empty list."""
    if not isinstance(payload, dict) or not payload:
        return []
    entry = payload.get("uniprot_entry") or {}
    resolved = entry.get("ac") or accession
    out: List[StructureRecord] = []
    for item in payload.get("structures") or []:
        rec = _to_record(resolved, item)
        if rec is not None:
            out.append(rec)
    return out


def summary(
    accession: str,
    *,
    provider: Optional[str] = None,
    exclude_provider: Optional[str] = None,
) -> List[StructureRecord]:
    """All models 3D-Beacons knows for ``accession``.

    Raises ``StructureFetchError`` on transport failure so the caller can tell
    that apart from a genuine empty result.
    """
    params: Dict[str, Any] = {}
    if provider:
        params["provider"] = provider
    if exclude_provider:
        params["exclude_provider"] = exclude_provider
    status, payload = http.get_json(
        SUMMARY_URL.format(qualifier=accession), params=params or None
    )
    if status == 404:
        return []
    if status >= 400:
        raise StructureFetchError(f"3D-Beacons returned HTTP {status} for {accession}")
    return parse_summary(accession, payload)


def summaries(accessions: Iterable[str]) -> Dict[str, List[StructureRecord]]:
    """Batch lookup, falling back to serial GETs if the batch route misbehaves.

    The POST route accepts an accession list; its response shape has varied
    across deployments, so anything unexpected degrades to per-accession GETs
    rather than silently returning nothing.
    """
    accs = [a for a in dict.fromkeys(accessions) if a]
    if not accs:
        return {}
    out: Dict[str, List[StructureRecord]] = {}
    try:
        status, payload = http.post_json(BATCH_SUMMARY_URL, {"accessions": accs})
        if status < 400 and isinstance(payload, list) and payload:
            for item in payload:
                if not isinstance(item, dict):
                    continue
                entry = item.get("uniprot_entry") or {}
                acc = entry.get("ac")
                if acc:
                    out[acc] = parse_summary(acc, item)
    except StructureFetchError:
        out = {}
    missing = [a for a in accs if a not in out]
    for acc in missing:
        try:
            out[acc] = summary(acc)
        except StructureFetchError:
            out[acc] = []
    return out


#: Providers whose "experimentally determined" models are not residue-resolved
#: atomic coordinates. SASBDB serves small-angle-scattering bead envelopes, which
#: 3D-Beacons still labels EXPERIMENTALLY DETERMINED with coverage 1.0 -- so an
#: unguarded coverage sort ranks five SAXS envelopes above 66 real titin crystal
#: structures. They are demoted, not dropped: they are legitimate data, just not
#: something you can map a glycosite onto.
LOW_INFORMATION_PROVIDERS = frozenset({"sasbdb", "ped"})


def _is_low_information(rec: StructureRecord) -> bool:
    return rec.provider.lower().replace("-", "").replace(" ", "") in LOW_INFORMATION_PROVIDERS


def _rank_key(rec: StructureRecord, prefer: str):
    """Sort key -- lower tuples sort first."""
    coverage = rec.coverage if rec.coverage is not None else 0.0
    plddt = rec.confidence_avg_local_score or 0.0
    resolution = rec.resolution if rec.resolution is not None else 99.0
    low = 1 if _is_low_information(rec) else 0

    if prefer == "experimental":
        tier = {EXPERIMENTAL: 0, TEMPLATE_BASED: 2, AB_INITIO: 3}.get(rec.model_category, 4)
        return (tier + low, -coverage, resolution)
    if prefer == "coverage":
        return (low, -coverage, 0 if rec.is_experimental else 1, resolution)
    # "predicted" (default): a full-length AlphaFold model is what showcase
    # figures want; experimental structures here are usually partial fragments.
    #
    # plddt is only a tiebreak WITHIN a tier, never across one: AlphaFold reports
    # pLDDT on 0-100 while SWISS-MODEL reports its own score on 0-1, so a global
    # sort on this field would rank every SWISS-MODEL entry last by accident.
    if rec.is_canonical_alphafold:
        tier = 0
    elif rec.is_alphafold:
        tier = 1          # checksum-keyed twin of the same sequence
    elif rec.model_category == TEMPLATE_BASED:
        tier = 2
    elif rec.is_experimental:
        tier = 3
    else:
        tier = 4
    return (tier + low, -plddt, -coverage)


def rank(records: List[StructureRecord], prefer: str = "predicted") -> List[StructureRecord]:
    """Order candidates. ``prefer`` is 'predicted' | 'experimental' | 'coverage'."""
    if prefer not in ("predicted", "experimental", "coverage"):
        raise ValueError(
            f"prefer must be 'predicted', 'experimental' or 'coverage', got {prefer!r}"
        )
    return sorted(records, key=lambda r: _rank_key(r, prefer))


def resolve(accession: str, *, prefer: str = "predicted") -> ResolveResult:
    """Resolve one accession to its ranked structure candidates."""
    try:
        records = summary(accession)
    except StructureFetchError as exc:
        return ResolveResult(
            requested=accession,
            accession=accession,
            status=FetchStatus.FETCH_FAILED,
            detail=str(exc),
        )
    resolved = records[0].accession if records else accession
    return ResolveResult(
        requested=accession,
        accession=resolved,
        records=rank(records, prefer),
        status=(
            FetchStatus.RESOLVED if records else FetchStatus.NO_MODEL_AVAILABLE
        ),
        detail="" if records else "3D-Beacons returned no models",
    )


def health() -> Dict[str, Any]:
    """Per-provider status. Providers can be down while the hub returns 200."""
    status, payload = http.get_json(HEALTH_URL)
    if status >= 400 or not isinstance(payload, dict):
        raise StructureFetchError(f"3D-Beacons health returned HTTP {status}")
    return payload
