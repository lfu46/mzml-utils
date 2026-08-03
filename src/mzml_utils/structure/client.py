"""The public entry points: resolve, fetch_structure, fetch_many.

``resolve`` is the function that replaces hand-built AlphaFold URLs.  It asks
3D-Beacons first (which federates PDBe, AlphaFold, SWISS-MODEL, AlphaFill and PED
and emits the correct AlphaFold version itself), then falls back to the
individual providers if the hub is unavailable or silent.

When nothing is found it does one more thing that a URL-builder cannot: it asks
UniProt whether the accession is even live.  A merged or deleted accession then
reports ``INACTIVE_ACCESSION`` with its successor, instead of being filed
alongside genuine "this protein has no structure" cases.
"""
from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Iterable, List, Optional

from . import alphafold, beacons, cache, pdbe, uniprot
from .records import (
    FetchResult,
    FetchStatus,
    ResolveResult,
    StructureFetchError,
    StructureRecord,
)


def _check_liveness(accession: str, result: ResolveResult) -> ResolveResult:
    """Turn an empty result into an INACTIVE_ACCESSION verdict where that is why."""
    try:
        state = uniprot.accession_state(accession)
    except StructureFetchError:
        return result  # cannot tell; leave the original verdict alone
    if state["active"]:
        return result
    replaced = state["replaced_by"]
    result.status = FetchStatus.INACTIVE_ACCESSION
    result.detail = (
        f"UniProt reports {accession} as {state['reason'] or 'INACTIVE'}"
        + (f"; replaced by {', '.join(replaced)}" if replaced else "")
    )
    return result


def resolve(
    accession: str,
    *,
    prefer: str = "predicted",
    fallback: bool = True,
    check_inactive: bool = True,
) -> ResolveResult:
    """Every structure available for ``accession``, ranked.

    ``prefer`` is ``'predicted'`` (default -- full-length AlphaFold first, which
    is what showcase figures want), ``'experimental'``, or ``'coverage'``.

    ``fallback`` queries AlphaFold and PDBe directly when 3D-Beacons is
    unreachable or returns nothing.  ``check_inactive`` spends one extra request
    on a miss to distinguish a dead accession from a protein with no model.

    A returned ``FETCH_FAILED`` status means the question is unanswered -- do not
    read it as "no structure exists".
    """
    result = beacons.resolve(accession, prefer=prefer)
    if result.records or not fallback:
        if not result.records and check_inactive and result.status.is_answered:
            return _check_liveness(accession, result)
        return result

    # 3D-Beacons was silent or unreachable -- ask the providers directly.
    records: List[StructureRecord] = []
    errors: List[str] = []
    for source in (alphafold.resolve, lambda a: _pdbe_resolve(a)):
        try:
            sub = source(accession)
        except StructureFetchError as exc:  # pragma: no cover - defensive
            errors.append(str(exc))
            continue
        if sub.status is FetchStatus.INVALID_ACCESSION:
            return sub
        if not sub.status.is_answered:
            errors.append(sub.detail)
        records.extend(sub.records)

    if records:
        return ResolveResult(
            requested=accession,
            accession=records[0].accession or accession,
            records=beacons.rank(records, prefer),
            status=FetchStatus.RESOLVED,
            detail="resolved via provider fallback (3D-Beacons returned nothing)",
        )

    # Nothing anywhere. Was the hub reachable at all?
    if not result.status.is_answered and errors:
        return ResolveResult(
            requested=accession,
            accession=accession,
            status=FetchStatus.FETCH_FAILED,
            detail="; ".join(filter(None, [result.detail] + errors)),
        )
    out = ResolveResult(
        requested=accession,
        accession=accession,
        status=FetchStatus.NO_MODEL_AVAILABLE,
        detail="no model at 3D-Beacons, AlphaFold or PDBe",
    )
    return _check_liveness(accession, out) if check_inactive else out


def _pdbe_resolve(accession: str) -> ResolveResult:
    try:
        recs = pdbe.best_structure_records(accession)
    except StructureFetchError as exc:
        return ResolveResult(
            requested=accession,
            accession=accession,
            status=FetchStatus.FETCH_FAILED,
            detail=str(exc),
        )
    return ResolveResult(
        requested=accession,
        accession=accession,
        records=recs,
        status=FetchStatus.RESOLVED if recs else FetchStatus.NO_MODEL_AVAILABLE,
    )


def fetch_structure(
    accession: str,
    *,
    cache_dir: Optional[os.PathLike] = None,
    prefer: str = "predicted",
    record: Optional[StructureRecord] = None,
    force: bool = False,
    check_inactive: bool = True,
) -> FetchResult:
    """Download one structure to the cache and report exactly what happened.

    Pass ``record`` to fetch a specific candidate from a previous
    :func:`resolve`; otherwise the best-ranked candidate is used.

    Never raises for an unreachable accession -- inspect ``.status``, and treat
    only ``status.is_answered`` results as evidence about coverage.
    """
    if record is None:
        res = resolve(accession, prefer=prefer, check_inactive=check_inactive)
        if not res.records:
            return FetchResult(
                accession=accession, status=res.status, detail=res.detail
            )
        record = res.best()

    path = cache.path_for(record, cache_dir)
    if not force and cache.is_cached(path):
        return FetchResult(
            accession=record.accession or accession,
            status=FetchStatus.CACHED,
            path=path,
            record=record,
        )

    from . import http  # local import keeps the module import graph shallow

    try:
        status, data = http.get_bytes(record.model_url)
    except StructureFetchError as exc:
        return FetchResult(
            accession=accession,
            status=FetchStatus.FETCH_FAILED,
            record=record,
            detail=str(exc),
        )
    if status == 404:
        return FetchResult(
            accession=accession,
            status=FetchStatus.NO_MODEL_AVAILABLE,
            record=record,
            detail=f"model file 404 at {record.model_url}",
        )
    if status >= 400 or not data:
        return FetchResult(
            accession=accession,
            status=FetchStatus.FETCH_FAILED,
            record=record,
            detail=f"HTTP {status} fetching {record.model_url}",
        )

    cache.write_atomic(path, data)
    cache.write_meta(path, record)
    return FetchResult(
        accession=record.accession or accession,
        status=FetchStatus.DOWNLOADED,
        path=path,
        record=record,
    )


def fetch_many(
    accessions: Iterable[str],
    *,
    cache_dir: Optional[os.PathLike] = None,
    prefer: str = "predicted",
    force: bool = False,
    max_workers: int = 4,
    check_inactive: bool = True,
) -> List[FetchResult]:
    """Fetch a batch. Always returns one result per input accession.

    Every accession gets a row whatever happens, so a batch summary can separate
    "no model exists" from "we failed to ask"::

        rows = [r.as_row() for r in fetch_many(accs)]
        unanswered = [r for r in rows if not r["is_answered"]]   # retry these
    """
    accs = [a for a in dict.fromkeys(accessions) if a]
    if not accs:
        return []

    def one(acc: str) -> FetchResult:
        try:
            return fetch_structure(
                acc,
                cache_dir=cache_dir,
                prefer=prefer,
                force=force,
                check_inactive=check_inactive,
            )
        except Exception as exc:  # never let one accession abort the batch
            return FetchResult(
                accession=acc, status=FetchStatus.FETCH_FAILED, detail=repr(exc)
            )

    if max_workers <= 1:
        return [one(a) for a in accs]
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        return list(pool.map(one, accs))


def cached_path(
    accession: str, *, cache_dir: Optional[os.PathLike] = None, prefer: str = "predicted"
) -> Optional[Path]:
    """Local path for an already-cached structure, without touching the network.

    ``None`` when it is not cached -- call :func:`fetch_structure` to get it.
    """
    root = cache.cache_root(cache_dir)
    if not root.is_dir():
        return None
    stem = accession.split("-")[0].upper()
    for provider_dir in sorted(root.iterdir()):
        if not provider_dir.is_dir():
            continue
        for path in sorted(provider_dir.glob(f"*{stem}*")):
            if path.name.endswith(".meta.json") or path.name.endswith(".tmp"):
                continue
            if cache.is_cached(path):
                return path
    return None
