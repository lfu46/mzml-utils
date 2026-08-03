"""On-disk cache for downloaded structure files.

Structure files are large, immutable for a given model version, and freely
re-downloadable, so the cache is a plain content directory with no expiry.  What
it does add is a **sidecar** recording which model version, provider and URL each
file actually came from.

That sidecar is the point.  A bare directory of ``P08571.cif`` files cannot
answer "which AlphaFold version is this, and when did we get it?" -- and once
v4 and v5 were deleted upstream, an undocumented local file became impossible to
place.  Every write records its provenance next to the data.

Layout::

    <root>/<provider>/<filename-from-url>
    <root>/<provider>/<filename-from-url>.meta.json

Root resolution: ``WU_LAB_STRUCTURE_CACHE`` env var, else
``~/.cache/wu-lab/structures``.
"""
from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional

from .records import StructureRecord

ENV_CACHE_ROOT = "WU_LAB_STRUCTURE_CACHE"
DEFAULT_CACHE_ROOT = Path.home() / ".cache" / "wu-lab" / "structures"


def cache_root(cache_dir: Optional[os.PathLike] = None) -> Path:
    """Resolve the cache root: explicit argument, then env var, then default."""
    if cache_dir is not None:
        return Path(cache_dir).expanduser()
    env = os.environ.get(ENV_CACHE_ROOT, "").strip()
    if env:
        return Path(env).expanduser()
    return DEFAULT_CACHE_ROOT


def path_for(
    record: StructureRecord, cache_dir: Optional[os.PathLike] = None
) -> Path:
    """Where this model's file belongs in the cache."""
    provider = (record.provider or "unknown").lower().replace(" ", "_")
    return cache_root(cache_dir) / provider / record.filename


def meta_path(path: Path) -> Path:
    """Sidecar path for a cached file."""
    return path.with_name(path.name + ".meta.json")


def write_atomic(path: Path, data: bytes) -> Path:
    """Write via a temp file in the same directory, then rename.

    Same pattern the existing fetchers use: a killed download can never leave a
    truncated file that a later run mistakes for a complete one.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_bytes(data)
    tmp.replace(path)
    return path


def write_meta(path: Path, record: StructureRecord, extra: Optional[Dict[str, Any]] = None) -> Path:
    """Record where a cached file came from."""
    meta: Dict[str, Any] = {
        "accession": record.accession,
        "provider": record.provider,
        "model_identifier": record.model_identifier,
        "model_version": record.model_version,
        "model_category": record.model_category,
        "model_url": record.model_url,
        "confidence_avg_local_score": record.confidence_avg_local_score,
        "fetched_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "bytes": path.stat().st_size if path.exists() else None,
    }
    if extra:
        meta.update(extra)
    mp = meta_path(path)
    tmp = mp.with_name(mp.name + ".tmp")
    tmp.write_text(json.dumps(meta, indent=2, sort_keys=True))
    tmp.replace(mp)
    return mp


def read_meta(path: Path) -> Optional[Dict[str, Any]]:
    """Provenance for a cached file, or ``None`` if it was written without one."""
    mp = meta_path(path)
    if not mp.is_file():
        return None
    try:
        return json.loads(mp.read_text())
    except (ValueError, OSError):
        return None


def is_cached(path: Path) -> bool:
    """True when a non-empty file exists at ``path``.

    Zero-byte files count as absent -- they are the signature of an interrupted
    write from a pre-atomic-write era, and re-fetching one is cheap.
    """
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False
