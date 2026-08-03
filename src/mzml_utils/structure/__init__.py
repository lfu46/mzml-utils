"""Protein structure retrieval -- the single place URLs are built.

**Never hand-write an AlphaFold, PDB or UniProt structure URL.**  Use
:func:`resolve` / :func:`fetch_structure` here instead.  Hand-built URLs are what
this package exists to retire, and they are already wrong in two ways:

* AlphaFold **v4 and v5 files have been deleted** -- any interpolated
  ``model_v4``/``model_v5`` is a 404 today.  v6 is current.
* AlphaFold v6 **abolished the ``-F2-``/``-F3-`` fragment scheme**.  Proteins over
  ~2700 aa are absent outright rather than split, so a fragment suffix never
  helps; isoform entries or an experimental structure are the route to them.

Typical use::

    from mzml_utils.structure import resolve, fetch_structure, fetch_many

    res = resolve("P08571")          # every provider, ranked
    res.best().model_url             # correct version, from the provider
    res.alphafold().confidence_avg_local_score

    hit = fetch_structure("P08571")  # download + cache
    if hit:
        parse(hit.path)

    rows = [r.as_row() for r in fetch_many(accessions)]

The status contract, which is the reason this package exists in this shape:
``FetchStatus.NO_MODEL_AVAILABLE`` (a provider answered, nothing exists) and
``FetchStatus.FETCH_FAILED`` (we never got an answer) are different states, and
``status.is_answered`` separates them.  Only answered results are evidence about
coverage.

Provider modules are importable directly when you need something specific:
:mod:`.beacons`, :mod:`.alphafold`, :mod:`.pdbe`, :mod:`.uniprot`,
:mod:`.interpro`.
"""
from __future__ import annotations

from . import alphafold, beacons, cache, http, interpro, pdbe, uniprot
from .cache import ENV_CACHE_ROOT, cache_root, read_meta
from .client import cached_path, fetch_many, fetch_structure, resolve
from .interpro import domains as interpro_domains
from .interpro import domains_at as interpro_domains_at
from .records import (
    ANSWERED_STATES,
    UNANSWERED_STATES,
    FetchResult,
    FetchStatus,
    ResolveResult,
    StructureError,
    StructureFetchError,
    StructureRecord,
)
from .uniprot import accession_state
from .uniprot import features as uniprot_features
from .uniprot import glycosylation_sites

__all__ = [
    # high-level API
    "resolve",
    "fetch_structure",
    "fetch_many",
    "cached_path",
    # annotation helpers
    "uniprot_features",
    "glycosylation_sites",
    "accession_state",
    "interpro_domains",
    "interpro_domains_at",
    # result types
    "StructureRecord",
    "ResolveResult",
    "FetchResult",
    "FetchStatus",
    "ANSWERED_STATES",
    "UNANSWERED_STATES",
    "StructureError",
    "StructureFetchError",
    # cache
    "cache_root",
    "read_meta",
    "ENV_CACHE_ROOT",
    # provider modules
    "beacons",
    "alphafold",
    "pdbe",
    "uniprot",
    "interpro",
    "cache",
    "http",
]
