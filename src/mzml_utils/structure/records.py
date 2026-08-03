"""Result types for structure retrieval.

The important thing in this module is ``FetchStatus``, and specifically that
``NO_MODEL_AVAILABLE`` and ``FETCH_FAILED`` are different states.

The fetchers this package replaces did not make that distinction: failures were
appended to an ``invalid`` list and surfaced only as a count, so "AlphaFold has no
model for this protein" and "the network was down while we asked" were recorded
identically.  Any coverage number computed from such a run is unfalsifiable -- you
cannot tell a real biological gap from a transient outage after the fact.

``is_answered`` is the predicate that keeps them apart.  A caller aggregating
coverage must filter on it; a caller retrying must filter on its negation.  Same
design as ``nglyco_retrieve``'s SEARCHED_STATES / UNSEARCHED_STATES split.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional


class StructureError(Exception):
    """Base class for structure-retrieval errors."""


class StructureFetchError(StructureError):
    """A request failed after the retry policy was exhausted.

    Raised only by the transport layer.  The high-level helpers catch it and
    return ``FetchStatus.FETCH_FAILED`` so a batch run is never aborted by one
    unreachable accession.
    """


class FetchStatus(str, Enum):
    """Outcome of a single accession's retrieval attempt."""

    RESOLVED = "resolved"                      # candidates found; nothing fetched yet
    CACHED = "cached"                          # already on disk; no network used
    DOWNLOADED = "downloaded"                  # retrieved successfully just now
    NO_MODEL_AVAILABLE = "no_model_available"  # provider answered: nothing exists
    INVALID_ACCESSION = "invalid_accession"    # malformed identifier (HTTP 400)
    INACTIVE_ACCESSION = "inactive_accession"  # UniProt says MERGED/DELETED/DEMERGED
    FETCH_FAILED = "fetch_failed"              # transport error after retries

    @property
    def is_answered(self) -> bool:
        """True when a provider actually answered the question.

        ``FETCH_FAILED`` is the only status that means "we still do not know".
        Absence of a model is evidence only when this is True.
        """
        return self is not FetchStatus.FETCH_FAILED

    @property
    def has_structure(self) -> bool:
        """True when a usable local file exists as a result of this attempt."""
        return self in (FetchStatus.CACHED, FetchStatus.DOWNLOADED)


#: Statuses that constitute an answer from a provider.
ANSWERED_STATES = frozenset(s for s in FetchStatus if s.is_answered)
#: Statuses that mean the question is still open and a retry is warranted.
UNANSWERED_STATES = frozenset(s for s in FetchStatus if not s.is_answered)


# Model categories as emitted by 3D-Beacons.
EXPERIMENTAL = "EXPERIMENTALLY DETERMINED"
TEMPLATE_BASED = "TEMPLATE-BASED"
AB_INITIO = "AB-INITIO"
ENSEMBLE = "CONFORMATIONAL ENSEMBLE"


@dataclass
class StructureRecord:
    """One candidate structure for an accession, from any provider.

    Field names follow the 3D-Beacons ``summary`` schema, which the other
    providers are normalised into so a caller never has to branch on provider.
    """

    accession: str
    provider: str
    model_identifier: str
    model_url: str
    model_category: str = ""
    model_format: str = ""
    model_page_url: str = ""
    #: Mean pLDDT for AlphaFold; None for experimental structures.
    confidence_avg_local_score: Optional[float] = None
    confidence_type: str = ""
    coverage: Optional[float] = None
    uniprot_start: Optional[int] = None
    uniprot_end: Optional[int] = None
    experimental_method: str = ""
    resolution: Optional[float] = None
    sequence_identity: Optional[float] = None
    #: Resolved AlphaFold model version, when the provider reported one.
    model_version: Optional[int] = None
    raw: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_experimental(self) -> bool:
        return self.model_category == EXPERIMENTAL

    @property
    def is_alphafold(self) -> bool:
        """True for AlphaFold DB models.

        Matched by prefix because 3D-Beacons reports the provider as
        ``"AlphaFold DB"`` while the AlphaFold API itself says ``"AlphaFold"``.
        An equality test here silently returned False for every federated
        record.  ``AlphaFill`` is a different provider and must not match.
        """
        return self.provider.lower().replace(" ", "").startswith("alphafold")

    @property
    def is_canonical_alphafold(self) -> bool:
        """True for the ``AF-{accession}-F1`` model, not a checksum-keyed twin.

        ``/api/prediction`` and 3D-Beacons both return extra models keyed by
        sequence checksum (``AF-0000000003510504``) alongside the real entry.
        """
        return self.is_alphafold and self.accession.upper() in self.model_identifier.upper()

    @property
    def filename(self) -> str:
        """Basename to cache this model under, taken from its URL."""
        name = self.model_url.split("?", 1)[0].rstrip("/").rsplit("/", 1)[-1]
        return name or f"{self.model_identifier}.cif"


@dataclass
class ResolveResult:
    """Every structure known for one accession, ranked."""

    requested: str
    accession: str                      # after any UniProt 303 redirect
    records: List[StructureRecord] = field(default_factory=list)
    status: FetchStatus = FetchStatus.NO_MODEL_AVAILABLE
    detail: str = ""

    def __bool__(self) -> bool:
        return bool(self.records)

    def __len__(self) -> int:
        return len(self.records)

    def __iter__(self):
        return iter(self.records)

    @property
    def redirected(self) -> bool:
        return self.accession != self.requested

    def best(self) -> Optional[StructureRecord]:
        return self.records[0] if self.records else None

    def alphafold(self) -> Optional[StructureRecord]:
        for r in self.records:
            if r.is_alphafold:
                return r
        return None

    def experimental(self) -> List[StructureRecord]:
        return [r for r in self.records if r.is_experimental]

    def by_provider(self, provider: str) -> List[StructureRecord]:
        p = provider.lower()
        return [r for r in self.records if r.provider.lower() == p]


@dataclass
class FetchResult:
    """Outcome of fetching one accession's structure file to disk."""

    accession: str
    status: FetchStatus
    path: Optional[Path] = None
    record: Optional[StructureRecord] = None
    detail: str = ""

    def __bool__(self) -> bool:
        return self.status.has_structure

    @property
    def is_answered(self) -> bool:
        return self.status.is_answered

    def as_row(self) -> Dict[str, Any]:
        """Flat dict for building a status table over a batch."""
        r = self.record
        return {
            "accession": self.accession,
            "status": self.status.value,
            "is_answered": self.status.is_answered,
            "path": str(self.path) if self.path else "",
            "provider": r.provider if r else "",
            "model_identifier": r.model_identifier if r else "",
            "model_version": r.model_version if r else None,
            "model_url": r.model_url if r else "",
            "mean_plddt": r.confidence_avg_local_score if r else None,
            "detail": self.detail,
        }
