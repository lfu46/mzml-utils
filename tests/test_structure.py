"""Tests for mzml_utils.structure.

Two layers:

* **Offline** -- payload shapes, ranking, cache behaviour.  These pin the
  decisions that live API responses forced, so a refactor cannot quietly undo
  them.
* **Network** (``-m network``) -- four anchors against the live services.  Each
  one corresponds to a way hand-written fetch code has actually broken.  They are
  deselected by default so the suite stays offline and fast.
"""
from __future__ import annotations

import json

import pytest

from mzml_utils.structure import (
    ANSWERED_STATES,
    UNANSWERED_STATES,
    FetchStatus,
    StructureRecord,
    alphafold,
    beacons,
    cache,
    interpro,
    uniprot,
)
from mzml_utils.structure.records import (
    AB_INITIO,
    EXPERIMENTAL,
    TEMPLATE_BASED,
    FetchResult,
)

# --------------------------------------------------------------------------
# The status contract
# --------------------------------------------------------------------------


def test_no_model_available_is_not_fetch_failed():
    """The distinction the whole package exists to preserve.

    A provider saying "nothing exists" is evidence. A transport failure is not.
    Collapsing them is what made the previous fetchers' coverage numbers
    unfalsifiable.
    """
    assert FetchStatus.NO_MODEL_AVAILABLE is not FetchStatus.FETCH_FAILED
    assert FetchStatus.NO_MODEL_AVAILABLE.is_answered is True
    assert FetchStatus.FETCH_FAILED.is_answered is False


def test_fetch_failed_is_the_only_unanswered_state():
    assert UNANSWERED_STATES == {FetchStatus.FETCH_FAILED}
    assert FetchStatus.NO_MODEL_AVAILABLE in ANSWERED_STATES
    assert FetchStatus.INVALID_ACCESSION in ANSWERED_STATES
    assert FetchStatus.INACTIVE_ACCESSION in ANSWERED_STATES


def test_only_cached_and_downloaded_yield_a_file():
    assert FetchStatus.CACHED.has_structure
    assert FetchStatus.DOWNLOADED.has_structure
    for s in (
        FetchStatus.RESOLVED,
        FetchStatus.NO_MODEL_AVAILABLE,
        FetchStatus.FETCH_FAILED,
        FetchStatus.INVALID_ACCESSION,
        FetchStatus.INACTIVE_ACCESSION,
    ):
        assert not s.has_structure


def test_batch_row_exposes_is_answered():
    row = FetchResult(accession="P1", status=FetchStatus.FETCH_FAILED).as_row()
    assert row["is_answered"] is False
    assert row["status"] == "fetch_failed"


# --------------------------------------------------------------------------
# Record semantics
# --------------------------------------------------------------------------


def _rec(**kw):
    base = dict(
        accession="P08571",
        provider="AlphaFold DB",
        model_identifier="AF-P08571-F1",
        model_url="https://alphafold.ebi.ac.uk/files/AF-P08571-F1-model_v6.cif",
    )
    base.update(kw)
    return StructureRecord(**base)


def test_is_alphafold_matches_the_federated_provider_string():
    """3D-Beacons says "AlphaFold DB"; the AlphaFold API says "AlphaFold".

    An equality test against "alphafold" returned False for every federated
    record, which silently dropped AlphaFold to the bottom of the ranking.
    """
    assert _rec(provider="AlphaFold DB").is_alphafold
    assert _rec(provider="AlphaFold").is_alphafold


def test_alphafill_is_not_alphafold():
    assert not _rec(provider="AlphaFill").is_alphafold


def test_canonical_alphafold_excludes_the_checksum_twin():
    assert _rec(model_identifier="AF-P08571-F1").is_canonical_alphafold
    assert not _rec(model_identifier="AF-0000000003510504").is_canonical_alphafold


def test_filename_comes_from_the_url_including_the_version():
    assert _rec().filename == "AF-P08571-F1-model_v6.cif"


def test_filename_strips_query_parameters():
    r = _rec(
        provider="SWISS-MODEL",
        model_url="https://swissmodel.expasy.org/3d-beacons/uniprot/P08571.cif?range=26-335",
    )
    assert r.filename == "P08571.cif"


# --------------------------------------------------------------------------
# 3D-Beacons parsing and ranking
# --------------------------------------------------------------------------

BEACONS_PAYLOAD = {
    "uniprot_entry": {"ac": "P08571", "id": "CD14_HUMAN", "sequence_length": 375},
    "structures": [
        {
            "summary": {
                "model_identifier": "P08571_26-335:4glp.1.A",
                "provider": "SWISS-MODEL",
                "model_category": TEMPLATE_BASED,
                "model_url": "https://swissmodel.expasy.org/3d-beacons/uniprot/P08571.cif",
                "confidence_avg_local_score": 0.742,
                "coverage": 0.827,
            }
        },
        {
            "summary": {
                "model_identifier": "4glp",
                "provider": "PDBe",
                "model_category": EXPERIMENTAL,
                "model_url": "https://www.ebi.ac.uk/pdbe/static/entry/4glp_updated.cif",
                "coverage": 0.871,
                "resolution": 4.002,
            }
        },
        {
            "summary": {
                "model_identifier": "AF-P08571-F1",
                "provider": "AlphaFold DB",
                "model_category": AB_INITIO,
                "model_url": "https://alphafold.ebi.ac.uk/files/AF-P08571-F1-model_v6.cif",
                "confidence_avg_local_score": 84.75,
                "coverage": 1.0,
            }
        },
        {
            "summary": {
                "model_identifier": "AF-0000000003510504",
                "provider": "AlphaFold DB",
                "model_category": AB_INITIO,
                "model_url": "https://alphafold.ebi.ac.uk/files/AF-0000000003510504-model_v1.cif",
                "confidence_avg_local_score": 82.25,
                "coverage": 1.0,
            }
        },
    ],
}


def test_parse_summary_reads_the_resolved_accession():
    recs = beacons.parse_summary("P08571", BEACONS_PAYLOAD)
    assert len(recs) == 4
    assert {r.accession for r in recs} == {"P08571"}


def test_empty_payload_is_no_models_not_an_error():
    """A missing accession returns {} with HTTP 200, not a 404."""
    assert beacons.parse_summary("P99999", {}) == []


def test_model_version_is_parsed_from_the_url():
    assert beacons.parse_model_version(
        "https://alphafold.ebi.ac.uk/files/AF-P08571-F1-model_v6.cif"
    ) == 6
    assert beacons.parse_model_version("https://example.org/thing.cif") is None


def test_predicted_ranking_puts_canonical_alphafold_first():
    ranked = beacons.rank(beacons.parse_summary("P08571", BEACONS_PAYLOAD), "predicted")
    assert ranked[0].model_identifier == "AF-P08571-F1"
    assert ranked[1].model_identifier == "AF-0000000003510504"


def test_experimental_ranking_puts_pdbe_first():
    ranked = beacons.rank(
        beacons.parse_summary("P08571", BEACONS_PAYLOAD), "experimental"
    )
    assert ranked[0].provider == "PDBe"


def test_sasbdb_envelopes_are_demoted_below_real_structures():
    """SAS bead models are tagged EXPERIMENTALLY DETERMINED with coverage 1.0.

    Against titin that put five SAXS envelopes above 66 crystal structures. They
    are not residue-resolved, so nothing can be mapped onto them.
    """
    records = [
        StructureRecord(
            accession="Q8WZ42",
            provider="SASBDB",
            model_identifier="6642",
            model_url="https://example.org/6642.pdb",
            model_category=EXPERIMENTAL,
            coverage=1.0,
        ),
        StructureRecord(
            accession="Q8WZ42",
            provider="PDBe",
            model_identifier="8g4l",
            model_url="https://example.org/8g4l.cif",
            model_category=EXPERIMENTAL,
            coverage=0.032,
            resolution=6.4,
        ),
    ]
    for prefer in ("predicted", "experimental", "coverage"):
        assert beacons.rank(records, prefer)[0].provider == "PDBe", prefer


def test_rank_rejects_an_unknown_preference():
    with pytest.raises(ValueError, match="prefer must be"):
        beacons.rank([], "whatever")


# --------------------------------------------------------------------------
# AlphaFold field handling
# --------------------------------------------------------------------------


def test_pick_prefers_the_new_field_name():
    assert alphafold._pick({"modelEntityId": "new", "entryId": "old"}, "modelEntityId") == "new"


def test_pick_falls_back_to_the_legacy_name():
    """The rename's sunset has passed but both sets are still served.

    Reading only one set breaks whenever the other disappears.
    """
    assert alphafold._pick({"entryId": "old"}, "modelEntityId") == "old"
    assert alphafold._pick({"uniprotStart": 26}, "sequenceStart") == 26


def test_pick_resolves_a_legacy_key_against_a_new_payload():
    assert alphafold._pick({"sequenceStart": 26}, "uniprotStart") == 26


def test_pick_returns_the_default_when_absent():
    assert alphafold._pick({}, "modelEntityId", "fallback") == "fallback"


def test_alphafold_record_carries_the_reported_version():
    rec = alphafold._to_record(
        {
            "uniprotAccession": "P08571",
            "modelEntityId": "AF-P08571-F1",
            "cifUrl": "https://alphafold.ebi.ac.uk/files/AF-P08571-F1-model_v6.cif",
            "globalMetricValue": 84.75,
            "latestVersion": 6,
        }
    )
    assert rec.model_version == 6
    assert rec.confidence_avg_local_score == 84.75
    assert rec.is_canonical_alphafold


# --------------------------------------------------------------------------
# UniProt features
# --------------------------------------------------------------------------

UNIPROT_ENTRY = {
    "primaryAccession": "P08571",
    "sequence": {"value": "MERASCLL"},
    "features": [
        {
            "type": "Glycosylation",
            "description": "N-linked (GlcNAc...) asparagine",
            "location": {
                "start": {"value": 37, "modifier": "EXACT"},
                "end": {"value": 37, "modifier": "EXACT"},
            },
        },
        {
            "type": "Glycosylation",
            "description": "O-linked (GalNAc...) threonine",
            "location": {
                "start": {"value": 100, "modifier": "EXACT"},
                "end": {"value": 100, "modifier": "EXACT"},
            },
        },
        {
            "type": "Glycosylation",
            "description": "N-linked (GlcNAc...) asparagine",
            "location": {
                "start": {"value": 200, "modifier": "UNKNOWN"},
                "end": {"value": 200, "modifier": "EXACT"},
            },
        },
        {
            "type": "Signal",
            "description": "",
            "location": {
                "start": {"value": 1, "modifier": "EXACT"},
                "end": {"value": 19, "modifier": "EXACT"},
            },
        },
    ],
    "uniProtKBCrossReferences": [
        {
            "database": "PDB",
            "id": "4GLP",
            "properties": [
                {"key": "Method", "value": "X-ray"},
                {"key": "Chains", "value": "A=26-335"},
            ],
        },
        {"database": "AlphaFoldDB", "id": "P08571", "properties": []},
    ],
}


def test_linkage_is_read_from_description_not_type():
    """Every site has type "Glycosylation"; only description says N or O."""
    n_sites = uniprot.glycosylation_sites("P08571", entry_json=UNIPROT_ENTRY)
    assert [s["start"] for s in n_sites] == [37, 200]
    o_sites = uniprot.glycosylation_sites("P08571", linkage="O", entry_json=UNIPROT_ENTRY)
    assert [s["start"] for s in o_sites] == [100]


def test_uncertain_positions_are_flagged_not_silently_used():
    sites = uniprot.glycosylation_sites("P08571", entry_json=UNIPROT_ENTRY)
    certain = {s["start"]: s["certain"] for s in sites}
    assert certain[37] is True
    assert certain[200] is False


def test_glycosylation_sites_rejects_a_bad_linkage():
    with pytest.raises(ValueError, match="linkage must be"):
        uniprot.glycosylation_sites("P08571", linkage="S", entry_json=UNIPROT_ENTRY)


def test_features_can_be_filtered_by_type():
    feats = uniprot.features("P08571", types=["Signal"], entry_json=UNIPROT_ENTRY)
    assert len(feats) == 1 and feats[0]["end"] == 19


def test_cross_reference_properties_are_flattened():
    xrefs = uniprot.cross_references("P08571", "PDB", entry_json=UNIPROT_ENTRY)
    assert xrefs[0]["properties"]["Chains"] == "A=26-335"
    assert uniprot.pdb_ids("P08571", entry_json=UNIPROT_ENTRY) == ["4GLP"]


# --------------------------------------------------------------------------
# InterPro fragments
# --------------------------------------------------------------------------

INTERPRO_RESULTS = [
    {
        "metadata": {
            "accession": "IPR000483",
            "name": "Cys-rich flank region",
            "type": "domain",
            "source_database": "interpro",
        },
        "proteins": [
            {
                "accession": "p08571",
                "protein_length": 375,
                "in_alphafold": True,
                "entry_protein_locations": [
                    {
                        "fragments": [
                            {"start": 30, "end": 60, "dc-status": "CONTINUOUS"},
                            {"start": 90, "end": 120, "dc-status": "CONTINUOUS"},
                        ]
                    }
                ],
            }
        ],
    }
]


def test_discontinuous_domains_keep_their_fragments(monkeypatch):
    monkeypatch.setattr(interpro, "entries", lambda *a, **k: INTERPRO_RESULTS)
    doms = interpro.domains("P08571")
    assert len(doms) == 1
    d = doms[0]
    assert d["discontinuous"] is True
    assert d["start"] == 30 and d["end"] == 120        # envelope
    assert [(f["start"], f["end"]) for f in d["fragments"]] == [(30, 60), (90, 120)]


def test_residue_in_the_gap_is_not_reported_as_in_the_domain(monkeypatch):
    """The envelope spans 30-120 but residue 75 sits in the gap.

    Collapsing fragments to min/max would claim the domain covers it.
    """
    monkeypatch.setattr(interpro, "entries", lambda *a, **k: INTERPRO_RESULTS)
    assert interpro.domains_at("P08571", 45) != []
    assert interpro.domains_at("P08571", 100) != []
    assert interpro.domains_at("P08571", 75) == []


# --------------------------------------------------------------------------
# Cache
# --------------------------------------------------------------------------


def test_cache_root_prefers_the_explicit_argument(tmp_path, monkeypatch):
    monkeypatch.setenv(cache.ENV_CACHE_ROOT, str(tmp_path / "from_env"))
    assert cache.cache_root(tmp_path / "explicit") == tmp_path / "explicit"
    assert cache.cache_root() == tmp_path / "from_env"


def test_cache_root_falls_back_to_the_default(monkeypatch):
    monkeypatch.delenv(cache.ENV_CACHE_ROOT, raising=False)
    assert cache.cache_root() == cache.DEFAULT_CACHE_ROOT


def test_atomic_write_leaves_no_temp_file(tmp_path):
    path = tmp_path / "sub" / "AF-P08571-F1-model_v6.cif"
    cache.write_atomic(path, b"data:loop_")
    assert path.read_bytes() == b"data:loop_"
    assert list(tmp_path.rglob("*.tmp")) == []


def test_sidecar_records_the_version_the_file_came_from(tmp_path):
    """A bare directory of .cif files cannot say which version it holds.

    Once v4 and v5 were deleted upstream, that became unanswerable after the
    fact -- so provenance is written next to the data.
    """
    path = tmp_path / "AF-P08571-F1-model_v6.cif"
    cache.write_atomic(path, b"x")
    cache.write_meta(path, _rec(model_version=6))
    meta = cache.read_meta(path)
    assert meta["model_version"] == 6
    assert meta["model_url"].endswith("model_v6.cif")
    assert meta["fetched_at"]


def test_zero_byte_files_count_as_absent(tmp_path):
    path = tmp_path / "empty.cif"
    path.write_bytes(b"")
    assert not cache.is_cached(path)


def test_read_meta_is_none_without_a_sidecar(tmp_path):
    path = tmp_path / "x.cif"
    path.write_bytes(b"x")
    assert cache.read_meta(path) is None


def test_path_for_separates_providers(tmp_path):
    af = cache.path_for(_rec(), tmp_path)
    pdb = cache.path_for(
        _rec(provider="PDBe", model_url="https://example.org/4glp_updated.cif"), tmp_path
    )
    assert af.parent.name == "alphafold_db"
    assert pdb.parent.name == "pdbe"


# --------------------------------------------------------------------------
# Network anchors -- run with: pytest -m network
# --------------------------------------------------------------------------


@pytest.mark.network
def test_alphafold_version_is_read_not_assumed():
    """v4 and v5 files are deleted upstream; a hard-coded version is a 404."""
    version = alphafold.latest_version("P08571")
    assert isinstance(version, int) and version >= 6


@pytest.mark.network
def test_resolve_returns_a_versioned_alphafold_url():
    from mzml_utils.structure import resolve

    rec = resolve("P08571").alphafold()
    assert rec is not None
    assert rec.model_version == rec.model_version and rec.model_version >= 6
    assert f"model_v{rec.model_version}.cif" in rec.model_url


@pytest.mark.network
def test_checksum_twin_is_filtered_out():
    """isUniProt does NOT discriminate -- both entries carry it. Entry id does."""
    ids = [e.get("entryId") or e.get("modelEntityId") for e in alphafold.prediction("P08571")]
    assert ids == ["AF-P08571-F1"]
    unfiltered = alphafold.prediction("P08571", canonical_only=False)
    assert len(unfiltered) > 1


@pytest.mark.network
def test_oversized_protein_is_answered_as_absent_not_failed():
    """Titin exceeds the AlphaFold size ceiling; v6 has no fragments for it.

    The critical part is that this is an *answer*, so it can be counted as a
    real coverage gap rather than a retryable error.
    """
    res = alphafold.resolve("Q8WZ42")
    assert res.status is FetchStatus.NO_MODEL_AVAILABLE
    assert res.status.is_answered


@pytest.mark.network
def test_federation_rescues_a_protein_alphafold_lacks():
    """The reason to resolve through 3D-Beacons rather than AlphaFold alone."""
    from mzml_utils.structure import resolve

    res = resolve("Q8WZ42")
    assert res.alphafold() is None
    assert res.records, "PDBe holds structures for titin"
    assert res.best().provider == "PDBe"


@pytest.mark.network
def test_merged_accession_is_reported_with_its_successor():
    state = uniprot.accession_state("Q9UEV2")
    assert state["active"] is False
    assert state["reason"] == "MERGED"
    assert "Q13683" in state["replaced_by"]


@pytest.mark.network
def test_merged_accession_resolves_through_the_redirect():
    from mzml_utils.structure import resolve

    res = resolve("Q9UEV2")
    assert res.redirected
    assert res.accession == "Q13683"


@pytest.mark.network
def test_malformed_identifier_is_distinct_from_absent():
    from mzml_utils.structure import resolve

    res = resolve("ZZZZZZ")
    assert res.status is FetchStatus.INVALID_ACCESSION
    assert res.status is not FetchStatus.NO_MODEL_AVAILABLE


@pytest.mark.network
def test_fetch_writes_a_file_and_its_provenance(tmp_path):
    from mzml_utils.structure import fetch_structure

    res = fetch_structure("P08571", cache_dir=tmp_path)
    assert res.status is FetchStatus.DOWNLOADED
    assert res.path.stat().st_size > 10_000
    assert res.path.read_bytes().lstrip().startswith(b"data_")
    meta = cache.read_meta(res.path)
    assert meta["model_version"] == res.record.model_version
    # second call must not re-download
    again = fetch_structure("P08571", cache_dir=tmp_path)
    assert again.status is FetchStatus.CACHED


@pytest.mark.network
def test_batch_returns_one_row_per_accession(tmp_path):
    from mzml_utils.structure import fetch_many

    accs = ["P08571", "ZZZZZZ"]
    rows = [r.as_row() for r in fetch_many(accs, cache_dir=tmp_path, max_workers=2)]
    assert len(rows) == len(accs)
    assert all("is_answered" in r for r in rows)
    assert json.dumps(rows)  # rows must be JSON-serialisable for status tables
