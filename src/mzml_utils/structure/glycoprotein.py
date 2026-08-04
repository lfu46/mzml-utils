"""Build a glycosylated protein model, and refuse to guess which glycan it is.

Two jobs, in this order:

1. **Turn a composition into a structure -- or refuse.**  Mass spectrometry
   reports a *composition* (``HexNAc(4)Hex(5)Fuc(1)NeuAc(2)``).  A 3D model needs
   a *structure*, and the mapping is one-to-many: GlycoShape's Asn library holds
   599 structures spanning only 270 distinct compositions.  ``N2H9`` has exactly
   one structure; ``N4H5F1A2`` has two, differing only in sialic-acid linkage;
   ``N5H6F1A1`` has thirteen.  :func:`resolve_glycan` raises
   :class:`AmbiguousGlycanError` rather than picking one, because a figure built
   from a silent guess is indistinguishable from a figure built from evidence.

   (Those counts are printed by the test suite, not typed from memory -- an
   earlier hand parse of the same library gave 216 and 16 and was wrong.)

2. **Build the model.**  :func:`build_glycoprotein` drives GlycoShape's Re-Glyco
   over an AlphaFold or PDB scaffold and returns the glycosylated coordinates
   plus the per-site steric report.

Reference: Ives, Singh, D'Andrea, Fogarty, Harbison, Satheesan, Tropea & Fadda,
*Nat Methods* 2024;21:2117-2127 (doi:10.1038/s41592-024-02464-7).  Data are CC
BY-NC-ND 4.0 for academic use.

This module never imports ``nglyco_*``: mzml-utils sits *below* those packages in
the install order.  Classification, plausibility and StrucGP decoding are
optional **injected callables** -- pass them and the report gets richer, omit
them and the structural ambiguity check still runs.
"""
from __future__ import annotations

import json
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from . import cache, http
from .records import StructureError

BASE = "https://glycoshape.org"
CONFIGURATIONS_URL = BASE + "/api/configurations"
SESSIONS_URL = BASE + "/api/sessions"
JOBS_URL = BASE + "/api/jobs"

CACHE_SUBDIR = "glycoshape"
LIBRARY_FILENAME = "library.json"

#: IUPAC-condensed base sugar -> composition key.  The first five keys are the
#: vocabulary the ``nglyco_*`` packages already speak; the rest exist so that
#: non-mammalian and GAG sugars parse into a composition that simply *cannot*
#: equal a human N-glycan's, instead of failing to parse.
#:
#: That distinction matters.  An entry that fails to parse would have to be
#: skipped, and a skipped entry silently shrinks the candidate set -- which is
#: exactly the under-reporting :class:`AmbiguousGlycanError` exists to prevent.
#: 174 of the library's 2061 entries contain these sugars.
_BASE_TO_COMP = {
    "Neu5Ac": "NeuAc", "Neu5Gc": "NeuGc",
    "GlcNAc": "HexNAc", "GalNAc": "HexNAc", "ManNAc": "HexNAc",
    "Man": "Hex", "Gal": "Hex", "Glc": "Hex", "Galf": "Hex", "Manf": "Hex",
    "Fuc": "Fuc",
    "GlcN": "HexN", "GalN": "HexN",
    "GlcA": "HexA", "IdoA": "HexA", "GalA": "HexA",
    "Xyl": "Pent", "Ara": "Pent", "Araf": "Pent",
    "Rha": "dHex", "Kdn": "Kdn",
}
_MOD_TO_COMP = {"S": "Sulfate", "P": "Phospho", "Me": "Methyl", "M": "Methyl",
                "Ac": "Acetyl", "Cho": "Choline"}

# Longest-first so 'GlcNAc' never matches as 'GlcN', and 'Araf' never as 'Ara'.
_BASES = sorted(_BASE_TO_COMP, key=len, reverse=True)
# A modifier run is a leading bare 'S' (GlcNS) and/or digit-led groups (6S, 2Me3M).
_TOKEN_RE = re.compile(
    r"(?:[DL]-)?"
    r"(?P<base>" + "|".join(re.escape(b) for b in _BASES) + r")"
    r"(?P<mod>S?(?:\d+(?:Me|Ac|S|P|M))*(?:Cho)?)"
)
_MOD_PART_RE = re.compile(r"\d*(Me|Ac|Cho|S|P|M)")
#: Linkage groups and branch brackets separate monosaccharides and carry none of
#: their own.  They are SPLIT ON, never deleted: deleting them would join
#: 'Gal' + 'Neu5Ac' into 'GalN' + 'eu5Ac' and invent a sugar that is not there.
_SEP_RE = re.compile(r"\([^)]*\)|[\[\]]")

_LINKAGE_RE = re.compile(r"\(([ab])(\d)-(\d)\)")


class GlycoShapeError(StructureError):
    """Re-Glyco or the GlycoShape database refused a request."""


class AmbiguousGlycanError(StructureError):
    """A composition matched more than one structure in the library.

    Carries every surviving candidate so the caller can name one, rather than
    reporting only that a choice is needed.
    """

    def __init__(self, composition: Dict[str, int], candidates: Sequence["GlycanCandidate"],
                 notes: Optional[Sequence[str]] = None):
        self.composition = dict(composition)
        self.candidates = list(candidates)
        self.notes = list(notes or [])
        super().__init__(self._render())

    def _render(self) -> str:
        comp = _comp_repr(self.composition)
        lines = [
            f"{comp} matches {len(self.candidates)} structures in the GlycoShape "
            "library; composition alone cannot choose between them.",
            "Name one with glytoucan=... or glycoshape_id=..., or supply a "
            "StrucGP structure_coding to narrow it.",
            "",
        ]
        for c in self.candidates:
            lines.append(f"  {c.glytoucan}  ({c.glycoshape_id})  {c.iupac}")
        differing = _differing_linkages(self.candidates)
        if differing:
            lines.append("")
            lines.append("They differ only in: " + ", ".join(sorted(differing)))
        lines.append("")
        lines.append(
            "Composition can never settle core (a1-6, FUT8) vs antennary "
            "(a1-3/4) fucose, nor a2-3 vs a2-6 sialylation."
        )
        lines.extend(self.notes)
        return "\n".join(lines)


class NoSuchGlycanError(StructureError):
    """No structure in the library carries this composition."""


@dataclass(frozen=True)
class GlycanCandidate:
    """One 3D structure from the GlycoShape database."""

    glycoshape_id: str
    glytoucan: str
    iupac: str
    mass: float
    composition: Dict[str, int] = field(default_factory=dict)

    @property
    def linkages(self) -> Tuple[str, ...]:
        """Every linkage descriptor in the IUPAC string, e.g. ``('a2-6', 'b1-4')``."""
        return tuple(sorted({f"{a}{i}-{j}" for a, i, j in _LINKAGE_RE.findall(self.iupac)}))

    @property
    def has_core_fucose(self) -> bool:
        """Fuc a1-6 on the reducing GlcNAc -- the FUT8 product."""
        return "Fuc(a1-6)" in self.iupac

    @property
    def antenna_count(self) -> int:
        """Antennary GlcNAc count: total HexNAc beyond the chitobiose core.

        Bisecting GlcNAc inflates this; it is a screening filter for narrowing
        candidates, not a structural claim.
        """
        return max(0, self.composition.get("HexNAc", 0) - 2)


@dataclass
class GlycanChoice:
    """A resolved, unambiguous glycan, with how confident that resolution is."""

    candidate: GlycanCandidate
    tier: int                       # 1 = structure-confirmed, 2 = composition-only
    glycan_type: str = ""           # from an injected classifier, when given
    plausibility_flags: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    @property
    def glytoucan(self) -> str:
        return self.candidate.glytoucan


def _comp_repr(comp: Dict[str, int]) -> str:
    """Short form for messages only.  Canonical short forms belong to nglyco_notation."""
    order = ("HexNAc", "Hex", "Fuc", "NeuAc", "NeuGc")
    short = {"HexNAc": "N", "Hex": "H", "Fuc": "F", "NeuAc": "A", "NeuGc": "G"}
    parts = [f"{short[k]}{comp[k]}" for k in order if comp.get(k)]
    extra = [f"{k}{v}" for k, v in sorted(comp.items()) if k not in order and v]
    return "".join(parts + extra) or "(empty)"


def _differing_linkages(candidates: Sequence[GlycanCandidate]) -> set:
    """Linkages present in some candidates but not all -- what actually differs."""
    if len(candidates) < 2:
        return set()
    sets = [set(c.linkages) for c in candidates]
    union, inter = set().union(*sets), set(sets[0]).intersection(*sets[1:])
    return union - inter


def iupac_to_composition(iupac: str) -> Dict[str, int]:
    """Count monosaccharides in an IUPAC-condensed glycan string.

    The existing parsers (``nglyco_notation.parse_glycan_comp``) read
    ``N4H5F1A2`` and ``HexNAc(4)Hex(5)``; neither reads IUPAC, so this is a new
    operation rather than a second copy of one.

    The whole string must be consumed.  An unrecognised sugar raises rather than
    being skipped, because a skipped sugar collapses two different glycans onto
    one composition and quietly makes an ambiguous choice look unambiguous.
    """
    if not iupac:
        return {}
    comp: Dict[str, int] = {}
    for segment in _SEP_RE.split(iupac):
        if not segment:
            continue
        m = _TOKEN_RE.fullmatch(segment)
        if not m:
            raise ValueError(
                f"unrecognised monosaccharide {segment!r} in {iupac!r}; add its "
                "base to _BASE_TO_COMP rather than letting it go uncounted"
            )
        key = _BASE_TO_COMP[m.group("base")]
        comp[key] = comp.get(key, 0) + 1
        for mod in _MOD_PART_RE.findall(m.group("mod") or ""):
            mkey = _MOD_TO_COMP[mod]
            comp[mkey] = comp.get(mkey, 0) + 1
    return comp


# --------------------------------------------------------------------------
# The glycan library
# --------------------------------------------------------------------------

def library_path(cache_dir: Optional[Any] = None) -> Path:
    return cache.cache_root(cache_dir) / CACHE_SUBDIR / LIBRARY_FILENAME


def fetch_library(*, cache_dir: Optional[Any] = None, force: bool = False) -> Dict[str, List[dict]]:
    """The GlycoShape structure library, keyed by attachment residue.

    Cached on disk so ambiguity checks are offline and reproducible: the answer
    to "how many structures share this composition" must not depend on whether
    the network was up.
    """
    path = library_path(cache_dir)
    if not force and cache.is_cached(path):
        try:
            return json.loads(path.read_text())["configurations"]
        except (ValueError, OSError, KeyError):
            pass  # fall through and refetch

    status, payload = http.get_json(CONFIGURATIONS_URL)
    if status != 200 or not isinstance(payload, dict) or "configurations" not in payload:
        raise GlycoShapeError(
            f"GlycoShape /api/configurations returned {status}; cannot build the "
            "glycan library, so ambiguity cannot be checked"
        )
    cache.write_atomic(path, json.dumps(payload).encode("utf-8"))
    return payload["configurations"]


def candidates_for(
    composition: Dict[str, int],
    *,
    residue: str = "ASN",
    library: Optional[Dict[str, List[dict]]] = None,
    cache_dir: Optional[Any] = None,
) -> List[GlycanCandidate]:
    """Every library structure whose composition equals ``composition``."""
    lib = library if library is not None else fetch_library(cache_dir=cache_dir)
    entries = lib.get(residue.upper())
    if entries is None:
        raise NoSuchGlycanError(
            f"no glycan library for residue {residue!r}; have {sorted(lib)}"
        )
    want = {k: v for k, v in composition.items() if v}
    out, unparsed = [], []
    for e in entries:
        try:
            comp = iupac_to_composition(e["iupac"])
        except ValueError as exc:
            unparsed.append(f"{e.get('glytoucan', e.get('ID'))}: {exc}")
            continue
        if {k: v for k, v in comp.items() if v} == want:
            out.append(
                GlycanCandidate(
                    glycoshape_id=e["ID"],
                    glytoucan=e["glytoucan"],
                    iupac=e["iupac"],
                    mass=float(e.get("mass", 0.0)),
                    composition=comp,
                )
            )
    if unparsed:
        # Fail closed.  Any unparsed entry means the candidate set is not
        # provably complete, and an incomplete set can turn a genuinely
        # ambiguous composition into an apparently unambiguous one.
        raise GlycoShapeError(
            f"{len(unparsed)} {residue} library entries did not parse, so the "
            f"candidate set for {_comp_repr(composition)} cannot be trusted:\n  "
            + "\n  ".join(unparsed[:5])
            + ("\n  ..." if len(unparsed) > 5 else "")
        )
    return out


def _find_by_id(ident: str, *, residue: str, library, cache_dir) -> Optional[GlycanCandidate]:
    lib = library if library is not None else fetch_library(cache_dir=cache_dir)
    for e in lib.get(residue.upper(), []):
        if ident in (e.get("glytoucan"), e.get("ID")):
            try:
                comp = iupac_to_composition(e["iupac"])
            except ValueError:
                comp = {}
            return GlycanCandidate(
                glycoshape_id=e["ID"], glytoucan=e["glytoucan"],
                iupac=e["iupac"], mass=float(e.get("mass", 0.0)), composition=comp,
            )
    return None


def resolve_glycan(
    composition: Optional[Dict[str, int]] = None,
    *,
    residue: str = "ASN",
    glytoucan: Optional[str] = None,
    glycoshape_id: Optional[str] = None,
    strucgp_code: Optional[str] = None,
    library: Optional[Dict[str, List[dict]]] = None,
    cache_dir: Optional[Any] = None,
    classify: Optional[Callable] = None,
    plausibility: Optional[Callable] = None,
    parse_strucgp: Optional[Callable] = None,
) -> GlycanChoice:
    """Pin a composition to exactly one 3D structure, or raise.

    Narrowing, in order:

    1. an explicit ``glytoucan=`` / ``glycoshape_id=`` -- unambiguous by
       construction, tier 1;
    2. ``strucgp_code=`` with ``parse_strucgp=`` -- real decoded topology
       (antenna count, core fucose, bisect) filters the candidates, tier 1;
    3. composition alone -- tier 2, and **raises** if more than one survives.

    ``classify`` and ``plausibility`` are the ``nglyco_*`` callables; inject them
    to get the glycan type and the R1-R13 flags in the result.
    """
    explicit = glytoucan or glycoshape_id
    if explicit:
        cand = _find_by_id(explicit, residue=residue, library=library, cache_dir=cache_dir)
        if cand is None:
            raise NoSuchGlycanError(
                f"{explicit!r} is not in the GlycoShape {residue} library"
            )
        return _finish(cand, tier=1, composition=cand.composition,
                       classify=classify, plausibility=plausibility,
                       strucgp_code=strucgp_code,
                       notes=["structure named explicitly"])

    if not composition:
        raise ValueError("give either a composition or glytoucan=/glycoshape_id=")

    cands = candidates_for(composition, residue=residue, library=library, cache_dir=cache_dir)
    if not cands:
        raise NoSuchGlycanError(
            f"no {residue} structure in the GlycoShape library has composition "
            f"{_comp_repr(composition)}; it cannot be rendered"
        )

    notes: List[str] = []
    tier = 2
    if strucgp_code and parse_strucgp is not None:
        decoded = parse_strucgp(strucgp_code)
        narrowed = [c for c in cands if _matches_topology(c, decoded)]
        if narrowed:
            if len(narrowed) < len(cands):
                notes.append(
                    f"StrucGP {strucgp_code} narrowed {len(cands)} candidates to {len(narrowed)}"
                )
            cands, tier = narrowed, 1
        else:
            notes.append(
                f"StrucGP {strucgp_code} matched none of the {len(cands)} candidates; "
                "falling back to composition only"
            )

    if len(cands) > 1:
        raise AmbiguousGlycanError(composition, cands, notes)

    return _finish(cands[0], tier=tier, composition=composition,
                   classify=classify, plausibility=plausibility,
                   strucgp_code=strucgp_code, notes=notes)


def _matches_topology(cand: GlycanCandidate, decoded: Dict[str, Any]) -> bool:
    """Does this library structure agree with a decoded StrucGP tree?"""
    want_antennae = decoded.get("antenna_count")
    if want_antennae is not None and cand.antenna_count != want_antennae:
        return False
    want_core_fuc = decoded.get("has_core_fuc")
    if want_core_fuc is not None and cand.has_core_fucose != bool(want_core_fuc):
        return False
    return True


def _finish(cand, *, tier, composition, classify, plausibility, strucgp_code, notes) -> GlycanChoice:
    gtype = ""
    flags: List[str] = []
    notes = list(notes)
    if classify is not None:
        result = classify(composition, strucgp_code) if strucgp_code else classify(composition)
        gtype = result[0] if isinstance(result, tuple) else result
    if plausibility is not None:
        flags = list(plausibility(composition, strucgp_code) if strucgp_code
                     else plausibility(composition))
        r11 = [f for f in flags if f.startswith("R11")]
        if r11:
            notes.append(
                "R11: 2 Fuc (292.12) and 1 NeuAc (291.10) differ by 1.02 Da -- "
                "confirm the composition before rendering it."
            )
    if tier == 2:
        notes.append(
            "tier 2 (composition-only): core vs antennary fucose and a2-3 vs "
            "a2-6 sialylation are chosen, not measured."
        )
    return GlycanChoice(candidate=cand, tier=tier, glycan_type=gtype,
                        plausibility_flags=flags, notes=notes)


# --------------------------------------------------------------------------
# Re-Glyco
# --------------------------------------------------------------------------

@dataclass
class SiteBuild:
    """What Re-Glyco actually did at one site."""

    residue: int
    glytoucan: str
    clash_free: Optional[bool] = None
    steric_score: Optional[float] = None
    phi: Optional[float] = None
    psi: Optional[float] = None


@dataclass
class GlycoproteinBuild:
    """A built glycoprotein model and the receipts for it."""

    accession: str
    path: Path
    source: str                                   # 'alphafold' | 'pdb' | uploaded
    sites: List[SiteBuild] = field(default_factory=list)
    uniprot_sites: List[dict] = field(default_factory=list)
    job_uuid: str = ""
    session_uuid: str = ""
    log: str = ""

    @property
    def all_clash_free(self) -> bool:
        return bool(self.sites) and all(s.clash_free for s in self.sites)


_SITE_LINE_RE = re.compile(
    r"site (?P<resi>\d+)_(?P<chain>\w+) Glycan (?P<gtc>\S+): steric score = "
    r"(?P<score>[\d.]+), phi = (?P<phi>[-\d.]+), psi = (?P<psi>[-\d.]+).*?"
    r"Clash Free = (?P<clash>True|False)"
)


def open_session(accession: str) -> dict:
    """Start a Re-Glyco session; the server fetches the scaffold itself.

    ``protID`` takes a UniProt accession (AlphaFold) or a PDB id.  Letting the
    server resolve it keeps us off the path this package exists to retire --
    hand-built AlphaFold URLs.
    """
    status, payload = http.post_json(SESSIONS_URL, {"protID": accession})
    if status not in (200, 201) or not isinstance(payload, dict):
        raise GlycoShapeError(f"Re-Glyco session for {accession!r} failed: HTTP {status}")
    if "session_uuid" not in payload:
        raise GlycoShapeError(f"Re-Glyco session for {accession!r} returned no session_uuid")
    return payload


def build_glycoprotein(
    accession: str,
    sites: Dict[int, str],
    *,
    chain: str = "A",
    population_size: int = 128,
    max_generations: int = 32,
    rotamer_scan: bool = True,
    cache_dir: Optional[Any] = None,
    timeout: float = 1800.0,
    poll_seconds: float = 15.0,
    session: Optional[dict] = None,
) -> GlycoproteinBuild:
    """Attach glycans at the named residues and return the built model.

    ``sites`` maps residue number to a **GlyTouCan accession** -- the value
    format the live API requires.  A GlycoShape ``GS…`` id or an IUPAC string is
    accepted by the endpoint's schema but crashes the worker, which is why
    :func:`resolve_glycan` hands back a candidate carrying its GlyTouCan.
    """
    sess = session or open_session(accession)
    sid = sess["session_uuid"]
    available = {int(s["residueID"]): s for s in sess.get("glycosylation", {}).get("available", [])}

    missing = [r for r in sites if r not in available]
    if missing:
        raise GlycoShapeError(
            f"{accession}: residues {missing} are not glycosylatable in this "
            f"structure (available Asn/Ser/Thr count = {len(available)})"
        )

    selected = {f"{r}_{available[r].get('residueChain', chain)}": gtc for r, gtc in sites.items()}
    request = {
        "session_uuid": sid,
        "jobType": "optimization",
        "selectedGlycans": selected,
        "parameters": {
            "populationSize": population_size,
            "maxGenerations": max_generations,
            "enableRotamerScan": rotamer_scan,
        },
        "output": {},
    }
    status, payload = http.post_json(JOBS_URL, request)
    if status not in (200, 201) or not isinstance(payload, dict) or "job_uuid" not in payload:
        raise GlycoShapeError(f"Re-Glyco job submission failed: HTTP {status}")
    jid = payload["job_uuid"]

    deadline = time.time() + timeout
    state: Dict[str, Any] = {}
    while time.time() < deadline:
        st, state = http.get_json(f"{JOBS_URL}/{jid}")
        if st != 200 or not isinstance(state, dict):
            raise GlycoShapeError(f"polling job {jid} returned HTTP {st}")
        if state.get("status") in ("completed", "failed", "error"):
            break
        time.sleep(poll_seconds)
    else:
        raise GlycoShapeError(f"Re-Glyco job {jid} did not finish within {timeout:.0f}s")

    results = state.get("results") or {}
    if state.get("status") != "completed":
        raise GlycoShapeError(
            f"Re-Glyco job {jid} {state.get('status')}: {results.get('error') or 'no detail'}"
        )

    st, blob = http.get_bytes(f"{JOBS_URL}/{jid}/files/output.pdb")
    if st != 200 or not blob:
        raise GlycoShapeError(f"job {jid} completed but output.pdb download returned HTTP {st}")

    out = cache.cache_root(cache_dir) / CACHE_SUBDIR / "models" / f"{accession}_{jid[:8]}.pdb"
    cache.write_atomic(out, blob)

    log = str(results.get("box_output") or "")
    built = [
        SiteBuild(
            residue=int(m.group("resi")),
            glytoucan=m.group("gtc"),
            clash_free=m.group("clash") == "True",
            steric_score=float(m.group("score")),
            phi=float(m.group("phi")),
            psi=float(m.group("psi")),
        )
        for m in _SITE_LINE_RE.finditer(log)
    ]
    return GlycoproteinBuild(
        accession=accession, path=out, source=sess.get("protein_info", {}).get("source", ""),
        sites=built, uniprot_sites=sess.get("glycosylation", {}).get("uniprot", []),
        job_uuid=jid, session_uuid=sid, log=log,
    )
