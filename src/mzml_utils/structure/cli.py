"""Command line for mzml_utils.structure.

    python3 -m mzml_utils.structure resolve P08571
    python3 -m mzml_utils.structure fetch P08571 Q9NR16 --out ./structures
    python3 -m mzml_utils.structure sites P08571 --linkage N
    python3 -m mzml_utils.structure domains P08571

``fetch`` prints one row per accession and exits non-zero only when something was
left **unanswered** -- an accession with genuinely no model is a result, not a
failure, and must not fail a pipeline step.
"""
from __future__ import annotations

import argparse
import json
import sys
from typing import Dict, List

from .client import fetch_many, resolve
from .interpro import domains as interpro_domains
from .records import FetchStatus
from .uniprot import accession_state, glycosylation_sites


def _read_accessions(values: List[str]) -> List[str]:
    """Accessions from the command line, or one per line from a file."""
    out: List[str] = []
    for v in values:
        if v.startswith("@"):
            with open(v[1:]) as fh:
                out.extend(line.strip() for line in fh if line.strip())
        else:
            out.append(v)
    return [a for a in dict.fromkeys(out) if a]


def _cmd_resolve(args) -> int:
    for acc in _read_accessions(args.accessions):
        res = resolve(acc, prefer=args.prefer)
        if args.json:
            print(json.dumps(
                {
                    "requested": res.requested,
                    "accession": res.accession,
                    "status": res.status.value,
                    "is_answered": res.status.is_answered,
                    "detail": res.detail,
                    "records": [
                        {
                            "provider": r.provider,
                            "model_identifier": r.model_identifier,
                            "model_category": r.model_category,
                            "model_version": r.model_version,
                            "coverage": r.coverage,
                            "resolution": r.resolution,
                            "confidence_avg_local_score": r.confidence_avg_local_score,
                            "model_url": r.model_url,
                        }
                        for r in res.records
                    ],
                }
            ))
            continue
        header = f"{res.requested}"
        if res.redirected:
            header += f" -> {res.accession}"
        print(f"{header}  [{res.status.value}]  {len(res)} model(s)")
        if res.detail:
            print(f"    {res.detail}")
        for r in res.records[: args.limit]:
            bits = [f"v{r.model_version}"] if r.model_version else []
            if r.confidence_avg_local_score is not None:
                bits.append(f"score={r.confidence_avg_local_score}")
            if r.coverage is not None:
                bits.append(f"cov={r.coverage}")
            if r.resolution is not None:
                bits.append(f"res={r.resolution}")
            print(f"    {r.provider:14s} {r.model_identifier:28s} {' '.join(bits)}")
            if args.urls:
                print(f"        {r.model_url}")
    return 0


def _cmd_fetch(args) -> int:
    results = fetch_many(
        _read_accessions(args.accessions),
        cache_dir=args.out,
        prefer=args.prefer,
        force=args.force,
        max_workers=args.workers,
    )
    unanswered = 0
    for res in results:
        row = res.as_row()
        if args.json:
            print(json.dumps(row))
        else:
            print(
                f"{row['accession']:12s} {row['status']:20s} "
                f"{row['path'] or row['detail']}"
            )
        if not res.status.is_answered:
            unanswered += 1

    have = sum(1 for r in results if r.status.has_structure)
    absent = sum(1 for r in results if r.status is FetchStatus.NO_MODEL_AVAILABLE)
    print(
        f"\n{have}/{len(results)} on disk | {absent} with no model | "
        f"{unanswered} unanswered",
        file=sys.stderr,
    )
    # Only unanswered accessions are a failure: "no model exists" is an answer.
    return 1 if unanswered else 0


def _cmd_sites(args) -> int:
    for acc in _read_accessions(args.accessions):
        for site in glycosylation_sites(acc, linkage=args.linkage):
            flag = "" if site["certain"] else "  (position uncertain)"
            print(f"{acc}\t{site['start']}\t{site['description']}{flag}")
    return 0


def _cmd_domains(args) -> int:
    for acc in _read_accessions(args.accessions):
        for d in interpro_domains(acc, db=args.db):
            span = f"{d['start']}-{d['end']}"
            if d["discontinuous"]:
                span += " (" + ",".join(
                    f"{f['start']}-{f['end']}" for f in d["fragments"]
                ) + ")"
            print(f"{acc}\t{d['entry_accession']}\t{d['type']}\t{span}\t{d['name']}")
    return 0


def _cmd_state(args) -> int:
    for acc in _read_accessions(args.accessions):
        print(json.dumps(accession_state(acc)))
    return 0


def _cmd_glycans(args) -> int:
    """Which 3D structures carry a composition -- and how many, which is the point."""
    from .glycoprotein import _comp_repr, candidates_for

    comp: Dict[str, int] = {}
    for part in args.composition.replace(",", " ").split():
        key, _, val = part.partition("=")
        comp[key] = int(val or 1)
    cands = candidates_for(comp, residue=args.residue)
    print(f"{_comp_repr(comp)} on {args.residue}: {len(cands)} structure(s)")
    for c in cands:
        print(f"  {c.glytoucan}\t{c.glycoshape_id}\t{c.iupac}")
    if len(cands) > 1:
        print(
            f"\nAmbiguous: composition alone cannot choose among {len(cands)}. "
            "Name one with its GlyTouCan id when rendering."
        )
    return 0 if len(cands) == 1 else 1


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python3 -m mzml_utils.structure",
        description="Protein structure retrieval (AlphaFold, PDBe, SWISS-MODEL, "
        "UniProt, InterPro) through one version-correct resolver.",
    )
    sub = p.add_subparsers(dest="command", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "accessions",
        nargs="+",
        help="UniProt accessions, or @file.txt to read one per line",
    )

    r = sub.add_parser("resolve", parents=[common], help="list available structures")
    r.add_argument("--prefer", default="predicted",
                   choices=["predicted", "experimental", "coverage"])
    r.add_argument("--limit", type=int, default=10, help="rows per accession")
    r.add_argument("--urls", action="store_true", help="print model URLs")
    r.add_argument("--json", action="store_true")
    r.set_defaults(func=_cmd_resolve)

    f = sub.add_parser("fetch", parents=[common], help="download structures to a cache")
    f.add_argument("--out", default=None,
                   help="cache directory (default: WU_LAB_STRUCTURE_CACHE or ~/.cache/wu-lab/structures)")
    f.add_argument("--prefer", default="predicted",
                   choices=["predicted", "experimental", "coverage"])
    f.add_argument("--force", action="store_true", help="re-download cached files")
    f.add_argument("--workers", type=int, default=4)
    f.add_argument("--json", action="store_true")
    f.set_defaults(func=_cmd_fetch)

    s = sub.add_parser("sites", parents=[common], help="annotated glycosylation sites")
    s.add_argument("--linkage", default="N", choices=["N", "O", "any"])
    s.set_defaults(func=_cmd_sites)

    d = sub.add_parser("domains", parents=[common], help="InterPro domains")
    d.add_argument("--db", default="interpro", help="interpro | all | pfam | ...")
    d.set_defaults(func=_cmd_domains)

    st = sub.add_parser("state", parents=[common], help="is this accession still live?")
    st.set_defaults(func=_cmd_state)

    g = sub.add_parser(
        "glycans",
        help="3D structures matching a composition (exit 1 if ambiguous)",
    )
    g.add_argument("composition",
                   help="e.g. 'HexNAc=4 Hex=5 Fuc=1 NeuAc=2'")
    g.add_argument("--residue", default="ASN",
                   help="attachment residue: ASN, SER, THR, TRP, HYP, PRO")
    g.set_defaults(func=_cmd_glycans)
    return p


def _cli(argv: List[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_cli())
