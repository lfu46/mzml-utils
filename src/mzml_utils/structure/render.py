"""Render a glycoprotein: an overview panel plus one zoom per glycosite.

**Label anchors are measured, never projected.**  Every render is done twice at
the same camera: a beauty pass, then a *detection pass* in which each glycosite
is painted a unique flat key colour with everything else hidden.  Anchors are
the centroids of those key colours in the finished image.

That is not a stylistic choice.  The predecessor of this module projected atom
coordinates through PyMOL's camera matrix and anchored labels on the Asn CA
while the glycan attaches at ND2 -- up to 78 px apart on a 2400 px frame, enough
that one CD14 leader line ended nearer a *different* site's marker than its own.
Measuring the marker in the rendered image makes anchor and marker the same
object, so that class of error cannot recur.  It also makes occlusion a
measurement: a site hidden behind the fold contributes no pixels.

The marker itself is the **glycosidic linkage** -- Asn ND2 plus C1 of the
reducing GlcNAc -- not the Cα, so the glycan reads as attached to its site
rather than floating near it.

PyMOL is a separate binary (not pip-installable), so it is driven by writing a
script and shelling out.  Import of this module never requires PyMOL; only
:func:`render_views` does.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

#: Near-black.  Chosen by measurement, not taste: against the five SNFG
#: monosaccharide colours plus a grey cartoon, an orange marker is a CIEDE2000
#: *blocker* under normal vision (dE 14.3 vs Fuc, threshold 15).  SNFG saturates
#: the hue circle, so lightness is the only axis left.
SITE_MARKER = "#1A1A1A"
PROTEIN_GREY = "#BFBFBF"

#: Key colours for the detection pass.  Widely separated in RGB so a centroid
#: cannot be contaminated by a neighbour; never seen by a reader.
_DET_KEYS = [
    (255, 0, 0), (0, 255, 0), (0, 0, 255), (255, 255, 0), (255, 0, 255),
    (0, 255, 255), (255, 128, 0), (128, 0, 255), (0, 128, 255), (128, 255, 0),
    (255, 0, 128), (0, 255, 128), (128, 128, 255), (255, 128, 128), (128, 255, 128),
]
_DET_TOL = 55


class RenderError(RuntimeError):
    """PyMOL is unavailable, or a render produced nothing usable."""


@dataclass
class LabelPlacement:
    """Where a label went, and how much structure is underneath it.

    ``ink`` is the measured fraction of non-white pixels in the grid cell the
    label landed on.  It is recorded rather than merely optimised so a panel
    that had nowhere clean to put a label says so.
    """

    residue: int
    xy: Tuple[float, float]
    anchor: Tuple[float, float]
    ink: float


@dataclass
class SiteView:
    """One glycosite as it appears in one rendered image.

    ``centroid`` comes from an *unoccluded* detection pass, so a label can be
    anchored on a site that is behind the fold.  ``pixels`` comes from an
    *occlusion-aware* pass, so ``visible`` reports what a reader can actually
    see.  Keeping them apart is what lets a hidden site be drawn with a dashed
    leader instead of being silently presented as facing the viewer.
    """

    residue: int
    centroid: Optional[Tuple[float, float]] = None   # position, occluded or not
    pixels: int = 0                                  # visible pixels only
    domain: str = ""

    @property
    def visible(self) -> bool:
        return self.pixels > 0

    @property
    def placeable(self) -> bool:
        return self.centroid is not None


@dataclass
class RenderedView:
    """A finished beauty render plus everything needed to label it."""

    name: str
    png: Path
    width: int
    height: int
    sites: Dict[int, SiteView] = field(default_factory=dict)
    focus: Optional[int] = None          # the site a zoom is centred on
    subtitle: str = ""

    def visible_sites(self) -> List[int]:
        return [r for r, s in self.sites.items() if s.visible]


def pymol_binary() -> str:
    """Locate PyMOL, or say so plainly."""
    exe = shutil.which("pymol") or "/opt/homebrew/bin/pymol"
    if not Path(exe).exists():
        raise RenderError(
            "PyMOL not found. It is not pip-installable; install the app or "
            "`brew install pymol`, then re-run."
        )
    return exe


def _rgb01(hex_color: str) -> List[float]:
    h = hex_color.lstrip("#")
    return [int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4)]


def map_glycan_chains(pdb: Path, chain: str = "A") -> Dict[int, List[str]]:
    """Which glycan chain hangs off which Asn, measured from the coordinates.

    Chain letters are assigned by the builder and do not follow residue order, so
    this is measured (nearest atom to the Asn ND2) rather than assumed.  A real
    N-glycosidic bond lands near 1.4 A; anything much longer is not a linkage.

    Needed because a glycan is a whole *chain*: selecting "the glycan at this
    site" by distance from the linkage catches only the reducing GlcNAc, which
    is how an earlier version greyed out the very glycan it was meant to show.
    """
    import numpy as np

    nd2: Dict[int, Any] = {}
    glyc: Dict[str, List[Any]] = {}
    for line in Path(pdb).read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        ch = line[21]
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        if line.startswith("ATOM") and ch == chain:
            if line[12:16].strip() == "ND2" and line[17:20].strip() == "ASN":
                nd2[int(line[22:26])] = np.array(xyz)
        else:
            glyc.setdefault(ch, []).append(xyz)

    out: Dict[int, List[str]] = {}
    for ch, coords in glyc.items():
        arr = np.asarray(coords)
        best, best_d = None, 1e9
        for resi, atom in nd2.items():
            d = float(np.linalg.norm(arr - atom, axis=1).min())
            if d < best_d:
                best, best_d = resi, d
        if best is not None and best_d < 2.0:
            out.setdefault(best, []).append(ch)
    return out


def _snfg_pymol_colors(snfg: Dict[str, Dict[str, str]]) -> Dict[str, List[float]]:
    """Map SNFG monosaccharide names onto the PDB residue codes PyMOL sees.

    ``snfg`` is the caller's SNFG table (``nglyco_figures.SNFG``) -- injected so
    the one authored definition of those colours stays the only one.
    """
    by_resn = {
        "NAG+NDG": snfg.get("GlcNAc"),
        "A2G+NGA": snfg.get("GalNAc") or snfg.get("GlcNAc"),
        "MAN+BMA": snfg.get("Man"),
        "GAL+GLA": snfg.get("Gal"),
        "FUC": snfg.get("Fuc"),
        "SIA": snfg.get("NeuAc"),
        "NGC": snfg.get("NeuGc") or snfg.get("NeuAc"),
        "XYS+XYP": snfg.get("Xyl"),
    }
    return {k: _rgb01(v["color"]) for k, v in by_resn.items() if v}


def _marker_selection(residue: int, chain: str = "A") -> str:
    """The glycosidic linkage: Asn ND2 and the reducing sugar's C1 next to it.

    Anchoring on the linkage rather than the Cα is what makes the label point at
    the place the glycan actually leaves the protein.
    """
    asn = f"(polymer.protein and chain {chain} and resi {residue} and name ND2+CG+OD1)"
    # PyMOL's `within` needs a bracketed left operand; without it the whole
    # selection is rejected as malformed.
    sugar = f"((not polymer.protein and name C1) within 2.5 of {asn})"
    return f"({asn} or {sugar})"


def _build_script(
    *,
    pdb: Path,
    sites: Sequence[int],
    chain: str,
    snfg: Dict[str, Dict[str, str]],
    views: Sequence[dict],
    out_dir: Path,
    px: int,
    domains: Optional[Dict[str, Tuple[str, List[Tuple[int, int]]]]] = None,
) -> str:
    """Emit a PyMOL script that renders every view twice (beauty + detection)."""
    lines = [
        "from pymol import cmd",
        "cmd.reinitialize()",
        f"cmd.load(r'{pdb}', 'gp')",
        "cmd.remove('solvent')",
        "cmd.bg_color('white')",
        "cmd.hide('everything', 'all')",
        f"cmd.set_color('gp_grey', {_rgb01(PROTEIN_GREY)})",
        f"cmd.set_color('gp_site', {_rgb01(SITE_MARKER)})",
    ]
    for i, rgb in enumerate(_DET_KEYS[:len(sites)]):
        lines.append(f"cmd.set_color('det{i}', {[c / 255.0 for c in rgb]})")

    lines += [
        "cmd.select('prot', 'polymer.protein')",
        "cmd.select('glyc', 'not polymer.protein')",
        "cmd.show('cartoon', 'prot')",
        "cmd.color('gp_grey', 'prot')",
        "cmd.set('cartoon_transparency', 0.25)",
        "cmd.show('sticks', 'glyc')",
        "cmd.set('stick_radius', 0.14)",
    ]
    for resn, rgb in _snfg_pymol_colors(snfg).items():
        lines.append(f"cmd.set_color('snfg_{resn.split('+')[0]}', {rgb})")
        lines.append(f"cmd.color('snfg_{resn.split('+')[0]}', 'glyc and resn {resn}')")

    if domains:
        for name, (color, spans) in domains.items():
            sel = " or ".join(f"resi {a}-{b}" for a, b in spans)
            safe = "dom_" + "".join(ch for ch in name if ch.isalnum())[:16]
            lines.append(f"cmd.set_color('{safe}', {_rgb01(color)})")
            lines.append(f"cmd.color('{safe}', 'prot and ({sel})')")

    lines += [
        "cmd.set('orthoscopic', 1)",
        "cmd.set('ray_shadows', 0)",
        "cmd.set('depth_cue', 0)",
        "cmd.set('specular', 0.2)",
        "cmd.set('ambient', 0.5)",
        "cmd.set('antialias', 2)",
        "cmd.set('ray_opaque_background', 1)",
        "cmd.set('sphere_quality', 3)",
    ]
    for i, r in enumerate(sites):
        lines.append(f"cmd.show('spheres', '{_marker_selection(r, chain)}')")
        lines.append(f"cmd.color('gp_site', '{_marker_selection(r, chain)}')")
    lines.append("cmd.set('sphere_scale', 0.42)")

    det_show = "\n".join(
        f"cmd.show('spheres', '{_marker_selection(r, chain)}')\n"
        f"cmd.color('det{i}', '{_marker_selection(r, chain)}')"
        for i, r in enumerate(sites)
    )

    for v in views:
        name, sel, buf = v["name"], v["selection"], v["buffer"]
        pre = v.get("pre", "")
        post = v.get("post", "")
        lines += [
            "# ---- %s" % name,
            pre,
            f"cmd.orient({sel!r})" if v.get("orient", True) else "",
            v.get("turn", ""),
            f"cmd.zoom({sel!r}, {buf})",
            "cmd.set('ray_trace_mode', 0)",
            f"cmd.png(r'{out_dir / (name + '.png')}', width={px}, height={px}, dpi=600, ray=1)",
            # --- detection pass 1: OCCLUSION-AWARE.  Same camera, flat key
            # colours on the markers, everything else forced to a neutral grey.
            # The grey matters: PyMOL's default object colour is green, which is
            # exactly one of the key colours, so an uncoloured cartoon would be
            # counted as a site marker.
            "cmd.hide('everything', 'all')",
            "cmd.set('ambient', 1.0)",
            "cmd.set('specular', 0.0)",
            "cmd.show('cartoon', 'prot')",
            "cmd.show('sticks', 'glyc')",
            "cmd.color('grey50', 'prot')",
            "cmd.color('grey50', 'glyc')",
            "cmd.set('cartoon_transparency', 0.0, 'prot')",
            det_show,
            "cmd.set('sphere_scale', 0.42)",
            f"cmd.png(r'{out_dir / (name + '_det.png')}', width={px}, height={px}, dpi=600, ray=1)",
            # --- detection pass 2: GEOMETRIC.  Markers only, nothing to hide
            # behind, so a site behind the fold still yields a position for its
            # label even though pass 1 correctly reported it as invisible.
            "cmd.hide('cartoon', 'all')",
            "cmd.hide('sticks', 'all')",
            f"cmd.png(r'{out_dir / (name + '_geo.png')}', width={px}, height={px}, dpi=600, ray=1)",
            # restore for the next view
            "cmd.set('ambient', 0.5)",
            "cmd.set('specular', 0.2)",
            "cmd.show('cartoon', 'prot')",
            "cmd.show('sticks', 'glyc')",
            "cmd.color('gp_grey', 'prot')",
            "cmd.set('cartoon_transparency', 0.25, 'prot')",
            post,
        ]
        for i, r in enumerate(sites):
            lines.append(f"cmd.show('spheres', '{_marker_selection(r, chain)}')")
            lines.append(f"cmd.color('gp_site', '{_marker_selection(r, chain)}')")

    lines.append("cmd.quit()")
    return "\n".join(x for x in lines if x)


def _detect(det_png: Path, sites: Sequence[int]) -> Dict[int, SiteView]:
    """Recover each site's pixel centroid from the flat-colour detection render."""
    import numpy as np
    from PIL import Image

    det = np.asarray(Image.open(det_png).convert("RGB")).astype(int)
    out: Dict[int, SiteView] = {}
    for i, resi in enumerate(sites):
        r, g, b = _DET_KEYS[i]
        mask = (
            (np.abs(det[..., 0] - r) < _DET_TOL)
            & (np.abs(det[..., 1] - g) < _DET_TOL)
            & (np.abs(det[..., 2] - b) < _DET_TOL)
        )
        ys, xs = np.where(mask)
        if len(xs) > 5:
            out[resi] = SiteView(resi, (float(xs.mean()), float(ys.mean())), int(len(xs)))
        else:
            out[resi] = SiteView(resi, None, int(len(xs)))
    return out


def choose_orientation(
    pdb: Path,
    sites: Sequence[int],
    *,
    chain: str = "A",
    px: int = 320,
    y_step: int = 30,
    x_angles: Sequence[int] = (-30, 0, 30),
    min_pixels: int = 8,
) -> Tuple[int, int, Dict[int, int]]:
    """Pick the view that leaves the fewest glycosites hidden.

    Scored on **pixels actually detected** in a cheap detection-only sweep, not
    on an atom-counting proxy: whether a marker is visible is a property of the
    finished image, so the image is what gets measured.

    Returns ``(turn_y, turn_x, pixels_per_site)`` for the winning view.  With
    sites distributed around a fold there may be no view that shows them all --
    the per-site counts are returned so the caller can say so rather than
    quietly present a view that hides one.
    """
    sites = list(sites)
    sweep_dir = Path(tempfile.mkdtemp(prefix="gp_sweep_"))
    try:
        candidates = [(ry, rx) for ry in range(0, 360, y_step) for rx in x_angles]
        det_show = "\n".join(
            f"cmd.show('spheres', '{_marker_selection(r, chain)}')\n"
            f"cmd.color('det{i}', '{_marker_selection(r, chain)}')"
            for i, r in enumerate(sites)
        )
        lines = [
            "from pymol import cmd",
            "cmd.reinitialize()",
            f"cmd.load(r'{pdb}', 'gp')",
            "cmd.remove('solvent')",
            "cmd.bg_color('white')",
            "cmd.hide('everything', 'all')",
        ]
        for i, rgb in enumerate(_DET_KEYS[:len(sites)]):
            lines.append(f"cmd.set_color('det{i}', {[c / 255.0 for c in rgb]})")
        lines += [
            "cmd.select('prot', 'polymer.protein')",
            "cmd.show('cartoon', 'prot')",   # shown so it can OCCLUDE the markers
            # PyMOL's default object colour is GREEN, which is exactly detection
            # key #1 -- an uncoloured cartoon gets counted as the second site's
            # marker.  Force a neutral grey, far from every key colour.
            "cmd.color('grey50', 'prot')",
            "cmd.color('grey50', 'not polymer.protein')",
            "cmd.show('sticks', 'not polymer.protein')",
            "cmd.set('orthoscopic', 1)",
            "cmd.set('ray_shadows', 0)",
            "cmd.set('ambient', 1.0)",
            "cmd.set('specular', 0.0)",
            "cmd.set('ray_trace_mode', 0)",
            det_show,
            "cmd.set('sphere_scale', 0.42)",
            "cmd.orient('all')",
            "base = cmd.get_view()",
        ]
        for k, (ry, rx) in enumerate(candidates):
            lines += [
                "cmd.set_view(base)",
                f"cmd.turn('y', {ry})",
                f"cmd.turn('x', {rx})",
                "cmd.zoom('all', 2.0)",
                f"cmd.png(r'{sweep_dir / f'sw{k}.png'}', width={px}, height={px}, dpi=72, ray=1)",
            ]
        lines.append("cmd.quit()")
        script = Path(sweep_dir) / "sweep.py"
        script.write_text("\n".join(lines))
        subprocess.run([pymol_binary(), "-cq", str(script)],
                       capture_output=True, text=True, timeout=1800)

        best, best_score = None, None
        for k, (ry, rx) in enumerate(candidates):
            png = sweep_dir / f"sw{k}.png"
            if not png.exists():
                continue
            found = {r: v.pixels for r, v in _detect(png, sites).items()}
            n_vis = sum(1 for v in found.values() if v >= min_pixels)
            score = (n_vis, min(found.values()) if found else 0)
            if best_score is None or score > best_score:
                best, best_score = (ry, rx, found), score
        if best is None:
            raise RenderError("orientation sweep produced no images")
        return best
    finally:
        shutil.rmtree(sweep_dir, ignore_errors=True)


def render_views(
    pdb: Path,
    sites: Sequence[int],
    *,
    out_dir: Path,
    snfg: Dict[str, Dict[str, str]],
    chain: str = "A",
    px: int = 2400,
    zooms: bool = True,
    domains: Optional[Dict[str, Tuple[str, List[Tuple[int, int]]]]] = None,
    keep_detection: bool = False,
    orient_search: bool = True,
) -> List[RenderedView]:
    """Render an overview plus one zoom per site, and measure every anchor.

    ``snfg`` is the caller's SNFG colour table; ``domains`` maps a display name
    to ``(hex_colour, [(start, end), ...])`` -- pass InterPro *fragments*, not
    envelopes, so a discontinuous domain does not colour the gap between its
    parts.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sites = list(sites)

    chain_map = map_glycan_chains(pdb, chain=chain)

    turn = ""
    if orient_search and len(sites) > 1:
        ry, rx, found = choose_orientation(pdb, sites, chain=chain)
        turn = f"cmd.turn('y', {ry})\ncmd.turn('x', {rx})"
    views: List[dict] = [
        {"name": "overview", "selection": "all", "buffer": 2.0, "turn": turn}
    ]
    if zooms:
        recolour = "\n".join(
            f"cmd.color('snfg_{resn.split('+')[0]}', 'glyc and resn {resn}')"
            for resn in _snfg_pymol_colors(snfg)
        )
        for r in sites:
            own = chain_map.get(r, [])
            if not own:
                raise RenderError(
                    f"no glycan chain is bonded to Asn{r} in {pdb.name}; a zoom "
                    "of an unglycosylated site would show nothing"
                )
            own_sel = "chain " + "+".join(own)
            focus = f"(({own_sel}) or (byres (polymer.protein within 12 of ({own_sel}))))"
            views.append({
                "name": f"site_N{r}",
                "selection": focus,
                "buffer": 2.5,
                # Fade every OTHER glycan: at this magnification a neighbouring
                # site's glycan is easily misread as part of this one.  Selected
                # by CHAIN -- a glycan is a whole chain, and selecting it by
                # distance from the linkage catches only the reducing GlcNAc.
                "pre": (
                    f"cmd.color('grey85', 'glyc and not ({own_sel})')\n"
                    f"cmd.set('stick_radius', 0.06, 'glyc and not ({own_sel})')\n"
                    f"cmd.set('stick_radius', 0.17, 'glyc and ({own_sel})')\n"
                    f"cmd.set('cartoon_transparency', 0.55, 'prot')\n"
                    f"cmd.set('cartoon_transparency', 0.15, 'prot and byres "
                    f"(prot within 10 of ({own_sel}))')"
                ),
                "post": (
                    "cmd.set('stick_radius', 0.14, 'glyc')\n"
                    "cmd.set('cartoon_transparency', 0.25, 'prot')\n"
                    + recolour
                ),
            })

    script = _build_script(pdb=pdb, sites=sites, chain=chain, snfg=snfg,
                           views=views, out_dir=out_dir, px=px, domains=domains)
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as fh:
        fh.write(script)
        script_path = Path(fh.name)
    try:
        proc = subprocess.run([pymol_binary(), "-cq", str(script_path)],
                              capture_output=True, text=True, timeout=1800)
    finally:
        script_path.unlink(missing_ok=True)

    rendered: List[RenderedView] = []
    for v in views:
        png = out_dir / f"{v['name']}.png"
        det = out_dir / f"{v['name']}_det.png"
        geo = out_dir / f"{v['name']}_geo.png"
        missing = [p.name for p in (png, det, geo) if not p.exists()]
        if missing:
            raise RenderError(
                f"PyMOL did not produce {missing} for {v['name']}; "
                f"stdout:\n{proc.stdout[-1500:]}\nstderr:\n{proc.stderr[-1500:]}"
            )
        from PIL import Image
        w, h = Image.open(png).size
        seen = _detect(det, sites)      # occlusion-aware -> visibility
        pos = _detect(geo, sites)       # unoccluded      -> position
        merged = {
            r: SiteView(residue=r,
                        centroid=pos[r].centroid,
                        pixels=seen[r].pixels)
            for r in sites
        }
        rv = RenderedView(name=v["name"], png=png, width=w, height=h, sites=merged)
        if v["name"].startswith("site_N"):
            rv.focus = int(v["name"].split("N")[1])
        rendered.append(rv)
        if not keep_detection:
            det.unlink(missing_ok=True)
            geo.unlink(missing_ok=True)
    return rendered


# --------------------------------------------------------------------------
# Label placement
# --------------------------------------------------------------------------

def occupancy_grid(png: Path, cells: int = 24, thresh: int = 250):
    """Coarse map of where the render has ink, for placing labels off it."""
    import numpy as np
    from PIL import Image

    grey = np.asarray(Image.open(png).convert("L"))
    h, w = grey.shape
    gy, gx = h // cells, w // cells
    grid = np.zeros((cells, cells), dtype=float)
    for i in range(cells):
        for j in range(cells):
            block = grey[i * gy:(i + 1) * gy, j * gx:(j + 1) * gx]
            grid[i, j] = float((block < thresh).mean()) if block.size else 1.0
    return grid


def place_labels(
    anchors: Dict[int, Tuple[float, float]],
    png: Path,
    *,
    cells: int = 24,
    max_ink: float = 0.02,
    min_sep_frac: float = 0.16,
    reserved: Sequence[Tuple[float, float, float, float]] = (),
) -> Dict[int, "LabelPlacement"]:
    """Put each label on background, near its anchor, and clear of other labels.

    Replaces fixed-radius placement, which put a label directly on top of the
    glycan it names.  Ink coverage in these renders is only 8-21% of frame, so a
    free cell almost always exists -- it simply has to be looked for.
    """
    import numpy as np

    grid = occupancy_grid(png, cells=cells)
    from PIL import Image
    w, h = Image.open(png).size
    cw, ch = w / cells, h / cells

    def _reserved(px_, py_) -> bool:
        """Areas the panel has already committed -- the subtitle, the cartoon."""
        return any(rx * w <= px_ <= (rx + rw) * w and ry * h <= py_ <= (ry + rh) * h
                   for rx, ry, rw, rh in reserved)

    free = [(i, j) for i in range(cells) for j in range(cells)
            if grid[i, j] <= max_ink
            and not _reserved((j + 0.5) * (w / cells), (i + 0.5) * (h / cells))]
    if not free:
        free = [(i, j) for i in range(cells) for j in range(cells)]

    placed: Dict[int, LabelPlacement] = {}
    taken: List[Tuple[float, float]] = []
    min_sep = min_sep_frac * max(w, h)
    cx, cy = w / 2.0, h / 2.0

    # Outermost anchors first: they have the fewest good options.
    order = sorted(anchors, key=lambda r: -((anchors[r][0] - cx) ** 2 + (anchors[r][1] - cy) ** 2))
    for resi in order:
        ax, ay = anchors[resi]
        # Prefer cells outward from the centre, close to the anchor, and away
        # from labels already placed.
        out_x, out_y = ax - cx, ay - cy
        norm = (out_x ** 2 + out_y ** 2) ** 0.5 or 1.0
        out_x, out_y = out_x / norm, out_y / norm

        best, best_score = None, None
        for i, j in free:
            px_, py_ = (j + 0.5) * cw, (i + 0.5) * ch
            dx, dy = px_ - ax, py_ - ay
            dist = (dx ** 2 + dy ** 2) ** 0.5
            if dist < 0.06 * max(w, h):
                continue
            outward = (dx * out_x + dy * out_y) / (dist or 1.0)
            clash = min((((px_ - tx) ** 2 + (py_ - ty) ** 2) ** 0.5 for tx, ty in taken),
                        default=1e9)
            if clash < min_sep:
                continue
            edge = min(px_, w - px_, py_, h - py_)
            if edge < 0.08 * max(w, h):
                continue
            score = dist / max(w, h) - 0.8 * outward
            if best_score is None or score < best_score:
                best, best_score = (px_, py_), score
        if best is None:
            # No free cell survived; place outward and RECORD the ink there, so
            # a crowded panel reports itself instead of looking clean.
            best = (float(np.clip(ax + out_x * 0.2 * w, 0.12 * w, 0.88 * w)),
                    float(np.clip(ay + out_y * 0.2 * h, 0.08 * h, 0.92 * h)))
        gi = min(cells - 1, max(0, int(best[1] / ch)))
        gj = min(cells - 1, max(0, int(best[0] / cw)))
        placed[resi] = LabelPlacement(residue=resi, xy=best,
                                      anchor=(ax, ay), ink=float(grid[gi, gj]))
        taken.append(best)
    return placed


# --------------------------------------------------------------------------
# Composition into a standard-compliant panel
# --------------------------------------------------------------------------

def _content_bbox(rgb, thresh: int = 250):
    import numpy as np
    grey = rgb.mean(axis=2)
    ys, xs = np.where(grey < thresh)
    if len(xs) == 0:
        return None
    return int(xs.min()), int(ys.min()), int(xs.max()), int(ys.max())


def compose_panel(
    view: RenderedView,
    *,
    out_stem,
    standard: Dict[str, Any],
    save,
    labels: Optional[Dict[int, str]] = None,
    subtitle: str = "",
    palette: Optional[Dict[str, str]] = None,
    inset=None,
    inset_rect: Tuple[float, float, float, float] = (0.60, 0.02, 0.38, 0.30),
):
    """Turn a rendered view into one gate-compliant panel file.

    ``standard`` is the lab figure standard **as data** and ``save`` is the
    project's saver.  Neither is imported here: the standard lives in
    wu-lab-skills and its reader lives in the consuming project, and mzml-utils
    sits below both.  Injecting them keeps exactly one copy of each.

    Everything is drawn into a **single axes** -- an inset axes would make
    ``n_axes`` 2 and trip the gate's ``combined-panels`` check.  ``inset`` is a
    callable ``(ax, (x0, y0, w, h)) -> None`` in axes fractions, used for the
    SNFG cartoon that gives monosaccharide identity a shape channel as well as a
    colour one.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from PIL import Image

    img = Image.open(view.png).convert("RGB")
    arr = np.asarray(img)
    anchors = {r: s.centroid for r, s in view.sites.items() if s.placeable}

    bbox = _content_bbox(arr)
    if bbox is None:
        raise RenderError(f"{view.png} is blank")
    pad = int(0.02 * img.size[0])
    x0 = max(0, bbox[0] - pad); y0 = max(0, bbox[1] - pad)
    x1 = min(img.size[0], bbox[2] + pad); y1 = min(img.size[1], bbox[3] + pad)
    img = img.crop((x0, y0, x1, y1))
    anchors = {r: (c[0] - x0, c[1] - y0) for r, c in anchors.items()}

    # Pad to square so the panel is EXACTLY the standard's size in both
    # directions.  Cropping alone leaves a non-square panel, which falls outside
    # the standard's size tolerance -- and a waiver is not ours to grant.
    cw, ch = img.size
    side = max(cw, ch)
    px, py = (side - cw) // 2, (side - ch) // 2
    square = Image.new("RGB", (side, side), (255, 255, 255))
    square.paste(img, (px, py))
    img = square
    anchors = {r: (c[0] + px, c[1] + py) for r, c in anchors.items()}

    # Keep labels off the subtitle strip and off the SNFG cartoon: both are
    # already-committed real estate, and the occupancy grid cannot see them
    # because neither is drawn yet.
    reserved: List[Tuple[float, float, float, float]] = []
    if subtitle:
        reserved.append((0.0, 0.0, 1.0, 0.10))
    if inset is not None:
        # inset_rect is axes fractions with y from the BOTTOM; `reserved` is
        # pixel fractions with y from the TOP.
        ix, iy, iw, ih = inset_rect
        reserved.append((ix, 1.0 - iy - ih, iw, ih))

    tmp = Path(str(out_stem) + "._compose.png")
    img.save(tmp)
    try:
        placed = place_labels(anchors, tmp, reserved=reserved) if anchors else {}
    finally:
        tmp.unlink(missing_ok=True)

    w_in = float(standard["width_in"])
    h_in = float(standard["height_in"])
    fs_max = float(standard["font_pt_max"])
    fs_min = float(standard["font_pt_min"])

    fig = plt.figure(figsize=(w_in, h_in), dpi=float(standard["dpi_save"]))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(img, extent=(0, side, side, 0), interpolation="antialiased")
    ax.set_xlim(0, side)
    ax.set_ylim(side, 0)
    ax.axis("off")

    for resi, lp in placed.items():
        lx, ly = lp.xy
        ax_, ay_ = anchors[resi]
        facing = view.sites[resi].visible
        ax.plot([ax_, lx], [ay_, ly], lw=0.6, color="#333333", zorder=3,
                solid_capstyle="round",
                linestyle="-" if facing else (0, (2, 1.6)))
        ax.text(lx, ly, (labels or {}).get(resi, f"N{resi}"),
                fontsize=fs_max,
                fontweight="bold" if facing else "normal",
                style="normal" if facing else "italic",
                color="#111111" if facing else "#555555",
                ha="center", va="center", zorder=4,
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.85))

    if subtitle:
        ax.text(0.5, 0.985, subtitle, transform=ax.transAxes, fontsize=fs_min,
                ha="center", va="top", color="#222222", zorder=4,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))

    if inset is not None:
        inset(ax, inset_rect)

    result = save(fig, str(out_stem), palette=palette) if palette is not None \
        else save(fig, str(out_stem))
    plt.close(fig)
    return {"saved": result, "labels": placed, "side": side}


# --------------------------------------------------------------------------
# The one call
# --------------------------------------------------------------------------

@dataclass
class GlycoproteinFigure:
    """Everything a run produced, including what it could not show."""

    accession: str
    model: Path
    panels: Dict[str, Any] = field(default_factory=dict)
    views: List[RenderedView] = field(default_factory=list)
    choices: Dict[int, Any] = field(default_factory=dict)
    domains: Dict[str, Any] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)

    def label_ink(self) -> Dict[str, Dict[int, float]]:
        """Structure ink measured under every label, per panel."""
        return {name: {r: lp.ink for r, lp in p["labels"].items()}
                for name, p in self.panels.items() if isinstance(p, dict)}

    def report(self) -> str:
        lines = [f"{self.accession}: {len(self.panels)} panels from {self.model.name}"]
        for v in self.views:
            hidden = [f"N{r}" for r, s in sorted(v.sites.items())
                      if s.placeable and not s.visible]
            off = [f"N{r}" for r, s in sorted(v.sites.items()) if not s.placeable]
            bits = []
            if hidden:
                bits.append("behind fold: " + ",".join(hidden))
            if off:
                bits.append("out of frame: " + ",".join(off))
            lines.append(f"  {v.name}: " + ("; ".join(bits) if bits else "all sites visible"))
        worst = [(n, r, v) for n, d in self.label_ink().items()
                 for r, v in d.items() if v > 0.02]
        for n, r, v in sorted(worst, key=lambda t: -t[2]):
            lines.append(f"  ! {n}: label N{r} sits on {100 * v:.0f}% structure ink")
        lines.extend("  ! " + w for w in self.warnings)
        return "\n".join(lines)


def render_glycoprotein(
    accession: str,
    sites: Dict[int, Any],
    *,
    out_dir,
    standard: Dict[str, Any],
    save,
    snfg: Dict[str, Dict[str, str]],
    domains: Any = False,
    domain_palette: Optional[Sequence[str]] = None,
    px: int = 2400,
    zooms: bool = True,
    palette: Optional[Dict[str, str]] = None,
    inset_factory=None,
    inset_rect: Tuple[float, float, float, float] = (0.60, 0.02, 0.38, 0.30),
    model: Optional[Path] = None,
    classify=None,
    plausibility=None,
) -> GlycoproteinFigure:
    """Overview panel plus one zoom per glycosite, ready for the figure gate.

    ``sites`` maps residue number to a GlyTouCan accession, a GlycoShape ``GS``
    id, or an already-resolved ``GlycanChoice``.  It does **not** accept a bare
    composition: use :func:`~.glycoprotein.resolve_glycan` first, so that an
    ambiguous composition raises where the choice is being made rather than
    silently becoming a picture.
    """
    from . import glycoprotein as gp
    from . import interpro

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    warnings: List[str] = []

    choices: Dict[int, gp.GlycanChoice] = {}
    for resi, value in sites.items():
        if isinstance(value, gp.GlycanChoice):
            choices[resi] = value
        elif isinstance(value, str):
            choices[resi] = gp.resolve_glycan(
                glytoucan=value if value.startswith("G") else None,
                glycoshape_id=None if value.startswith("G") else value,
                classify=classify, plausibility=plausibility,
            )
        else:
            raise TypeError(
                f"site {resi}: expected a GlyTouCan id, a GlycoShape id or a "
                f"GlycanChoice, got {type(value).__name__}. A bare composition "
                "must go through resolve_glycan() so ambiguity is raised there."
            )
        for note in choices[resi].notes:
            warnings.append(f"N{resi}: {note}")

    if model is None:
        build = gp.build_glycoprotein(
            accession, {r: c.glytoucan for r, c in choices.items()}
        )
        model = build.path
        for s in build.sites:
            if not s.clash_free:
                warnings.append(
                    f"N{s.residue}: Re-Glyco could not place {s.glytoucan} without a clash"
                )
    model = Path(model)

    dom_spec: Optional[Dict[str, Tuple[str, List[Tuple[int, int]]]]] = None
    dom_rows: Dict[str, Any] = {}
    if domains:
        if isinstance(domains, dict):
            dom_spec = domains
        else:
            rows = interpro.domains(accession, types=["domain"])
            if not rows:
                warnings.append(f"InterPro reports no domains for {accession}")
            if rows and not domain_palette:
                raise ValueError(
                    "domains=True needs domain_palette=[...] -- colours are the "
                    "project's to choose, not this module's to invent"
                )
            dom_spec = {}
            for i, row in enumerate(rows):
                # Fragments, not the envelope: a discontinuous domain must not
                # colour the gap between its parts.
                spans = [(int(f["start"]), int(f["end"]))
                         for f in (row.get("fragments") or [])] or \
                        [(int(row["start"]), int(row["end"]))]
                name = row.get("name") or row.get("accession") or f"domain{i}"
                dom_spec[name] = (domain_palette[i % len(domain_palette)], spans)
                dom_rows[name] = row

    views = render_views(model, sorted(choices), out_dir=out_dir / "_raw",
                         snfg=snfg, px=px, zooms=zooms, domains=dom_spec)

    panels: Dict[str, Any] = {}
    for v in views:
        sub = ""
        if v.focus is not None:
            ch = choices[v.focus]
            sub = gp._comp_repr(ch.candidate.composition)
            if ch.glycan_type:
                sub += f"  ({ch.glycan_type.replace('_', ' ')})"
        ins = None
        if inset_factory is not None and v.focus is not None:
            ins = inset_factory(choices[v.focus])
        stem = out_dir / f"{accession}_{v.name}"
        panels[v.name] = compose_panel(
            v, out_stem=stem, standard=standard, save=save,
            subtitle=sub, palette=palette, inset=ins, inset_rect=inset_rect,
        )
        print(f"wrote {stem}.png")

    return GlycoproteinFigure(accession=accession, model=model, panels=panels,
                              views=views, choices=choices, domains=dom_rows,
                              warnings=warnings)
