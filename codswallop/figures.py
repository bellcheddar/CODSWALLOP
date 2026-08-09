"""Figures for the dossier, drawn server-side and embedded.

The dossier fetches nothing. That is the whole point of it: a file that still opens from a
mailbox in ten years with no network, no CDN and no CODSWALLOP. So every picture in it is
either inline SVG generated here, or a raster embedded as a `data:` URI.

The app's own panels draw the same things in the browser from the same numbers. These are
deliberately NOT a port of that code: the browser versions are interactive and read CSS
custom properties, and a document that has neither needs its own, simpler drawing with the
colours written in.
"""

from __future__ import annotations

import base64
import logging
from typing import Optional

from . import db, http

logger = logging.getLogger(__name__)

# Paper palette, matching the dossier's own CSS. Written in rather than read from a variable
# because there is no stylesheet in the document to read one from.
INK = "#1D2430"
DIM = "#55627A"
MUTE = "#7C8AA3"
LINE = "#D3C7B0"
BRASS = "#8A6D40"
ACCENT = "#1F5FD8"
SURFACE = "#F1ECE0"

# How many ligand depictions are embedded. Each is a real fetch and perhaps 4 kB of base64 in
# the document, and a family can hold hundreds of components.
MAX_THUMBS = 12


def data_uri(url: str, mime: Optional[str] = None) -> Optional[str]:
    """Fetch a URL once and return it as a `data:` URI, cached.

    Cached because the dossier is regenerated on every request and these are the only network
    calls in it: without this, a family with twelve ligands would make twelve requests to the
    RCSB every time a link checker touched the page.
    """
    def fetch():
        raw = http.get_bytes(url)
        if not raw:
            return None
        kind = mime or ("image/svg+xml" if url.endswith(".svg") else "image/jpeg")
        return f"data:{kind};base64," + base64.b64encode(raw).decode("ascii")

    try:
        return db.cached(("data_uri", url), fetch)
    except Exception:                           # noqa: BLE001 - a figure is never fatal
        logger.warning("could not embed %s", url, exc_info=True)
        return None


def ligand_thumbs(fam: dict, limit: int = MAX_THUMBS) -> list:
    """The RCSB's own 2D depiction for each component somebody meant to be there.

    Ions and waters are skipped, as they are in the app: a picture of a sodium is a circle
    with Na in it and tells a reader nothing they did not already know from the label.
    """
    from . import ligands as ligand_engine
    out = []
    for c in ((fam.get("ligands") or {}).get("components") or []):
        if c.get("klass") not in ligand_engine.COUNTS_AS_BOUND:
            continue
        cid = (c.get("id") or "").upper()
        if not cid:
            continue
        uri = data_uri(f"https://cdn.rcsb.org/images/ccd/unlabeled/{cid[0]}/{cid}.svg")
        if uri:
            out.append({**c, "uri": uri})
        if len(out) >= limit:
            break
    return out


def structure_image(fam: dict) -> Optional[dict]:
    """A rendered picture of the reference structure.

    The app shows a live Mol* viewport here and a document cannot: Mol* is 5 MB and fetches
    the coordinates itself. The RCSB already renders every assembly, so the still is embedded
    instead. It is the same structure the family is superposed onto, so the picture and the
    superposition agree about what "the reference" means.
    """
    pdb_id = (fam.get("map") or {}).get("reference")
    if not pdb_id:
        return None
    uri = data_uri(f"https://cdn.rcsb.org/images/structures/{pdb_id.lower()}_assembly-1.jpeg")
    return {"pdb_id": pdb_id.upper(), "uri": uri} if uri else None


def _esc(text) -> str:
    return (str(text or "").replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def domain_svg(fam: dict, width: int = 860) -> str:
    """Domain architecture: every source's domains on one seed axis.

    One row per source (Pfam, CATH, SCOP, InterPro) rather than all of them merged, because
    they disagree about boundaries and merging them would invent a consensus none of them
    stated.
    """
    doms = (fam.get("domains") or {}).get("domains") or []
    length = fam.get("seed_length") or 0
    if not doms or not length:
        return ""

    by_source: dict = {}
    for d in doms:
        if d.get("start") and d.get("end"):
            by_source.setdefault(d.get("source") or "?", []).append(d)

    pad, row_h, label_w = 8, 22, 62
    span = width - label_w - pad * 2
    height = pad * 2 + row_h * len(by_source) + 16
    parts = [f'<svg viewBox="0 0 {width} {height}" width="100%" height="{height}" '
             'role="img" aria-label="Domain architecture on the seed sequence">']

    def x(pos):
        return label_w + pad + span * (pos - 1) / max(1, length - 1)

    for i, (source, items) in enumerate(sorted(by_source.items())):
        y = pad + i * row_h
        parts.append(f'<text x="{pad}" y="{y + 13}" font-size="10" fill="{DIM}" '
                     f'font-family="ui-sans-serif,system-ui">{_esc(source)}</text>')
        parts.append(f'<line x1="{x(1):.1f}" y1="{y + 9:.1f}" x2="{x(length):.1f}" '
                     f'y2="{y + 9:.1f}" stroke="{LINE}" stroke-width="2"/>')
        for d in items:
            x0, x1 = x(d["start"]), x(min(d["end"], length))
            parts.append(
                f'<rect x="{x0:.1f}" y="{y:.1f}" width="{max(2, x1 - x0):.1f}" height="18" '
                f'rx="3" fill="{SURFACE}" stroke="{BRASS}" stroke-width="1.2"/>')
            if x1 - x0 > 44:
                label = d.get("name") or d.get("id") or ""
                parts.append(
                    f'<text x="{(x0 + x1) / 2:.1f}" y="{y + 13:.1f}" font-size="9.5" '
                    f'fill="{INK}" text-anchor="middle" '
                    f'font-family="ui-sans-serif,system-ui">{_esc(label)[:26]}</text>')

    base = pad + row_h * len(by_source) + 4
    for tick in (1, length // 2, length):
        parts.append(f'<line x1="{x(tick):.1f}" y1="{base:.1f}" x2="{x(tick):.1f}" '
                     f'y2="{base + 4:.1f}" stroke="{MUTE}"/>')
        parts.append(f'<text x="{x(tick):.1f}" y="{base + 14:.1f}" font-size="9" '
                     f'fill="{MUTE}" text-anchor="middle" '
                     f'font-family="ui-monospace,Menlo,monospace">{tick}</text>')
    parts.append("</svg>")
    return "".join(parts)


def coverage_svg(fam: dict, width: int = 860, height: int = 92) -> str:
    """How many of the family's constructs contain each residue of the seed.

    A filled step profile, not a line: the value is a count per residue, and a line between
    two residues implies a value in the gap where there is no gap.
    """
    cov = (fam.get("stats") or {}).get("coverage") or {}
    depth = cov.get("depth") or []
    if not depth:
        return ""
    n, peak = len(depth), max(depth) or 1
    pad = 8
    span, plot = width - pad * 2, height - 26

    pts = [f"{pad},{pad + plot}"]
    for i, d in enumerate(depth):
        px = pad + span * i / max(1, n - 1)
        py = pad + plot - plot * (d / peak)
        pts.append(f"{px:.1f},{py:.1f}")
    pts.append(f"{pad + span},{pad + plot}")

    parts = [f'<svg viewBox="0 0 {width} {height}" width="100%" height="{height}" '
             'role="img" aria-label="Construct coverage along the seed sequence">',
             f'<polygon points="{" ".join(pts)}" fill="{SURFACE}" stroke="{BRASS}" '
             'stroke-width="1.2"/>',
             f'<line x1="{pad}" y1="{pad + plot}" x2="{pad + span}" y2="{pad + plot}" '
             f'stroke="{LINE}"/>']
    for tick in (1, n // 2, n):
        px = pad + span * (tick - 1) / max(1, n - 1)
        parts.append(f'<text x="{px:.1f}" y="{height - 10}" font-size="9" fill="{MUTE}" '
                     f'text-anchor="middle" font-family="ui-monospace,Menlo,monospace">'
                     f'{tick}</text>')
    parts.append(f'<text x="{pad}" y="{pad + 9}" font-size="9" fill="{MUTE}" '
                 f'font-family="ui-sans-serif,system-ui">{peak} constructs</text>')
    parts.append("</svg>")
    return "".join(parts)
