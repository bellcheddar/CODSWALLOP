"""Functional sites on the seed, grounded in what the family actually resolved.

Two sources, deliberately unequal in standing and labelled as such on the panel:

* **UniProt curated features.** Somebody read a paper and recorded that this residue is
  phosphorylated on this protein. Evidence about the real molecule.
* **PROSITE signatures.** A pattern or profile matched this sequence. Evidence about a
  string, which may or may not be evidence about the protein.

The distinction is not pedantry. Scanning chicken lysozyme with the defaults returns a PKC
phosphorylation site and two N-myristoylation sites, on a secreted protein that has neither:
they are PROSITE's own high-probability patterns, which match by chance in almost any
sequence. `skip=1` drops them and leaves the lysozyme-like domain profile and the glycosyl
hydrolase family 22 signature, which are the two things actually worth reporting. Six hits
become two, and the two are right.

What makes this a family panel rather than a sequence-annotation panel is the last step: a
site is reported with how many of the family's constructs contain it and how often it is
actually resolved. A predicted site nobody has ever built into a construct is a different
thing from one that fifty structures resolve, and the difference is invisible on any single
entry's page.
"""

from __future__ import annotations

import logging
import urllib.parse
from typing import Optional

from . import config, db, http

logger = logging.getLogger(__name__)

# Part of every cache key here. See rcsb.PARSE_VERSION for what forgetting it costs.
MOTIF_VERSION = 1

SCANPROSITE_URL = "https://prosite.expasy.org/cgi-bin/prosite/scanprosite/PSScan.cgi"
INTERPRO_PROSITE_URL = "https://www.ebi.ac.uk/interpro/api/entry/prosite/{ac}"

# How the UniProt feature keys are presented, and in what order. Curated observations lead,
# because they are the ones a reader should weight.
FEATURE_LABELS = [
    ("active_site", "Active site", "curated"),
    ("binding_site", "Binding site", "curated"),
    ("modified", "Modified residue", "curated"),
    ("glycosylation", "Glycosylation", "curated"),
    ("lipidation", "Lipidation", "curated"),
    ("disulphide", "Disulphide bond", "curated"),
    ("transmembrane", "Transmembrane", "curated"),
    ("signal_peptide", "Signal peptide", "curated"),
    ("motif", "Short linear motif", "curated"),
    ("site", "Site of interest", "curated"),
]

# Single-residue feature kinds: `uniprot._features_from` expands these to bare positions
# rather than to spans, so they need reading differently.
POSITION_KINDS = {"active_site", "binding_site"}


def scan_prosite(sequence: str, include_frequent: bool = False) -> list[dict]:
    """PROSITE signatures on one sequence, without the ones that match everything.

    `include_frequent` maps to ScanProsite's `skip`, inverted. It defaults off, and that
    default is the whole value of this function: see the module docstring.
    """
    if not sequence or len(sequence) < 12:
        return []

    def fetch():
        body = http.post_form(
            SCANPROSITE_URL,
            {"seq": sequence, "output": "json", "skip": "0" if include_frequent else "1"},
        )
        if not body:
            return []
        out = []
        for m in body.get("matchset") or []:
            ac, start, stop = m.get("signature_ac"), m.get("start"), m.get("stop")
            if not ac or start is None or stop is None:
                continue
            out.append({
                "accession": ac,
                "start": int(start),
                "end": int(stop),
                "score": m.get("score"),
                # A PROSITE accession beginning PS5 is a profile or a rule, which is scored
                # and far more specific than a PS00-series regular expression. Worth saying,
                # because a reader deciding whether to believe a hit wants to know which.
                "kind": "profile" if ac.upper().startswith("PS5") else "pattern",
            })
        return out

    return db.cached(("prosite", MOTIF_VERSION, sequence, include_frequent), fetch) or []


def prosite_name(accession: str) -> Optional[dict]:
    """The human-readable name for a PROSITE accession, via InterPro.

    ScanProsite returns `signature_id: null` for every hit, so the accessions arrive
    anonymous and a panel listing PS00128 and PS51348 tells the reader nothing at all.
    """
    def fetch():
        # InterPro does not carry every PROSITE signature: it has the PS00-series patterns
        # but answered for the PS51348 profile with an HTML page rather than a 404, so
        # `get_json` raised and one unnamed hit took the whole panel down through the
        # optional-panel handler. An accession with no name is not an error, it is a hit
        # that gets displayed by its accession.
        try:
            rec = http.get_json(
                INTERPRO_PROSITE_URL.format(ac=urllib.parse.quote(accession)) + "/")
        except Exception:                       # noqa: BLE001
            logger.info("no InterPro record for %s", accession)
            return None
        meta = ((rec or {}).get("metadata") or {})
        if not meta:
            return None
        name = meta.get("name") or {}
        return {"name": (name.get("name") if isinstance(name, dict) else name)
                        or meta.get("accession"),
                "short": name.get("short") if isinstance(name, dict) else None,
                "type": meta.get("type")}

    return db.cached(("prosite_name", MOTIF_VERSION, accession.upper()), fetch)


def _ground(start: int, end: int, depth: list, seen: list, max_depth: int) -> dict:
    """How much of the family actually contains this span, and how much resolves it.

    This is the step that turns a sequence annotation into a family observation. Averaged
    over the span rather than taken at its first residue, because a transmembrane helix that
    is half-resolved is a different statement from one that is not resolved at all.
    """
    lo, hi = max(1, start), min(end, len(depth))
    if hi < lo or not max_depth:
        return {"in_constructs": None, "resolved": None, "residues": 0}
    span = range(lo - 1, hi)
    d = [depth[i] for i in span]
    s = [seen[i] for i in span]
    return {
        "residues": hi - lo + 1,
        "in_constructs": round(100.0 * (sum(d) / len(d)) / max_depth, 1),
        # `seen` counts chains that resolved the residue; expressed against the constructs
        # that contained it, so "resolved" means "resolved when present" and a site nobody
        # cloned does not read as disordered.
        "resolved": round(100.0 * sum(s) / sum(d), 1) if sum(d) else None,
    }


def build(fam: dict, features: Optional[dict] = None) -> dict:
    """Assemble the motifs panel for one family."""
    seq = fam.get("seed_sequence") or ""
    cov = (fam.get("stats") or {}).get("coverage") or {}
    depth = cov.get("depth") or []
    seen = cov.get("seen") or []
    max_depth = cov.get("max_depth") or 0

    def ground(a, b):
        return _ground(a, b, depth, seen, max_depth) if depth else {
            "in_constructs": None, "resolved": None, "residues": max(0, b - a + 1)}

    rows: list[dict] = []

    # ---- curated, first ---------------------------------------------------------------
    feats = features or {}
    for key, label, standing in FEATURE_LABELS:
        vals = feats.get(key) or []
        if key in POSITION_KINDS:
            for pos in sorted(set(vals)):
                rows.append({"source": "UniProt", "standing": standing, "kind": label,
                             "name": f"{label} {pos}", "start": pos, "end": pos,
                             "description": "", **ground(pos, pos)})
        else:
            for f in vals:
                a, b = f["start"], f["end"]
                # A disulphide is two cysteines, so it is grounded on the pair rather than
                # on the sequence lying between them, which is not part of the bond.
                g = (ground(a, a) if key == "disulphide" else ground(a, b))
                if key == "disulphide":
                    gb = ground(b, b)
                    for fld in ("in_constructs", "resolved"):
                        if g.get(fld) is not None and gb.get(fld) is not None:
                            g[fld] = round(min(g[fld], gb[fld]), 1)
                    g["residues"] = 2
                rows.append({"source": "UniProt", "standing": standing, "kind": label,
                             "name": (f.get("description")
                                      or (f"{label} {a}\u2013{b}" if key == "disulphide"
                                          else label)),
                             "start": a, "end": b,
                             "description": f.get("description") or "", **g})

    # ---- predicted, second and labelled ----------------------------------------------
    for hit in scan_prosite(seq):
        meta = prosite_name(hit["accession"]) or {}
        rows.append({
            "source": "PROSITE", "standing": "predicted",
            "kind": "Profile" if hit["kind"] == "profile" else "Pattern",
            "name": meta.get("name") or hit["accession"],
            "accession": hit["accession"], "score": hit.get("score"),
            "start": hit["start"], "end": hit["end"], "description": meta.get("short") or "",
            **ground(hit["start"], hit["end"]),
        })

    rows.sort(key=lambda r: (r["standing"] != "curated", r["start"], r["end"]))
    return {
        "n": len(rows),
        "rows": rows[:400],
        "n_curated": sum(1 for r in rows if r["standing"] == "curated"),
        "n_predicted": sum(1 for r in rows if r["standing"] == "predicted"),
        "seed_length": len(seq),
        "grounded": bool(depth),
    }
