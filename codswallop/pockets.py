"""Family-specific reference numbering: KLIFS for kinases, GPCRdb for receptors.

Both give a family a vocabulary its own residue numbers do not. "The gatekeeper" and "3.50"
mean the same thing in every member of the family, where "Thr315" and "Arg131" mean it in
exactly one, and a page that exists to compare 1,988 structures should be speaking the
former.

The two are used differently because the two databases are shaped differently:

* **GPCRdb** returns per-residue generic numbering keyed on the **UniProt sequence number**,
  so it lands directly in seed coordinates alongside conservation, the coverage census and
  the PLIP hot residues. Nothing has to be inferred.
* **KLIFS** numbers its 85-residue pocket against each deposited structure's own author
  numbering, which is not the seed's. Rather than guess a mapping, this reports the thing
  KLIFS is uniquely able to say about a *family*: which conformation each of its structures
  was caught in. DFG-in against DFG-out across 158 ABL1 structures is a statement no single
  entry page can make, and it needs no coordinate mapping to be true.
"""

from __future__ import annotations

import logging
import re
from collections import Counter
from typing import Optional

from . import db, http

logger = logging.getLogger(__name__)

# 2: GPCRdb names arrive as HTML and are now stripped. Bumped rather than deleted by hand,
# because cache keys are hashed, so a DELETE ... LIKE '%gpcrdb%' matches nothing at all and
# looks exactly like a cache that was already clear.
POCKET_VERSION = 2

KLIFS_KINASE_URL = "https://klifs.net/api/kinase_ID"
KLIFS_STRUCTURES_URL = "https://klifs.net/api/structures_list"
GPCRDB_PROTEIN_URL = "https://gpcrdb.org/services/protein/accession/{acc}/"
GPCRDB_RESIDUES_URL = "https://gpcrdb.org/services/residues/{entry}/"

# The micro-switches, by generic number rather than by residue number, which is the whole
# point of a generic numbering scheme. Ranges are inclusive.
GPCR_SWITCHES = [
    ("DRY", "3.49", "3.51", "The ionic lock at the cytoplasmic end of TM3: its arginine "
                            "holds the receptor inactive."),
    ("CWxP", "6.47", "6.50", "The rotamer toggle on TM6, directly under the orthosteric "
                             "pocket."),
    ("NPxxY", "7.49", "7.53", "TM7's activation switch, which repacks against TM6 when the "
                              "receptor turns on."),
    ("PIF", "3.40", "3.40", "The connector region between the pocket and the "
                            "cytoplasmic switches."),
]


def _num(generic: Optional[str]) -> Optional[float]:
    """A Ballesteros-Weinstein number as a float, for range tests.

    GPCRdb decorates the display form ("3.50x50" for the structure-based scheme, and a bare
    "3.50" otherwise), so the decoration is stripped rather than parsed.
    """
    if not generic:
        return None
    head = str(generic).split("x")[0].strip()
    try:
        return float(head)
    except ValueError:
        return None


# --------------------------------------------------------------------------------------
# Kinases
# --------------------------------------------------------------------------------------
def klifs_kinase(accession: str) -> Optional[dict]:
    """The KLIFS record for an accession, or None if it is not a kinase."""
    def fetch():
        try:
            rows = http.get_json(KLIFS_KINASE_URL,
                                 params={"kinase_name": accession, "species": "HUMAN"})
        except Exception:                       # noqa: BLE001
            return None
        # KLIFS matches on name, so an accession that is not a kinase returns an error body
        # or an empty list rather than a 404. Confirm by the accession it hands back.
        for r in rows or []:
            if (r.get("uniprot") or "").upper() == accession.upper():
                return {"kinase_id": r.get("kinase_ID"), "name": r.get("name"),
                        "family": r.get("family"), "group": r.get("group"),
                        "full_name": r.get("full_name"), "pocket": r.get("pocket")}
        return None

    return db.cached(("klifs_kinase", POCKET_VERSION, accession.upper()), fetch)


def klifs_conformations(kinase_id: int, pdb_ids: set) -> dict:
    """How the family's own structures were caught: DFG-in or out, alpha-C in or out.

    Restricted to the entries this family actually contains, so a family filtered to one
    organism does not report the conformations of a different species' structures.
    """
    def fetch():
        try:
            return http.get_json(KLIFS_STRUCTURES_URL, params={"kinase_ID": kinase_id}) or []
        except Exception:                       # noqa: BLE001
            return []

    rows = db.cached(("klifs_structures", POCKET_VERSION, kinase_id), fetch) or []
    mine = [r for r in rows if (r.get("pdb") or "").upper() in pdb_ids]
    if not mine:
        return {"n": 0, "dfg": [], "achelix": [], "n_klifs_total": len(rows)}

    def tally(field):
        c = Counter((r.get(field) or "unknown").lower() for r in mine)
        total = sum(c.values())
        return [{"state": k, "n": v, "pct": round(100.0 * v / total, 1)}
                for k, v in c.most_common()]

    return {"n": len(mine), "dfg": tally("DFG"), "achelix": tally("aC_helix"),
            "n_klifs_total": len(rows)}


# --------------------------------------------------------------------------------------
# GPCRs
# --------------------------------------------------------------------------------------
def gpcrdb_protein(accession: str) -> Optional[dict]:
    """The GPCRdb record for an accession, or None if it is not a receptor it knows."""
    def fetch():
        try:
            r = http.get_json(GPCRDB_PROTEIN_URL.format(acc=accession.upper()))
        except Exception:                       # noqa: BLE001
            return None
        if not r or not r.get("entry_name"):
            return None
        # GPCRdb names carry markup ("A<sub>2A</sub> receptor"). It is somebody else's HTML
        # arriving over the wire, so it is stripped here rather than trusted downstream: the
        # renderer escapes, which would print the tags, and not escaping would be worse.
        return {"entry_name": r["entry_name"],
                "name": re.sub(r"<[^>]+>", "", r.get("name") or "").strip(),
                "family": r.get("family"), "scheme": r.get("residue_numbering_scheme")}

    return db.cached(("gpcrdb_protein", POCKET_VERSION, accession.upper()), fetch)


def gpcrdb_residues(entry_name: str) -> list[dict]:
    """Per-residue segment and generic number, keyed on the UniProt sequence number."""
    def fetch():
        try:
            rows = http.get_json(GPCRDB_RESIDUES_URL.format(entry=entry_name)) or []
        except Exception:                       # noqa: BLE001
            return []
        return [{"pos": r.get("sequence_number"), "aa": r.get("amino_acid"),
                 "segment": r.get("protein_segment"),
                 "generic": r.get("display_generic_number")}
                for r in rows if r.get("sequence_number")]

    return db.cached(("gpcrdb_residues", POCKET_VERSION, entry_name), fetch) or []


def _segments(residues: list[dict]) -> list[dict]:
    """Collapse the per-residue list into contiguous named segments (TM1, ICL2, ...)."""
    out: list[dict] = []
    for r in residues:
        seg = r.get("segment")
        if not seg:
            continue
        if out and out[-1]["name"] == seg and r["pos"] == out[-1]["end"] + 1:
            out[-1]["end"] = r["pos"]
        else:
            out.append({"name": seg, "start": r["pos"], "end": r["pos"]})
    return out


def _switches(residues: list[dict]) -> list[dict]:
    """Locate the micro-switches by generic number, and report what this receptor has there.

    Reported with the residues actually present rather than as a name alone: a receptor whose
    "DRY" is a DRF is exactly the kind of thing worth seeing, and naming the motif without
    showing the sequence would hide it.
    """
    by_generic = [(r, _num(r.get("generic"))) for r in residues]
    out = []
    for name, lo, hi, why in GPCR_SWITCHES:
        a, b = float(lo), float(hi)
        hits = [r for r, g in by_generic if g is not None and a - 1e-9 <= g <= b + 1e-9]
        if not hits:
            continue
        out.append({
            "name": name, "why": why,
            "generic": f"{lo}–{hi}" if lo != hi else lo,
            "start": hits[0]["pos"], "end": hits[-1]["pos"],
            "sequence": "".join(h.get("aa") or "?" for h in hits),
        })
    return out


# --------------------------------------------------------------------------------------
def build(fam: dict) -> dict:
    """Whichever reference numbering applies to this family, or none.

    A family that is neither a kinase nor a receptor gets an empty result, and the panel says
    so rather than showing an empty kinase pocket for a lysozyme.
    """
    seed = (fam.get("seed") or "").upper()
    if fam.get("kind") != "uniprot" or not seed:
        return {"kind": None}

    kin = klifs_kinase(seed)
    if kin and kin.get("kinase_id"):
        pdb_ids = {(m.get("pdb_id") or "").upper() for m in fam.get("members") or []}
        return {"kind": "kinase", "kinase": kin,
                "conformations": klifs_conformations(kin["kinase_id"], pdb_ids),
                "source": "KLIFS", "source_url": "https://klifs.net/"}

    rec = gpcrdb_protein(seed)
    if rec:
        residues = gpcrdb_residues(rec["entry_name"])
        return {"kind": "gpcr", "receptor": rec,
                "segments": _segments(residues), "switches": _switches(residues),
                "n_numbered": sum(1 for r in residues if r.get("generic")),
                "source": "GPCRdb", "source_url": "https://gpcrdb.org/"}

    return {"kind": None}
