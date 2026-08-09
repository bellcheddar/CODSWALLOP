"""Secondary structure and sheet topology for a family, in the seed's own coordinates.

Built in-house rather than scraped from PDBsum, whose diagrams are copyrighted.

Two things make this a family panel rather than a picture of one entry:

* Elements land on **seed coordinates**, the same axis as the conservation track, the
  coverage census, the domain ribbon and the motifs tab. A strand here is the same strand
  there, and the reader never has to translate between an entry's author numbering and the
  canonical sequence.
* The **sheet pairing** comes from DSSP's bridge partners, so the diagram shows which strands
  actually hydrogen-bond to which. That is the part a linear track cannot say and the part
  that makes a topology diagram worth drawing at all.

DSSP proper is used when `mkdssp` is on the PATH, and biotite's P-SEA implementation
otherwise. Which one produced the artefact is recorded in it: they do not always agree at the
ends of elements, and a reader comparing two families deserves to know they were measured the
same way.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from Bio.Align import PairwiseAligner, substitution_matrices

from . import embed, http, topology_io
from .topology_io import TOPOLOGY_DIR, VERSION

logger = logging.getLogger(__name__)

# DSSP's eight states, collapsed to the three a topology cartoon can draw. G (3-10) and I
# (pi) are helices and are drawn as such: splitting them out would put a second kind of
# cylinder on the diagram for a distinction the reader did not ask about.
HELIX = set("HGI")
STRAND = set("EB")

# Shorter than this and it is a wobble in the backbone rather than an element worth drawing.
MIN_HELIX = 4
MIN_STRAND = 2

_ALIGNER = PairwiseAligner()
_ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
_ALIGNER.open_gap_score = -11
_ALIGNER.extend_gap_score = -1
_ALIGNER.mode = "global"
# Free end gaps: a deposited chain is a fragment of the canonical, and charging it for the
# residues it never contained pushes the alignment off register. The combined setter, as in
# constructs.py: the separate target_/query_ properties are deprecated in Biopython 1.85.
_ALIGNER.end_gap_score = 0.0

_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M",
}


def dssp_available() -> Optional[str]:
    """The mkdssp binary, if there is one. CCP4 ships it and is where Marc's lives."""
    for candidate in ("mkdssp", "dssp"):
        found = shutil.which(candidate)
        if found:
            return found
    ccp4 = Path("/Applications/ccp4-9/bin/mkdssp")
    return str(ccp4) if ccp4.exists() else None


def _parse_dssp(text: str, chain: Optional[str]) -> list[dict]:
    """DSSP's fixed-column table into per-residue records.

    Fixed columns, not split(): the fields butt against each other on a residue with a
    four-digit number and a chain break line has almost nothing in it at all.
    """
    lines = text.splitlines()
    try:
        start = next(i for i, l in enumerate(lines) if l.startswith("  #  RESIDUE")) + 1
    except StopIteration:
        return []
    out = []
    for line in lines[start:]:
        if len(line) < 38 or line[13:14] == "!":      # chain break
            continue
        ch = line[11:12].strip()
        if chain and ch != chain:
            continue
        try:
            resnum = int(line[5:10])
        except ValueError:
            continue
        code = line[16:17]
        out.append({
            "res": resnum,
            "chain": ch,
            "aa": line[13:14],
            "sse": code if code.strip() else "-",
            # Bridge partners are DSSP's own sequential index, not the residue number, so
            # they are resolved against the '#' column after the whole table is read.
            "idx": int(line[0:5]) if line[0:5].strip() else None,
            "bp1": int(line[25:29]) if line[25:29].strip() else 0,
            "bp2": int(line[29:33]) if line[29:33].strip() else 0,
        })
    return out


def _sse_from_biotite(pdb_id: str, chain: Optional[str]) -> list[dict]:
    """P-SEA secondary structure, for a workstation with no DSSP.

    Three states only, and no bridge partners: the sheet pairing is simply absent rather
    than guessed at, and the artefact says which method produced it.
    """
    import biotite.structure as struc
    import biotite.structure.io.pdbx as pdbx

    path = embed._cif_path(pdb_id)
    if not path.exists():
        return []
    arr = pdbx.get_structure(pdbx.CIFFile.read(str(path)), model=1)
    arr = arr[struc.filter_amino_acids(arr)]
    if chain:
        arr = arr[arr.chain_id == chain]
    if arr.array_length() == 0:
        return []
    sse = struc.annotate_sse(arr)
    res_ids = struc.get_residues(arr)[0]
    names = struc.get_residues(arr)[1]
    out = []
    for i, r in enumerate(res_ids):
        if i >= len(sse):
            break
        code = {"a": "H", "b": "E", "c": "-"}.get(str(sse[i]), "-")
        out.append({"res": int(r), "chain": chain or "", "sse": code,
                    "aa": _THREE_TO_ONE.get(str(names[i]).upper(), "X"),
                    "idx": i + 1, "bp1": 0, "bp2": 0})
    return out


def secondary_structure(pdb_id: str, chain: Optional[str] = None) -> tuple:
    """(records, method) for one structure, DSSP where available and P-SEA otherwise."""
    binary = dssp_available()
    if binary:
        path = embed._cif_path(pdb_id)
        if not path.exists():
            embed.ca_trace(pdb_id, chain)      # downloads it as a side effect
        if path.exists():
            with tempfile.NamedTemporaryFile(suffix=".dssp", delete=False) as tmp:
                out_path = tmp.name
            try:
                proc = subprocess.run([binary, "--output-format", "dssp",
                                       str(path), out_path],
                                      capture_output=True, text=True, timeout=180)
                if proc.returncode == 0:
                    text = Path(out_path).read_text()
                    recs = _parse_dssp(text, chain)
                    if recs:
                        return recs, "DSSP"
                logger.info("mkdssp failed on %s: %s", pdb_id, (proc.stderr or "")[:200])
            except (OSError, subprocess.SubprocessError):
                logger.warning("mkdssp could not run on %s", pdb_id, exc_info=True)
            finally:
                Path(out_path).unlink(missing_ok=True)
    try:
        return _sse_from_biotite(pdb_id, chain), "P-SEA"
    except Exception:                               # noqa: BLE001
        logger.warning("no secondary structure for %s", pdb_id, exc_info=True)
        return [], "none"


def map_to_seed(records: list[dict], seed_sequence: str) -> dict:
    """Structure residue number -> seed position, by aligning the observed sequence.

    Aligned rather than offset. A deposited entry may be numbered on the mature protein, on
    the construct from 1, or on the canonical, and an assumed offset is right often enough
    to look like it works and wrong often enough to matter.
    """
    obs = "".join(r["aa"] if r["aa"].isalpha() else "X" for r in records)
    if not obs or not seed_sequence:
        return {}
    # Lowercase in DSSP marks a cysteine in a disulphide bond; the aligner wants the residue.
    obs = obs.upper()
    # Through the construct engine's own sanitiser, not a second copy of it. BLOSUM62 has no
    # U and Biopython refuses the whole alignment rather than skipping the letter, so one
    # selenocysteine in serum albumin cost that family its entire topology artefact. Fixing
    # it in constructs.py alone left this aligner, which has the same matrix and the same
    # problem, still failing.
    from .constructs import _alignable
    aln = _ALIGNER.align(_alignable(seed_sequence), _alignable(obs))[0]
    mapping = {}
    for (s_start, s_end), (o_start, o_end) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(s_end - s_start):
            rec = records[o_start + k]
            mapping[rec["res"]] = int(s_start + k) + 1        # seed positions are 1-based
    return mapping


def _elements(records: list[dict], mapping: dict) -> list[dict]:
    """Collapse per-residue codes into drawable elements, in seed coordinates."""
    out: list[dict] = []
    run: Optional[dict] = None

    def close(run):
        if not run:
            return
        length = run["end"] - run["start"] + 1
        floor = MIN_HELIX if run["kind"] == "helix" else MIN_STRAND
        if length >= floor:
            out.append(run)

    for r in records:
        pos = mapping.get(r["res"])
        kind = "helix" if r["sse"] in HELIX else ("strand" if r["sse"] in STRAND else None)
        if pos is None or kind is None:
            close(run)
            run = None
            continue
        # A gap in the seed numbering ends the element: two halves either side of an
        # unmodelled loop are two elements, not one long one spanning residues nobody saw.
        if run and run["kind"] == kind and pos == run["end"] + 1:
            run["end"] = pos
            run["residues"].append(r["idx"])
        else:
            close(run)
            run = {"kind": kind, "start": pos, "end": pos, "residues": [r["idx"]]}
    close(run)

    for i, e in enumerate(out):
        e["id"] = i
        e["length"] = e["end"] - e["start"] + 1
    return out


def _pairings(records: list[dict], elements: list[dict]) -> list[dict]:
    """Which strands hydrogen-bond to which, from DSSP's bridge partners.

    This is the part a linear track cannot say. Reported as element pairs with the number of
    bridges supporting each, so a single spurious bridge does not draw the same line as a
    ten-residue pairing.
    """
    by_idx = {}
    for e in elements:
        if e["kind"] != "strand":
            continue
        for idx in e["residues"]:
            by_idx[idx] = e["id"]

    counts: dict = {}
    for r in records:
        a = by_idx.get(r["idx"])
        if a is None:
            continue
        for bp in (r.get("bp1"), r.get("bp2")):
            b = by_idx.get(bp) if bp else None
            if b is None or b == a:
                continue
            key = (min(a, b), max(a, b))
            counts[key] = counts.get(key, 0) + 1
    return [{"a": a, "b": b, "bridges": n}
            for (a, b), n in sorted(counts.items(), key=lambda kv: -kv[1]) if n >= 2]


PDBE_TOPOLOGY_URL = "https://www.ebi.ac.uk/pdbe/api/topology/entry/{pdb_id}"


def pdbe_diagram(pdb_id: str, chain: Optional[str]) -> Optional[dict]:
    """PDBe's own 2D topology layout for one chain: strands, helices, coils and termini.

    Not drawn here. Laying out a fold in two dimensions is a real algorithm with a
    literature behind it, and the PDBe already runs one over the whole archive: this asks
    for the answer rather than inventing a worse one. What comes back is geometry, in the
    form of SVG coordinate paths, so the drawing is ours and the science is theirs.

    Independently checked against the DSSP run above: PDBe reports 3 strands for lysozyme
    1AKI and so does DSSP, which is the agreement worth having between two assignments made
    by different programs.
    """
    from . import db

    # Cached by hand rather than through `db.cached`, to keep two different answers apart.
    # "The PDBe has no topology for this entry" is a result and is worth remembering; "the
    # request failed" is not, and caching it poisons the entry permanently. That is exactly
    # what happened here: one transient error on the first call, and every later build read
    # back a cached None and reported the diagram as unavailable for a structure the PDBe
    # was serving perfectly well.
    # Its own key namespace. `db.cached` wraps every value in a one-element list so that a
    # cached None is not refetched forever, and reading its entries directly hands back the
    # wrapper: the first attempt at this returned [None] and the diagram came out as a list.
    # Two caching schemes must not share a key.
    key = db.cache_key("pdbe_topo_raw", VERSION, pdb_id.upper(), chain or "")
    hit = db.cache_get(key)
    if isinstance(hit, dict):
        return None if hit.get("none") else hit

    def fetch():
        try:
            body = http.get_json(PDBE_TOPOLOGY_URL.format(pdb_id=pdb_id.lower()))
        except Exception:                       # noqa: BLE001
            logger.warning("PDBe topology request failed for %s; not caching", pdb_id)
            raise
        entry = (body or {}).get(pdb_id.lower()) or {}
        # Keyed by entity id then chain. The chain we superposed on is the one to draw, and
        # anything else would be a diagram of a different molecule in the same crystal.
        for _entity, chains in entry.items():
            if not isinstance(chains, dict):
                continue
            got = chains.get(chain) if chain else None
            if got is None and chains:
                got = next(iter(chains.values()))
            if not got:
                continue
            return {
                "extents": got.get("extents"),
                "strands": [{"start": x.get("start"), "stop": x.get("stop"),
                             "path": x.get("path")} for x in got.get("strands") or []],
                "helices": [{"start": x.get("start"), "stop": x.get("stop"),
                             "path": x.get("path"), "major": x.get("majoraxis"),
                             "minor": x.get("minoraxis")} for x in got.get("helices") or []],
                "coils": [{"path": x.get("path")} for x in got.get("coils") or []],
                "terms": [{"type": x.get("type"), "resnum": x.get("resnum"),
                           "path": x.get("path")} for x in got.get("terms") or []],
            }
        return {"none": True}          # the PDBe genuinely has nothing for this chain

    try:
        got = fetch()
    except Exception:                           # noqa: BLE001
        # Logged with the traceback, not swallowed to a one-line note: a diagram that is
        # quietly absent is indistinguishable from one the PDBe does not have, and the two
        # want different responses.
        logger.warning("PDBe topology unavailable for %s %s", pdb_id, chain, exc_info=True)
        return None                             # transient: try again next build
    db.cache_put(key, got)
    return None if got.get("none") else got


def build(fam: dict) -> Optional[dict]:
    """The topology artefact for one family, computed from its reference structure."""
    slug = fam["slug"]
    seed = fam.get("seed_sequence") or ""
    if not seed:
        return None

    # The same structure the superposition uses, so the diagram and the 3D view agree about
    # which entry is the family's representative. Falls back to the best-resolution member
    # when no embedding has been built.
    art = None
    try:
        from . import embed_io
        art = embed_io.load(slug)
    except Exception:                               # noqa: BLE001
        art = None
    pdb_id = (art or {}).get("reference")
    chain = None
    if pdb_id:
        for r in (art or {}).get("representatives") or []:
            if r.get("pdb_id") == pdb_id:
                chain = r.get("chain")
                break
    if not pdb_id:
        members = [m for m in fam.get("members") or [] if m.get("resolution")]
        if not members:
            return None
        best = min(members, key=lambda m: m["resolution"])
        pdb_id, chain = best["pdb_id"], (best.get("chains") or [None])[0]

    records, method = secondary_structure(pdb_id, chain)
    if not records:
        return None

    mapping = map_to_seed(records, seed)
    if not mapping:
        return None
    elements = _elements(records, mapping)
    if not elements:
        return None
    pairs = _pairings(records, elements) if method == "DSSP" else []

    for e in elements:
        e.pop("residues", None)                     # DSSP indices: internal, not shipped

    return {
        "version": VERSION,
        "slug": slug,
        "reference": pdb_id,
        # The 2D fold layout, from the PDBe. Absent is a real answer: not every entry has
        # one, and drawing something else in its place would be worse than saying so.
        "diagram": pdbe_diagram(pdb_id, chain),
        "chain": chain,
        "method": method,
        "seed_length": len(seed),
        "elements": elements,
        "pairings": pairs,
        "n_helices": sum(1 for e in elements if e["kind"] == "helix"),
        "n_strands": sum(1 for e in elements if e["kind"] == "strand"),
        "coverage": round(100.0 * sum(e["length"] for e in elements) / len(seed), 1),
    }


def write(fam: dict) -> Optional[dict]:
    """Build and persist, returning the artefact."""
    art = build(fam)
    if not art:
        return None
    TOPOLOGY_DIR.mkdir(parents=True, exist_ok=True)
    path = topology_io.artefact_path(fam["slug"])
    tmp = path.with_suffix(".json.part")
    import json
    tmp.write_text(json.dumps(art))
    tmp.replace(path)
    return art
