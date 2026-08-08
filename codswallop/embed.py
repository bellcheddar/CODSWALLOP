"""The structural embedding: pairwise TM-scores, and the map positions they produce.

**This module runs on a workstation, never on the droplet.** It downloads mmCIF files and
does real numerical work, and the droplet has two cores shared with eight apps and 19 GB
free. The web app never imports it: it reads the JSON artefact this writes, and falls back
to the sequence-identity placeholder when there is none. `tmtools` and `biotite` are
therefore in requirements-dev.txt, not requirements.txt. (numpy is on the droplet either
way, as a biopython dependency; it is the mmCIF parsing and TM-align that are not.)

Two decisions worth stating, because both bound the cost:

* **One structure per distinct construct, not per entry.** A family is mostly the same
  construct solved repeatedly (carbonic anhydrase II: 1,490 entities, 232 distinct
  sequences), and superposing two crystal forms of an identical construct measures
  crystallography, not biology. Every member inherits the position of its construct's
  representative, so the map still shows every entry.
* **Capped, and the cap is on pairs.** All-versus-all is quadratic: 232 representatives is
  27,000 alignments, which is minutes; 800 would be 320,000, which is not. Past the cap the
  representatives are chosen by entity count, so the constructs most people actually made
  are the ones that get placed.
"""

from __future__ import annotations

import json
import logging
import math
import time
from pathlib import Path
from typing import Optional

import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure
from tmtools import tm_align

from . import config, http

logger = logging.getLogger(__name__)

# Above this many representatives the pair count stops being worth the wall clock.
MAX_REPRESENTATIVES = 260
# Chains longer than this are truncated for the alignment. TM-align is O(n^2) per pair too,
# and a 1,200-residue chain against 260 partners is an hour on its own.
MAX_RESIDUES = 700

EMBED_DIR = config.DATA_DIR / "embeddings"
STRUCT_DIR = config.DATA_DIR / "structures"

# The 20 standard residues, as biotite reports them, plus the selenomethionine that appears
# wherever anyone phased by SAD.
_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M",
}


# VERSION lives in embed_io and is imported, not redeclared. Two constants that must
# agree will eventually not: bumping it here for the superposition transforms while
# embed_io still said 1 would have silently rejected every artefact, and bumping
# embed_io alone silently accepted stale ones. Either way the map falls back to the
# placeholder with no error anywhere.
from .embed_io import VERSION, artefact_path, load  # noqa: F401


# --------------------------------------------------------------------------------------
# Structure fetching and parsing
# --------------------------------------------------------------------------------------
def _cif_path(pdb_id: str) -> Path:
    return STRUCT_DIR / f"{pdb_id.upper()}.cif"


def ca_trace(pdb_id: str, chain: Optional[str] = None) -> Optional[tuple]:
    """Alpha-carbon coordinates and the one-letter sequence for one chain.

    Returns (coords, sequence) or None. Chains are truncated at MAX_RESIDUES because
    TM-align is quadratic per pair as well as across the family.
    """
    path = _cif_path(pdb_id)
    if not path.exists():
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
        if http.download(url, path) is None:
            return None
    try:
        struct = get_structure(CIFFile.read(str(path)), model=1)
    except Exception:
        logger.warning("could not parse %s", pdb_id, exc_info=True)
        return None

    mask = struct.atom_name == "CA"
    if chain:
        mask &= struct.chain_id == chain
    ca = struct[mask]
    if ca.array_length() == 0:
        return None
    # First chain present, if the caller did not name one.
    if not chain:
        first = ca.chain_id[0]
        ca = ca[ca.chain_id == first]

    seq = "".join(_THREE_TO_ONE.get(r, "X") for r in ca.res_name)
    coords = np.asarray(ca.coord, dtype=np.float64)
    if len(seq) > MAX_RESIDUES:
        seq, coords = seq[:MAX_RESIDUES], coords[:MAX_RESIDUES]
    if len(seq) < 20:
        return None
    return coords, seq


# --------------------------------------------------------------------------------------
# The matrix, and the embedding
# --------------------------------------------------------------------------------------
def transforms_to_reference(traces: list[tuple], ref: int = 0) -> list[Optional[dict]]:
    """The rigid-body transform that puts each structure onto the reference.

    This falls out of the same TM-align run that builds the matrix, so superposition costs
    one extra row of alignments rather than a second pipeline. tmtools returns `u` (a 3x3
    rotation) and `t` (a translation) that map chain 1 onto chain 2, so aligning each
    structure AS chain 1 against the reference AS chain 2 gives the transform in the
    direction the viewer wants: everything onto the reference's frame.
    """
    out: list[Optional[dict]] = []
    cr, sr = traces[ref]
    for i, (ci, si) in enumerate(traces):
        if i == ref:
            out.append({"u": np.eye(3).tolist(), "t": [0.0, 0.0, 0.0], "tm": 1.0})
            continue
        try:
            r = tm_align(ci, cr, si, sr)
            out.append({
                "u": np.asarray(r.u, dtype=float).tolist(),
                "t": np.asarray(r.t, dtype=float).tolist(),
                "tm": round(float(max(r.tm_norm_chain1, r.tm_norm_chain2)), 3),
            })
        except Exception:
            # A structure that will not align is left without a transform rather than
            # given the identity, which would silently drop it in the wrong place.
            out.append(None)
    return out


def pairwise_tm(traces: list[tuple]) -> np.ndarray:
    """Symmetric TM-score matrix.

    TM-align is directional: normalising by chain 1 or chain 2 gives different numbers when
    the chains differ in length, which in this archive they routinely do (a truncated
    construct against a full-length one). The larger of the two is taken, which is the
    convention for asking "are these the same fold" rather than "is one contained in the
    other".
    """
    n = len(traces)
    tm = np.eye(n)
    for i in range(n):
        ci, si = traces[i]
        for j in range(i + 1, n):
            cj, sj = traces[j]
            try:
                r = tm_align(ci, cj, si, sj)
                score = max(r.tm_norm_chain1, r.tm_norm_chain2)
            except Exception:
                score = 0.0
            tm[i, j] = tm[j, i] = score
    return tm


def embed(tm: np.ndarray) -> np.ndarray:
    """Classical multidimensional scaling of the TM-score matrix, into two dimensions.

    Distance is 1 - TM, which is the usual reading: TM 1.0 is the same structure, and the
    0.5 mark is the conventional same-fold threshold, so a distance of 0.5 is the fold
    boundary and the map has a meaningful scale rather than an arbitrary one.

    Classical MDS by eigendecomposition rather than a scipy or sklearn dependency: the
    matrix is a few hundred square, this is twenty lines, and it keeps the pipeline's
    dependency list short.
    """
    n = tm.shape[0]
    if n < 3:
        return np.zeros((n, 2))
    d = 1.0 - tm
    np.fill_diagonal(d, 0.0)
    d2 = d ** 2
    # Double centring: B = -1/2 J d^2 J, with J = I - 1/n
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ d2 @ j
    vals, vecs = np.linalg.eigh(b)
    order = np.argsort(vals)[::-1][:2]
    vals, vecs = vals[order], vecs[:, order]
    # Negative eigenvalues mean the distances are not Euclidean, which 1-TM is not. Clamped
    # rather than failing: the two leading components still carry the structure, and this is
    # a map to look at, not a metric to compute on.
    coords = vecs * np.sqrt(np.maximum(vals, 0))
    span = np.abs(coords).max() or 1.0
    return coords / span


AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"


def alphafold_trace(accession: str, span: Optional[tuple] = None) -> Optional[tuple]:
    """Download the AlphaFold model for an accession and return its CA trace.

    The file URL comes from the API, never constructed: the model version is in the filename
    and it moves (the DB was on v4 when this was written and serves v6 now), so a built URL
    is a 404 waiting to happen.
    """
    from . import http

    try:
        rows = http.get_json(AFDB_API.format(accession=accession.upper()))
    except Exception:
        logger.warning("AlphaFold DB lookup failed for %s", accession, exc_info=True)
        return None
    if not rows or not rows[0].get("cifUrl"):
        return None

    url = rows[0]["cifUrl"]
    dest = STRUCT_DIR / f"AF-{accession.upper()}.cif"
    if http.download(url, dest) is None:
        return None
    try:
        struct = get_structure(CIFFile.read(str(dest)), model=1)
    except Exception:
        logger.warning("could not parse the AlphaFold model for %s", accession, exc_info=True)
        return None

    ca = struct[struct.atom_name == "CA"]
    if ca.array_length() == 0:
        return None
    seq = "".join(_THREE_TO_ONE.get(r, "X") for r in ca.res_name)
    coords = np.asarray(ca.coord, dtype=np.float64)

    # Trim to the part of the protein the reference structure actually covers, NOT to the
    # first MAX_RESIDUES. An AlphaFold model is full-length and a crystal structure rarely
    # is: EGFR's model is 1,210 residues and its reference covers seed residues 716-974, so
    # taking the first 700 handed TM-align the extracellular domain to superpose onto the
    # kinase domain. It scored 0.275, which reads as a bad prediction rather than a bad slice.
    #
    # The span is only applied when the caller could establish that seed coordinates ARE
    # canonical coordinates, which is not automatic. A family seeded from a PDB entity is
    # numbered against that construct: lysozyme's seed is the 129-residue MATURE protein
    # while P00698's canonical is 147 residues including an 18-residue signal peptide, so
    # slicing the model at seed [1,129] takes the signal peptide plus the wrong end and
    # dropped its TM from 0.993 to 0.854. Where the frames cannot be reconciled, no slice:
    # TM-align handles partial overlap perfectly well on its own.
    if span:
        beg, end = span
        lo, hi = max(0, int(beg) - 1), min(len(seq), int(end))
        if hi - lo >= 20:
            seq, coords = seq[lo:hi], coords[lo:hi]
    if len(seq) > MAX_RESIDUES:
        seq, coords = seq[:MAX_RESIDUES], coords[:MAX_RESIDUES]
    return (coords, seq, url) if len(seq) >= 20 else None


def alphafold_transform(accession: str, reference: tuple,
                        span: Optional[tuple] = None) -> Optional[dict]:
    """Align the AlphaFold model onto the reference structure.

    This is why the overlay can be a real superposition rather than two models sitting in
    different frames: the alignment is the same TM-align call used for everything else, and
    it belongs in the pipeline where TM-align lives rather than in a browser that has no
    way to do it.
    """
    got = alphafold_trace(accession, span)
    if not got:
        return None
    coords, seq, url = got
    cr, sr = reference
    try:
        r = tm_align(coords, cr, seq, sr)
    except Exception:
        logger.warning("could not align the AlphaFold model for %s", accession, exc_info=True)
        return None
    return {
        "accession": accession.upper(),
        "url": url,
        "u": np.asarray(r.u, dtype=float).tolist(),
        "t": np.asarray(r.t, dtype=float).tolist(),
        # How well the prediction matches the experiment, which is worth stating on the
        # panel: an overlay of something that does not superpose is a misleading picture.
        "tm": round(float(max(r.tm_norm_chain1, r.tm_norm_chain2)), 3),
        "length": len(seq),
        "span": list(span) if span else None,
    }


def choose_representatives(constructs: list, seed_acc: str, budget: int) -> list:
    """Which constructs get a node on the map, half of them reserved for the subject.

    Ranking purely by how many entities use a construct hands the map to whichever member of
    the superfamily has been deposited most. Of the 80 representatives an ABL1 search
    produced, two were ABL1, and both of those were its 63-residue SH3 domain: the reader
    asked about ABL1 and got a map of EGFR, superposed onto an EGFR structure. Reserving half
    the budget keeps the family context the app exists to show while guaranteeing the subject
    is on its own map, and in enough copies to anchor the frame.

    Families whose seed is already dominant are unaffected: the reserved half fills with the
    same constructs the open half would have taken.
    """
    ranked = sorted(range(len(constructs)), key=lambda i: -constructs[i]["n_entities"])
    mine = lambda i: (constructs[i].get("uniprot") or "").upper() == seed_acc   # noqa: E731
    own = [i for i in ranked if mine(i)]
    if not seed_acc or not own:
        return [constructs[i] for i in ranked[:budget]]

    # Usage still ranks the map. The quota only swaps the subject in when usage alone would
    # have left it under-represented, so a family whose seed is dominant keeps exactly the
    # nodes it had: partitioning the budget instead of topping it up displaced a construct
    # used by 92 entities with one used by a single entity.
    chosen = ranked[:budget]
    quota = min(len(own), max(1, budget // 2))
    have = [i for i in chosen if mine(i)]
    need = quota - len(have)
    if need > 0:
        seen = set(chosen)
        extras = [i for i in own if i not in seen][:need]
        others = [i for i in chosen if not mine(i)]
        drop = set(others[len(others) - len(extras):])
        chosen = [i for i in chosen if i not in drop] + extras
    return [constructs[i] for i in sorted(chosen, key=lambda i: -constructs[i]["n_entities"])]


def reference_index(reps: list, centrality, seed_acc: str) -> int:
    """Which representative everything else is superposed onto.

    Centrality (most similar to everything else) rather than best resolution, because
    superposing a family onto an outlier makes every other structure look wrong.

    Restricted first to the protein the reader actually asked for. Families are assembled at
    30% identity, so a search for ABL1 legitimately returns most of the tyrosine kinases,
    and the most central structure of that superfamily was an EGFR entry (4L7S): the page
    superposed ABL1's family onto a different protein, with no indication it had done so.
    Centrality still decides between the seed's own structures.
    """
    own = [i for i, r in enumerate(reps) if (r.get("uniprot") or "").upper() == seed_acc]
    if seed_acc and own:
        return int(max(own, key=lambda i: centrality[i]))
    return int(np.argmax(centrality))


def build(fam: dict, max_representatives: int = MAX_REPRESENTATIVES,
          progress=None) -> Optional[dict]:
    """Compute the embedding for one assembled family and return the artefact."""
    EMBED_DIR.mkdir(parents=True, exist_ok=True)
    STRUCT_DIR.mkdir(parents=True, exist_ok=True)

    constructs = fam.get("constructs") or []
    if not constructs:
        return None
    seed_acc = (fam.get("seed") or "").upper() if fam.get("kind") == "uniprot" else ""
    chosen = choose_representatives(constructs, seed_acc, max_representatives)

    by_pdb = {}
    for m in fam["members"]:
        by_pdb.setdefault(m["pdb_id"], m)

    reps, traces = [], []
    for i, c in enumerate(chosen):
        pdb_id = c.get("best_pdb_id")
        if not pdb_id:
            continue
        member = by_pdb.get(pdb_id)
        chain = (member.get("chains") or [None])[0] if member else None
        if progress:
            progress("fetch", i + 1, len(chosen), pdb_id)
        trace = ca_trace(pdb_id, chain)
        if trace is None:
            continue
        reps.append({"seq_id": c["seq_id"], "pdb_id": pdb_id, "chain": chain,
                     "n_entities": c["n_entities"],
                     "uniprot": (member or {}).get("uniprot")})
        traces.append(trace)

    if len(reps) < 3:
        logger.warning("only %d usable structures for %s; not embedding", len(reps), fam["slug"])
        return None

    if progress:
        progress("align", 0, len(reps) * (len(reps) - 1) // 2, "")
    t0 = time.time()
    tm = pairwise_tm(traces)
    coords = embed(tm)
    # The reference for superposition: the representative most similar to everything else,
    # not simply the best-resolution one. Superposing a family onto an outlier makes every
    # other structure look wrong.
    #
    ref = reference_index(reps, tm.sum(axis=1), seed_acc)
    transforms = transforms_to_reference(traces, ref)

    # The AlphaFold model, aligned onto the same reference so the viewer can superpose it
    # like any other structure rather than dropping it in its own frame.
    from collections import Counter as _Counter
    accs = _Counter(m["uniprot"] for m in fam["members"] if m.get("uniprot"))
    af = None
    if accs:
        if progress:
            progress("fetch", len(chosen), len(chosen), "AlphaFold")
        # The seed span of the reference structure, so the model is trimmed to the same
        # region rather than to its own first 700 residues -- but ONLY when the family's
        # seed is that same UniProt sequence, because otherwise the two are not in the same
        # coordinate frame and the slice lands in the wrong place.
        # The seed's own accession, never the family's modal one. The modal accession is a
        # popularity contest the subject frequently loses: an ABL1 family at 30% identity is
        # mostly EGFR by entry count, so it fetched EGFR's model, and every A2A receptor
        # structure carries a BRIL fusion, so A2A fetched the model of E. coli cytochrome
        # b562. Both rendered as a confident superposition of the wrong protein.
        accession = seed_acc or accs.most_common(1)[0][0]
        ref_pdb = reps[ref]["pdb_id"]
        ref_member = by_pdb.get(ref_pdb) or {}
        span = None
        seed_is_canonical = bool(seed_acc) and seed_acc == accession.upper()
        if seed_is_canonical and ref_member.get("query_beg") and ref_member.get("query_end"):
            span = (ref_member["query_beg"], ref_member["query_end"])
        af = alphafold_transform(accession, traces[ref], span)
    elapsed = time.time() - t0

    artefact = {
        "version": VERSION,
        "slug": fam["slug"],
        "built_at": int(time.time()),
        "n_representatives": len(reps),
        "n_pairs": len(reps) * (len(reps) - 1) // 2,
        "seconds": round(elapsed, 1),
        "reference": reps[ref]["pdb_id"],
        "alphafold": af,
        "representatives": [
            {**r, "x": round(float(coords[i, 0]), 4), "y": round(float(coords[i, 1]), 4),
             "transform": transforms[i]}
            for i, r in enumerate(reps)
        ],
        # The matrix itself, rounded: the Structures panel draws it as a clustered heatmap,
        # and two decimals is well inside TM-align's own reproducibility.
        "tm": [[round(float(v), 2) for v in row] for row in tm],
        "median_tm": round(float(np.median(tm[np.triu_indices(len(reps), 1)])), 3),
    }
    artefact_path(fam["slug"]).write_text(json.dumps(artefact, separators=(",", ":")))
    return artefact
