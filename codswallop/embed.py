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

VERSION = 2

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


from .embed_io import artefact_path, load  # noqa: F401  (re-exported for the CLI)


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


def build(fam: dict, max_representatives: int = MAX_REPRESENTATIVES,
          progress=None) -> Optional[dict]:
    """Compute the embedding for one assembled family and return the artefact."""
    EMBED_DIR.mkdir(parents=True, exist_ok=True)
    STRUCT_DIR.mkdir(parents=True, exist_ok=True)

    constructs = fam.get("constructs") or []
    if not constructs:
        return None
    chosen = sorted(constructs, key=lambda c: -c["n_entities"])[:max_representatives]

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
                     "n_entities": c["n_entities"]})
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
    ref = int(np.argmax(tm.sum(axis=1)))
    transforms = transforms_to_reference(traces, ref)
    elapsed = time.time() - t0

    artefact = {
        "version": VERSION,
        "slug": fam["slug"],
        "built_at": int(time.time()),
        "n_representatives": len(reps),
        "n_pairs": len(reps) * (len(reps) - 1) // 2,
        "seconds": round(elapsed, 1),
        "reference": reps[ref]["pdb_id"],
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
