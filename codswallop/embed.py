"""The structural embedding: pairwise TM-scores, and the map positions they produce.

**This module runs on the droplet now.** It used to say the opposite, and the reason it
could not was memory: parsing a whole deposited assembly to extract one chain's alpha
carbons peaked at 6.3 GB on a box with 2.2 GB free and no swap. Fetching the single chain
from the Model Server instead peaks at about 100 MB, so the droplet runs it comfortably and
two families at once barely move the free-memory figure. `tmtools` and `biotite` live in
requirements-compute.txt and are installed there.

The WEB process still never imports this. It reads the JSON artefact through `embed_io`,
which imports nothing beyond the standard library, and falls back to the sequence-identity
placeholder when there is none: a web worker that has never imported biotite cannot be
killed by parsing a structure.

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


# One chain from the RCSB's Model Server instead of the whole deposited file.
#
# This is what makes the embedding runnable anywhere but a workstation. A large assembly is
# hundreds of megabytes of text and biotite holds all of it in memory to parse it: 8GLV is a
# 453 MB file that peaks at 6.3 GB RSS, and yields 426 alpha carbons. gemmi is twelve times
# faster on the same file and still needs about six gigabytes, because the file really is
# that big; the parser was never the problem. Asking for the one chain we want returns 533 kB
# instead, an 850-fold reduction, and takes peak memory into the tens of megabytes.
#
# The trade is latency: the Model Server assembles the response, so a request takes ten to
# twenty seconds against under one for the static file. That is the right trade for a machine
# with 2 GB of RAM and no swap, where the alternative is not "slower" but "killed".
MODEL_SERVER = "https://models.rcsb.org/v1/{pdb}/atoms"
# Only fetch the whole file when it is small enough to be safe to parse. Measured against the
# structures already cached here: the median is about 1 MB and the tail runs to 453 MB.
MAX_WHOLE_FILE_MB = 12


def _chain_cif(pdb_id: str, chain: str) -> Optional[Path]:
    """One chain, cached, from the Model Server. None if it cannot be had."""
    path = STRUCT_DIR / f"{pdb_id.upper()}_{chain}.cif"
    if path.exists():
        return path
    url = MODEL_SERVER.format(pdb=pdb_id.lower())
    params = {"auth_asym_id": chain, "encoding": "cif", "copy_all_categories": "false"}
    if http.download(url, path, params=params) is None:
        return None
    # A Model Server miss is a 200 with a near-empty body rather than a 404, so size is the
    # only signal that the chain was not there.
    if path.stat().st_size < 2048:
        path.unlink(missing_ok=True)
        return None
    return path


def structure_path(pdb_id: str, chain: Optional[str] = None) -> Optional[Path]:
    """A local mmCIF for this structure, preferring the single chain.

    The one place anything else should ask for a file. `ca_trace` used to download the whole
    entry as a side effect and topology quietly depended on that: when the fetch moved to the
    Model Server, `_cif_path` stopped existing and topology got no records at all, with no
    error, because a missing file was already a legitimate "no layout for this structure".
    A side effect that another module relies on is not an interface.
    """
    if chain:
        got = _chain_cif(pdb_id, chain)
        if got is not None:
            return got
    path = _cif_path(pdb_id)
    if not path.exists():
        if http.download(f"https://files.rcsb.org/download/{pdb_id.upper()}.cif",
                         path) is None:
            return None
    if path.stat().st_size / 1e6 > MAX_WHOLE_FILE_MB:
        logger.warning("skipping %s: %.0f MB whole-file parse", pdb_id,
                       path.stat().st_size / 1e6)
        return None
    return path


def ca_trace(pdb_id: str, chain: Optional[str] = None) -> Optional[tuple]:
    """Alpha-carbon coordinates and the one-letter sequence for one chain.

    Returns (coords, sequence) or None. Chains are truncated at MAX_RESIDUES because
    TM-align is quadratic per pair as well as across the family.
    """
    path = structure_path(pdb_id, chain)
    if path is None:
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


def pairwise_tm(traces: list[tuple], progress=None) -> np.ndarray:
    """Symmetric TM-score matrix.

    TM-align is directional: normalising by chain 1 or chain 2 gives different numbers when
    the chains differ in length, which in this archive they routinely do (a truncated
    construct against a full-length one). The larger of the two is taken, which is the
    convention for asking "are these the same fold" rather than "is one contained in the
    other".
    """
    n = len(traces)
    tm = np.eye(n)
    total = n * (n - 1) // 2
    done = 0
    # Reported from inside the loop, not once before it. This is the longest silent stretch
    # in the whole pipeline: spike is 32 minutes of it, and an unattended supervisor watching
    # for silence cannot tell that from a wedged socket. Announcing the alignment before
    # starting it says nothing about whether it is still going.
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
            done += 1
            # Every ten pairs, not every two hundred. The consumer throttles by time, and
            # it is the only one that can: a pair of 60-residue chains aligns in a
            # millisecond and a pair of 766-residue ones takes a second and a half, so any
            # fixed pair count is either noise on small families or five minutes of silence
            # on large ones. Reporting is a function call; throttling is the caller's job.
            if progress and done % 10 == 0:
                progress("align", done, total, "")
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
        # The construct's accession, which the family assigned, and not the member's first
        # UniProt cross-reference. The order RCSB lists cross-references in is not
        # meaningful: 4EIY, the classic A2A-BRIL structure, lists BRIL first, so 24 of A2A's
        # 36 representatives did not look like A2A and the reference was chosen from 12
        # candidates instead of 36. This is the same trap the construct diff already fell
        # into with 2RH1 and T4 lysozyme.
        reps.append({"seq_id": c["seq_id"], "pdb_id": pdb_id, "chain": chain,
                     "n_entities": c["n_entities"],
                     "uniprot": c.get("uniprot") or (member or {}).get("uniprot")})
        traces.append(trace)

    if len(reps) < 3:
        logger.warning("only %d usable structures for %s; not embedding", len(reps), fam["slug"])
        return None

    if progress:
        progress("align", 0, len(reps) * (len(reps) - 1) // 2, "")
    t0 = time.time()
    tm = pairwise_tm(traces, progress=progress)
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
