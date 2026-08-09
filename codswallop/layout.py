"""Placeholder layout for the hero map.

**This is scaffolding with a known replacement date.** Phase 3 computes a pairwise TM-score
matrix for the family and positions every node in a real structural embedding. Until then,
the map has to show something, and the honest something is the one similarity measure Phase
1 actually has: percent identity to the seed, which the RCSB sequence search hands us per
member for free.

So: distance from the centre is sequence distance from the seed, and the angular sector is
the source organism. That makes orthologue groups legible (the thing a family map is for)
without pretending to structural knowledge the app does not yet have. The panel says so on
the face of it, and the Structures pip in the rail carries the phase number.
"""

from __future__ import annotations

import math
from typing import Optional

# How many organism clusters get their own labelled sector before the tail is pooled.
MAX_CLUSTERS = 8
# Total edges drawn. Past this the constellation stops reading as structure and starts
# reading as a hairball, and the SVG gets slow on a phone.
MAX_EDGES = 600
NEIGHBOURS = 2

# ---------------------------------------------------------------------------------------
# Stacking.
#
# Members sharing a construct share an MDS coordinate exactly, and used to be fanned out
# around it by a golden-angle spiral of radius 0.012*sqrt(k). That is a sunflower packing,
# and past a couple of dozen members it draws its own parastichies: the visible "arms" in a
# dense family were phyllotaxis, not structure. Worse, the radius was never bounded, so on
# ABL1 (137 members on one point) the decoration spanned 0.281 against a median
# representative-to-representative distance of 0.104. The artefact was 2.7x the signal, and
# two nodes on opposite rims of one disc are the same sequence.
#
# So: past STACK_MAX a point is drawn as a single node scaled by how many entries it holds,
# and expands on click. The fan survives below that, where it does what it was meant to do.
STACK_MAX = 12
# Even an expanded stack is bounded, as a fraction of the distance to the nearest other
# representative: a fan must never reach far enough to be mistaken for a neighbouring group.
FAN_FRACTION = 0.35
# ... but never collapsed to nothing, or a tight cluster's members become unpickable.
FAN_MIN = 0.02

# Where the dendrogram is cut for the cluster labels printed on the field. The same value as
# `DEFAULTS.cut` in family.js, and the same average linkage as `TmHeatmap.cluster`, so the
# labels on the map name the same groups the Structures heatmap shows and the cluster filter
# selects. Two clusterings of one matrix that disagree would be worse than none.
CLUSTER_CUT = 0.18
# Clusters below this share of the family are left unlabelled: the field has room for a
# handful of names, and a singleton's label costs more legibility than it buys.
LABEL_MIN_SHARE = 0.02
MAX_LABELS = 6


def _short_organism(name: Optional[str]) -> str:
    """`Escherichia coli K-12` -> `E. coli`. Cluster labels are printed on the field, so
    they have to fit on it."""
    if not name:
        return "Unknown"
    parts = name.split()
    if len(parts) >= 2 and parts[0][:1].isupper():
        return f"{parts[0][0]}. {parts[1]}"
    return parts[0]


def _golden_jitter(i: int) -> float:
    """Deterministic spread in [0, 1). The golden ratio keeps successive values far apart,
    so members of a cluster fill their sector evenly instead of clumping the way a random
    jitter does at small counts. Deterministic so a family lays out identically every time:
    a map that reshuffles on reload is a map nobody trusts."""
    return (i * 0.6180339887498949) % 1.0


def third_axis(tm: list) -> Optional[list]:
    """The third principal coordinate of the TM matrix, so the map can be rotated.

    Computed here rather than stored in the artefact, because the artefact already ships the
    whole matrix for the heatmap: the third axis is a twenty-line eigendecomposition of
    something already on hand, and deriving it costs less than a pipeline version bump would
    have cost in rebuilds. Same classical MDS as `embed.embed`, one component further along.

    numpy only, and only here: it is already installed everywhere as a biopython dependency,
    which is not true of biotite or tmtools.
    """
    try:
        import numpy as np
    except ImportError:                         # pragma: no cover - numpy is a hard dep
        return None
    n = len(tm or [])
    if n < 4:
        return None
    try:
        m = np.asarray(tm, dtype=float)
        if m.shape != (n, n):
            return None
        d = 1.0 - m
        np.fill_diagonal(d, 0.0)
        j = np.eye(n) - np.ones((n, n)) / n
        b = -0.5 * j @ (d ** 2) @ j
        vals, vecs = np.linalg.eigh(b)
        order = np.argsort(vals)[::-1][:3]
        if len(order) < 3 or vals[order[2]] <= 0:
            # A family that is genuinely flat has no third axis to show, and inventing one
            # from a negative eigenvalue would be drawing noise as structure.
            return None
        coords = vecs[:, order] * np.sqrt(np.maximum(vals[order], 0))
        # Scaled with x and y, not independently: the point of the third axis is that it is
        # the same kind of distance, so normalising it alone would exaggerate a shallow one.
        scale = float(np.max(np.abs(coords[:, :2]))) or 1.0
        return [round(float(v / scale), 4) for v in coords[:, 2]]
    except Exception:                           # noqa: BLE001
        return None


def _average_linkage(tm: list, cut: float) -> Optional[list]:
    """Cluster the representatives, returning one cluster index per representative.

    Average linkage on 1 - TM, cut at `cut`. Deliberately the same algorithm and the same
    default height as `TmHeatmap.cluster`/`cut` in the browser, because the map's labels and
    the heatmap's groups should name the same things: a reader who cuts the matrix into four
    groups and then counts five labelled clumps on the map has been told the app cannot keep
    its story straight.

    "Should", not "will". Checked against scipy's `linkage(method="average")` over all 71
    built families: 69 partitions identical, and the two that differ (ribonuclease, spike)
    agree exactly once ties are broken by a 1e-9 jitter. The matrix ships rounded to 2 dp for
    the heatmap, which leaves tie groups of up to 344 equal pairs, and any two
    implementations then merge in whatever order they happen to scan. Both trees are valid
    average linkage; a handful of borderline representatives can fall either side. The counts
    printed on the labels are exact for the grouping shown.

    Single linkage would be the cheaper choice and is useless here: every member of a family
    is the same fold, so a single-linkage tree chains straight through and puts 244 of
    hemoglobin's 244 representatives in one group.

    numpy, guarded, exactly as `third_axis`: it arrives as a biopython dependency rather than
    being requested, and a droplet without it should lose the labels, not the page.
    """
    try:
        import numpy as np
    except ImportError:                         # pragma: no cover - numpy is a hard dep
        return None
    n = len(tm or [])
    if n == 0:
        return None
    if n == 1:
        return [0]
    try:
        m = np.asarray(tm, dtype=float)
        if m.shape != (n, n):
            return None
        d = 1.0 - m
        np.fill_diagonal(d, 0.0)
        # The stored matrix is rounded to 2 dp per row, so it is very nearly but not exactly
        # symmetric. Averaging costs nothing and stops the merge order depending on which
        # triangle a value was read from.
        d = (d + d.T) / 2.0
    except Exception:                           # noqa: BLE001
        return None

    inf = float("inf")
    work = d.copy()
    np.fill_diagonal(work, inf)
    active = np.ones(n, dtype=bool)
    size = np.ones(n, dtype=float)
    leaves = [[i] for i in range(n)]

    for _ in range(n - 1):
        sub = np.where(active[:, None] & active[None, :], work, inf)
        flat = int(np.argmin(sub))
        i, j = flat // n, flat % n
        height = float(sub[i, j])
        # Average linkage is monotonic, so the first merge above the cut height is also the
        # last: everything remaining is further apart than this.
        if not (height == height) or height == inf or height > cut:
            break
        si, sj = float(size[i]), float(size[j])
        row = (work[i] * si + work[j] * sj) / (si + sj)
        work[i, :] = row
        work[:, i] = row
        work[i, i] = inf
        size[i] = si + sj
        leaves[i] = leaves[i] + leaves[j]
        active[j] = False

    out = [0] * n
    order = sorted((i for i in range(n) if active[i]),
                   key=lambda i: -len(leaves[i]))
    for c, i in enumerate(order):
        for leaf in leaves[i]:
            out[leaf] = c
    return out


def _tidy_description(text: Optional[str]) -> str:
    """A deposited description, made fit to print on the field.

    Entries from the 1990s are shouted in full caps (`PROTO-ONCOGENE TYROSINE-PROTEIN KINASE
    ABL`) and sit next to modern mixed-case ones in the same cluster, so the raw mode of the
    two is decided by typography rather than by content.
    """
    t = (text or "").strip()
    if not t:
        return ""
    letters = [c for c in t if c.isalpha()]
    if letters and all(c.isupper() for c in letters):
        t = t.title()
    for prefix in ("Isoform Short Of ", "Isoform Long Of ",
                   "Isoform Short of ", "Isoform Long of "):
        if t.startswith(prefix):
            t = t[len(prefix):]
            break
    return t


# Words a truncated label must not end on. Pfam's human-readable names are frequently long
# conjunctions ("Protein tyrosine and serine/threonine kinase"), and cutting one at a word
# boundary alone left "Protein tyrosine and…", which reads as a mistake rather than as an
# abbreviation.
_TRAILING = {"and", "or", "of", "the", "a", "an", "in", "to", "with", "for", "from", "on"}


def _shorten(text: str, limit: int = 40) -> str:
    if len(text) <= limit:
        return text
    cut = text[:limit].rsplit(" ", 1)[0]
    if len(cut) < limit // 2:
        cut = text[:limit]
    words = cut.split(" ")
    while len(words) > 1 and words[-1].lower().strip(",;:-") in _TRAILING:
        words.pop()
    return " ".join(words).rstrip(" ,;:-") + "…"


# A cluster is named after one protein only when that protein actually owns it. Below this
# the cluster is a mixture and gets named after what its members share instead.
DOMINANT_ACCESSION = 0.5
DOMINANT_DOMAIN = 0.6


def _domains_of(member: dict, annotations: Optional[dict]) -> list:
    """The member's Pfam list, whether or not the payload has been compacted yet.

    `_compact` lifts the Pfam and InterPro lists off every member into one family-level
    lookup and leaves an `annot_id` behind, because near-identical domain lists on 2,000
    members was the heaviest field in the payload. The map is built after that has happened,
    so reading `m["pfam"]` here finds nothing on every real request: checked against the live
    server, 0 of 2,000 members carried it, which would have left the domain rule below dead
    and every mixed cluster named after whichever protein held a plurality.
    """
    direct = member.get("pfam")
    if direct:
        return direct
    if annotations:
        entry = annotations.get(member.get("annot_id"))
        if entry:
            return entry.get("pfam") or []
    return []


def _label_clusters(groups: list[list[dict]], annotations: Optional[dict] = None) -> list[str]:
    """Name each cluster after whatever is actually true of it.

    Counted over entities rather than over representatives: a representative can stand for a
    hundred entries, and naming a cluster after whichever construct happened to be picked
    would let one PDB entry outvote a hundred.

    Three rules, in order, because a family search at 30 % identity does not return one
    protein. It returns a superfamily, and the biggest structural cluster is usually a
    mixture with no majority member at all:

    1. One UniProt accession holds half the cluster -> name it after that protein. Keyed on
       the accession and not on the deposited description, which splits the vote: ABL1 is
       variously "Tyrosine-protein kinase ABL1", "Proto-oncogene tyrosine-protein kinase
       ABL1" and "ABL TYROSINE KINASE", so counting free text lets a contaminant with one
       consistent spelling out-poll the family's own subject.

    2. No majority protein, but a shared Pfam domain -> name the domain. ABL1's largest
       cluster is 1,932 entities of which the top accession (EGFR) is 17 %; calling it
       "Epidermal growth factor receptor" states something false. 95 % of it carries the
       protein tyrosine kinase domain, which is both true and the reason those structures
       cluster.

    3. Neither -> the plurality description, which at that point is the best on offer.

    Where two clusters still take the same name, the median construct length disambiguates
    them, because that is usually what separates them: ABL1 splits into a ~290 aa kinase
    domain and a ~62 aa SH3/SH2 fragment set, both truthfully "Tyrosine-protein kinase ABL1".
    """
    import collections

    base, lengths = [], []
    for group in groups:
        n = len(group) or 1
        accessions = collections.Counter()
        by_accession: dict = {}
        domains = collections.Counter()
        descriptions = collections.Counter()
        lens = []
        for m in group:
            desc = _tidy_description(m.get("description"))
            if desc:
                descriptions[desc] += 1
            acc = m.get("uniprot")
            if acc:
                accessions[acc] += 1
                if desc:
                    by_accession.setdefault(acc, collections.Counter())[desc] += 1
            # Counted once per member however many copies of the domain it carries, so the
            # share is "how much of this cluster has it" rather than a domain census.
            for dom in {(d or {}).get("name") for d in _domains_of(m, annotations)}:
                if dom:
                    domains[dom] += 1
            if m.get("seq_length"):
                lens.append(int(m["seq_length"]))

        label = ""
        if accessions:
            acc, count = accessions.most_common(1)[0]
            if count / n >= DOMINANT_ACCESSION:
                named = by_accession.get(acc)
                label = named.most_common(1)[0][0] if named else ""
        if not label and domains:
            dom, count = domains.most_common(1)[0]
            if count / n >= DOMINANT_DOMAIN:
                label = _domain_name(dom)
        if not label and descriptions:
            label = descriptions.most_common(1)[0][0]
        if not label:
            orgs = collections.Counter(_short_organism(m.get("organism")) for m in group)
            label = orgs.most_common(1)[0][0] if orgs else "Group"

        base.append(label)
        lens.sort()
        lengths.append(lens[len(lens) // 2] if lens else None)

    clashes = {name for name, c in collections.Counter(base).items() if c > 1}
    out = []
    for label, median_len in zip(base, lengths):
        if label in clashes and median_len:
            out.append(_shorten(label, 30) + ", " + str(median_len) + " aa")
        else:
            out.append(_shorten(label))

    # Length does not always separate them. Alpha-synuclein's six clusters are six
    # conformations of one 140-residue protein, so all six came back "Alpha-synuclein,
    # 140 aa" and the field printed the same label six times. The groups are real and the
    # nodes stay; only the largest keeps the name, because a label repeated is a label that
    # identifies nothing. `groups` arrives largest-first, so the one that keeps it is the one
    # worth naming.
    seen: set = set()
    for i, label in enumerate(out):
        if label in seen:
            out[i] = ""
        else:
            seen.add(label)
    return out


def _domain_name(name: str) -> str:
    """`Protein kinase domain (Pkinase)` -> `Protein kinase domain`.

    Pfam carries the accession's short name in brackets after the human one. It is the more
    precise of the two and the less readable, and the field has room for one.
    """
    text = (name or "").strip()
    if text.endswith(")") and "(" in text:
        head = text[:text.rindex("(")].strip()
        if head:
            return head
    return text


def _nearest_neighbour(points: list[tuple]) -> dict:
    """Distance from each distinct representative point to the closest other one.

    This is what bounds the fan: a group of identical constructs may be spread far enough to
    be individually pickable, and no further than a share of the way to the nearest thing it
    could be confused with.
    """
    out: dict = {}
    for i, a in enumerate(points):
        best = float("inf")
        for j, b in enumerate(points):
            if i == j:
                continue
            dist = math.hypot(a[0] - b[0], a[1] - b[1])
            if dist < best:
                best = dist
        out[a] = best if best < float("inf") else 1.0
    return out


def _from_embedding(members: list[dict], embedding: dict,
                    annotations: Optional[dict] = None) -> dict:
    """Position every member from the pairwise TM-score matrix.

    Representatives carry their MDS coordinates directly. Every other member inherits the
    position of its own construct's representative, which is exact: they are the same
    sequence, so they would land on the same point anyway.

    Members whose construct was not among the representatives (the cap is on pairs, and it
    bites on a family with hundreds of distinct constructs) inherit from the representative
    closest in identity to the seed. That is an approximation and it is labelled as one on
    the panel, rather than being left to look like a measurement.
    """
    reps = embedding["representatives"]
    # Attach the third coordinate before anything reads a representative, so the nodes and
    # the "z" flag below agree about whether this family has one.
    zs = third_axis(embedding.get("tm") or [])
    if zs and len(zs) == len(reps):
        for r, z in zip(reps, zs):
            r["z"] = z
    by_seq = {r["seq_id"]: r for r in reps}
    # Representatives sorted by identity, for the nearest-identity fallback.
    ident_of = {}
    for m in members:
        r = by_seq.get(m.get("seq_id"))
        if r and m.get("identity") is not None:
            ident_of.setdefault(r["seq_id"], m["identity"])
    ranked = sorted(((ident_of.get(r["seq_id"], 0.0), r) for r in reps), key=lambda t: t[0])

    # ---- pass 1: which representative does each member sit on? -------------------------
    # Separated from placement because the fan needs to know how many members share a point
    # before it can decide how far to spread the first of them.
    #
    # The fallback used to be a plain `min(ranked, key=|identity - target|)`. In these
    # families hundreds of representatives sit at 100 % identity, so the tie fell to sort
    # order and every ambiguous member landed on one representative: 126 of ABL1's, 194 of
    # spike's, which is most of what made the largest discs. The tie carries no information,
    # so it is spread round-robin over the representatives that tie rather than concentrated
    # into a group that looks like a finding. It is still a guess, and it is now flagged as
    # one on the node.
    tie_turn: dict[float, int] = {}
    assigned, approximated = [], 0
    for m in members:
        r = by_seq.get(m.get("seq_id"))
        approx = r is None
        if approx:
            approximated += 1
            target = m.get("identity")
            if target is None or not ranked:
                r = reps[0]
            else:
                best = min(abs(ident - target) for ident, _ in ranked)
                tied = [rep for ident, rep in ranked if abs(ident - target) <= best + 1e-9]
                turn = tie_turn.get(best, 0)
                tie_turn[best] = turn + 1
                r = tied[turn % len(tied)]
        assigned.append((m, r, approx))

    per_point: dict[tuple, int] = {}
    for _, r, _ in assigned:
        key = (r["x"], r["y"])
        per_point[key] = per_point.get(key, 0) + 1
    near = _nearest_neighbour(list(per_point))

    # ---- clusters, from the same tree the Structures heatmap cuts ----------------------
    cluster_of = _average_linkage(embedding.get("tm") or [], CLUSTER_CUT)
    rep_cluster = {}
    if cluster_of and len(cluster_of) == len(reps):
        rep_cluster = {reps[i]["seq_id"]: c for i, c in enumerate(cluster_of)}
    grouped: dict[int, list] = {}
    for m, r, _ in assigned:
        grouped.setdefault(rep_cluster.get(r["seq_id"], 0), []).append(m)
    # Which clusters earn a label is settled BEFORE they are named, so that clashes are
    # resolved only among the names a reader will actually see. Deciding it the other way
    # round let a 46-entity cluster that never gets drawn force the 1,827-entity one next to
    # it to carry a disambiguating suffix against a label that was not on the field.
    floor = max(1, LABEL_MIN_SHARE * len(assigned))
    keys = [k for k in sorted(grouped, key=lambda k: -len(grouped[k])) if len(grouped[k]) >= floor]
    keys = keys[:MAX_LABELS]
    names = {}
    if rep_cluster and keys:
        for key, name in zip(keys, _label_clusters([grouped[k] for k in keys], annotations)):
            # A cluster whose name duplicated a larger one's comes back blank, and an unnamed
            # label is not a label: it would print as a bare count floating on the field.
            if name:
                names[key] = name

    # ---- pass 2: place ------------------------------------------------------------------
    nodes = []
    seen: dict[tuple, int] = {}
    for m, r, approx in assigned:
        key = (r["x"], r["y"])
        k = seen.get(key, 0)
        seen[key] = k + 1
        total = per_point[key]
        # Identical constructs share a point exactly, so spread them just enough to be
        # individually pickable without implying a difference that is not there. Bounded by
        # a share of the distance to the nearest other representative, so the fan can never
        # reach far enough to be read as a group of its own: unbounded, ABL1's largest fan
        # was 2.7x the median distance between genuinely distinct structures.
        limit = max(FAN_MIN, FAN_FRACTION * near.get(key, 1.0))
        radius = min(0.012 * math.sqrt(k), limit)
        angle = 2 * math.pi * _golden_jitter(k)
        ci = rep_cluster.get(r["seq_id"], 0)
        node = {
            "id": m["entity_id"], "pdb_id": m["pdb_id"],
            "cluster": ci, "cluster_name": names.get(ci, ""),
            "x": round(r["x"] + radius * math.cos(angle), 4),
            "y": round(r["y"] + radius * math.sin(angle), 4),
            # The point this node shares, so the renderer can draw a crowded one as a single
            # node and expand it on demand, and `k` so it can pick which one to draw.
            # Keyed on the coordinate rather than on the representative's seq_id, so that the
            # group the renderer assembles and the `sn` count computed here are guaranteed to
            # be counting the same thing even if two representatives land on one point.
            "stack": "%.4f,%.4f" % (r["x"], r["y"]),
            "si": k,
            "sn": total,
        }
        # Flagged rather than quietly drawn like a measurement. 65 % of ABL1's nodes are at
        # an inferred position, and the map gave no way to tell which.
        if approx:
            node["approx"] = True
        if r.get("z") is not None:
            node["z"] = round(r["z"] + radius * math.sin(angle * 1.7), 4)
        nodes.append(node)

    # ---- cluster labels, positioned off the group they name ----------------------------
    clusters = []
    if names:
        by_cluster: dict[int, list] = {}
        for n in nodes:
            by_cluster.setdefault(n["cluster"], []).append(n)
        for ci, group in by_cluster.items():
            if ci not in names:
                continue
            gx = sum(n["x"] for n in group) / len(group)
            gy = sum(n["y"] for n in group) / len(group)
            spread = max(math.hypot(n["x"] - gx, n["y"] - gy) for n in group)
            # Pushed out along the line from the map's centre, so a label sits clear of its
            # own nodes instead of on top of the densest part of them.
            norm = math.hypot(gx, gy) or 1.0
            out = spread + 0.08
            entry = {
                "index": ci, "name": names.get(ci, ""), "count": len(group),
                "x": round(gx + out * (gx / norm), 4),
                "y": round(gy + out * (gy / norm), 4),
            }
            zs_group = [n["z"] for n in group if n.get("z") is not None]
            if zs_group:
                entry["z"] = round(sum(zs_group) / len(zs_group), 4)
            clusters.append(entry)
        clusters.sort(key=lambda c: -c["count"])
        clusters = clusters[:MAX_LABELS]

    # Edges between representatives above the conventional same-fold threshold. This is what
    # the plan asked for all along: "edges between entries above the identity threshold",
    # with TM-score standing in for identity now that there is one.
    tm = embedding.get("tm") or []
    edges = []
    for i in range(len(reps)):
        for j in range(i + 1, len(reps)):
            if i < len(tm) and j < len(tm[i]) and tm[i][j] >= 0.5:
                edges.append({"a": reps[i]["seq_id"], "b": reps[j]["seq_id"]})
    edges = edges[:MAX_EDGES]

    return {
        "nodes": nodes,
        # Edges reference construct ids, not entity ids, so the renderer resolves them
        # through the representative each node came from.
        "edges": [],
        "clusters": clusters,
        # Where the labels came from, so the panel can say so and the Structures heatmap can
        # be moved to the same height rather than the two quietly disagreeing.
        "cluster_cut": CLUSTER_CUT,
        # Past this many entries on one point the renderer draws a single node and expands it
        # on click. Sent rather than hard-coded in the JS so the two cannot drift.
        "stack_max": STACK_MAX,
        "placeholder": False,
        "embedded": True,
        # Whether the map can be turned. A flat family has no third axis worth drawing, and
        # the control says so rather than offering a rotation that does nothing.
        "three_d": bool(zs and len(zs) == len(reps)),
        # Shipped so the Structures panel can draw the matrix and cluster it in the browser,
        # where the cut height can be a live control rather than a decision baked into the
        # artefact. ~36 kB for an 80-representative family.
        "tm": tm,
        "representatives": reps,
        # The structure everything else was superposed onto, so the panel can
        # name it rather than guessing at the first representative.
        "reference": embedding.get("reference"),
        # The AlphaFold model and the transform that puts it on the reference.
        "alphafold": embedding.get("alphafold"),
        "n_representatives": len(reps),
        "n_pairs": embedding.get("n_pairs"),
        "median_tm": embedding.get("median_tm"),
        "approximated": approximated,
    }


def compute(members: list[dict], identity_floor: int = 30,
            embedding: Optional[dict] = None,
            annotations: Optional[dict] = None) -> dict:
    """Position every member, returning nodes, edges and cluster labels.

    Coordinates are in a -1..1 box; the client scales them to whatever the panel is. Doing
    it here rather than in JavaScript means the table, the exports and the map all agree on
    one set of numbers, and a headless dossier export (Phase 4) gets the same picture.
    """
    if not members:
        return {"nodes": [], "edges": [], "clusters": []}

    # A real structural embedding retires everything below. The placeholder exists only
    # because the matrix is expensive; once it has been computed for a family, the map means
    # what section 1.2 of the plan always said it meant.
    if embedding and embedding.get("representatives"):
        return _from_embedding(members, embedding, annotations)

    # ---- cluster by source organism ---------------------------------------------------
    counts: dict[str, int] = {}
    for m in members:
        counts[_short_organism(m.get("organism"))] = counts.get(_short_organism(m.get("organism")), 0) + 1
    ordered = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    named = [name for name, _ in ordered[:MAX_CLUSTERS]]
    # The tail only earns a sector of its own if there is a tail.
    has_other = len(ordered) > MAX_CLUSTERS
    sectors = named + (["Other"] if has_other else [])
    index = {name: i for i, name in enumerate(sectors)}

    # ---- sector widths ------------------------------------------------------------------
    # Proportional to sqrt(count), not to count. A family is typically dominated by one
    # organism (1,423 of carbonic anhydrase II's 1,490 entities are human), and splitting the
    # circle by raw count gives that cluster 95 % of it and leaves the others as slivers too
    # thin to see or label. The square root keeps the dominant cluster visibly dominant while
    # the minor organisms stay legible, which is the thing the sectors exist to show.
    sizes = {name: (counts.get(name, 0) if name != "Other"
                    else sum(c for _, c in ordered[MAX_CLUSTERS:])) for name in sectors}
    weights = [max(1.0, math.sqrt(sizes[name])) for name in sectors]
    total_w = sum(weights)
    bounds = []
    acc = 0.0
    for w in weights:
        start = 2 * math.pi * acc / total_w
        acc += w
        bounds.append((start, 2 * math.pi * acc / total_w))

    # ---- radial axis: identity RANK, not raw identity ----------------------------------
    # A family's identity distribution is savagely skewed: 1,423 of carbonic anhydrase II's
    # 1,490 entities sit above 99 %, with a thin tail running down to 31 %. Mapping identity
    # linearly to radius therefore piles almost the whole family onto one circle and scatters
    # the tail into concentric ripples, which is what the first version of this drew.
    #
    # Ranking instead spends the full radius on whatever spread the family actually has.
    # Reading outward still means "less like the seed", which is the property the axis is
    # there to carry; it is the spacing, not the ordering, that stops being linear. The panel
    # note says so, and Phase 3's real embedding retires the question.
    ordered_ids = sorted(
        members,
        key=lambda m: (-(m.get("identity") if m.get("identity") is not None else -1), m["entity_id"]),
    )
    denom = max(1, len(ordered_ids) - 1)
    rank_of = {m["entity_id"]: i / denom for i, m in enumerate(ordered_ids)}

    # ---- place ------------------------------------------------------------------------
    per_cluster: dict[int, int] = {}
    nodes = []
    for m in members:
        org = _short_organism(m.get("organism"))
        ci = index.get(org, len(sectors) - 1 if has_other else 0)
        n = per_cluster.get(ci, 0)
        per_cluster[ci] = n + 1
        lo, hi = bounds[ci]

        # Radial position: most-like-the-seed near the centre, least-like at the rim.
        # sqrt of the rank, so nodes spread over the disc by AREA rather than by radius:
        # without it the inner rings are sparse and the outer ones jammed, because a ring's
        # circumference grows with its radius.
        r = 0.13 + 0.87 * math.sqrt(rank_of[m["entity_id"]])

        # A small deterministic wobble on top, so members sharing a rank neighbourhood do not
        # line up into visible arcs.
        r = max(0.05, min(1.0, r + 0.04 * (_golden_jitter(n * 3 + 1) - 0.5)))

        # Angular position: filled evenly across the cluster's own sector, inset from its
        # edges so neighbouring clusters stay visually separate. A second irrational
        # multiplier, so angle and radius do not march in step and stripe the sector.
        theta = lo + (hi - lo) * (0.06 + 0.88 * _golden_jitter(n))
        nodes.append({
            "id": m["entity_id"],
            "pdb_id": m["pdb_id"],
            "cluster": ci,
            "cluster_name": sectors[ci],
            "x": round(r * math.cos(theta), 4),
            "y": round(r * math.sin(theta), 4),
        })

    # ---- edges: nearest neighbours within a cluster -----------------------------------
    # Genuinely a placeholder for the identity-threshold edges of the finished map, which
    # need the pairwise matrix Phase 3 computes. Within-cluster nearest neighbours at least
    # draw the real thing an organism group is: a set of repeat solutions of one protein.
    by_cluster: dict[int, list[dict]] = {}
    for nd in nodes:
        by_cluster.setdefault(nd["cluster"], []).append(nd)

    edges: list[dict] = []
    seen: set[tuple[str, str]] = set()
    for group in by_cluster.values():
        # Cap the per-cluster cost: this is O(n^2) in the cluster, and a family can have a
        # thousand entries from one organism.
        pool = group[:200]
        for a in pool:
            # key= on the distance alone: two nodes at an identical distance are common
            # (a family solved repeatedly at one resolution) and tuple ordering would fall
            # through to comparing the node dicts themselves, which raises.
            dists = sorted(
                (((a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2, b) for b in pool if b is not a),
                key=lambda pair: pair[0],
            )
            for _, b in dists[:NEIGHBOURS]:
                key = tuple(sorted((a["id"], b["id"])))
                if key in seen:
                    continue
                seen.add(key)
                edges.append({"a": a["id"], "b": b["id"]})
        if len(edges) >= MAX_EDGES:
            break

    clusters = [
        {"index": i, "name": name, "count": sizes[name],
         # Where to print the label: the middle of the sector, out past the rim.
         "x": round(1.15 * math.cos((bounds[i][0] + bounds[i][1]) / 2), 4),
         "y": round(1.15 * math.sin((bounds[i][0] + bounds[i][1]) / 2), 4)}
        for i, name in enumerate(sectors)
    ]

    return {"nodes": nodes, "edges": edges[:MAX_EDGES], "clusters": clusters,
            "placeholder": True}
