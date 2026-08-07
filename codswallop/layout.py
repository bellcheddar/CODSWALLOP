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


def compute(members: list[dict], identity_floor: int = 30) -> dict:
    """Position every member, returning nodes, edges and cluster labels.

    Coordinates are in a -1..1 box; the client scales them to whatever the panel is. Doing
    it here rather than in JavaScript means the table, the exports and the map all agree on
    one set of numbers, and a headless dossier export (Phase 4) gets the same picture.
    """
    if not members:
        return {"nodes": [], "edges": [], "clusters": []}

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
