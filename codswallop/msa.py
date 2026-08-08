"""Family alignment: per-column conservation and the sequence logo.

A **star alignment** against the seed, not a progressive multiple alignment. Every distinct
deposited sequence is aligned pairwise to the seed and its residues are read off into the
seed's coordinate frame. That is a real and well-understood approximation, and it is the
right one here for three reasons:

* Every other panel already speaks seed coordinates, so a column here is the same column as
  in the coverage census, the domain ribbon and the disorder profile. A progressive MSA would
  introduce its own gapped frame that nothing else shares.
* The pairwise aligner is already loaded and each alignment is ~2 ms, so a 794-construct
  family costs under two seconds and needs no external binary on the droplet.
* Insertions relative to the seed genuinely have nowhere to go in a seed-framed view. They
  are counted and reported per column rather than silently dropped, so a column with many
  insertions is visible as such.

What it is not: a phylogenetically weighted alignment. Counts are weighted by how many
entities used each construct, so conservation reflects what was deposited, not what evolved.
The panel says so.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Optional

from .constructs import _ALIGNER

PARSE_VERSION = 6

# See the calibration note on `conserved` below before changing this.
CONSERVED_CUT = 0.97

# The 20, plus the gap. Anything else (X, B, Z, U) is counted as unknown and excluded from
# the entropy, rather than being treated as a 21st residue that is conserved wherever it
# happens to appear.
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_SET = set(AA)

# Kyte-Doolittle-ish grouping, used only to colour the logo. Chemistry, not conservation.
GROUPS = {
    **{c: "hydrophobic" for c in "AVLIMFWC"},
    **{c: "polar" for c in "STNQGY"},
    **{c: "basic" for c in "KRH"},
    **{c: "acidic" for c in "DE"},
    **{c: "special" for c in "P"},
}


def _engineered_position(col: dict, max_depth: int) -> Optional[dict]:
    """A position somebody deliberately mutated, or None.

    Keyed on the fraction carrying something OTHER than the wild-type residue, and
    deliberately NOT on whether the wild type is the most common residue. Requiring the wild
    type to top the column excludes exactly the positions worth surfacing: for a heavily
    studied hotspot the mutant can outnumber it. p53's R273 and carbonic anhydrase II's H64
    both vanished from the list for precisely that reason, and they are the two positions a
    structural biologist would name first in those families.

    The depth floor keeps out the N-terminal columns, where truncated constructs starting at
    different points collect a spread of residues that is cloning noise, not design.
    """
    top = col.get("top") or []
    if not top or col["depth"] < 0.5 * max_depth:
        return None
    seed_res = col["seed"]
    seed_frac = next((t["f"] for t in top if t["aa"] == seed_res), 0.0)
    others = [t for t in top if t["aa"] != seed_res and t["f"] >= 0.005]
    # The wild type must still be a real presence (this is a mutant series, not a different
    # protein) and something else must be a real presence too.
    if not (0.3 <= seed_frac <= 0.995) or not others:
        return None
    return {
        "pos": col["pos"], "seed": seed_res, "conservation": col["conservation"],
        "depth": col["depth"], "wild_type_fraction": round(seed_frac, 3),
        "substituted": round(1 - seed_frac, 3), "variants": others,
    }


def build(seed: str, sequences: dict[str, str], weights: dict[str, int],
          max_columns: int = 4000) -> Optional[dict]:
    """Align every distinct construct to the seed and summarise each column.

    `sequences` maps an id to a sequence; `weights` maps the same id to how many entities
    used it. Returns None when there is nothing to align.
    """
    if not seed or not sequences:
        return None
    n = len(seed)
    if n > max_columns:
        # A seed longer than this is a multi-domain giant where a per-residue logo is not
        # readable anyway, and the payload would be larger than the rest of the family.
        return {"length": n, "too_long": True, "columns": []}

    counts: list[Counter] = [Counter() for _ in range(n + 1)]
    insertions = [0] * (n + 1)
    total_weight = 0

    for sid, seq in sequences.items():
        if not seq:
            continue
        w = max(1, weights.get(sid, 1))
        total_weight += w
        try:
            aln = _ALIGNER.align(seed, seq)[0]
        except Exception:
            continue
        blocks = list(zip(aln.aligned[0], aln.aligned[1]))
        for (cs, ce), (qs, qe) in blocks:
            for off in range(int(ce) - int(cs)):
                counts[int(cs) + off + 1][seq[int(qs) + off]] += w
        # Residues the construct has that the seed does not: counted against the column they
        # follow, so a heavily-inserted position is visible rather than invisible.
        for i in range(len(blocks) - 1):
            (_, c_end), (_, q_end) = blocks[i]
            (c_start, _), (q_start, _) = blocks[i + 1]
            gap = int(q_start) - int(q_end)
            if gap > 0 and 1 <= int(c_end) <= n:
                insertions[int(c_end)] += w

    columns = []
    for pos in range(1, n + 1):
        c = counts[pos]
        observed = sum(v for k, v in c.items() if k in AA_SET)
        if not observed:
            columns.append({"pos": pos, "seed": seed[pos - 1], "depth": 0,
                            "conservation": None, "top": [], "insertions": insertions[pos]})
            continue
        # Shannon entropy over the 20 residues, normalised so 1.0 is a single residue
        # everywhere and 0.0 is a uniform mixture of all twenty.
        h = 0.0
        top = []
        for res, k in c.most_common():
            if res not in AA_SET:
                continue
            p = k / observed
            h -= p * math.log(p, 2)
            if len(top) < 4:
                top.append({"aa": res, "f": round(p, 3), "group": GROUPS.get(res, "other")})
        conservation = round(1 - h / math.log(20, 2), 3)
        columns.append({
            "pos": pos, "seed": seed[pos - 1], "depth": observed,
            "conservation": conservation, "top": top, "insertions": insertions[pos],
        })

    scored = [c["conservation"] for c in columns if c["conservation"] is not None]
    max_depth = max((c["depth"] for c in columns), default=0)
    return {
        "length": n,
        "too_long": False,
        "n_sequences": len(sequences),
        "total_weight": total_weight,
        "columns": columns,
        "mean_conservation": round(sum(scored) / len(scored), 3) if scored else None,
        # The positions a reader should look at first. The cut is well below 1.0, and that is
        # the whole point: a family's most functionally important residues are almost never
        # perfectly invariant, because those are exactly the positions somebody deliberately
        # mutated.
        #
        # Calibrated against four residues whose importance is not in question. Normalised
        # entropy, not top-residue frequency, which is a different and higher number:
        #
        #   carbonic anhydrase II  H94  0.979   (zinc-coordinating; the archive holds H94A)
        #                          H96  0.991
        #   p53                    R175 0.991   (DNA-contacting; R175H is a cancer hotspot)
        #                          R248 0.982
        #
        # 0.999 found none of them. 0.99 found half. 0.97 finds all four and still excludes
        # the bulk of the protein, whose mean sits near 0.88.
        "conserved": [c["pos"] for c in columns
                      if c["conservation"] is not None and c["conservation"] >= CONSERVED_CUT][:300],
        # Near-invariant but not quite: conserved everywhere except where it was engineered.
        # For a well-studied family this is the most informative list on the panel.
        "conserved_with_exceptions": [
            {"pos": c["pos"], "seed": c["seed"], "conservation": c["conservation"],
             "variants": [t for t in c["top"][1:] if t["f"] > 0]}
            for c in columns
            if c["conservation"] is not None and CONSERVED_CUT <= c["conservation"] < 1.0
            and len(c["top"]) > 1
        ][:60],
        # And the opposite: the positions that vary most.
        "variable": [c["pos"] for c in sorted(
            (c for c in columns if c["conservation"] is not None),
            key=lambda c: c["conservation"])[:40]],
        # The list this panel exists for, and the one neither of the two above finds.
        #
        # A position where the wild-type residue still dominates but a real minority of
        # constructs carry something else is a position somebody deliberately mutated. Those
        # are the most interesting residues in any well-studied family, and they sit in the
        # gap between "conserved" and "most variable": carbonic anhydrase II's proton-shuttle
        # H64 scores 0.943 and p53's R273 scores 0.875, so the conserved cut excludes both,
        # while neither is variable enough to reach the top-40 variable list. The variants
        # they carry are exactly the ones you would name from memory: H64 to A and Y, R273
        # to H and C.
        # Ranked by how much of the family carries a substitution, not by how many DIFFERENT
        # substitutions appear, and only over columns most constructs actually contain.
        # Both conditions are load-bearing. Ranking by variant count put M1, S2, H3 and H4 at
        # the top of carbonic anhydrase II's list: N-terminal columns collect many different
        # residues because truncated constructs start at different points, which is cloning
        # noise rather than a designed mutation. Requiring real depth removes them, because
        # those same columns are present in only a fraction of the family.
        "engineered": sorted(
            filter(None, (_engineered_position(c, max_depth) for c in columns)),
            key=lambda x: -x["substituted"])[:150],   # 60 was arbitrary and cut p53 R273 (0.084) just below the line
    }
