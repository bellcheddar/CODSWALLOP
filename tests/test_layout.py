"""The constellation layout. Every case here is one the map got visibly wrong.

The theme running through this file: the map is the one panel a reader believes without
checking, because a picture of a cloud of points looks like a measurement whatever produced
it. Three of the four bugs below drew something with no data in it at all.
"""
import math

from codswallop import layout


def _reps(coords, tm=None):
    return [{"seq_id": "seq%d" % i, "pdb_id": "P%03d" % i, "x": x, "y": y}
            for i, (x, y) in enumerate(coords)]


def _members(counts, seq_ids=None, **extra):
    """`counts[i]` entities all carrying construct `i`."""
    out = []
    for i, n in enumerate(counts):
        for k in range(n):
            out.append({
                "entity_id": "E%d_%d" % (i, k), "pdb_id": "P%03d" % i,
                "seq_id": (seq_ids or {}).get(i, "seq%d" % i),
                "identity": 100.0, **extra,
            })
    return out


def _embedding(coords, tm=None):
    n = len(coords)
    if tm is None:
        tm = [[1.0 if i == j else 0.9 for j in range(n)] for i in range(n)]
    return {"representatives": _reps(coords), "tm": tm, "reference": "P000"}


# ---------------------------------------------------------------------------------------
# The fan.
# ---------------------------------------------------------------------------------------

def test_a_crowded_point_never_fans_further_than_its_nearest_neighbour():
    """The regression this exists for.

    The spread was `0.012 * sqrt(k)` with nothing bounding it, so a construct solved 137
    times drew a disc 0.281 across while the median distance between two genuinely distinct
    representatives in that family was 0.104. The decoration was 2.7x the signal, and two
    nodes on opposite rims of one disc were the same sequence.
    """
    emb = _embedding([(0.0, 0.0), (0.30, 0.0)])
    fam = layout.compute(_members([400, 1]), 30, embedding=emb)
    on_first = [n for n in fam["nodes"] if n["pdb_id"] == "P000"]
    spread = max(math.hypot(n["x"], n["y"]) for n in on_first)
    # The tolerance is the coordinate rounding, not slack in the bound: positions ship at
    # 4 dp, so the bound holds to about 1e-4 and no tighter.
    assert spread <= layout.FAN_FRACTION * 0.30 + 1e-3, "the fan reached its neighbour"
    # And it is genuinely bounded rather than never having been reached: unbounded, 400
    # members would have spread 0.012 * sqrt(399) = 0.24.
    assert spread < 0.24


def test_the_fan_is_not_collapsed_to_nothing_in_a_tight_cluster():
    """Bounding it purely as a fraction of the nearest neighbour makes every member of a
    tight family land on one pixel, which trades a false spread for an unpickable map."""
    emb = _embedding([(0.0, 0.0), (0.002, 0.0)])
    fam = layout.compute(_members([40, 1]), 30, embedding=emb)
    on_first = [n for n in fam["nodes"] if n["pdb_id"] == "P000"]
    assert max(math.hypot(n["x"], n["y"]) for n in on_first) >= layout.FAN_MIN * 0.9


def test_stack_metadata_counts_the_point_not_the_construct():
    emb = _embedding([(0.0, 0.0), (0.5, 0.5)])
    fam = layout.compute(_members([7, 3]), 30, embedding=emb)
    first = [n for n in fam["nodes"] if n["pdb_id"] == "P000"]
    assert {n["sn"] for n in first} == {7}
    assert sorted(n["si"] for n in first) == list(range(7))
    # One node per stack sits exactly on the representative, so a collapsed stack is drawn
    # where the measurement is rather than at a fanned-out offset.
    origin = [n for n in first if n["si"] == 0][0]
    assert (origin["x"], origin["y"]) == (0.0, 0.0)


def test_every_node_on_one_point_shares_a_stack_key():
    emb = _embedding([(0.25, -0.25), (0.5, 0.5)])
    fam = layout.compute(_members([5, 5]), 30, embedding=emb)
    keys = {n["stack"] for n in fam["nodes"] if n["pdb_id"] == "P000"}
    assert len(keys) == 1


# ---------------------------------------------------------------------------------------
# Inferred positions.
# ---------------------------------------------------------------------------------------

def test_a_member_with_no_representative_is_flagged_as_inferred():
    emb = _embedding([(0.0, 0.0), (0.5, 0.0)])
    members = _members([2, 2]) + _members([3], seq_ids={0: "not-a-representative"})
    fam = layout.compute(members, 30, embedding=emb)
    assert fam["approximated"] == 3
    assert sum(1 for n in fam["nodes"] if n.get("approx")) == 3
    assert all(not n.get("approx") for n in fam["nodes"] if n["pdb_id"] in {"P001"})


def test_tied_identities_are_spread_rather_than_piled_onto_one_representative():
    """The fallback placed a member at the representative closest in identity to the seed.

    In a real family hundreds of representatives sit at 100 %, so `min()` broke the tie by
    sort order and every ambiguous member landed on the same one: 126 of ABL1's, 194 of
    spike's. That is most of what made the largest discs, and it looked like a finding.
    """
    emb = _embedding([(0.0, 0.0), (0.4, 0.0), (0.8, 0.0), (1.2, 0.0)])
    members = _members([1, 1, 1, 1]) + _members([40], seq_ids={0: "unknown"})
    fam = layout.compute(members, 30, embedding=emb)
    landed = {}
    for n in fam["nodes"]:
        if n.get("approx"):
            landed[n["stack"]] = landed.get(n["stack"], 0) + 1
    assert len(landed) > 1, "every tied member went to one representative"
    assert max(landed.values()) <= 40 * 0.6


# ---------------------------------------------------------------------------------------
# Clustering and labels.
# ---------------------------------------------------------------------------------------

def test_average_linkage_separates_two_groups_single_linkage_would_chain():
    # Two tight groups, far apart. Every member of a family is the same fold, so a
    # single-linkage tree chains straight through and returns one cluster for everything.
    tm = [[1.0] * 6 for _ in range(6)]
    for i in range(6):
        for j in range(6):
            tm[i][j] = 1.0 if i == j else (0.95 if (i < 3) == (j < 3) else 0.30)
    labels = layout._average_linkage(tm, layout.CLUSTER_CUT)
    assert labels is not None
    assert len({labels[0], labels[1], labels[2]}) == 1
    assert len({labels[3], labels[4], labels[5]}) == 1
    assert labels[0] != labels[3]


def test_a_mixed_cluster_is_named_for_its_shared_domain_not_its_plurality_protein():
    """ABL1's largest cluster is 1,932 entities whose commonest accession is EGFR at 17 %.

    Naming it "Epidermal growth factor receptor" states something false about 83 % of it.
    95 % of it carries the protein tyrosine kinase domain, which is both true and the reason
    those structures cluster together in the first place.
    """
    kinase = [{"id": "PF07714", "name": "Protein tyrosine kinase (PK_Tyr_Ser-Thr)"}]
    # The real shares from ABL1's largest cluster: the top accession is 17 % and the tail is
    # long, so no protein owns it.
    group = []
    for accession, description, n in [("P00533", "Epidermal growth factor receptor", 17),
                                      ("Q06187", "Tyrosine-protein kinase BTK", 13),
                                      ("P08581", "Hepatocyte growth factor receptor", 12),
                                      ("P43405", "Tyrosine-protein kinase SYK", 9),
                                      ("P12931", "Proto-oncogene tyrosine-protein kinase Src", 9)]:
        group += [{"description": description, "uniprot": accession,
                   "pfam": kinase, "seq_length": 300}] * n
    # The remaining 40 %, one structure each, which is what a superfamily search returns.
    group += [{"description": "Kinase %d" % i, "uniprot": "X%05d" % i,
               "pfam": kinase, "seq_length": 300} for i in range(40)]
    assert layout._label_clusters([group]) == ["Protein tyrosine kinase"]


def test_domains_are_read_through_the_annotation_lookup_after_compaction():
    """The bug this file exists to stop happening twice.

    `_compact` lifts every member's Pfam list into one family-level lookup and leaves an
    `annot_id`, and the map is built after that has run. So the domain rule read `m["pfam"]`,
    found nothing on all 2,000 members of every real request, and fell through to naming the
    cluster after its plurality protein: the exact mislabelling it was written to prevent.
    The unit tests all passed, because they built their members by hand with `pfam` on them.

    So this one uses the shape the layout is actually handed, and asserts the label a live
    request produces.
    """
    annotations = {"a1": {"pfam": [{"id": "PF07714",
                                    "name": "Protein tyrosine kinase (PK_Tyr_Ser-Thr)"}],
                          "interpro": []}}
    group = [{"description": "Epidermal growth factor receptor", "uniprot": "P%05d" % i,
              "annot_id": "a1", "seq_length": 300} for i in range(30)]
    group += [{"description": "Kinase %d" % i, "uniprot": "Q%05d" % i,
               "annot_id": "a1", "seq_length": 300} for i in range(70)]

    assert layout._label_clusters([group]) == ["Epidermal growth factor receptor"], \
        "without the lookup the plurality name is all that is available"
    assert layout._label_clusters([group], annotations) == ["Protein tyrosine kinase"]


def test_a_cluster_one_protein_owns_is_named_for_that_protein():
    group = [{"description": "Tyrosine-protein kinase ABL1", "uniprot": "P00519",
              "pfam": [{"id": "PF00018", "name": "SH3 domain (SH3_1)"}],
              "seq_length": 62}] * 9
    group += [{"description": "Something else", "uniprot": "Q99999", "pfam": [],
               "seq_length": 62}]
    assert layout._label_clusters([group]) == ["Tyrosine-protein kinase ABL1"]


def test_the_deposited_description_does_not_split_a_protein_s_own_vote():
    """One accession, spelled three ways, against a contaminant spelled one way. Counting
    free text lets the contaminant win a cluster it holds a third of."""
    group = (
        [{"description": "Tyrosine-protein kinase ABL1", "uniprot": "P00519", "pfam": [],
          "seq_length": 290}] * 12
        + [{"description": "Proto-oncogene tyrosine-protein kinase ABL1",
            "uniprot": "P00519", "pfam": [], "seq_length": 290}] * 11
        + [{"description": "ABL TYROSINE KINASE", "uniprot": "P00519", "pfam": [],
            "seq_length": 290}] * 10
        + [{"description": "Insulin receptor", "uniprot": "P06213", "pfam": [],
            "seq_length": 290}] * 15
    )
    assert layout._label_clusters([group])[0].startswith("Tyrosine-protein kinase ABL1")


def test_two_clusters_with_one_name_are_separated_by_construct_length():
    def group(length):
        return [{"description": "Tyrosine-protein kinase ABL1", "uniprot": "P00519",
                 "pfam": [], "seq_length": length}] * 10
    out = layout._label_clusters([group(290), group(62)])
    assert out[0] != out[1]
    assert "290 aa" in out[0] and "62 aa" in out[1]


def test_a_name_that_still_duplicates_is_dropped_rather_than_printed_twice():
    """Alpha-synuclein's clusters are six conformations of one 140-residue protein, so the
    length cannot separate them either and the field printed one label six times."""
    def group():
        return [{"description": "Alpha-synuclein", "uniprot": "P37840", "pfam": [],
                 "seq_length": 140}] * 10
    out = layout._label_clusters([group(), group(), group()])
    assert out[0].startswith("Alpha-synuclein")
    assert out[1] == "" and out[2] == "", "the same name printed more than once"


def test_a_truncated_label_does_not_end_on_a_conjunction():
    assert not layout._shorten(
        "Protein tyrosine and serine/threonine kinase", 34).endswith("and…")


def test_domain_names_drop_the_pfam_short_code():
    assert layout._domain_name("Protein kinase domain (Pkinase)") == "Protein kinase domain"
    assert layout._domain_name("Bromodomain") == "Bromodomain"


def test_a_cluster_label_counts_every_entity_it_covers():
    """The count printed next to a name is the whole group, not the representatives in it.

    Only labelled clusters appear in `clusters` (a name that duplicated a larger one's is
    dropped), so the total across labels need not be the whole family; what must hold is that
    each label's number is exactly the nodes it names.
    """
    tm = [[1.0 if i == j else (0.95 if (i < 2) == (j < 2) else 0.2) for j in range(4)]
          for i in range(4)]
    emb = _embedding([(0.0, 0.0), (0.05, 0.0), (1.0, 0.0), (1.05, 0.0)], tm=tm)
    members = (_members([30, 30], description="Kinase A", uniprot="A", pfam=[])
               + [dict(m, description="Kinase B", uniprot="B")
                  for m in _members([0, 0, 20, 20])])
    fam = layout.compute(members, 30, embedding=emb)
    assert fam["clusters"], "a two-group matrix produced no labels at all"
    for c in fam["clusters"]:
        covered = sum(1 for n in fam["nodes"] if n["cluster"] == c["index"])
        assert c["count"] == covered


def test_the_placeholder_layout_is_untouched_by_any_of_this():
    members = [{"entity_id": "E%d" % i, "pdb_id": "P%d" % i, "identity": 90 - i,
                "organism": "Homo sapiens" if i % 2 else "Escherichia coli"}
               for i in range(20)]
    fam = layout.compute(members, 30)
    assert fam["placeholder"] is True
    assert fam["clusters"] and all("name" in c for c in fam["clusters"])
    assert all("stack" not in n for n in fam["nodes"])
