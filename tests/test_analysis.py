"""Ligand classification, crystallisation parsing, the coverage census and the map layout."""
import pytest
from codswallop import crystals, family, layout, ligands


# ---- ligand classification ------------------------------------------------------------
@pytest.mark.parametrize("comp,expected", [
    ("HOH", "water"), ("ZN", "ion"), ("SO4", "ion"),
    ("GOL", "cryoprotectant"), ("EDO", "cryoprotectant"), ("PEG", "cryoprotectant"),
    ("TRS", "buffer"), ("EPE", "buffer"),
    ("BCN", "buffer"),          # bicine: read as a ligand until carbonic anhydrase II showed it
    ("HEM", "cofactor"), ("ATP", "cofactor"), ("NAD", "cofactor"),
    ("OLC", "lipid/detergent"),
    ("STI", "ligand"),          # imatinib: an actual designed ligand
])
def test_component_classes(comp, expected):
    assert ligands.classify(comp) == expected


def test_only_designed_components_count_as_ligand_bound():
    """The whole point of the classification: without it, "ligand-bound" means "was frozen"."""
    assert "ligand" in ligands.COUNTS_AS_BOUND
    assert "cofactor" in ligands.COUNTS_AS_BOUND
    for junk in ("cryoprotectant", "buffer", "ion", "water", "solvent"):
        assert junk not in ligands.COUNTS_AS_BOUND


# ---- crystallisation parsing ----------------------------------------------------------
def test_parses_a_real_lipidic_cubic_phase_condition():
    p = crystals.parse({
        "method": "LIPIDIC CUBIC PHASE", "ph": 6.5, "temp_k": 293.0,
        "details": "100 mM MES pH 6.2-6.7, 40-100 mM ammonium phosphate dibasic, "
                   "18-24% PEG400, LIPIDIC CUBIC PHASE, temperature 293K",
    })
    assert p["parsed"]
    assert "PEG" in p["precipitants"]
    assert "Ammonium phosphate" in p["precipitants"]
    assert "MES" in p["buffers"]
    assert p["method"] == "Lipidic cubic phase"
    assert p["temp_c"] == 19.9


def test_sized_peg_does_not_also_report_unspecified_peg():
    p = crystals.parse({"details": "20% PEG 3350, 0.1M Tris pH 8.5", "ph": 8.5, "temp_k": None})
    assert p["precipitants"].count("PEG") == 1
    assert "PEG (unspecified)" not in p["precipitants"]


def test_bis_tris_propane_is_not_also_reported_as_bis_tris():
    p = crystals.parse({"details": "0.1 M Bis-Tris propane pH 7.0", "ph": 7.0, "temp_k": None})
    assert p["buffers"] == ["Bis-Tris propane"]


def test_unrecognised_condition_is_flagged_not_silently_empty():
    p = crystals.parse({"details": "grown from an in-house screen", "ph": None, "temp_k": None})
    assert p["parsed"] is False
    assert p["details"], "the verbatim text must survive so a reader can see it"


# ---- the coverage census --------------------------------------------------------------
def _member(beg, end, unobserved=None):
    return {"query_beg": beg, "query_end": end, "unobserved_seed": unobserved}


def test_census_reports_a_ratio_not_a_saturating_absolute():
    """One full-length construct makes every "never" statistic read zero.

    This is the regression that mattered: with 267 p53 constructs something resolved almost
    every residue at least once, so "resolved in no construct" read 0 % for a protein a third
    of which is disordered. The ratio has to discriminate where the absolute cannot.
    """
    fam = {"seed_length": 100}
    members = [_member(1, 100, [(50, 60)]) for _ in range(19)]
    members.append(_member(1, 100, []))          # the one entry that resolved the loop
    c = family.coverage_census(fam, members)

    assert c["never_seen"] == 0, "the absolute saturates, which is why it is not the headline"
    assert c["rarely_resolved"] == 11, "the ratio still finds the eleven-residue loop"
    assert c["disorder_runs"][0]["start"] == 50
    assert c["disorder_runs"][0]["end"] == 60


def test_census_says_not_measured_rather_than_zero_when_density_is_absent():
    c = family.coverage_census({"seed_length": 50}, [_member(1, 50)])
    assert c["seen"] is None
    assert c["rarely_resolved_pct"] is None, "absent data must not render as perfect order"


# ---- map layout -----------------------------------------------------------------------
def _node(i, identity, organism):
    return {"entity_id": f"E{i}", "pdb_id": f"P{i:03d}", "identity": identity,
            "organism": organism, "resolution": 2.0}


def test_dominant_cluster_gets_more_room_but_not_all_of_it():
    members = [_node(i, 99.9, "Homo sapiens") for i in range(400)]
    members += [_node(1000 + i, 60.0, "Mus musculus") for i in range(4)]
    out = layout.compute(members, 30)
    spans = {c["name"]: c["count"] for c in out["clusters"]}
    assert spans["H. sapiens"] == 400 and spans["M. musculus"] == 4
    # sqrt weighting: the big cluster dominates without crushing the small one to a sliver.
    assert len(out["nodes"]) == 404


def test_near_identical_family_still_spreads_across_the_radius():
    """Rank, not raw identity. A linear map piles a 95 %-identical family onto one circle."""
    members = [_node(i, 99.9 if i else 100.0, "Homo sapiens") for i in range(200)]
    out = layout.compute(members, 30)
    radii = sorted((n["x"] ** 2 + n["y"] ** 2) ** 0.5 for n in out["nodes"])
    assert radii[-1] - radii[0] > 0.5, "the radial axis must use its range, not collapse"


def test_layout_is_deterministic():
    members = [_node(i, 100 - i * 0.1, "Homo sapiens") for i in range(50)]
    a = layout.compute(members, 30)
    b = layout.compute(members, 30)
    assert a["nodes"] == b["nodes"], "a map that reshuffles on reload is a map nobody trusts"


def test_plain_tris_is_still_found_when_it_is_genuinely_there():
    p = crystals.parse({"details": "0.1 M Tris pH 8.0, 2 M ammonium sulfate",
                        "ph": 8.0, "temp_k": None})
    assert "Tris" in p["buffers"]


def test_both_buffers_survive_when_a_condition_names_both():
    p = crystals.parse({"details": "0.1 M Bis-Tris pH 6.5 and 50 mM Tris pH 7.5",
                        "ph": 6.5, "temp_k": None})
    assert set(p["buffers"]) >= {"Bis-Tris", "Tris"}


# ---- conservation ----------------------------------------------------------------------
def test_conservation_finds_a_position_that_is_conserved_except_where_engineered():
    """The cut has to be below 1.0, or it finds none of the interesting positions.

    A family's functionally important residues are almost never perfectly invariant: those
    are precisely the positions somebody made a mutant of. Carbonic anhydrase II's
    zinc-coordinating H94 sits at 0.991 because the archive holds H94A knockouts.
    """
    from codswallop import msa
    seed = "MKHAWGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKYDPSL"
    # Weighted the way a real family is: the wild type is most of the archive and the
    # knockout is a handful of entries. Carbonic anhydrase II's H94 sits at 0.991 with about
    # 0.9 % variation, so the fixture has to be that lopsided for the threshold to mean
    # anything: an even split of 50 to 1 is 2 % and lands below the cut.
    seqs = {f"s{i}": seed for i in range(20)}
    weights = {k: 60 for k in seqs}
    seqs["mut"] = seed[:2] + "A" + seed[3:]      # one deliberate knockout at position 3
    weights["mut"] = 8                            # ~0.7 % of the family, like a real knockout

    out = msa.build(seed, seqs, weights)
    col = out["columns"][2]
    assert col["seed"] == "H"
    assert msa.CONSERVED_CUT <= col["conservation"] < 1.0
    assert 3 in out["conserved"]
    flagged = [c["pos"] for c in out["conserved_with_exceptions"]]
    assert 3 in flagged, "conserved-except-where-engineered is the informative list"


def test_conservation_is_weighted_by_how_many_entities_used_each_construct():
    from codswallop import msa
    seed = "ACDEFGHIKLMNPQRSTVWY"
    out = msa.build(seed, {"a": seed, "b": "W" + seed[1:]}, {"a": 99, "b": 1})
    col = out["columns"][0]
    assert col["top"][0]["aa"] == "A" and col["top"][0]["f"] > 0.9


def test_deliberately_mutated_positions_get_their_own_list():
    """The gap between "conserved" and "most variable" is where the interesting residues are.

    Carbonic anhydrase II's proton-shuttle H64 scores 0.943 and p53's R273 scores 0.875:
    below the conserved cut, but nowhere near variable enough to reach the top-40 variable
    list. Neither list found them, and they are the two positions in those families that a
    structural biologist would name first.
    """
    from codswallop import msa
    seed = "MKHAWGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKYDPSL"
    seqs = {f"wt{i}": seed for i in range(10)}
    weights = {k: 10 for k in seqs}
    for i, sub in enumerate("AY"):               # the classic knockout pair
        seqs[f"m{i}"] = seed[:2] + sub + seed[3:]
        weights[f"m{i}"] = 6
    # A truncated construct, so the N-terminal columns are shallow the way they really are.
    seqs["trunc"] = seed[8:]
    weights["trunc"] = 40

    out = msa.build(seed, seqs, weights)
    hits = {x["pos"]: x for x in out["engineered"]}
    assert 3 in hits, "a position with real minority substitution is an engineered position"
    assert {v["aa"] for v in hits[3]["variants"]} == {"A", "Y"}
    assert 3 not in out["conserved"], "and it is deliberately NOT in the conserved list"


# ---- the structural embedding -----------------------------------------------------------
def test_embedding_positions_members_by_their_construct():
    """Members of the same construct land on the same point; the rest are approximated.

    The regression: the map used to be computed before _compact put seq_id on members, so
    nothing matched a representative and all 1,688 lysozyme entries were reported as
    approximated while looking like a measurement.
    """
    emb = {
        "representatives": [
            {"seq_id": "a", "pdb_id": "1AAA", "x": -0.8, "y": 0.0},
            {"seq_id": "b", "pdb_id": "1BBB", "x": 0.8, "y": 0.0},
        ],
        "tm": [[1.0, 0.9], [0.9, 1.0]], "n_pairs": 1, "median_tm": 0.9,
    }
    members = [
        {"entity_id": "E1", "pdb_id": "1AAA", "seq_id": "a", "identity": 100.0},
        {"entity_id": "E2", "pdb_id": "1AAC", "seq_id": "a", "identity": 100.0},
        {"entity_id": "E3", "pdb_id": "1BBB", "seq_id": "b", "identity": 60.0},
        {"entity_id": "E4", "pdb_id": "1ZZZ", "seq_id": "z", "identity": 61.0},
    ]
    out = layout.compute(members, 30, embedding=emb)
    assert out["embedded"] is True
    assert out["approximated"] == 1, "only the unrepresented construct is approximated"
    by_id = {n["id"]: n for n in out["nodes"]}
    # The two entries sharing construct "a" sit together, well away from "b".
    assert abs(by_id["E1"]["x"] - by_id["E2"]["x"]) < 0.1
    assert by_id["E1"]["x"] < 0 and by_id["E3"]["x"] > 0
    # The unrepresented one inherits from the nearest representative by identity.
    assert by_id["E4"]["x"] > 0


def test_no_embedding_falls_back_to_the_placeholder():
    members = [{"entity_id": f"E{i}", "pdb_id": f"P{i:03d}", "identity": 90.0 - i,
                "organism": "Homo sapiens", "resolution": 2.0} for i in range(10)]
    out = layout.compute(members, 30, embedding=None)
    assert out.get("embedded") is not True
    assert out.get("placeholder") is True


def test_mds_puts_dissimilar_structures_further_apart():
    import numpy as np
    from codswallop.embed import embed as mds
    # Two tight pairs, far from each other.
    tm = np.array([[1.0, 0.95, 0.2, 0.2],
                   [0.95, 1.0, 0.2, 0.2],
                   [0.2, 0.2, 1.0, 0.95],
                   [0.2, 0.2, 0.95, 1.0]])
    c = mds(tm)
    within = np.linalg.norm(c[0] - c[1])
    between = np.linalg.norm(c[0] - c[2])
    assert between > within * 3, "the embedding must separate the two folds"


def test_artefact_versions_have_exactly_one_definition():
    """A version constant declared in two places will eventually disagree.

    It already did: bumping embed.VERSION to 2 for the superposition transforms while
    embed_io still said 1 made every previously computed artefact fail its version check, so
    three families silently reverted to the placeholder map with no error anywhere.
    """
    import inspect
    from codswallop import contacts_io, embed, embed_io
    assert embed.VERSION is embed_io.VERSION
    for mod in (embed, contacts_io):
        src = inspect.getsource(mod)
        assert src.count("\nVERSION = ") <= 1, f"{mod.__name__} redeclares VERSION"
