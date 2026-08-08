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
