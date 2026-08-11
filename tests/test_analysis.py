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


def test_artefact_survey_flags_a_wrong_version_as_stale(tmp_path, monkeypatch):
    """Present-but-old must count as stale, not as present."""
    import json, threading
    from codswallop import artefacts, config, db, embed_io

    monkeypatch.setattr(config, "DB_PATH", tmp_path / "t.db")
    monkeypatch.setattr(config, "DATA_DIR", tmp_path / "data")
    monkeypatch.setattr(embed_io, "EMBED_DIR", tmp_path / "data" / "embeddings")
    db._local = threading.local()
    db.init()
    db.connect().execute(
        "INSERT INTO family(slug, query, kind, identity_threshold, built_at) "
        "VALUES ('fam', 'P00918', 'uniprot', 30, 0)")
    db.connect().commit()

    (tmp_path / "data" / "embeddings").mkdir(parents=True)
    (tmp_path / "data" / "embeddings" / "fam.json").write_text(
        json.dumps({"version": embed_io.VERSION - 1}))

    st = artefacts.status("fam")
    assert st["embedding"]["present"] is True
    assert st["embedding"]["current"] is False, "an old version must not read as present"
    assert [s["slug"] for s in artefacts.stale("embedding")] == ["fam"]


def test_a_metal_seen_once_is_not_promoted_to_a_cofactor():
    """A heavy-atom derivative is not a catalytic centre.

    Mercury was promoted in carbonic anhydrase II on the strength of a single entry: soaked
    in for phasing and coordinated by whatever was nearby. A catalytic metal is present
    every time somebody solves the protein.
    """
    comps = [{"id": "HG", "klass": "ion"}, {"id": "ZN", "klass": "ion"}]
    contacts = {"metal_coordination": {
        "HG": {"entries": 1, "residues": [{"pos": p, "restype": "CYS", "n": 1}
                                          for p in (10, 20, 30)]},
        "ZN": {"entries": 40, "residues": [{"pos": p, "restype": "HIS", "n": 30}
                                           for p in (94, 96, 119)]},
    }}
    conservation = {10: 0.99, 20: 0.99, 30: 0.99, 94: 0.98, 96: 0.99, 119: 0.99}
    out = ligands.promote_catalytic_metals(comps, contacts, conservation)
    by = {c["id"]: c for c in out}
    assert by["ZN"]["klass"] == "cofactor"
    assert by["HG"]["klass"] == "ion", "one entry is not evidence"


def test_water_does_not_count_as_a_coordinating_residue():
    comps = [{"id": "ZN", "klass": "ion"}]
    contacts = {"metal_coordination": {"ZN": {
        "entries": 40,
        "residues": [{"pos": 94, "restype": "HIS", "n": 30},
                     {"pos": 421, "restype": "HOH", "n": 30},
                     {"pos": 422, "restype": "HOH", "n": 30}],
    }}}
    out = ligands.promote_catalytic_metals(comps, contacts, {94: 0.99, 421: 0.99, 422: 0.99})
    assert out[0]["klass"] == "ion", "two waters and a histidine is not a metal site"


def test_alphafold_span_is_only_applied_in_the_canonical_frame():
    """Seed coordinates are not canonical coordinates unless the seed IS the canonical.

    Lysozyme's family is seeded from a PDB entity: its seed is the 129-residue mature
    protein, while P00698's canonical is 147 residues including an 18-residue signal
    peptide. Slicing the AlphaFold model at seed [1,129] therefore takes the signal peptide
    and the wrong end, and dropped the superposition from TM 0.993 to 0.854.
    """
    import inspect
    from codswallop import embed
    src = inspect.getsource(embed.build)
    assert "seed_is_canonical" in src
    # The guard must require BOTH that the family was seeded from a UniProt accession and
    # that it is the same accession the model belongs to.
    assert 'fam.get("kind") == "uniprot"' in src
    assert "== accession.upper()" in src


# Families are assembled at 30% identity, so the family of a well-studied protein is really
# its superfamily, and the subject is routinely outvoted inside its own family. Both of these
# shipped a confident superposition of the wrong protein onto the page.
def test_reference_prefers_the_protein_that_was_searched_for():
    from codswallop import embed_io  # noqa: F401  (keeps the droplet import boundary honest)
    from codswallop.embed import reference_index
    reps = [{"pdb_id": "4L7S", "uniprot": "P00533"},   # EGFR: most central of the kinases
            {"pdb_id": "2HYY", "uniprot": "P00519"},   # ABL1: the subject, less central
            {"pdb_id": "3POZ", "uniprot": "P00533"}]
    centrality = [9.0, 4.0, 8.0]
    assert reference_index(reps, centrality, "P00519") == 1
    # Among the subject's own structures, centrality still decides.
    reps.append({"pdb_id": "1IEP", "uniprot": "P00519"})
    assert reference_index(reps, centrality + [6.0], "P00519") == 3


def test_reference_falls_back_to_centrality_when_the_seed_is_not_an_accession():
    from codswallop.embed import reference_index
    reps = [{"pdb_id": "1AKI", "uniprot": "P00698"}, {"pdb_id": "2LZM", "uniprot": "P00720"}]
    assert reference_index(reps, [3.0, 7.0], "") == 1
    # A UniProt seed with no structure of its own among the representatives must not crash.
    assert reference_index(reps, [3.0, 7.0], "Q99999") == 1


def test_representatives_reserve_half_the_map_for_the_subject():
    from codswallop.embed import choose_representatives
    # A superfamily where the subject is heavily outnumbered, as ABL1 is by EGFR.
    cs = ([{"uniprot": "P00533", "n_entities": 100 - i} for i in range(50)] +
          [{"uniprot": "P00519", "n_entities": 5 - (i % 5)} for i in range(20)])
    got = choose_representatives(cs, "P00519", 10)
    assert len(got) == 10
    assert sum(1 for c in got if c["uniprot"] == "P00519") == 5, \
        "the protein the reader searched for must be on its own map"


def test_a_dominant_seed_keeps_every_slot_it_would_have_had():
    from codswallop.embed import choose_representatives
    cs = ([{"uniprot": "P00698", "n_entities": 100 - i} for i in range(20)] +
          [{"uniprot": "P00720", "n_entities": 1}])
    got = choose_representatives(cs, "P00698", 10)
    assert [c["uniprot"] for c in got] == ["P00698"] * 10
    assert [c["n_entities"] for c in got] == list(range(100, 90, -1))


def test_representatives_fall_back_when_the_seed_has_no_constructs():
    from codswallop.embed import choose_representatives
    cs = [{"uniprot": "P00533", "n_entities": 9}, {"uniprot": "P00533", "n_entities": 4}]
    assert len(choose_representatives(cs, "Q99999", 5)) == 2
    assert len(choose_representatives(cs, "", 5)) == 2


def test_representative_accession_agrees_with_the_construct_table():
    """The two halves of the fix must count the subject the same way.

    choose_representatives reads the construct's family-assigned accession while the
    representative record was storing the member's first UniProt cross-reference. RCSB's
    ordering is not meaningful, so A2A's own reference was picked from 12 candidates rather
    than 36: a quota that selects on one definition and a reference that filters on another
    silently narrows to their intersection.
    """
    import inspect
    from codswallop import embed
    src = inspect.getsource(embed.build)
    marker = 'reps.append('
    block = src[src.index(marker):src.index(marker) + 400]
    assert 'c.get("uniprot")' in block, \
        "representatives must carry the construct's accession, not the member's first xref"


# ---- assembly and oligomeric state ---------------------------------------------------
def _entry(pdb, count, details, prov, alts=None, area=None, ifres=None):
    return {"pdb_id": pdb, "assembly": {
        "count": count, "details": details, "provenance": prov, "method": "PISA",
        "buried_area": area, "interface_residues": ifres,
        "n_assemblies": 1 + len(alts or []), "ambiguous": bool(alts),
        "alternatives": alts or []}}


def test_oligomeric_states_group_on_the_count_not_the_wording():
    """`oligomeric_details` is free text whose capitalisation is not consistent across the
    archive: lysozyme reported "trimeric" on 64 entries and "Trimeric" on 3, which listed
    one oligomeric state twice in the same table as though they were different."""
    from codswallop.family import build_assemblies
    a = build_assemblies([_entry("1AAA", 3, "trimeric", "author_defined_assembly"),
                          _entry("2BBB", 3, "Trimeric", "author_defined_assembly"),
                          _entry("3CCC", 3, "  TRIMERIC ", "author_defined_assembly")])
    assert [s["count"] for s in a["states"]] == [3]
    assert a["states"][0]["entries"] == 3
    assert a["states"][0]["details"] == "trimeric"


def test_provenance_is_three_counts_and_never_an_agreement_rate():
    """author_defined means PISA returned nothing or never ran, not that it disagreed.
    Scoring it as disagreement would invent a conflict on 758 of 1,686 lysozyme entries."""
    from codswallop.family import build_assemblies
    a = build_assemblies([
        _entry("1AAA", 1, "monomeric", "author_and_software_defined_assembly"),
        _entry("2BBB", 1, "monomeric", "author_defined_assembly"),
        _entry("3CCC", 2, "dimeric", "software_defined_assembly"),
    ])
    assert a["provenance"] == {"both": 1, "author": 1, "software": 1}
    assert "disagree" not in str(a).lower()


def test_only_an_entry_that_contradicts_itself_counts_as_ambiguous():
    from codswallop.family import build_assemblies
    a = build_assemblies([
        _entry("1AAA", 1, "monomeric", "author_defined_assembly"),
        _entry("2BBB", 2, "dimeric", "author_defined_assembly", alts=[1]),
    ])
    assert a["n_ambiguous"] == 1
    assert a["ambiguous"][0]["pdb_id"] == "2BBB"
    assert a["ambiguous"][0]["alternatives"] == [1]


def test_interface_area_is_reported_as_quartiles_not_a_range():
    """Buried area scales with the whole assembly, so one 60-mer sets a maximum three orders
    of magnitude above the median (lysozyme: median 1,476, maximum 615,514) and a min-to-max
    range describes that one entry rather than the family."""
    from codswallop.family import build_assemblies
    rows = [_entry(f"{i}XXX", 2, "dimeric", "author_defined_assembly", area=1000.0 + i,
                   ifres=50) for i in range(10)]
    rows.append(_entry("BIGX", 60, "60-meric", "author_defined_assembly",
                       area=615514.0, ifres=9000))
    a = build_assemblies(rows)
    assert a["interfaces"]["median_area"] < 2000, "one huge assembly must not move the median"
    assert "max_area" not in a["interfaces"]
    assert a["interfaces"]["q3_area"] < 2000


def test_a_family_with_no_assembly_annotation_returns_a_complete_shape():
    """The empty path must carry every key the populated one does, or the panel raises on
    the family it was meant to degrade gracefully for."""
    from codswallop.family import build_assemblies
    empty = build_assemblies([{"pdb_id": "1AAA"}])
    full = build_assemblies([_entry("1AAA", 1, "monomeric", "author_defined_assembly")])
    assert set(empty) == set(full)


def test_hot_residues_carry_their_own_type_breakdown():
    """"Residue 199 makes 47 contacts" is a ranking. "37 of them are hydrogen bonds and 6
    are metal coordination" is what says whether it is the catalytic centre or a wall of the
    pocket, and it can only be built while the raw rows still exist."""
    import inspect
    from codswallop import contacts
    src = inspect.getsource(contacts.build)
    assert '"types": res_types[pos].most_common()' in src
    assert '"ligands": res_ligands[pos].most_common(12)' in src
    # Counted over distinct entries, not contact rows: one entry with forty contacts is not
    # forty entries.
    assert '"entries": len(res_entries[pos])' in src


# ---- drugs ---------------------------------------------------------------------------
def test_a_drug_is_filed_under_every_class_it_belongs_to():
    """Taking the first ATC code filed acetazolamide, the best-known carbonic anhydrase
    inhibitor there is, under "genito-urinary": G01AE10 happens to sort before S01EC01."""
    from codswallop.drugs import therapeutic_classes
    got = therapeutic_classes(["G01AE10", "S01EC01"])
    assert "Genito-urinary & sex hormones" in got and "Sensory organs" in got


def test_anti_infective_and_antineoplastic_split_at_the_second_level():
    """"Anti-infective" hides the difference between an antibacterial and an antiviral."""
    from codswallop.drugs import therapeutic_classes
    assert therapeutic_classes(["J01CA04"]) == ["Antibacterial"]
    assert therapeutic_classes(["J05AP01"]) == ["Antiviral"]
    assert therapeutic_classes(["L01EA01"]) == ["Oncology"]
    assert therapeutic_classes(["L04AA01"]) == ["Immunology (suppressant)"]
    assert therapeutic_classes([]) == ["Unclassified"]


def test_an_isoform_number_is_not_a_substring_match():
    """"Carbonic anhydrase 2" and "Carbonic anhydrase 12" are different proteins, and a
    naive substring test makes 2 match 12."""
    from codswallop.drugs import on_target
    fam = {"name": "Carbonic anhydrase 2"}
    assert on_target([{"name": "Carbonic anhydrase 12"}], fam) == "annotated", \
        "across a family assembled at 30% identity, an isoform is on target"
    assert on_target([{"name": "GTPase HRas"}], {"name": "GTPase KRas"}) is None


def test_a_described_tie_is_not_reported_as_an_annotated_one():
    """The text match is the weaker signal and has to stay labelled as such: DrugBank's
    target lists lag, and KRAS has approved inhibitors with no annotation at all."""
    from codswallop.drugs import on_target
    fam = {"name": "Carbonic anhydrase 2"}
    assert on_target([], fam, "inhibits carbonic anhydrase in the proximal tubule") == "described"
    assert on_target([], fam, "an antibiotic that binds the ribosome") is None


def test_the_stage_ladder_prefers_the_strongest_evidence():
    from codswallop.drugs import _stage
    assert _stage(["approved", "investigational"], {"phase": "Phase 1", "n_studies": 3})[0] \
        == "Approved", "an approved drug is approved whatever its trials say"
    assert _stage(["investigational"], {"phase": "Phase 3", "n_studies": 9})[0] == "Phase 3"
    assert _stage(["experimental"], None)[0] == "Preclinical"
    assert _stage([], None)[0] == "Unknown"


def test_a_drug_untied_to_the_protein_is_not_a_drug_for_this_pocket():
    """Isopropyl alcohol is an approved topical antiseptic and a cryoprotectant in four KRAS
    crystals. Both are true; only one is about the protein. Listing every DrugBank-linked
    component answers "what happens to be in DrugBank" under a heading that claims to answer
    "what drugs this target", and every drug the KRAS family produced was untied."""
    from codswallop import drugs

    fam = {"name": "GTPase KRas",
           "ligands": {"components": [
               {"id": "IPA", "klass": "ligand", "count": 4, "name": "ISOPROPYL ALCOHOL"}]}}
    ann = {"IPA": {"drugbank_id": "DB02325", "generic": "Isopropyl alcohol",
                   "groups": ["approved"], "atc": ["D08AX05"], "brands": [],
                   "targets": [{"name": "Tumor necrosis factor"}],
                   "indication": "topical antiseptic", "mechanism": "denatures protein"}}
    real_annotate = drugs.annotate
    drugs.annotate = lambda ids: ann
    try:
        out = drugs.build(fam)
    finally:
        drugs.annotate = real_annotate

    assert out["n"] == 0, "an untied component must not be listed as a drug for this target"
    assert out["n_untied"] == 1, "but it is counted and named, not silently dropped"
    assert out["untied"][0]["generic"] == "Isopropyl alcohol"


def test_small_solvents_are_cryoprotectants_not_ligands():
    """They arrive from the cryo or the mother liquor, never from biology."""
    from codswallop.ligands import classify
    for cid in ("IPA", "EOH", "DIO", "TBU"):
        assert classify(cid, "") == "cryoprotectant", f"{cid} should not read as a ligand"


def test_the_request_path_never_fetches_drug_enrichment():
    """One third-party call per drug, and ABL1 has 65 drug-like components: doing them
    inside a page load put the request past two minutes on a two-core droplet. The page
    reads what is cached; the warm is what fills the cache."""
    import inspect
    from codswallop import drugs
    sig = inspect.signature(drugs.build)
    assert sig.parameters["fetch_missing"].default is False
    for fn in (drugs.highest_trial_phase, drugs.fda_record):
        src = inspect.getsource(fn)
        assert "if not fetch_missing:" in src and "cache_get" in src, \
            f"{fn.__name__} must be able to answer from the cache alone"


def test_the_third_map_axis_comes_from_the_matrix_already_shipped():
    """Derived from the TM matrix the artefact already carries for the heatmap, rather than
    stored: a new field would have meant a pipeline version bump and a full rebuild for a
    twenty-line eigendecomposition of something already on hand."""
    from codswallop.layout import third_axis
    # Four points on a line: the third principal coordinate carries nothing, and inventing
    # one from a non-positive eigenvalue would draw noise as structure.
    flat = [[1.0, 0.9, 0.8, 0.7], [0.9, 1.0, 0.9, 0.8],
            [0.8, 0.9, 1.0, 0.9], [0.7, 0.8, 0.9, 1.0]]
    assert third_axis(flat) is None

    # A tetrahedron is genuinely three-dimensional and must produce one.
    tetra = [[1.0, 0.5, 0.5, 0.5], [0.5, 1.0, 0.5, 0.5],
             [0.5, 0.5, 1.0, 0.5], [0.5, 0.5, 0.5, 1.0]]
    zs = third_axis(tetra)
    assert zs is not None and len(zs) == 4
    assert max(abs(z) for z in zs) > 0, "a real third axis is not all zeros"


def test_too_few_points_have_no_third_axis():
    from codswallop.layout import third_axis
    assert third_axis([[1.0, 0.5], [0.5, 1.0]]) is None
    assert third_axis([]) is None


def test_a_huge_assembly_is_not_parsed_whole_when_a_chain_will_do(monkeypatch, tmp_path):
    """8GLV is a 453 MB deposited file that yields 426 alpha carbons, and biotite holds all
    of it in memory to get them: 6.3 GB peak RSS, on a droplet with 2.2 GB and no swap.

    Asking the Model Server for the one chain returns 533 kB and peaks at 102 MB. So the
    chain route must be tried FIRST whenever a chain is known, not as a fallback.
    """
    from codswallop import embed

    calls = []

    def fake_download(url, dest, skip_if_exists=True, params=None):
        calls.append(url)
        from pathlib import Path
        Path(dest).parent.mkdir(parents=True, exist_ok=True)
        Path(dest).write_bytes(b"x" * 4096)
        return dest

    monkeypatch.setattr(embed.http, "download", fake_download)
    monkeypatch.setattr(embed, "STRUCT_DIR", tmp_path)
    monkeypatch.setattr(embed, "get_structure", lambda *a, **k: (_ for _ in ()).throw(ValueError()))
    monkeypatch.setattr(embed, "CIFFile", type("F", (), {"read": staticmethod(lambda p: None)}))

    embed.ca_trace("8GLV", "0A")
    assert calls and "models.rcsb.org" in calls[0], \
        "the whole deposited file was fetched even though the chain was known"


def test_a_whole_file_too_big_to_parse_safely_is_refused(monkeypatch, tmp_path):
    """Better to lose one structure from the family than to be killed by the OOM reaper
    half way through it, which on a box with no swap takes the other apps down too."""
    from codswallop import embed

    big = tmp_path / "8GLV.cif"
    big.write_bytes(b"x" * int((embed.MAX_WHOLE_FILE_MB + 1) * 1e6))
    monkeypatch.setattr(embed, "STRUCT_DIR", tmp_path)
    monkeypatch.setattr(embed, "_chain_cif", lambda *a: None)
    assert embed.ca_trace("8GLV", "0A") is None


def test_a_structure_too_big_for_contacts_is_counted_apart_from_a_failure(monkeypatch, tmp_path):
    """PLIP needs the WHOLE structure, so the single-chain fetch that made the embedding
    runnable on a small box does not apply and the 6.3 GB worst case returns.

    A skip and a failure must not be added together: a conversion that broke is a bug to
    chase, a structure too large to hold in memory is a decision this made, and the panel
    says so rather than implying those complexes had no interactions.
    """
    from codswallop import contacts

    monkeypatch.setattr(contacts, "STRUCT_DIR", tmp_path)
    monkeypatch.setattr(contacts, "_load_cif2plip", lambda: None)
    # download() returns None when it aborts past the cap.
    import codswallop.http as http
    monkeypatch.setattr(http, "download", lambda *a, **k: None)
    assert contacts.contacts_for("8GLV", tmp_path) == "too_big"


def test_an_oversized_download_is_aborted_rather_than_finished_and_refused(monkeypatch, tmp_path):
    """A HEAD on files.rcsb.org times out, so the size cannot be known in advance. Streaming
    until the cap is passed costs at most the cap; finishing a 453 MB transfer and then
    declining to parse it costs 453 MB of bandwidth and disk on a box with neither spare."""
    from codswallop import http

    class FakeResp:
        status_code = 200
        def raise_for_status(self): pass
        def iter_content(self, chunk_size=0):
            for _ in range(100):
                yield b"x" * 1000
        def close(self): pass

    monkeypatch.setattr(http, "get", lambda *a, **k: FakeResp())
    dest = tmp_path / "big.cif"
    assert http.download("http://x/big.cif", dest, max_bytes=10_000) is None
    assert not dest.exists(), "a refused download must leave nothing behind"
