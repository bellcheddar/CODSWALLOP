"""Names and people. Both panels had a confound that made a plausible wrong answer."""
from codswallop import provenance

SEED = {
    "accession": "P00698",
    "name": "Lysozyme C",
    "alt_names": ["1,4-beta-N-acetylmuramidase C", "Allergen Gal d IV"],
    "short_names": [],
    "genes": ["LYZ"],
}


def _m(desc, acc="P00698", **extra):
    return {"description": desc, "uniprot": acc, **extra}


# ---------------------------------------------------------------------------------------
# Names.
# ---------------------------------------------------------------------------------------

def test_case_and_punctuation_variants_collapse_but_real_differences_do_not():
    out = provenance.build_names(
        [_m("Lysozyme C"), _m("LYSOZYME C"), _m("lysozyme c"), _m("Lysozyme")], SEED)
    names = {r["name"]: r["n"] for r in out["rows"]}
    assert names["Lysozyme C"] == 3, "case variants are the same name"
    assert names["Lysozyme"] == 1, "'Lysozyme' is a different name from 'Lysozyme C'"


def test_a_relatives_correct_name_is_not_counted_as_a_misnaming_of_the_seed():
    """The confound this panel shipped with.

    A family is assembled at 30 % identity, so it legitimately holds other proteins.
    Counting description strings across all of it reported ALPHA-LACTALBUMIN and HUMAN
    LYSOZYME as unrecognised names for hen lysozyme, when they are the correct names of
    different proteins. Only entities carrying the seed's accession can say what the seed is
    called.
    """
    members = [_m("Lysozyme C")] * 8 + [
        _m("ALPHA-LACTALBUMIN", "P00711"), _m("HUMAN LYSOZYME", "P61626")]
    out = provenance.build_names(members, SEED)
    assert out["total"] == 8, "only the seed's own entities count towards its names"
    assert out["unrecognised"] == 0
    assert out["n_other"] == 2
    assert {o["accession"] for o in out["others"]} == {"P00711", "P61626"}


def test_the_same_name_used_by_two_proteins_is_reported_as_a_collision():
    """The direction that actually breaks a literature search: within this one family,
    'Lysozyme C' is the deposited description of both chicken P00698 and human P61626."""
    members = [_m("Lysozyme C")] * 5 + [_m("Lysozyme C", "P61626")] * 2
    out = provenance.build_names(members, SEED)
    hit = [c for c in out["collisions"] if c["name"] == "Lysozyme C"]
    assert hit, "a name shared by two accessions was not reported"
    assert {a["accession"] for a in hit[0]["accessions"]} == {"P00698", "P61626"}


def test_a_name_used_by_one_protein_is_not_a_collision():
    out = provenance.build_names([_m("Lysozyme C")] * 4, SEED)
    assert out["collisions"] == []


def test_alternative_and_gene_names_are_recognised_not_just_the_recommended_one():
    out = provenance.build_names(
        [_m("Lysozyme C"), _m("Allergen Gal d IV"), _m("LYZ")], SEED)
    kinds = {r["name"]: r["kind"] for r in out["rows"]}
    assert kinds["Lysozyme C"] == "recommended"
    assert kinds["Allergen Gal d IV"] == "alternative"
    assert kinds["LYZ"] == "gene"
    assert out["unrecognised"] == 0


def test_an_isoform_qualifier_is_not_a_made_up_name():
    """UniProt's own convention. Untreated it made a third of KRAS's entities read as
    deposited under a name nobody recognises."""
    seed = {"accession": "P01116", "name": "GTPase KRas", "alt_names": [],
            "short_names": [], "genes": ["KRAS"]}
    members = ([_m("GTPase KRas", "P01116")] * 3
               + [_m("Isoform 2B of GTPase KRas", "P01116")] * 2
               + [_m("GTPase KRas, N-terminally processed", "P01116")])
    out = provenance.build_names(members, seed)
    kinds = {r["name"]: r["kind"] for r in out["rows"]}
    assert kinds["Isoform 2B of GTPase KRas"] == "variant"
    assert kinds["GTPase KRas, N-terminally processed"] == "variant"
    assert out["unrecognised"] == 0


def test_a_genuinely_invented_name_is_still_reported_as_one():
    out = provenance.build_names([_m("Lysozyme C"), _m("HEN EGG WHITE LYSOZYME")], SEED)
    kinds = {r["name"]: r["kind"] for r in out["rows"]}
    assert kinds["HEN EGG WHITE LYSOZYME"] is None
    assert out["unrecognised"] == 1


def test_the_recognised_count_does_not_move_when_the_tail_is_pooled():
    """The headline figure is counted over every spelling, not over the rows that survived
    the display cut, or pooling the tail would quietly change the statistic."""
    members = [_m("Lysozyme C")] * 400 + [_m("Odd name %d" % i) for i in range(40)]
    out = provenance.build_names(members, SEED)
    assert out["pooled"] > 0, "this fixture is meant to exercise pooling"
    assert out["recognised"] == 400
    assert out["recognised"] + out["unrecognised"] == out["total"]


def test_no_seed_accession_still_lists_what_the_archive_says():
    out = provenance.build_names([_m("Lysozyme C", None)] * 3,
                                 {"name": None, "alt_names": [], "genes": []})
    assert out["total"] == 3 and out["rows"][0]["kind"] is None


# ---------------------------------------------------------------------------------------
# People.
# ---------------------------------------------------------------------------------------

def _e(pi, first="Smith, A.", year="2020-01-01", doi="10.1/a", **extra):
    return {"deposit_date": year, "method": "X-ray",
            "citation": {"authors": [first, pi], "doi": doi}, **extra}


def test_groups_are_counted_by_entry_not_by_paper():
    """One paper routinely covers a series of depositions. Counting papers would say a group
    that solved thirty structures and published them together did less work than one that
    published three papers about one structure each."""
    entries = [_e("Yutani, K.", doi="10.1/one") for _ in range(30)]
    entries += [_e("Other, B.", doi="10.1/%d" % i) for i in range(3)]
    out = provenance.build_people(entries)
    top = out["rows"][0]
    assert top["pi"] == "Yutani, K."
    assert top["entries"] == 30 and top["papers"] == 1


def test_entries_with_no_citation_are_counted_in_the_total_but_in_no_group():
    entries = [_e("Yutani, K.") for _ in range(4)] + [
        {"deposit_date": "1990-01-01", "citation": None}]
    out = provenance.build_people(entries)
    assert out["n_entries"] == 5 and out["n_cited"] == 4 and out["n_uncited"] == 1
    assert sum(r["entries"] for r in out["rows"]) == 4


def test_the_timeline_counts_every_entry_including_uncited_ones():
    entries = [_e("A, A.", year="2001-05-05"), _e("B, B.", year="2001-06-06"),
               {"deposit_date": "1999-01-01", "citation": None}]
    out = provenance.build_people(entries)
    counts = {t["year"]: t["n"] for t in out["timeline"]}
    assert counts == {1999: 1, 2001: 2}
    assert out["span"] == [1999, 2001]


def test_a_single_author_paper_does_not_credit_them_as_their_own_first_author():
    out = provenance.build_people([{"deposit_date": "2020-01-01",
                                    "citation": {"authors": ["Solo, S."], "doi": "10.1/x"}}])
    assert out["rows"][0]["pi"] == "Solo, S."
    assert out["rows"][0]["top_first_author"] is None


def test_concentration_is_reported_against_every_entry_not_only_the_cited_ones():
    """Otherwise a family where half the entries have no citation reports the top group as
    twice the share of the archive it actually holds."""
    entries = [_e("Big, B.") for _ in range(5)] + [
        {"deposit_date": "1990-01-01", "citation": None} for _ in range(5)]
    out = provenance.build_people(entries)
    assert out["top_share"] == 0.5


def test_an_empty_family_does_not_raise():
    out = provenance.build_people([])
    assert out["n_groups"] == 0 and out["timeline"] == [] and out["span"] is None


# ---------------------------------------------------------------------------------------
# Curated sites on a seed that is not the canonical sequence.
# ---------------------------------------------------------------------------------------

def test_uniprot_features_land_on_a_mature_seed_by_alignment(monkeypatch):
    """A PDB-seeded family used to get no curated sites at all: `lysozyme-1aki-1` showed no
    active site and no disulphides while `lysozyme-c-p00698`, the same protein, showed both.

    The guard existed for a real reason. UniProt numbers on the canonical, and a PDB-seeded
    family's seed is typically the mature protein: hen lysozyme is 129 residues against the
    canonical's 147, so every feature is eighteen out. Shifting by a fixed 18 would be the
    assumed-offset mistake that has already cost this codebase two panels, so it aligns.
    """
    from codswallop import family, uniprot

    canonical = "MRSLLILVLCFLPLAALG" + "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
    mature = canonical[18:]
    fam = {"kind": "pdb_id", "seed": "1AKI_1", "seed_sequence": mature,
           "members": [{"entity_id": "1AKI_1", "uniprot": "P00698"}]}
    feats = {
        "active_site": [53, 70],
        "signal_peptide": [{"start": 1, "end": 18, "description": ""}],
        "disulphide": [{"start": 24, "end": 145, "description": ""}],
    }
    monkeypatch.setattr(uniprot, "entry_with_features",
                        lambda acc: ({"sequence": canonical}, feats))

    out = family._features_on_seed(fam)
    # Glu35 and Asp52 are hen lysozyme's catalytic pair in mature numbering.
    assert out["active_site"] == [35, 52]
    assert mature[34] == "E" and mature[51] == "D"
    # Cys6-Cys127 is the mature form of the canonical's 24-145.
    assert out["disulphide"] == [{"start": 6, "end": 127, "description": ""}]
    # The signal peptide is entirely outside the mature protein and is dropped whole rather
    # than clamped to a shorter peptide that does not exist.
    assert out["signal_peptide"] == []
    assert all(isinstance(p, int) for p in out["active_site"]), \
        "numpy integers are not JSON-serialisable and would fail at render time"


def test_features_are_not_moved_onto_a_seed_that_is_a_different_protein(monkeypatch):
    from codswallop import family, uniprot
    monkeypatch.setattr(uniprot, "entry_with_features",
                        lambda acc: ({"sequence": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"},
                                     {"active_site": [5]}))
    fam = {"kind": "pdb_id", "seed": "9XYZ_1",
           "seed_sequence": "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW",
           "members": [{"entity_id": "9XYZ_1", "uniprot": "P00000"}]}
    assert family._features_on_seed(fam) == {}
