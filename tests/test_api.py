"""The v1 summary. A contract other projects build against, so its shape is the test."""
from codswallop import api

FAM = {
    "slug": "lysozyme-1aki-1", "query": "1AKI", "kind": "pdb_id", "seed": "1AKI_1",
    "name": "LYSOZYME", "organism": "Gallus gallus", "seed_length": 129,
    "seed_sequence": "KVFGRCELAA",
    "stats": {"entries": 1687, "entities": 1688, "constructs": 281, "organisms": 26,
              "holo_entries": 324, "best_resolution": 0.65, "median_resolution": 1.74,
              "identity_min": 35.6, "identity_max": 100.0, "methods": {"X-ray": 1638}},
    "ligands": {"n": 277, "components": [
        {"id": "NAG", "name": "N-acetylglucosamine", "klass": "ligand", "count": 40,
         "best_resolution": 1.1},
        {"id": "GOL", "name": "Glycerol", "klass": "cryoprotectant", "count": 300,
         "best_resolution": 0.9},
    ]},
    "contacts": {"hot_residues": [{"pos": 35, "restype": "GLU", "contacts": 40,
                                   "entries": 12}]},
    "topology": {"elements": [1], "n_strands": 3, "n_helices": 6, "method": "DSSP",
                 "reference": "1IOT"},
    "crystals": {"n_parsed": 1414, "precipitants": [{"name": "Sodium chloride"}],
                 "buffers": [{"name": "Sodium acetate"}]},
    "assemblies": {"states": [{"count": 1, "details": "monomeric", "fraction": 90.8}]},
    "map": {"embedded": True, "reference": "1IOT"},
    "provenance": {
        "names": {"accession": "P00698", "n": 7, "unrecognised": 235, "n_collisions": 5,
                  "official": [{"name": "Lysozyme C", "kind": "recommended"}]},
        "people": {"n_groups": 380, "rows": [{"pi": "Yutani, K."}], "top_share": 0.083,
                   "span": [1975, 2026]},
    },
    "pfam": [{"id": "PF00062", "name": "C-type lysozyme"}], "interpro": [],
    "built_at": 1,
}


def test_the_summary_is_small_and_carries_its_schema_version():
    import json
    out = api.summarise(FAM, base_url="https://example.org")
    assert out["schema_version"] == api.SCHEMA_VERSION
    assert out["url"] == "https://example.org/f/lysozyme-1aki-1"
    # The point of a separate endpoint is that it is cheap to poll.
    assert len(json.dumps(out)) < 20000


def test_every_panel_is_optional_and_a_bare_family_does_not_raise():
    """Each panel in this app degrades on its own, so a consumer must get a family with
    empty fields rather than a 500 because PLIP has not run."""
    out = api.summarise({"slug": "x", "stats": {}})
    assert out["slug"] == "x"
    assert out["ligands"] == [] and out["hot_residues"] == []
    assert out["topology"] is None and out["assembly"] is None
    assert out["artefacts"] == {"embedding": False, "contacts": False, "topology": False}


def test_only_components_somebody_meant_to_be_there_are_listed():
    """A summary that put glycerol above the inhibitor would actively mislead a consumer
    looking for chemistry."""
    out = api.summarise(FAM)
    assert [l["id"] for l in out["ligands"]] == ["NAG"]


def test_topology_counts_are_taken_from_the_build_not_recounted():
    out = api.summarise(FAM)
    assert out["topology"] == {"strands": 3, "helices": 6, "method": "DSSP",
                               "reference": "1IOT"}


def test_a_contact_position_past_the_end_of_the_seed_withholds_the_whole_set():
    """`contacts.py` adds the query offset to PLIP's author residue number, which already
    follows the canonical numbering, so the offset is counted twice: JAK1's hot residues
    come out at 1,340 and 2,110 on a 1,154-residue seed. Carbonic anhydrase has query_beg 1,
    so its offset is zero and its numbers are right, which is why the panel validated on it.

    The in-range positions of a broken family are wrong by the same offset and merely happen
    to land inside the sequence, so the set is withheld rather than filtered.
    """
    broken = dict(FAM, seed_length=1154, contacts={"hot_residues": [
        {"pos": 1340, "restype": "LEU", "contacts": 40, "entries": 3},
        {"pos": 400, "restype": "ALA", "contacts": 10, "entries": 2},
    ]})
    assert api.contacts_positions_look_sane(broken) is False
    out = api.summarise(broken)
    assert out["hot_residues"] == []
    assert out["warnings"] and "seed coordinates" in out["warnings"][0]
    # The artefact still exists; it is the coordinates that are not to be believed.
    assert out["artefacts"]["contacts"] is True


def test_a_family_whose_positions_all_fit_is_published_with_no_warning():
    out = api.summarise(FAM)
    assert out["warnings"] == []
    assert out["hot_residues"][0] == {"seed_position": 35, "residue": "GLU",
                                      "contacts": 40, "entries": 12}


def test_the_seed_accession_comes_from_the_reconciled_names_not_a_vote():
    out = api.summarise(FAM)
    assert out["seed"]["accession"] == "P00698"
    assert out["names"]["recommended"] == "Lysozyme C"
