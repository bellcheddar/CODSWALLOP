"""Functional sites. Every case here is one that shipped a wrong or useless answer."""
from codswallop import motifs, uniprot


def _fam(depth, seen, max_depth, seq="A" * 147):
    return {"seed_sequence": seq,
            "stats": {"coverage": {"depth": depth, "seen": seen, "max_depth": max_depth}}}


def test_a_disulphide_is_two_cysteines_and_not_the_span_between_them():
    """UniProt gives a disulphide's two paired cysteines as start and end. Expanding that to
    a range turned lysozyme's 24-145 bond into 122 consecutive "disulphide bonds", one per
    residue, and a 10-row panel into a 128-row one that was wrong about all but four."""
    rec = {"features": [
        {"type": "Disulfide bond",
         "location": {"start": {"value": 24}, "end": {"value": 145}}}]}
    out = uniprot._features_from(rec)
    assert out["disulphide"] == [{"start": 24, "end": 145, "description": ""}]


def test_grounding_expresses_resolved_against_what_was_cloned():
    """Otherwise a site nobody has ever built into a construct reads as disordered, which is
    the opposite of what the data says: it says nothing at all about it."""
    # 10 residues; the first 5 are in every construct and resolved, the last 5 in none.
    depth = [100] * 5 + [0] * 5
    seen = [90] * 5 + [0] * 5
    g = motifs._ground(1, 5, depth, seen, 100)
    assert g["in_constructs"] == 100.0 and g["resolved"] == 90.0
    g2 = motifs._ground(6, 10, depth, seen, 100)
    assert g2["in_constructs"] == 0.0
    assert g2["resolved"] is None, "never cloned is unknown, not unresolved"


def test_grounding_averages_over_the_span():
    """A transmembrane helix half of which is resolved is a different statement from one
    that is not resolved at all, so a span cannot be judged at its first residue."""
    depth = [100] * 10
    seen = [100] * 5 + [0] * 5
    assert motifs._ground(1, 10, depth, seen, 100)["resolved"] == 50.0


def test_a_site_off_the_end_of_the_seed_does_not_crash_or_lie():
    g = motifs._ground(200, 260, [10] * 10, [10] * 10, 10)
    assert g["in_constructs"] is None and g["resolved"] is None


def test_curated_sorts_before_predicted():
    """The panel's whole argument is that the two are different claims, so the ordering has
    to carry it even for a reader who never gets past the first screen."""
    fam = _fam([10] * 147, [10] * 147, 10)
    feats = {"active_site": [53], "signal_peptide": [{"start": 1, "end": 18}]}
    out = motifs.build(fam, feats)          # no network: PROSITE only runs on a real seq
    standings = [r["standing"] for r in out["rows"]]
    assert standings == sorted(standings, key=lambda s: s != "curated")
    assert out["n_curated"] == 2 and out["n_predicted"] == 0


def test_the_high_probability_prosite_patterns_are_off_by_default():
    """Scanning a secreted lysozyme with them returns a PKC phosphorylation site and two
    N-myristoylation sites, none of which happen to it. Six hits become two."""
    import inspect
    src = inspect.getsource(motifs.scan_prosite)
    assert 'include_frequent: bool = False' in src
    assert '"0" if include_frequent else "1"' in src


def test_uniprot_descriptions_are_tidied_not_invented():
    """A transmembrane helix arrives as `Helical; Name=1`, which is flat-file encoding and
    was being printed to the reader verbatim. The number survives; the encoding does not."""
    assert motifs.tidy("Helical; Name=1", "Transmembrane") == "Transmembrane helix 1"
    assert motifs.tidy("Helical", "Transmembrane") == "Transmembrane helix"
    assert motifs.tidy("", "Transmembrane") == "Transmembrane"
    # Anything that is already prose is left exactly alone.
    assert motifs.tidy("N-linked (GlcNAc...) asparagine", "Glycosylation") == \
        "N-linked (GlcNAc...) asparagine"


# ---- family-specific reference numbering ---------------------------------------------
def test_generic_numbers_survive_the_structure_based_suffix():
    """GPCRdb writes 3.50 as "3.50x50" under the structure-based scheme. Parsing that as a
    float raises, and a switch whose number will not parse is a switch silently missing."""
    from codswallop import pockets
    assert pockets._num("3.50x50") == 3.50
    assert pockets._num("7.49") == 7.49
    assert pockets._num(None) is None and pockets._num("N-term") is None


def test_switches_are_located_by_generic_number_not_by_sequence():
    """The point of a generic numbering scheme is that 3.50 is 3.50 in every receptor, so a
    receptor whose DRY is a DRF still has its switch found and shown as DRF."""
    from codswallop import pockets
    residues = [
        {"pos": 101, "aa": "D", "segment": "TM3", "generic": "3.49"},
        {"pos": 102, "aa": "R", "segment": "TM3", "generic": "3.50"},
        {"pos": 103, "aa": "F", "segment": "TM3", "generic": "3.51"},
    ]
    got = [s for s in pockets._switches(residues) if s["name"] == "DRY"]
    assert got and got[0]["sequence"] == "DRF"
    assert (got[0]["start"], got[0]["end"]) == (101, 103)


def test_segments_collapse_only_when_contiguous():
    """A segment that reappears after a gap is two spans, not one running through the gap."""
    from codswallop import pockets
    res = [{"pos": 1, "segment": "TM1"}, {"pos": 2, "segment": "TM1"},
           {"pos": 9, "segment": "TM1"}]
    assert pockets._segments(res) == [{"name": "TM1", "start": 1, "end": 2},
                                      {"name": "TM1", "start": 9, "end": 9}]


def test_a_family_that_is_neither_gets_neither():
    """The tab must not show an empty kinase pocket for a lysozyme."""
    from codswallop import pockets
    assert pockets.build({"kind": "pdb_id", "seed": "1AKI"})["kind"] is None
    assert pockets.build({"kind": "uniprot", "seed": ""})["kind"] is None
