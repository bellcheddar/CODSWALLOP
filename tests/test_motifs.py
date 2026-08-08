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
