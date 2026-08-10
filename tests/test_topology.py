"""Fold topology. The parsing and the mapping, which are where the wrong answers live."""
from codswallop import topology

# Two residues of real DSSP output, with the header line the parser looks for. Fixed columns,
# so this is copied verbatim rather than reconstructed: the whole point is the column offsets.
DSSP = (
    "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N\n"
    "    1    1 A M              0   0  142      0\n"
    "    2    2 A T  E     -a   51   0A  35     48\n"
    "    3    3 A E  E     -a   52   0A  79     -2\n"
)


def test_dssp_is_parsed_by_column_not_by_split():
    """The fields butt against each other on a four-digit residue number, and a chain-break
    line has almost nothing in it, so whitespace splitting drifts silently."""
    recs = topology._parse_dssp(DSSP, None)
    assert [r["res"] for r in recs] == [1, 2, 3]
    assert [r["sse"] for r in recs] == ["-", "E", "E"]
    assert [r["aa"] for r in recs] == ["M", "T", "E"]
    assert recs[1]["bp1"] == 51 and recs[2]["bp1"] == 52


def test_a_chain_break_line_is_skipped():
    text = DSSP + "    4        !              0   0    0      0\n"
    assert len(topology._parse_dssp(text, None)) == 3


def test_elements_stop_at_a_gap_in_the_seed():
    """Two halves either side of an unmodelled loop are two elements, not one long one
    spanning residues nobody saw."""
    recs = [{"res": i, "sse": "H", "idx": i, "aa": "A", "bp1": 0, "bp2": 0}
            for i in range(1, 13)]
    mapping = {i: i for i in range(1, 7)}          # 7-12 unobserved in the seed
    mapping.update({i: i + 20 for i in range(7, 13)})
    els = topology._elements(recs, mapping)
    assert len(els) == 2
    assert (els[0]["start"], els[0]["end"]) == (1, 6)
    assert (els[1]["start"], els[1]["end"]) == (27, 32)


def test_short_runs_are_not_drawn():
    """A three-residue 'helix' is a wobble in the backbone, not an element."""
    recs = [{"res": i, "sse": "H", "idx": i, "aa": "A", "bp1": 0, "bp2": 0} for i in (1, 2, 3)]
    assert topology._elements(recs, {1: 1, 2: 2, 3: 3}) == []


def test_pairings_need_more_than_one_bridge():
    """One bridge is noise and would draw the same line as a ten-residue pairing."""
    recs = [
        {"res": 1, "sse": "E", "idx": 1, "aa": "A", "bp1": 10, "bp2": 0},
        {"res": 2, "sse": "E", "idx": 2, "aa": "A", "bp1": 0, "bp2": 0},
        {"res": 10, "sse": "E", "idx": 10, "aa": "A", "bp1": 1, "bp2": 0},
        {"res": 11, "sse": "E", "idx": 11, "aa": "A", "bp1": 0, "bp2": 0},
    ]
    els = [{"id": 0, "kind": "strand", "residues": [1, 2]},
           {"id": 1, "kind": "strand", "residues": [10, 11]}]
    assert topology._pairings(recs, els) == [{"a": 0, "b": 1, "bridges": 2}]


def test_mapping_is_by_alignment_not_by_offset():
    """An entry may be numbered on the mature protein, on the construct from 1, or on the
    canonical. An assumed offset is right often enough to look like it works."""
    seed = "MKVLAAGIVGLNLQ"
    # The structure starts at the seed's fifth residue and numbers it 101.
    recs = [{"res": 100 + i, "aa": c, "sse": "-", "idx": i, "bp1": 0, "bp2": 0}
            for i, c in enumerate(seed[4:], start=1)]
    mapping = topology.map_to_seed(recs, seed)
    assert mapping[101] == 5, "the first observed residue is the seed's fifth"
    assert mapping[110] == 14


def test_the_topology_aligner_survives_an_unknown_residue_too():
    """The same BLOSUM62 gap, in the second of the two aligners. Fixing it in constructs.py
    alone left serum albumin failing here instead: one selenocysteine, and the family lost
    its whole topology artefact while its embedding and contacts built perfectly."""
    seed = "MKWVTFISLLFLFSSAYSRGVFRR"
    recs = [{"res": i, "aa": c, "sse": "-", "idx": i, "bp1": 0, "bp2": 0}
            for i, c in enumerate(seed.replace("F", "U", 1), start=1)]
    mapping = topology.map_to_seed(recs, seed)
    assert mapping, "an unknown residue must not empty the mapping"
    assert len(mapping) >= len(seed) - 2


def test_topology_asks_for_a_path_rather_than_relying_on_a_download_side_effect(monkeypatch):
    """`ca_trace` used to fetch the whole entry and topology depended on that happening.

    When the fetch moved to the Model Server the whole file stopped appearing, `_cif_path`
    no longer existed, and topology returned no records at all with no error, because a
    missing file was already a legitimate "no 2D layout for this structure". A side effect
    another module relies on is not an interface, so there is one now.
    """
    from codswallop import embed, topology

    asked = []
    monkeypatch.setattr(embed, "structure_path",
                        lambda pdb, chain=None: asked.append((pdb, chain)) or None)
    assert topology._sse_from_biotite("1ABC", "A") == []
    assert asked == [("1ABC", "A")], "topology no longer goes through structure_path"


def test_structure_path_prefers_the_chain(monkeypatch, tmp_path):
    from codswallop import embed
    monkeypatch.setattr(embed, "_chain_cif", lambda p, c: tmp_path / f"{p}_{c}.cif")
    got = embed.structure_path("8GLV", "0A")
    assert got.name == "8GLV_0A.cif", "the whole assembly was chosen over the chain"
