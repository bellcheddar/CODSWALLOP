"""The construct diff engine. Every case here is one that shipped a wrong answer."""
from codswallop import constructs

# A short stand-in for a target protein, and the real T4 lysozyme sequence, which is what
# a GPCR crystallisation construct carries in place of intracellular loop 3.
TARGET = ("MGQPGNGSAFLLAPNRSHAPDHDVTQQRDEVWVVGMGIVMSLIVLAIVFGNVLVITAIAKFERLQTVTNYFITSLA"
          "CADLVMGLAVVPFGAAHILMKMWTFGNFWCEFWTSIDVLCVTASIETLCVIAVDRYFAITSPFKYQSLLTKNKARV"
          "IILMVWIVSGLTSFLPIQMHWYRATHQEAINCYANETCCDFFTNQAYAIASSIVSFYVPLVIMVFVYSRVFQEAKR")
T4L = constructs.FUSIONS[[f[0] for f in constructs.FUSIONS].index("T4 lysozyme")][1]


def test_his_tag_and_protease_site_are_named():
    d = constructs.diff(TARGET, "MGSSHHHHHHSSGLVPRGSH" + TARGET)
    assert "His6" in d["tags"]
    assert "Thrombin" in d["proteases"]


def test_truncation_is_reported_as_a_span_not_as_mutations():
    d = constructs.diff(TARGET, TARGET[30:200])
    assert d["canonical_span"] == [31, 200]
    assert d["mutation_count"] == 0, "a truncation must not read as a wall of substitutions"


def test_point_mutation_is_labelled_with_its_position():
    mutant = TARGET[:49] + "W" + TARGET[50:]
    d = constructs.diff(TARGET, mutant)
    assert [m["label"] for m in d["mutations"]] == [f"{TARGET[49]}50W"]


def test_internal_fusion_is_found_and_does_not_shred_the_alignment():
    """A GPCR with T4 lysozyme replacing an internal loop.

    The regression this guards: a 164-residue insertion costs a global aligner about what
    mis-aligning those residues costs, so it used to carve T4 lysozyme into small fake indels
    spread across the receptor, report no fusion, and invent ~20 point mutations on the way.
    """
    chimera = TARGET[:150] + T4L + TARGET[150:]
    d = constructs.diff(TARGET, chimera)
    assert "T4 lysozyme" in d["fusions"]
    assert "T4 lysozyme" in d["internal_fusions"]
    assert d["mutation_count"] == 0, "the partner must be excised, not aligned to the target"
    site = d["fusion_sites"][0]
    assert site["where"] == "internal"
    assert abs(site["after"] - 150) <= 2


def test_terminal_fusion_is_placed_at_the_terminus():
    d = constructs.diff(TARGET, T4L + TARGET)
    assert "T4 lysozyme" in d["fusions"]
    assert d["fusion_sites"][0]["where"] == "N-terminal"
    assert "T4 lysozyme" not in d.get("internal_fusions", [])


def test_active_site_mutation_is_described_not_diagnosed():
    mutant = TARGET[:63] + "A" + TARGET[64:]
    d = constructs.diff(TARGET, mutant, {"active_site": [64]})
    classes = d["mutations"][0]["classes"]
    assert any("active site" in c for c in classes)
    # It must never claim the enzyme was inactivated: that is a result, not an annotation.
    assert not any("inactiv" in c.lower() for c in classes)


def test_identical_sequence_reports_nothing():
    d = constructs.diff(TARGET, TARGET)
    assert constructs.summarise(d) == "matches the canonical sequence"
    assert not d["is_engineered"]
