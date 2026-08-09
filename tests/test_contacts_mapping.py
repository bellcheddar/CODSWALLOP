"""Putting PLIP's residue numbers into seed coordinates.

This shipped wrong for 52 of 71 built families and looked right on the one it was validated
against, which is the whole reason these tests exist.
"""
from codswallop import contacts

# Carbonic anhydrase II, whose author numbering follows the canonical: this is the family
# the Contacts panel was validated on, and the family the bug could not show up in.
CA2 = ("MSHHWGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKYDPSLKPLSVSYDQATSLRILNNGHAFNVEFDDSQDKAVLK"
       "GGPLDGTYRLIQFHFHWGSLDGQGSEHTVDKKKYAAELHLVHWNTKYGDFGKAVQQPDGLAVLGIFLKVGSAKPGLQKVV"
       "DVLDSIKTKGKSADFTNFDPRGLLPESLDYWTYPGSLTTPPLLECVTWIVLKEPISVSSEQVLKFRKLNFNGEGEPEELM"
       "VDNWRPAQPLKNRQIKASFK")
# An unrelated chain of the sort that turns up co-crystallised.
PARTNER = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR"


def _chain(seq, start=1):
    return [(start + i, aa) for i, aa in enumerate(seq)]


def test_canonical_numbering_maps_to_itself():
    """The zinc triad and the gatekeeper threonines, which is what the panel is known for."""
    m = contacts._chain_mappings({"A": _chain(CA2)}, CA2, only={"A"})
    for pos in (94, 96, 119, 199, 200):
        assert m[("A", pos)] == pos


def test_a_construct_numbered_on_the_canonical_is_not_offset_twice():
    """The regression.

    `seed_pos = resnr + (query_beg - 1)` added the seed offset to a number that already
    carried it. JAK1 (query_beg 879, seed 1,154 residues) reported hot residues at 1,340 and
    2,110: positions that do not exist. Every position must land inside the sequence.
    """
    frag = CA2[100:200]
    m = contacts._chain_mappings({"A": _chain(frag, start=879)}, CA2, only={"A"})
    assert m, "the fragment did not map at all"
    assert all(1 <= v <= len(CA2) for v in m.values())
    # And it maps where the residues actually are, not merely somewhere in range.
    assert m[("A", 879)] == 101
    assert m[("A", 879 + 99)] == 200


def test_a_construct_numbered_from_one_also_lands_correctly():
    """The other numbering convention in the archive: a construct numbered from 1 whatever
    the canonical says. An assumed offset cannot be right for both this and the case above."""
    frag = CA2[100:200]
    m = contacts._chain_mappings({"A": _chain(frag, start=1)}, CA2, only={"A"})
    assert m[("A", 1)] == 101 and m[("A", 100)] == 200


def test_a_partner_chain_is_not_filed_against_the_seeds_residues():
    """Which chains are the family member is READ from the entity record, not inferred.

    Inference does not work here: this partner aligns to carbonic anhydrase at 37.8 %
    identity over 74 columns, which is inside the twilight zone and above any threshold that
    still admits a legitimate 35 %-identity orthologue.
    """
    obs = {"A": _chain(CA2), "B": _chain(PARTNER)}
    m = contacts._chain_mappings(obs, CA2, only={"A"})
    assert not [k for k in m if k[0] == "B"]
    assert len([k for k in m if k[0] == "A"]) == len(CA2)


def test_every_copy_of_the_member_is_mapped_not_just_the_first():
    """A contact made by chain B of a homodimer is a contact this family made."""
    obs = {"A": _chain(CA2), "B": _chain(CA2)}
    m = contacts._chain_mappings(obs, CA2, only={"A", "B"})
    assert m[("A", 94)] == 94 and m[("B", 94)] == 94


def test_no_seed_sequence_maps_nothing_rather_than_guessing():
    assert contacts._chain_mappings({"A": _chain(CA2)}, "", only={"A"}) == {}


def test_an_unmappable_residue_is_dropped_rather_than_placed():
    """A residue outside the aligned region has no seed position, and a contact placed on
    the wrong residue is worse than a contact nobody counted."""
    frag = CA2[100:140]
    m = contacts._chain_mappings({"A": _chain(frag, start=1)}, CA2, only={"A"})
    assert ("A", 500) not in m
