"""A small, stable family summary for other programs.

`/api/family/<slug>` already exists and returns the whole internal payload: several
megabytes, shaped by whatever the front end needed that week, and no kind of contract. It is
fine for the page that ships with it and it is not something another project should build
against, because every change to a panel would be a breaking change to them.

So this is the deliberately narrow version. Field names are a promise; `schema_version` is
bumped if one of them has to change meaning, and fields are added rather than repurposed.
Everything is nullable: each panel in this app is independently optional, and a consumer
should get a family with an empty `hot_residues` rather than a 500 because PLIP has not run.

Kept to a few kilobytes so it is cheap to poll: the heavy things (every entry, the full
TM-score matrix, per-residue tracks) stay on the internal endpoint.
"""

from __future__ import annotations

from typing import Optional

SCHEMA_VERSION = 1

# How many rows each list-shaped field carries. A summary that returned all 239 of
# lysozyme's ligands would not be a summary.
TOP_N = 10


def _round(value, places=2):
    return round(value, places) if isinstance(value, (int, float)) else None


def _seed(fam: dict) -> dict:
    stats = fam.get("stats") or {}
    return {
        "kind": fam.get("kind"),
        "id": fam.get("seed"),
        "accession": ((fam.get("provenance") or {}).get("names") or {}).get("accession"),
        "name": fam.get("name"),
        "organism": fam.get("organism"),
        "length": fam.get("seed_length"),
        "sequence": fam.get("seed_sequence"),
        "identity_range": [stats.get("identity_min"), stats.get("identity_max")],
    }


def _counts(fam: dict) -> dict:
    stats = fam.get("stats") or {}
    return {
        "entries": stats.get("entries"),
        "entities": stats.get("entities"),
        "constructs": stats.get("constructs"),
        "organisms": stats.get("organisms"),
        "ligand_bound_entries": stats.get("holo_entries"),
        "distinct_ligands": (fam.get("ligands") or {}).get("n"),
    }


def _ligands(fam: dict) -> list:
    """Only components somebody meant to be there, which is what the Ligands panel decides.

    A summary listing glycerol and sulphate above the inhibitor would be actively
    misleading to a consumer that is looking for chemistry.
    """
    from . import ligands as ligand_engine
    out = []
    for c in ((fam.get("ligands") or {}).get("components") or []):
        if c.get("klass") not in ligand_engine.COUNTS_AS_BOUND:
            continue
        out.append({
            "id": c.get("id"), "name": c.get("name"), "class": c.get("klass"),
            "entries": c.get("count"), "best_resolution": _round(c.get("best_resolution")),
        })
        if len(out) >= TOP_N:
            break
    return out


def contacts_positions_look_sane(fam: dict) -> bool:
    """Whether this family's contact positions can be believed as seed coordinates.

    `contacts.py` maps PLIP's residue numbers with `seed_pos = resnr + (query_beg - 1)`.
    PLIP reports the AUTHOR residue number, which in a well-annotated entry already follows
    the canonical numbering, so adding the query offset counts it twice. On JAK1
    (query_beg 879, seed 1,154 residues) the hot residues come out at 1,340 and 2,110; on
    carbonic anhydrase query_beg is 1, the offset is zero, and the numbers are the correct
    His94/His96/Thr199/Thr200. 52 of 71 built families are affected.

    A position past the end of the seed is proof the mapping is wrong. It is not proof of
    the converse, so a family that passes this is only *not obviously* broken: one that
    fails is certainly broken and its numbers must not be published as seed coordinates.
    """
    length = fam.get("seed_length") or 0
    rows = (fam.get("contacts") or {}).get("hot_residues") or []
    if not length or not rows:
        return bool(rows)
    return not any((r.get("pos") or 0) > length for r in rows)


def _hot_residues(fam: dict) -> list:
    c = fam.get("contacts") or {}
    if not contacts_positions_look_sane(fam):
        # Withheld rather than filtered: the in-range positions of a family whose mapping is
        # broken are wrong by the same offset as the rest and merely happen to land inside
        # the sequence. Publishing those would be worse than publishing none.
        return []
    out = []
    for r in (c.get("hot_residues") or [])[:TOP_N]:
        out.append({
            "seed_position": r.get("pos"),
            "residue": r.get("restype"),
            "contacts": r.get("contacts"),
            "entries": r.get("entries"),
        })
    return out


def _topology(fam: dict) -> Optional[dict]:
    t = fam.get("topology") or {}
    if not t.get("elements"):
        return None
    # n_strands and n_helices are already counted by the topology build; recounting the
    # element list here would be a second definition of the same number, free to disagree.
    return {"strands": t.get("n_strands"), "helices": t.get("n_helices"),
            "method": t.get("method"), "reference": t.get("reference")}


def summarise(fam: dict, base_url: str = "") -> dict:
    """A decorated family, reduced to the parts another program is likely to want."""
    stats = fam.get("stats") or {}
    crystals = fam.get("crystals") or {}
    assemblies = fam.get("assemblies") or {}
    names = (fam.get("provenance") or {}).get("names") or {}
    people = (fam.get("provenance") or {}).get("people") or {}
    fam_map = fam.get("map") or {}
    top_state = (assemblies.get("states") or [None])[0]

    return {
        "schema_version": SCHEMA_VERSION,
        "slug": fam.get("slug"),
        "url": f"{base_url}/f/{fam.get('slug')}" if base_url else None,
        "query": fam.get("query"),
        "built_at": fam.get("built_at"),
        "seed": _seed(fam),
        "counts": _counts(fam),
        "resolution": {
            "best": _round(stats.get("best_resolution")),
            "median": _round(stats.get("median_resolution")),
        },
        "methods": stats.get("methods") or {},
        "domains": {
            "pfam": fam.get("pfam") or [],
            "interpro": fam.get("interpro") or [],
        },
        # The structure everything else was superposed onto, which is the one a consumer
        # should start from if it wants a single representative of the family.
        "representative": fam_map.get("reference"),
        "ligands": _ligands(fam),
        "hot_residues": _hot_residues(fam),
        "topology": _topology(fam),
        "crystallisation": {
            "parsed": crystals.get("n_parsed"),
            "top_precipitants": [p.get("name") for p in (crystals.get("precipitants") or [])[:5]],
            "top_buffers": [b.get("name") for b in (crystals.get("buffers") or [])[:5]],
        },
        "assembly": ({"state": top_state.get("details") or f"{top_state.get('count')}-mer",
                      "chains": top_state.get("count"),
                      "fraction": _round((top_state.get("fraction") or 0) / 100, 4)}
                     if top_state else None),
        "names": {
            "recommended": next((o["name"] for o in (names.get("official") or [])
                                 if o.get("kind") == "recommended"), None),
            "deposited_spellings": names.get("n"),
            "unrecognised_entities": names.get("unrecognised"),
            "name_collisions": names.get("n_collisions"),
        },
        "people": {
            "groups": people.get("n_groups"),
            "top_group": (people.get("rows") or [{}])[0].get("pi"),
            "top_group_share": people.get("top_share"),
            "active": people.get("span"),
        },
        # Which of the workstation-built artefacts this family actually has, so a consumer
        # can tell "no interactions" from "interactions not computed yet" rather than
        # reading an empty list as a finding.
        "artefacts": {
            "embedding": bool(fam_map.get("embedded")),
            "contacts": bool(fam.get("contacts")),
            "topology": bool(fam.get("topology")),
        },
        # Known problems with THIS family's data, so a consumer can tell an empty field that
        # means "nothing found" from one that means "we do not trust what we have".
        "warnings": ([] if contacts_positions_look_sane(fam) else [
            "contact positions for this family are not reliable seed coordinates "
            "(see codswallop.api.contacts_positions_look_sane); hot_residues withheld"
        ]),
    }
