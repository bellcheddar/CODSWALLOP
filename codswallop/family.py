"""Assemble a family: seed in, complete inventory out.

The build runs once, at the *lowest* identity threshold the app offers, and stores every
member with its own identity to the seed. The threshold slider, the orthologue toggle and
the chimera toggle then filter that stored set in the browser, instantly, with no rebuild.
Membership is therefore one cached artefact rather than one per slider position, and moving
the slider costs nothing: the panels all re-filter from the same selection state.
"""

from __future__ import annotations

import hashlib
import json
import re
import time
from collections import Counter
from typing import Optional

from . import config, db, layout, rcsb


def slugify(text: str, fallback: str = "family") -> str:
    s = re.sub(r"[^a-z0-9]+", "-", (text or "").lower()).strip("-")
    return s[:48].strip("-") or fallback


def slug_for(seed: dict) -> str:
    """A stable, readable URL for a family.

    Keyed on the seed alone, deliberately: the filters are view state, not identity, so
    `/f/lysozyme-c-p00698` is one page whose controls change what it shows rather than a
    different page per slider position.
    """
    name = slugify(seed.get("name") or "")
    ident = seed.get("seed") or ("seq-" + hashlib.sha1(seed["sequence"].encode()).hexdigest()[:8])
    ident = slugify(ident, "seed")
    return f"{name}-{ident}".strip("-") if name else ident


# --------------------------------------------------------------------------------------
# Build
# --------------------------------------------------------------------------------------
def build(seed: dict, query: str) -> dict:
    """Search, fetch, classify and lay out. Returns the family as stored."""
    hits, total = rcsb.sequence_search(seed["sequence"], config.IDENTITY_MIN)
    if not hits:
        raise ValueError("No structures in the PDB match that sequence at 30 % identity.")

    truncated = total > len(hits)
    by_id = {h["id"]: h for h in hits}

    raw = rcsb.fetch_entities(list(by_id.keys()))
    parsed = [rcsb.parse_entity(r) for r in raw]

    seed_organism = (seed.get("organism") or "").strip().lower()
    seed_len = seed.get("length") or len(seed["sequence"])

    entities: list[dict] = []
    entries: dict[str, dict] = {}
    for p in parsed:
        ent, entry = p["entity"], p["entry"]
        hit = by_id.get(ent["entity_id"], {})

        ent["identity"] = hit.get("identity")
        ent["aligned_length"] = hit.get("aligned_length")
        ent["query_beg"] = hit.get("query_beg")
        ent["query_end"] = hit.get("query_end")

        # An entity from a different source organism than the seed. With a pasted sequence
        # there is no seed organism, so nothing can be called an orthologue and the filter
        # is disabled in the UI rather than silently excluding everything.
        ent["is_orthologue"] = bool(
            seed_organism and (ent.get("organism") or "").strip().lower() != seed_organism
        )

        # Blunt Phase 1 fusion/chimera test: the deposited construct is much longer than the
        # part of it that aligns to the seed. Every common fusion partner (MBP, GST, SUMO,
        # BRIL, T4 lysozyme) clears the threshold comfortably. Phase 2's construct diff
        # engine replaces this with a real alignment that names the partner.
        seq_len = ent.get("seq_length") or 0
        aligned = ent.get("aligned_length") or 0
        ent["is_fusion"] = bool(seq_len and aligned and seq_len - aligned >= config.FUSION_EXCESS_RESIDUES)

        entities.append(ent)
        entries.setdefault(entry["pdb_id"], entry)

    fam = {
        "slug": slug_for(seed),
        "query": query,
        "kind": seed["kind"],
        "seed": seed.get("seed"),
        "name": seed.get("name") or f"{seed_len}-residue query sequence",
        "organism": seed.get("organism"),
        "seed_sequence": seed["sequence"],
        "seed_length": seed_len,
        "pfam": seed.get("pfam") or [],
        "interpro": seed.get("interpro") or [],
        "identity_threshold": config.IDENTITY_DEFAULT,
        "total_hits": total,
        "truncated": truncated,
    }
    db.save_family(fam, list(entries.values()), entities)
    return fam


def get_or_build(seed: dict, query: str, force: bool = False) -> dict:
    """Return the assembled family, building it only if there is no fresh copy."""
    slug = slug_for(seed)
    if force or not db.family_fresh(slug):
        build(seed, query)
    fam = db.load_family(slug)
    if fam is None:                       # a build that wrote nothing: treat as a hard fail
        raise ValueError("The family could not be assembled.")
    return decorate(fam)


# --------------------------------------------------------------------------------------
# Derived views
# --------------------------------------------------------------------------------------
def decorate(fam: dict) -> dict:
    """Attach the derived pieces the UI needs: per-member view rows, stats and the map."""
    entries = {e["pdb_id"]: e for e in fam["entries"]}

    members = []
    for ent in fam["entities"]:
        entry = entries.get(ent["pdb_id"], {})
        members.append({
            **ent,
            "title": entry.get("title"),
            "method": entry.get("method"),
            "resolution": entry.get("resolution"),
            "r_work": entry.get("r_work"),
            "r_free": entry.get("r_free"),
            "space_group": entry.get("space_group"),
            "cell": entry.get("cell"),
            "deposit_date": entry.get("deposit_date"),
            "release_date": entry.get("release_date"),
            "chain_count": entry.get("chain_count"),
            "assembly_count": entry.get("assembly_count"),
            "ligands": entry.get("ligands") or [],
            "citation": entry.get("citation"),
            # "Holo" here means the entry has at least one non-polymer component that is not
            # obviously crystallisation chemistry. Phase 3 classifies these properly, ligand
            # by ligand, against the full buffer/cryoprotectant/ion list; this is the
            # cheap version that gets the amber halo roughly right on day one.
            "has_ligand": any(
                (lig.get("id") or "").upper() not in _NON_LIGANDS
                for lig in (entry.get("ligands") or [])
            ),
        })

    members.sort(key=lambda m: (m.get("resolution") is None, m.get("resolution") or 999, m["entity_id"]))
    fam["members"] = members
    fam["stats"] = summarise(fam, members)
    # Lay out against the identity range actually present, not the nominal 30 % floor. On a
    # family large enough to hit the cap, every surviving member can sit above 95 % identity,
    # and spreading those across the full radius is the difference between a readable map and
    # every node stacked on the centre point.
    identities = [m["identity"] for m in members if m.get("identity") is not None]
    floor = min(identities) if identities else config.IDENTITY_MIN
    fam["map"] = layout.compute(members, min(floor, config.IDENTITY_MAX - 1))

    _compact(fam, members)
    return fam


def _compact(fam: dict, members: list[dict]) -> None:
    """Shrink the payload before it goes over the wire, in place.

    Two easy wins, both worth taking: a 1,490-entry family was shipping a 6.3 MB page.

    1. `entries` and `entities` are the two halves that `members` already merges, so all
       three going out means sending everything three times. The client only ever reads
       `members`, so the sources are dropped once the derived views are computed.
    2. Sequences are deduplicated into a lookup table. Carbonic anhydrase II has 1,490
       entities but only 232 distinct deposited sequences, and that ratio is the norm rather
       than the exception: a family is mostly the same construct solved again.
    """
    sequences: dict[str, str] = {}
    for m in members:
        seq = m.pop("sequence", None)
        if not seq:
            m["seq_id"] = None
            continue
        key = hashlib.sha1(seq.encode()).hexdigest()[:10]
        sequences.setdefault(key, seq)
        m["seq_id"] = key

    # Primary citations repeat too: one paper typically covers a whole series of entries.
    citations: dict[str, dict] = {}
    for m in members:
        cit = m.pop("citation", None)
        if not cit:
            m["cite_id"] = None
            continue
        key = (cit.get("doi") or cit.get("title") or "")[:120].lower()
        key = hashlib.sha1(key.encode()).hexdigest()[:10] if key else None
        if key:
            citations.setdefault(key, cit)
        m["cite_id"] = key

    # Chemical components: the member keeps the CCD ids it contains, and the names and
    # formulae move to one family-level dictionary. ZN appearing in 900 entries was shipping
    # its name 900 times.
    components: dict[str, dict] = {}
    for m in members:
        ids = []
        for lig in m.pop("ligands", None) or []:
            cid = (lig.get("id") or "").upper()
            if not cid:
                continue
            components.setdefault(cid, {"name": lig.get("name"), "formula": lig.get("formula")})
            ids.append(cid)
        m["ligands"] = ids

    # Domain annotations: every member of a family carries near-identical Pfam and InterPro
    # lists, so the distinct *sets* go in a lookup and members hold a key. This was the
    # single heaviest field in the payload.
    annotations: dict[str, dict] = {}
    for m in members:
        pf, ip = m.pop("pfam", None) or [], m.pop("interpro", None) or []
        if not pf and not ip:
            m["annot_id"] = None
            continue
        blob = json.dumps({"pfam": pf, "interpro": ip}, sort_keys=True, separators=(",", ":"))
        key = hashlib.sha1(blob.encode()).hexdigest()[:10]
        annotations.setdefault(key, {"pfam": pf, "interpro": ip})
        m["annot_id"] = key

    fam["sequences"] = sequences
    fam["citations"] = citations
    fam["components"] = components
    fam["annotations"] = annotations
    fam.pop("entries", None)
    fam.pop("entities", None)


# Components that are almost always crystallisation chemistry rather than a ligand anyone
# meant to be there: waters, common ions, buffers, cryoprotectants and precipitants. Phase 3
# turns this into a proper classification panel with a per-component verdict; here it exists
# only so that "has a ligand" does not mean "was frozen in glycerol".
_NON_LIGANDS = {
    "HOH", "DOD", "SO4", "PO4", "CL", "NA", "K", "MG", "CA", "ZN", "MN", "FE", "NI", "CO",
    "CU", "CD", "HG", "BR", "IOD", "ACT", "EDO", "GOL", "PEG", "PG4", "PGE", "1PE", "2PE",
    "MPD", "DMS", "TRS", "MES", "EPE", "IMD", "FMT", "CIT", "TAR", "MLA", "ACY", "NH4",
    "CAC", "BME", "DTT", "TCE", "AZI", "NO3", "CO3", "F", "LI", "CS", "RB", "SR", "BA",
}


def summarise(fam: dict, members: list[dict]) -> dict:
    """The header stat strip, plus the breakdowns the filter chips are built from."""
    entries = fam["entries"]
    resolutions = sorted(e["resolution"] for e in entries if e.get("resolution") is not None)
    dates = sorted(e["release_date"] for e in entries if e.get("release_date"))

    methods = Counter(e.get("method") or "Unknown" for e in entries)
    organisms = Counter(m.get("organism") or "Unknown" for m in members)
    hosts = Counter(m.get("host_organism") or "Unknown" for m in members)
    space_groups = Counter(e["space_group"] for e in entries if e.get("space_group"))

    ligand_counts: Counter = Counter()
    ligand_names: dict[str, str] = {}
    for e in entries:
        for lig in e.get("ligands") or []:
            cid = (lig.get("id") or "").upper()
            if cid and cid not in _NON_LIGANDS:
                ligand_counts[cid] += 1
                ligand_names.setdefault(cid, lig.get("name") or cid)

    # Distinct deposited sequences: the honest Phase 1 answer to "how many different things
    # did people actually make?". Phase 2 refines it into named constructs (tag, cleavage
    # scar, fusion partner, point mutation), which is where it stops being a count and
    # starts being advice.
    constructs = len({m["sequence"] for m in members if m.get("sequence")})

    # The identity range the family actually spans. The threshold slider bounds itself to
    # this rather than to the nominal 30-100 %, so it never offers a range that would filter
    # out everything or nothing across most of its travel.
    identities = sorted(m["identity"] for m in members if m.get("identity") is not None)

    return {
        "identity_min": identities[0] if identities else None,
        "identity_max": identities[-1] if identities else None,
        "entries": len(entries),
        "entities": len(members),
        "constructs": constructs,
        "organisms": len([o for o in organisms if o != "Unknown"]),
        "ligands": len(ligand_counts),
        "holo_entries": sum(1 for m in members if m.get("has_ligand")),
        "best_resolution": resolutions[0] if resolutions else None,
        "median_resolution": resolutions[len(resolutions) // 2] if resolutions else None,
        "first_release": dates[0] if dates else None,
        "latest_release": dates[-1] if dates else None,
        "coverage": coverage_census(fam, members),
        "methods": methods.most_common(),
        "top_organisms": organisms.most_common(12),
        "top_hosts": hosts.most_common(8),
        "top_ligands": [
            {"id": cid, "name": ligand_names[cid], "count": n}
            for cid, n in ligand_counts.most_common(20)
        ],
        "top_space_groups": space_groups.most_common(10),
        "fusions": sum(1 for m in members if m.get("is_fusion")),
        "orthologues": sum(1 for m in members if m.get("is_orthologue")),
    }


# A residue covered by fewer than this fraction of the family's constructs is "thin": it
# exists in somebody's construct, but hardly anybody's.
THIN_FRACTION = 0.05


def coverage_census(fam: dict, members: list[dict]) -> dict:
    """How many of the family's constructs contain each residue of the seed.

    A **depth** profile, not a binary union, and the difference matters. The union answers
    "has anyone ever put this residue in a construct?", which for any well-studied protein is
    yes everywhere: 140 of the 2,000 spike entities carry the full 1,273 residues, so the
    union reports perfect coverage for a protein whose cytoplasmic tail is essentially never
    studied. The depth profile says the useful thing instead: that the tail appears in 7 % of
    constructs and the receptor-binding domain in nearly all of them.

    Free from the alignment spans the sequence search already returned, so it needs no
    coordinates. Phase 2 layers the harder question on top of the same axis: of the residues
    that were *in* a construct, which has nobody ever seen density for. The two are different
    questions and the UI labels them as such rather than letting one stand in for the other.
    """
    n = fam.get("seed_length") or 0
    spans = [(m["query_beg"], m["query_end"]) for m in members
             if m.get("query_beg") and m.get("query_end")]
    if not n or not spans:
        return {"length": n, "depth": [], "uncovered": 0, "pct": 0.0, "thin": 0,
                "max_depth": 0, "median_coverage": None, "gaps": []}

    # Difference array: O(members + length) rather than O(members x length), which matters
    # when a 1,273-residue seed meets 2,000 members.
    delta = [0] * (n + 2)
    for beg, end in spans:
        delta[max(1, beg)] += 1
        delta[min(n, end) + 1] -= 1
    depth, running = [0] * (n + 1), 0
    for i in range(1, n + 1):
        running += delta[i]
        depth[i] = running

    total = len(spans)
    thin_cut = max(1, int(total * THIN_FRACTION))
    uncovered = sum(1 for i in range(1, n + 1) if depth[i] == 0)
    thin = sum(1 for i in range(1, n + 1) if 0 < depth[i] < thin_cut)

    # Runs of thin-or-absent coverage: the regions a construct designer should know about.
    gaps, start = [], None
    for i in range(1, n + 2):
        poor = i <= n and depth[i] < thin_cut
        if poor and start is None:
            start = i
        elif not poor and start is not None:
            gaps.append({"start": start, "end": i - 1,
                         "depth": min(depth[start:i]), "length": i - start})
            start = None

    coverages = sorted((min(n, e) - max(1, b) + 1) / n for b, e in spans)
    return {
        "length": n,
        # One value per residue. A 1,273-residue seed is 1,273 small integers: a few kB of
        # JSON, and the header sparkline and Phase 2's stacked plot both read it directly.
        "depth": depth[1:],
        "max_depth": max(depth[1:]),
        "uncovered": uncovered,
        "pct": round(100 * uncovered / n, 1),
        "thin": thin,
        "thin_pct": round(100 * thin / n, 1),
        "thin_cut": thin_cut,
        # The median construct covers this fraction of the seed. For a multi-domain protein
        # this is the number that says "everybody solves one piece of it".
        "median_coverage": round(100 * coverages[len(coverages) // 2], 1),
        # Longest first: a 40-residue region almost nobody expresses is the headline, and a
        # scatter of single thin residues at the termini is not.
        "gaps": sorted(gaps, key=lambda g: -g["length"])[:12],
    }
