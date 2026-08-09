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
import logging
import re
import time
from collections import Counter
from typing import Optional

from . import (config, constructs as construct_engine, crystals as crystal_engine, db,
               instances, layout, ligands as ligand_engine, msa as msa_engine, provenance,
               rcsb, uniprot)

# How many distinct UniProt references a family will fetch canonicals for. A family is
# dominated by a handful of accessions (carbonic anhydrase II: 1,249 of 1,490 entities are
# P00918, and 12 accessions cover almost all the rest), so this cap costs nothing real while
# stopping a diverse family from making a hundred UniProt calls.
MAX_REFERENCES = 40

logger = logging.getLogger(__name__)


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
    # Density BEFORE the statistics, not after: the coverage census reads `unobserved_seed`
    # off each member, and summarising first would compute the whole census against members
    # that do not have it yet and silently report "not measured" for every family.
    _attach_density(fam, members)
    fam["stats"] = summarise(fam, members)
    # Lay out against the identity range actually present, not the nominal 30 % floor. On a
    # family large enough to hit the cap, every surviving member can sit above 95 % identity,
    # and spreading those across the full radius is the difference between a readable map and
    # every node stacked on the centre point.

    # Crystallisation and validation summaries read the entry rows, and _compact drops those
    # once the derived views exist. Computed before it, not after.
    #
    # Each enrichment is optional. The family inventory is the page; a panel that could not
    # be built because an upstream API had a bad afternoon should cost the reader that panel,
    # not the whole family. Only the enrichments that make their own network calls need this
    # (ligands fetches chem_comp detail); the pure-parse ones are wrapped for symmetry so a
    # future change cannot quietly make one of them fatal.
    fam["crystals"] = _optional("crystals", lambda: crystal_engine.summarise(fam["entries"]),
                                {"n": 0, "n_parsed": 0})
    fam["quality"] = _optional("quality", lambda: build_quality(fam["entries"]),
                               {"n": 0, "rows": []})
    # Here, with the other entry-level panels, and not further down beside the orthologue
    # matrix: `_compact` pops `fam["entries"]` to shrink the payload, so anything reading it
    # after that point silently receives an empty list and reports a family with no
    # assemblies at all. That is the fourth bug in this function to come from its order.
    # Motifs need the coverage census, so this runs after `_attach_density` has put `depth`
    # and `seen` on the stats: without them every site reports an unknown share of the
    # family rather than being grounded in it, which is the whole point of the panel.
    fam["motifs"] = _optional("motifs", lambda: build_motifs(fam),
                              {"n": 0, "rows": [], "n_curated": 0, "n_predicted": 0})
    # KLIFS or GPCRdb, when the family is one of theirs. Separate from the motifs build so a
    # kinase whose reference numbering is unavailable still gets its sequence sites.
    fam["pocket"] = _optional("pocket", lambda: _build_pocket(fam), {"kind": None})
    fam["assemblies"] = _optional("assemblies", lambda: build_assemblies(fam["entries"]),
                                  {"n": 0, "states": [], "provenance": {}, "ambiguous": [],
                                   "n_ambiguous": 0, "interfaces": None})
    fam["ligands"] = _optional("ligands", lambda: ligand_engine.summarise(fam["entries"]),
                               {"components": [], "by_class": {}, "n": 0})
    # Also here, and for the same reason as assemblies: it reads the citation and deposit
    # date off `fam["entries"]`, which `_compact` is about to pop.
    fam["provenance"] = _optional(
        "provenance", lambda: provenance.build(fam, members, _seed_names(fam, members)),
        {"names": {"n": 0, "rows": []}, "people": {"n_groups": 0, "rows": []}})

    # The classification supersedes the Phase 1 exclusion list that decided the amber halo.
    # "Ligand-bound" now means a component somebody put there on purpose (a ligand or a
    # cofactor), not merely a component: with the old list, "how many structures are
    # ligand-bound" was close to "how many were frozen in glycerol".
    _klass = {c["id"]: c["klass"] for c in fam["ligands"]["components"]}
    for m in members:
        m["has_ligand_heuristic"] = m["has_ligand"]
        # This runs BEFORE _compact, so a member's ligands are still full dicts rather than
        # the bare id list the client eventually receives.
        m["has_ligand"] = any(
            _klass.get(((lig.get("id") if isinstance(lig, dict) else lig) or "").upper())
            in ligand_engine.COUNTS_AS_BOUND
            for lig in (m.get("ligands") or [])
        ) if _klass else m["has_ligand"]
    fam["stats"]["holo_entries_heuristic"] = fam["stats"]["holo_entries"]
    fam["stats"]["holo_entries"] = sum(1 for m in members if m["has_ligand"])

    _compact(fam, members)

    # Constructs need the deduplicated sequence table and the seq_id back-references that
    # _compact creates, so this runs after it. Cached on the exact set of sequences, because
    # the diff is ~0.4 s for a 232-construct family and nothing about it changes between
    # requests: a warm page must stay warm.
    seed_acc = (fam.get("seed") or "") if fam.get("kind") == "uniprot" else ""
    fam["constructs"] = db.cached(
        ("constructs", construct_engine.ENGINE_VERSION, seed_acc, sorted(fam["sequences"])),
        lambda: build_constructs(members, fam["sequences"], seed_acc),
    )
    _apply_constructs(fam, members)

    # Conservation, in the seed's own coordinate frame so a column here is the same column
    # as in the coverage census and the domain ribbon. Weighted by how many entities used
    # each construct, and cached with the constructs because it costs the same alignments.
    weights = {c["seq_id"]: c["n_entities"] for c in fam["constructs"]}
    fam["msa"] = _optional("msa", lambda: db.cached(
        ("msa", msa_engine.PARSE_VERSION, fam["slug"], sorted(fam["sequences"])),
        lambda: msa_engine.build(fam.get("seed_sequence") or "", fam["sequences"], weights),
    ), None)
    # The map goes LAST, after _compact has put seq_id on every member and the constructs
    # exist. Placed any earlier it silently falls back to the placeholder even when a real
    # embedding is present, because the representatives are keyed on seq_id and nothing
    # matches: the first version reported all 1,688 members as approximated.
    #
    # Three separate bugs in this function have now come from the same cause. The order here
    # is a dependency chain, not a list: density -> stats -> compaction -> constructs -> map.
    identities = [m["identity"] for m in members if m.get("identity") is not None]
    floor = min(identities) if identities else config.IDENTITY_MIN
    # The structural embedding when a workstation has computed one for this family, the
    # sequence-identity placeholder when it has not. Read through embed_io, which imports
    # nothing beyond the standard library, so the droplet never needs biotite or tmtools.
    # `annotations` explicitly: _compact has already lifted every member's Pfam list into
    # that lookup by the time the map is built, so the cluster labeller reading m["pfam"]
    # found nothing on every real request and silently fell through to naming each cluster
    # after whichever protein held a plurality of it.
    fam["map"] = layout.compute(members, min(floor, config.IDENTITY_MAX - 1),
                                embedding=_embedding_for(fam["slug"]),
                                annotations=fam.get("annotations"))

    fam["contacts"] = _optional("contacts", lambda: _contacts_for(fam["slug"]), None)
    fam["topology"] = _optional("topology", lambda: _topology_for(fam["slug"]), None)
    # After ligands: it reads the classified component list rather than the raw one, so a
    # buffer never gets a DrugBank lookup.
    fam["drugs"] = _optional("drugs", lambda: _drugs_for(fam),
                             {"n": 0, "classes": [], "drugs": []})

    # Ask a workstation for what this machine cannot build. Only the droplet ever needs
    # this, but it costs one indexed upsert and guessing which machine we are on from
    # inside the request path would be worse than doing it unconditionally.
    _queue_missing_artefacts(fam, len(members))

    # After the msa, because the promotion needs conservation as well as contacts.
    _optional("metals", lambda: _promote_metals(fam, members), None)

    fam["domains"] = _optional("domains", lambda: build_domains(fam, members),
                               {"sources": [], "domains": []})
    fam["orthologues"] = _optional(
        "orthologues", lambda: build_orthologue_matrix(members, fam.get("seed_length") or 0), [])
    fam.pop("domains_by_entity", None)      # per-chain detail: aggregated, not shipped
    return fam


def _seed_names(fam: dict, members: list[dict]) -> dict:
    """The UniProt record for the family's OWN seed, for the nomenclature panel.

    The accession is taken from the seed and never from a vote among the members. A family
    is assembled at 30 % identity, so the commonest accession in it is frequently a different
    protein: that popularity contest has already sent ABL1's AlphaFold overlay to EGFR and
    A2A's to a cytochrome, and here it would reconcile lysozyme's names against whatever
    relative happened to outnumber it.

    A UniProt-seeded family carries the accession directly. A PDB- or entity-seeded one
    carries an entity id, so the accession is read off that exact entity: `1AKI_1` gives
    P00698, which is a lookup rather than a guess.
    """
    accession = None
    if fam.get("kind") == "uniprot":
        accession = fam.get("seed")
    else:
        target = (fam.get("seed") or "").upper()
        for m in members:
            if (m.get("entity_id") or "").upper() == target and m.get("uniprot"):
                accession = m["uniprot"]
                break
    if not accession:
        # No canonical record to reconcile against. The panel still lists what the archive
        # actually says, and simply cannot mark any of it official.
        return {"name": fam.get("name"), "alt_names": [], "short_names": [], "genes": []}
    rec = uniprot.entry(accession) or {}
    return {
        "accession": accession,
        "name": rec.get("name") or fam.get("name"),
        "alt_names": rec.get("alt_names") or [],
        "short_names": rec.get("short_names") or [],
        "genes": rec.get("genes") or [],
    }


def _attach_density(fam: dict, members: list[dict]) -> None:
    """Fetch per-chain unobserved regions and map them onto seed coordinates.

    The API reports unobserved residues in the entity's own SEQRES numbering. Everything else
    in this app speaks seed coordinates, and the two differ by exactly the offset the sequence
    search already measured: an entity aligning from its residue `subject_beg` to the seed's
    `query_beg`. Mapping through that offset is what lets one census cover a family whose
    members have different tags, truncations and numbering.
    """
    ids = instances.representative_ids(members)
    if not ids:
        return
    try:
        data = instances.fetch(ids)
    except Exception:
        # Density is an enhancement, not the page. A family that still renders without it is
        # better than one that 502s because a second API had a bad afternoon.
        return

    domains_by_entity: dict[str, list] = {}
    for m in members:
        chains = m.get("chains") or []
        if not chains:
            continue
        rec = data.get(f"{m['pdb_id']}.{chains[0]}")
        if not rec:
            continue

        # SEQRES -> seed. The aligned block starts at seed `query_beg` and at construct
        # position `q_start`; RCSB gives the latter as 1-based, and the search hit carries
        # the former. Without the offset every truncated construct would report its disorder
        # shifted by the length of its own truncation.
        beg = m.get("query_beg")
        if beg is None:
            continue
        offset = beg - 1                      # seed position of SEQRES residue 1, 0-based
        m["unobserved_seed"] = [
            (u_start + offset, u_end + offset) for u_start, u_end in rec["unobserved"]
        ]
        m["chain_sampled"] = rec.get("chain")
        if rec["domains"]:
            domains_by_entity[m["entity_id"]] = rec["domains"]

    fam["domains_by_entity"] = domains_by_entity


def _embedding_for(slug: str) -> Optional[dict]:
    """Read the structural embedding artefact, if a workstation has produced one."""
    try:
        from .embed_io import load
        return load(slug)
    except Exception:
        logger.warning("could not read the embedding for %s", slug, exc_info=True)
        return None


def _drugs_for(fam: dict) -> dict:
    from . import drugs as drug_engine
    return drug_engine.build(fam)


def _topology_for(slug: str) -> Optional[dict]:
    """Read the topology artefact, if a workstation has produced one."""
    try:
        from .topology_io import load
        return load(slug)
    except Exception:                       # noqa: BLE001
        logger.warning("could not read the topology for %s", slug, exc_info=True)
        return None


def _queue_missing_artefacts(fam: dict, n_entries: int) -> None:
    """Record a family that is being served in its degraded form.

    The map falls back to a sequence-identity placeholder and the contacts panel disappears
    entirely, both without an error, because the artefacts they need are built on a
    workstation and rsynced across. Nothing outside a pre-warmed list ever asked for them,
    so a family a reader assembles is stuck that way until somebody happens to notice.
    """
    try:
        from . import artefacts
        st = artefacts.status(fam["slug"])
        want = [k for k in ("embedding", "contacts") if not st[k]["current"]]
        if not want:
            return
        db.request_artefact(
            slug=fam["slug"],
            query=fam.get("query") or fam.get("seed") or fam["slug"],
            kind="both" if len(want) == 2 else want[0],
            name=fam.get("name"),
            n_entries=n_entries,
        )
    except Exception:
        logger.warning("could not queue artefacts for %s", fam.get("slug"), exc_info=True)


def _promote_metals(fam: dict, members: list[dict]) -> None:
    """Let the interaction data decide whether this family's metals are catalytic.

    Retires the caveat the Ligands panel used to print. A metal coordinated by several
    conserved residues is the catalytic centre and counts as ligand-bound; one with a couple
    of surface contacts is crystallisation chemistry and does not. The evidence is per
    family, which is the only level at which the question has an answer.
    """
    lig = fam.get("ligands")
    contacts = fam.get("contacts")
    if not lig or not lig.get("components") or not contacts:
        return

    conservation = {}
    msa = fam.get("msa")
    if msa and msa.get("columns"):
        conservation = {c["pos"]: c["conservation"] for c in msa["columns"]
                        if c["conservation"] is not None}

    before = {c["id"]: c["klass"] for c in lig["components"]}
    ligand_engine.promote_catalytic_metals(lig["components"], contacts, conservation)
    promoted = [c["id"] for c in lig["components"]
                if c.get("promoted") and before.get(c["id"]) != c["klass"]]
    if not promoted:
        return

    # Recount what "ligand-bound" means now that the metals have moved.
    lig["by_class"] = {}
    for c in lig["components"]:
        lig["by_class"][c["klass"]] = lig["by_class"].get(c["klass"], 0) + 1
    lig["promoted_metals"] = promoted

    klass = {c["id"]: c["klass"] for c in lig["components"]}
    for m in members:
        m["has_ligand"] = any(
            klass.get((lid or "").upper()) in ligand_engine.COUNTS_AS_BOUND
            for lid in (m.get("ligands") or [])
        )
    fam["stats"]["holo_entries"] = sum(1 for m in members if m["has_ligand"])


def _contacts_for(slug: str) -> Optional[dict]:
    """The PLIP fingerprint, if a workstation has produced one. Imports nothing heavy."""
    from .contacts_io import load
    return load(slug)


def _optional(name: str, build, fallback):
    """Build one enrichment, or fall back rather than take the whole family down with it."""
    try:
        return build()
    except Exception:                       # noqa: BLE001 - deliberately broad
        logger.warning("optional panel %r failed; serving the family without it",
                       name, exc_info=True)
        return fallback


def _apply_constructs(fam: dict, members: list[dict]) -> None:
    """Push each construct's verdict back onto the members that use it, and correct the
    Phase 1 statistics it supersedes."""
    by_seq = {c["seq_id"]: c for c in fam["constructs"]}
    real_fusions = 0
    engineered = 0
    tagged = 0
    for m in members:
        c = by_seq.get(m.get("seq_id"))
        if not c:
            continue
        m["construct_summary"] = c["summary"]
        m["tags"] = c["tags"]
        m["fusions"] = c["fusions"]
        m["engineered"] = c["engineered"]
        m["mutation_count"] = c["mutation_count"]
        # The real answer replaces the Phase 1 length heuristic. Keep the old verdict under
        # its own name rather than silently overwriting it: the two disagreeing is the
        # interesting case, and it is how the heuristic's error rate was measured.
        m["is_fusion_heuristic"] = m["is_fusion"]
        m["is_fusion"] = bool(c["fusions"])
        if c["fusions"]:
            real_fusions += 1
        if c["engineered"]:
            engineered += 1
        if c["tags"]:
            tagged += 1

    s = fam["stats"]
    s["fusions_heuristic"] = s["fusions"]
    s["fusions"] = real_fusions
    s["engineered"] = engineered
    s["tagged"] = tagged
    s["constructs_engineered"] = sum(1 for c in fam["constructs"] if c["engineered"])
    s["constructs_unreferenced"] = sum(1 for c in fam["constructs"] if c["engineered"] is None)


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


def build_constructs(members: list[dict], sequences: dict[str, str],
                     seed_acc: str = "") -> list[dict]:
    """One row per unique deposited sequence, diffed against its own reference.

    The diff is run against **each construct's own UniProt canonical**, not against the
    family seed. Diffing everything against the seed looks reasonable and is wrong: a
    paralogue sitting at 45 % identity produces 150 "mutations", which is arithmetically
    true and useless as construct advice. The first version did exactly that and reported
    carbonic anhydrase I entries as hundred-mutation variants of carbonic anhydrase II.

    Entities with no UniProt cross-reference keep a row and say so, rather than being
    silently dropped or diffed against something they are not.
    """
    from collections import Counter

    by_seq: dict[str, list[dict]] = {}
    for m in members:
        if m.get("seq_id"):
            by_seq.setdefault(m["seq_id"], []).append(m)

    # Which references are worth fetching, most-used first. Counted over *every* accession an
    # entity carries, not just its first, so a family's own protein still dominates the count
    # when half its chimeras happen to list the fusion partner first.
    acc_counts: Counter = Counter()
    for m in members:
        for a in (m.get("uniprot_ids") or ([m["uniprot"]] if m.get("uniprot") else [])):
            acc_counts[a] += 1
    wanted = [a for a, _ in acc_counts.most_common(MAX_REFERENCES)]
    # The family's own protein: the reference a chimera must be diffed against. Getting it
    # from the family rather than from the entity is the whole point of being family-centric.
    #
    # The seed's accession outranks the modal one, because a modal vote loses the subject in
    # exactly the families where it matters most. Families are assembled at 30% identity, so
    # a search for ABL1 returns most of the tyrosine kinases: EGFR carried 332 entities to
    # ABL1's 85, and every ABL1 construct was being diffed against EGFR's canonical sequence.
    # A2A won its own vote by 189 to BRIL's 172, which is a coin flip, not a safety margin.
    dominant = acc_counts.most_common(1)[0][0] if acc_counts else None
    if seed_acc and acc_counts.get(seed_acc):
        dominant = seed_acc
    if seed_acc and seed_acc not in wanted and acc_counts.get(seed_acc):
        wanted.append(seed_acc)      # its sequence must be fetched even if rarely carried
    # Fetched concurrently: on a diverse family this was 36 sequential UniProt round trips
    # before anything else could start.
    from .http import parallel_map
    references: dict[str, dict] = {}
    for acc, (rec, feats) in zip(wanted, parallel_map(uniprot.entry_with_features, wanted)):
        if rec and rec.get("sequence"):
            references[acc] = {"sequence": rec["sequence"], "name": rec.get("name"),
                               "organism": rec.get("organism"), "length": rec.get("length"),
                               "features": feats or {}}

    rows: list[dict] = []
    for seq_id, users in by_seq.items():
        seq = sequences.get(seq_id)
        if not seq:
            continue
        # Which of this construct's accessions to diff against.
        #
        # The family's own protein wins whenever the construct carries it. A chimera is
        # cross-referenced to everything it is made of, and the order RCSB lists them in is
        # not meaningful: 2RH1, the beta-2 adrenergic receptor-T4 lysozyme fusion, lists
        # T4 lysozyme (164 aa) first. Taking the head of that list made the diff run
        # backwards, reporting the 500-residue chimera as a T4 lysozyme with the entire
        # receptor hanging off it as two unexplained overhangs, and no fusion at all.
        # Against the receptor, the same construct reads as a receptor with T4 lysozyme
        # replacing intracellular loop 3, which is what it is.
        candidates: list[str] = []
        for m in users:
            candidates.extend(m.get("uniprot_ids") or
                              ([m["uniprot"]] if m.get("uniprot") else []))
        acc = None
        if dominant and dominant in candidates:
            acc = dominant
        elif candidates:
            acc = Counter(candidates).most_common(1)[0][0]
        ref = references.get(acc)

        entry = {
            "seq_id": seq_id,
            "length": len(seq),
            "n_entities": len(users),
            "n_entries": len({m["pdb_id"] for m in users}),
            "pdb_ids": sorted({m["pdb_id"] for m in users})[:60],
            "entity_ids": [m["entity_id"] for m in users],
            "uniprot": acc,
            "reference_name": (ref or {}).get("name"),
            "organism": Counter(m.get("organism") for m in users).most_common(1)[0][0],
            "description": Counter(m.get("description") for m in users).most_common(1)[0][0],
        }

        resolutions = [m["resolution"] for m in users if m.get("resolution") is not None]
        entry["best_resolution"] = min(resolutions) if resolutions else None
        entry["median_resolution"] = (sorted(resolutions)[len(resolutions) // 2]
                                      if resolutions else None)
        entry["best_pdb_id"] = min(
            (m for m in users if m.get("resolution") is not None),
            key=lambda m: m["resolution"], default={}).get("pdb_id")
        entry["holo"] = sum(1 for m in users if m.get("has_ligand"))

        if ref:
            d = construct_engine.diff(ref["sequence"], seq, ref["features"])
            entry["diff"] = d
            entry["summary"] = construct_engine.summarise(d)
            entry["engineered"] = bool(d.get("is_engineered"))
            entry["tags"] = d.get("tags") or []
            entry["fusions"] = d.get("fusions") or []
            entry["proteases"] = d.get("proteases") or []
            entry["mutation_count"] = d.get("mutation_count") or 0
        else:
            entry["diff"] = None
            entry["summary"] = ("no UniProt reference for this entity, so it cannot be "
                                "diffed against a canonical sequence")
            entry["engineered"] = None
            entry["tags"], entry["fusions"], entry["proteases"] = [], [], []
            entry["mutation_count"] = 0

        rows.append(entry)

    # Most-used constructs first: "what did most people make" is the question this answers.
    rows.sort(key=lambda r: (-r["n_entities"], r["length"]))
    return rows


def build_quality(entries: list[dict]) -> dict:
    """A blunt triage of which entries to trust.

    Thresholds are the ones the wwPDB validation report itself flags on, not invented here:
    a clashscore over 20 and RSRZ outliers over 5 % are the standard "worse than most of the
    archive at this resolution" marks. Reported as a traffic light per entry plus the
    family's distribution, because a single entry's clashscore means little without knowing
    what the rest of the family managed.
    """
    rows = []
    for e in entries:
        v = e.get("validation")
        if not v:
            continue
        flags = []
        if (v.get("clashscore") or 0) > 20:
            flags.append("clashscore")
        if (v.get("rsrz_outliers") or 0) > 5:
            flags.append("RSRZ outliers")
        if (v.get("rota_outliers") or 0) > 5:
            flags.append("rotamer outliers")
        if (v.get("rama_outliers") or 0) > 0.5:
            flags.append("Ramachandran outliers")
        gap = None
        if e.get("r_free") is not None and e.get("r_work") is not None:
            gap = round(e["r_free"] - e["r_work"], 3)
            # A wide R-free/R-work gap is the classic overfitting tell.
            if gap > 0.07:
                flags.append("R-free gap")
        rows.append({
            "pdb_id": e["pdb_id"], "resolution": e.get("resolution"),
            "clashscore": v.get("clashscore"), "rsrz": v.get("rsrz_outliers"),
            "rama": v.get("rama_outliers"), "rota": v.get("rota_outliers"),
            "eds_r": v.get("eds_r"), "completeness": v.get("completeness"),
            "r_gap": gap, "has_sf": e.get("has_sf"),
            "flags": flags,
            "verdict": "poor" if len(flags) >= 2 else ("check" if flags else "ok"),
        })

    rows.sort(key=lambda r: (-len(r["flags"]), r.get("resolution") or 99))
    def med(key):
        vals = sorted(r[key] for r in rows if r.get(key) is not None)
        return vals[len(vals) // 2] if vals else None

    return {
        "n": len(rows),
        "ok": sum(1 for r in rows if r["verdict"] == "ok"),
        "check": sum(1 for r in rows if r["verdict"] == "check"),
        "poor": sum(1 for r in rows if r["verdict"] == "poor"),
        "with_sf": sum(1 for r in rows if r["has_sf"]),
        "median_clashscore": med("clashscore"),
        "median_rsrz": med("rsrz"),
        "median_r_gap": med("r_gap"),
        "rows": rows[:400],
    }


def build_domains(fam: dict, members: list[dict]) -> dict:
    """Consensus domain architecture, in seed coordinates.

    Domains are assigned per deposited chain, so a family holds hundreds of assignments of
    the same domain at slightly different boundaries. Reporting all of them is noise;
    reporting one is a guess. This reports each distinct domain once, with the *median*
    boundary across every chain that carries it and the number of chains that agreed, so a
    reader can see both where the domain is and how firmly the databases agree it is there.
    """
    by_entity = fam.get("domains_by_entity") or {}
    if not by_entity:
        return {"sources": [], "domains": []}

    offsets = {m["entity_id"]: (m.get("query_beg") or 1) - 1 for m in members}
    grouped: dict[tuple, dict] = {}
    for entity_id, doms in by_entity.items():
        off = offsets.get(entity_id, 0)
        for d in doms:
            key = (d["source"], d["id"] or d["name"])
            g = grouped.setdefault(key, {
                "source": d["source"], "id": d["id"], "name": d["name"],
                "starts": [], "ends": [], "n_chains": 0,
            })
            g["n_chains"] += 1
            for beg, end in d["spans"]:
                g["starts"].append(beg + off)
                g["ends"].append(end + off)

    seed_len = fam.get("seed_length") or 0
    # A domain has to be seen on more than one chain, and on a non-trivial share of them, to
    # be reported. Without this the list fills with domains belonging to *other proteins* in
    # the complexes: p53's family picked up "Green fluorescent protein", "Annexin" and
    # "Blc2-like", each on a single chain, because a partner chain in one entry carries them
    # and the seed-coordinate offset that is correct for p53 is meaningless for them.
    support = max(2, int(0.01 * max(1, len(members))))
    out = []
    dropped = 0
    for g in grouped.values():
        if not g["starts"] or g["n_chains"] < support:
            dropped += 1
            continue
        starts, ends = sorted(g["starts"]), sorted(g["ends"])
        start = starts[len(starts) // 2]
        end = ends[len(ends) // 2]
        # Independent medians of the two boundaries can invert on a domain whose assignments
        # disagree wildly (p53 produced an "Annexin" at 418-393). An inverted or
        # out-of-range span is evidence the mapping does not apply to this domain at all,
        # so it is dropped rather than clamped into something plausible-looking.
        if end <= start or start < 1 or (seed_len and start > seed_len):
            dropped += 1
            continue
        out.append({
            "source": g["source"], "id": g["id"], "name": g["name"],
            "start": start, "end": min(seed_len, end) if seed_len else end,
            "n_chains": g["n_chains"],
        })

    out.sort(key=lambda d: (d["source"], d["start"], -d["n_chains"]))
    sources = sorted({d["source"] for d in out})
    return {"sources": sources, "domains": out, "support": support,
            "dropped": dropped}


def _build_pocket(fam: dict) -> dict:
    from . import pockets
    return pockets.build(fam)


def build_motifs(fam: dict) -> dict:
    """Functional sites on the seed, cached on the seed rather than on the family.

    The scan is a property of the sequence, so two families seeded on the same protein at
    different identity thresholds share it. The grounding is not, so only the scan is cached.
    """
    from . import motifs as motif_engine

    seed = (fam.get("seed") or "").upper()
    feats = {}
    if fam.get("kind") == "uniprot" and seed:
        try:
            feats = uniprot.features(seed) or {}
        except Exception:                       # noqa: BLE001
            logger.warning("no UniProt features for %s", seed, exc_info=True)
    else:
        # A family seeded on a PDB entry rather than an accession used to get no curated
        # sites at all: `lysozyme-1aki-1` showed no active site and no disulphides while
        # `lysozyme-c-p00698`, the same protein, showed both. The guard was there for a real
        # reason, since UniProt's positions are on the CANONICAL sequence and a PDB-seeded
        # family's seed is whatever that entity was, typically the mature protein: hen
        # lysozyme is 129 residues against the canonical's 147, so the features are eighteen
        # out. Skipping them was safe and lost the panel; shifting them by a fixed 18 would
        # be the assumed-offset mistake that has already cost this codebase two panels.
        feats = _features_on_seed(fam)
    return motif_engine.build(fam, feats)


def _features_on_seed(fam: dict) -> dict:
    """UniProt's curated features, moved onto a seed that is not the canonical sequence.

    Aligned, never offset. The seed may be the mature protein, a single domain, or a
    construct numbered from 1, and each of those is a different shift; an alignment finds
    whichever it is, and finds nothing when the two are not the same protein.
    """
    members = fam.get("members") or []
    names = _seed_names(fam, members)
    accession = names.get("accession")
    seed_seq = fam.get("seed_sequence") or ""
    if not accession or not seed_seq:
        return {}
    try:
        rec, feats = uniprot.entry_with_features(accession)
    except Exception:                           # noqa: BLE001
        logger.warning("no UniProt features for %s", accession, exc_info=True)
        return {}
    canonical = (rec or {}).get("sequence") or ""
    if not canonical or not feats:
        return {}

    from .constructs import _ALIGNER, _alignable
    try:
        aln = _ALIGNER.align(_alignable(canonical), _alignable(seed_seq))[0]
    except Exception:                           # noqa: BLE001
        return {}
    # canonical position (1-based) -> seed position (1-based)
    mapping, matches, pairs = {}, 0, 0
    for (c0, c1), (s0, s1) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(c1 - c0):
            # int(), because Biopython's aligned blocks are numpy integers and a numpy int
            # is not JSON-serialisable: the panel would build correctly and then take the
            # whole family down at render time.
            mapping[int(c0) + k + 1] = int(s0) + k + 1
            pairs += 1
            if canonical[c0 + k] == seed_seq[s0 + k]:
                matches += 1
    # The same agreement floor the conservation colouring uses, and for the same reason: a
    # seed that is not this protein must come back with nothing rather than with sites moved
    # onto whichever residues the aligner happened to pair them with.
    if not pairs or matches / pairs < 0.6:
        return {}

    out: dict = {}
    for kind, positions in feats.items():
        moved = []
        for f in positions or []:
            # Two shapes, both from `uniprot._features_from`: a single-residue kind is a bare
            # integer, a ranged one is a dict carrying start, end and a description.
            if isinstance(f, dict):
                beg, end = mapping.get(f.get("start")), mapping.get(f.get("end"))
                # A feature only partly inside the seed is dropped whole. Half a signal
                # peptide is not a shorter signal peptide, and a disulphide with one of its
                # cysteines missing is not a bond.
                if beg and end:
                    moved.append({**f, "start": beg, "end": end})
            else:
                q = mapping.get(f)
                if q:
                    moved.append(q)
        out[kind] = moved
    return out


def build_assemblies(entries: list[dict]) -> dict:
    """What oligomeric state this family crystallises in, and how well attested it is.

    The family-level question is not "what is the assembly of 4HHB", which the RCSB already
    answers on its own page, but whether a protein deposited 1,400 times has ever been seen
    as anything other than the state everyone reports.

    Provenance is reported as three separate counts and not as an agreement rate, because
    `author_defined_assembly` says only that the depositor stated one: PISA may have returned
    nothing or never run. Treating it as disagreement would invent a conflict on 68 of 150
    thrombin entries. The ambiguity that is real, an entry whose own assemblies disagree
    about the count, is 4 of that 150, and it is reported as a list rather than a rate.
    """
    rows = [e["assembly"] for e in entries if e.get("assembly")]
    if not rows:
        return {"n": 0, "states": [], "provenance": {}, "ambiguous": [],
                "n_ambiguous": 0, "interfaces": None}

    # Keyed on the chain count, never on the wording. `oligomeric_details` is free text whose
    # capitalisation is not consistent across the archive, so lysozyme reported "trimeric" on
    # 64 entries and "Trimeric" on 3 as two different oligomeric states sitting in the same
    # table. The count is the fact; the wording is a label for it.
    states: Counter = Counter()
    labels: dict[int, Counter] = {}
    for a in rows:
        n = a.get("count")
        if not n:
            continue
        states[n] += 1
        if a.get("details"):
            labels.setdefault(n, Counter())[a["details"].strip().lower()] += 1

    prov = Counter(a.get("provenance") for a in rows if a.get("provenance"))

    ambiguous = [
        {"pdb_id": e["pdb_id"], "count": e["assembly"]["count"],
         "alternatives": e["assembly"].get("alternatives") or []}
        for e in entries
        if e.get("assembly") and e["assembly"].get("ambiguous")
    ]
    ambiguous.sort(key=lambda r: r["pdb_id"])

    areas = sorted(a["buried_area"] for a in rows if a.get("buried_area"))
    ifres = sorted(a["interface_residues"] for a in rows if a.get("interface_residues"))
    interfaces = None
    if areas:
        # Quartiles rather than min and max. Buried area scales with the whole assembly, so
        # one 60-mer sets a maximum three orders of magnitude above the median (lysozyme:
        # median 1,476 A^2, maximum 615,514) and a min-to-max range says nothing about the
        # family it is supposed to describe.
        def q(xs, f):
            return round(xs[min(int(f * len(xs)), len(xs) - 1)], 1)
        interfaces = {
            "n": len(areas),
            "median_area": q(areas, 0.5),
            "q1_area": q(areas, 0.25),
            "q3_area": q(areas, 0.75),
            "median_residues": ifres[len(ifres) // 2] if ifres else None,
        }

    total = sum(states.values())
    return {
        "n": len(rows),
        "states": [
            {"count": c, "entries": n,
             "details": labels[c].most_common(1)[0][0] if labels.get(c) else f"{c}-meric",
             "fraction": round(100.0 * n / total, 1) if total else 0.0}
            for c, n in states.most_common(12)
        ],
        "provenance": {
            "both": prov.get("author_and_software_defined_assembly", 0),
            "author": prov.get("author_defined_assembly", 0),
            "software": prov.get("software_defined_assembly", 0),
        },
        "ambiguous": ambiguous[:40],
        "n_ambiguous": len(ambiguous),
        "interfaces": interfaces,
    }


def build_orthologue_matrix(members: list[dict], seed_length: int) -> list[dict]:
    """Which organisms have structures, at what coverage and what quality.

    The question this answers is the one asked at the start of a project: somebody else has
    probably solved this in a more tractable organism, and if their construct covers the
    region you care about, that is where to start.
    """
    from collections import defaultdict
    rows: dict[str, dict] = defaultdict(
        lambda: {"entities": 0, "entries": set(), "resolutions": [], "covered": set(),
                 "accessions": set(), "holo": 0})
    for m in members:
        org = m.get("organism") or "Unknown"
        r = rows[org]
        r["entities"] += 1
        r["entries"].add(m["pdb_id"])
        if m.get("resolution") is not None:
            r["resolutions"].append(m["resolution"])
        if m.get("uniprot"):
            r["accessions"].add(m["uniprot"])
        if m.get("has_ligand"):
            r["holo"] += 1
        beg, end = m.get("query_beg"), m.get("query_end")
        if beg and end:
            r["covered"].update(range(max(1, beg), min(seed_length or end, end) + 1))

    out = []
    for org, r in rows.items():
        res = sorted(r["resolutions"])
        out.append({
            "organism": org,
            "entities": r["entities"],
            "entries": len(r["entries"]),
            "best_resolution": res[0] if res else None,
            "median_resolution": res[len(res) // 2] if res else None,
            "coverage_pct": round(100 * len(r["covered"]) / seed_length, 1) if seed_length else None,
            "accessions": sorted(r["accessions"])[:4],
            "holo": r["holo"],
        })
    out.sort(key=lambda x: -x["entities"])
    return out


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

# A residue resolved in under this fraction of the constructs that contained it is
# "rarely resolved": present in the crystal, absent from the density.
RARELY_RESOLVED = 0.25


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

    # ---- the density layer, where it is available -----------------------------------
    # How many constructs actually resolved each residue, as opposed to merely containing
    # it. Kept as a separate array on the same axis rather than folded into `depth`: they
    # are different questions and a reader who conflates them will design the wrong
    # construct. Seen depth can never exceed construct depth, and the gap between the two
    # curves is exactly the disorder.
    seen = [0] * (n + 1)
    have_density = False
    for m in members:
        beg, end = m.get("query_beg"), m.get("query_end")
        unobs = m.get("unobserved_seed")
        if not beg or not end or unobs is None:
            continue
        have_density = True
        blocked = set()
        for u_start, u_end in unobs:
            blocked.update(range(max(1, u_start), min(n, u_end) + 1))
        # Counted per residue directly. A difference array would be faster but the
        # unobserved set fragments the span into arbitrarily many pieces, and the
        # bookkeeping to express that as deltas is where an off-by-one lives.
        for i in range(max(1, beg), min(n, end) + 1):
            if i not in blocked:
                seen[i] += 1

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
    # "Resolved in no construct at all" is as strict a bar as the binary union was, and fails
    # the same way: with 267 p53 constructs, something resolved almost every residue once, so
    # the figure reads 0 % for a protein a third of which is famously disordered.
    #
    # The informative measure is the RATIO: of the constructs that contained this residue,
    # what fraction actually resolved it. p53 residue 63 sits in 14 constructs and is resolved
    # in 2; spike 1265 sits in 140 and is resolved in 6. That is the disorder, and it is what
    # the flagged figure in the header reports.
    never_seen = disorder = rarely = None
    worst = []
    if have_density:
        never_seen = sum(1 for i in range(1, n + 1) if depth[i] > 0 and seen[i] == 0)
        ratios = [(seen[i] / depth[i]) if depth[i] else None for i in range(n + 1)]
        rarely = sum(1 for i in range(1, n + 1)
                     if ratios[i] is not None and ratios[i] < RARELY_RESOLVED)
        disorder = [round(r, 3) if r is not None else None for r in ratios[1:]]
        # The runs a construct designer should know about: stretches resolved in under a
        # quarter of the constructs that contained them, longest first.
        run_start = None
        for i in range(1, n + 2):
            poor = i <= n and ratios[i] is not None and ratios[i] < RARELY_RESOLVED
            if poor and run_start is None:
                run_start = i
            elif not poor and run_start is not None:
                worst.append({"start": run_start, "end": i - 1, "length": i - run_start,
                              "resolved_fraction": round(
                                  min(ratios[j] for j in range(run_start, i)), 2)})
                run_start = None
        worst.sort(key=lambda w: -w["length"])
        worst = worst[:12]

    return {
        "length": n,
        # One value per residue. A 1,273-residue seed is 1,273 small integers: a few kB of
        # JSON, and the header sparkline and the stacked plot both read it directly.
        "depth": depth[1:],
        # The second curve: how many constructs actually RESOLVED each residue. None when no
        # per-chain data was fetched, so the UI can say "not measured" rather than draw a
        # flat zero and let a reader read it as total disorder.
        "seen": seen[1:] if have_density else None,
        "disorder": disorder,
        "never_seen": never_seen,
        "never_seen_pct": round(100 * never_seen / n, 1) if never_seen is not None else None,
        "rarely_resolved": rarely,
        "rarely_resolved_pct": round(100 * rarely / n, 1) if rarely is not None else None,
        "rarely_cut": RARELY_RESOLVED,
        "disorder_runs": worst,
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
