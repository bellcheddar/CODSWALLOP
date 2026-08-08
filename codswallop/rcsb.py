"""RCSB Search API v2 and Data GraphQL clients.

Every call routes through `db.cached`, and every search is parsed down to the fields we
actually use *before* it is cached: a verbose sequence search carries both aligned sequences
for every hit, which for a large family is megabytes of text we would otherwise store, read
and parse forever to get at one identity number.
"""

from __future__ import annotations

import re
from typing import Any, Iterable, Optional

from . import config, db, http

# --------------------------------------------------------------------------------------
# Search API v2
# --------------------------------------------------------------------------------------
# RCSB caps a single page; 1000 is comfortably inside it and keeps each cached response
# small enough that one flaky request is cheap to retry.
_PAGE = 1000

# Bumped whenever the *shape* of a parsed/cached result changes. It is part of every cache
# key, so old entries are orphaned rather than served.
#
# This is not bookkeeping for its own sake. Because searches are cached in parsed form, a
# new field added to the parser is simply absent from every cached row, and the code that
# consumes it fails silently: adding the alignment spans without bumping this left the
# "residues in no construct" figure reading 100 % and the fusion count reading 0, both of
# which look like plausible answers rather than missing inputs.
PARSE_VERSION = 5


def _paged_search(query: dict, return_type: str, limit: int, verbose: bool = False) -> tuple[list[dict], int]:
    """Run a search, following pagination up to `limit` hits.

    Returns (hits, total_count). Each hit is `{"id": ..., "score": ...}`, plus the sequence
    service's alignment numbers when `verbose` is set.
    """
    hits: list[dict] = []
    total = 0
    start = 0
    while start < limit:
        options: dict[str, Any] = {
            "paginate": {"start": start, "rows": min(_PAGE, limit - start)},
            "results_content_type": ["experimental"],
        }
        if verbose:
            options["results_verbosity"] = "verbose"
        body = http.post_search(
            config.RCSB_SEARCH_URL,
            {"query": query, "return_type": return_type, "request_options": options},
        )
        # 204: the search ran and matched nothing. A real answer, not an error.
        if not body:
            break
        total = body.get("total_count", 0)
        page = body.get("result_set") or []
        if not page:
            break
        for r in page:
            hit = {"id": r["identifier"], "score": r.get("score")}
            if verbose:
                hit.update(_alignment_of(r))
            hits.append(hit)
        start += len(page)
        if start >= total:
            break
    return hits, total


def _alignment_of(result: dict) -> dict:
    """Pull the best alignment out of a verbose sequence-service result.

    A hit can carry several match_contexts (repeats, multi-domain matches). We take the one
    with the highest bitscore, which is the alignment the hit was actually ranked on.
    """
    best = None
    for service in result.get("services") or []:
        if service.get("service_type") != "sequence":
            continue
        for node in service.get("nodes") or []:
            for ctx in node.get("match_context") or []:
                if best is None or (ctx.get("bitscore") or 0) > (best.get("bitscore") or 0):
                    best = ctx
    if not best:
        return {}
    return {
        # sequence_identity comes back as a fraction; percent is what the UI and the
        # threshold slider speak.
        "identity": round(100 * best["sequence_identity"], 1) if best.get("sequence_identity") is not None else None,
        "aligned_length": best.get("alignment_length"),
        "evalue": best.get("evalue"),
        "subject_length": best.get("subject_length"),
        # The span of the *seed* this member aligns to. Kept per hit rather than reduced to
        # a coverage percentage, because the union of these spans across the whole family is
        # what answers "which residues of this protein has nobody ever put in a construct?".
        "query_beg": best.get("query_beg"),
        "query_end": best.get("query_end"),
        "query_length": best.get("query_length"),
        "query_coverage": (
            round(100 * (best["query_end"] - best["query_beg"] + 1) / best["query_length"], 1)
            if best.get("query_length") else None
        ),
    }


def sequence_search(sequence: str, identity_pct: int, limit: Optional[int] = None) -> tuple[list[dict], int]:
    """Every polymer entity in the PDB above `identity_pct` identity to `sequence`.

    This is the single definition of family membership in CODSWALLOP: whatever the user
    typed, the resolver turns it into a seed sequence and this call turns that seed into a
    family. One consequence worth stating plainly in the UI: every member carries a real
    identity number to the same seed, so the threshold slider means something exact rather
    than "whatever the annotation happened to say".
    """
    limit = limit or config.MAX_FAMILY_ENTITIES
    query = {
        "type": "terminal",
        "service": "sequence",
        "parameters": {
            "sequence_type": "protein",
            "value": sequence,
            "identity_cutoff": identity_pct / 100.0,
            "evalue_cutoff": 1.0,
        },
    }
    return db.cached(
        ("seqsearch", PARSE_VERSION, sequence, identity_pct, limit),
        lambda: _paged_search(query, "polymer_entity", limit, verbose=True),
    )


def entities_by_uniprot(accession: str, limit: int = 50) -> list[str]:
    """Polymer entities whose reference sequence is this UniProt accession."""
    query = {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers"
                         ".reference_sequence_identifiers.database_accession",
            "operator": "exact_match",
            "value": accession,
        },
    }
    hits, _ = db.cached(
        ("uniprot_entities", PARSE_VERSION, accession, limit),
        lambda: _paged_search(query, "polymer_entity", limit),
    )
    return [h["id"] for h in hits]


def entities_by_annotation(annotation_id: str, limit: int = 50) -> list[str]:
    """Polymer entities carrying a given Pfam or InterPro annotation."""
    query = {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": "rcsb_polymer_entity_annotation.annotation_id",
            "operator": "exact_match",
            "value": annotation_id,
        },
    }
    hits, total = db.cached(
        ("annotation_entities", PARSE_VERSION, annotation_id, limit),
        lambda: _paged_search(query, "polymer_entity", limit),
    )
    return [h["id"] for h in hits]


def entities_by_text(text: str, limit: int = 25) -> list[str]:
    """Polymer entities matching a free-text query, best first."""
    query = {"type": "terminal", "service": "full_text", "parameters": {"value": text}}
    hits, _ = db.cached(
        ("text_entities", PARSE_VERSION, text, limit),
        lambda: _paged_search(query, "polymer_entity", limit),
    )
    return [h["id"] for h in hits]


def entities_of_entry(pdb_id: str) -> list[str]:
    """Polymer entity ids belonging to one PDB entry."""
    query = {
        "type": "terminal",
        "service": "text",
        "parameters": {"attribute": "rcsb_entry_container_identifiers.entry_id",
                       "operator": "exact_match", "value": pdb_id.upper()},
    }
    hits, _ = db.cached(
        ("entry_entities", PARSE_VERSION, pdb_id.upper()),
        lambda: _paged_search(query, "polymer_entity", 50),
    )
    return [h["id"] for h in hits]


# --------------------------------------------------------------------------------------
# Data API (GraphQL)
# --------------------------------------------------------------------------------------
# Everything Phase 1 promises per entry, in one round trip per 50 entities: identity,
# description, sequence, organism, expression host, UniProt cross-reference, domain
# annotations, and the whole entry-level block (method, resolution, R factors, cell,
# dates, assemblies, ligands, primary citation).
_ENTITY_QUERY = """
query($ids: [String!]!) {
  polymer_entities(entity_ids: $ids) {
    rcsb_id
    rcsb_polymer_entity { pdbx_description formula_weight }
    entity_poly { pdbx_seq_one_letter_code_can rcsb_sample_sequence_length type }
    rcsb_entity_source_organism { ncbi_scientific_name ncbi_taxonomy_id }
    rcsb_entity_host_organism { ncbi_scientific_name }
    rcsb_polymer_entity_container_identifiers {
      auth_asym_ids
      reference_sequence_identifiers { database_accession database_name }
    }
    rcsb_polymer_entity_annotation { type annotation_id name }
    entry {
      rcsb_id
      struct { title }
      exptl { method }
      rcsb_entry_info {
        resolution_combined experimental_method polymer_entity_count
        deposited_polymer_entity_instance_count
      }
      rcsb_accession_info { deposit_date initial_release_date has_released_experimental_data }
      refine { ls_R_factor_R_work ls_R_factor_R_free ls_R_factor_obs ls_R_factor_all }
      symmetry { space_group_name_H_M }
      cell { length_a length_b length_c angle_alpha angle_beta angle_gamma }
      rcsb_entry_container_identifiers { assembly_ids }
      exptl_crystal_grow { method pH temp pdbx_details }
      pdbx_vrpt_summary_geometry {
        clashscore percent_ramachandran_outliers percent_rotamer_outliers
      }
      pdbx_vrpt_summary_diffraction {
        percent_RSRZ_outliers EDS_R data_completeness
      }
      nonpolymer_entities { nonpolymer_comp { chem_comp { id name formula } } }
      citation {
        id title journal_abbrev year pdbx_database_id_DOI rcsb_authors
        journal_volume page_first page_last rcsb_is_primary
      }
    }
  }
}
"""


def _chunks(items: list, n: int) -> Iterable[list]:
    for i in range(0, len(items), n):
        yield items[i:i + n]


def fetch_entities(entity_ids: list[str]) -> list[dict]:
    """Batch-fetch full metadata for polymer entities, in GRAPHQL_BATCH-sized requests.

    Each batch is cached independently, so two families that overlap (a mutant series and
    its parent, say) share every batch they have in common.
    """
    out: list[dict] = []
    for batch in _chunks(entity_ids, config.GRAPHQL_BATCH):
        # Sort within the batch so the same set of ids in a different order is one cache
        # entry, not two.
        key_ids = sorted(batch)
        data = db.cached(
            ("entities", PARSE_VERSION, key_ids),
            lambda ids=key_ids: http.graphql(config.RCSB_GRAPHQL_URL, _ENTITY_QUERY, {"ids": ids}),
        )
        out.extend([e for e in (data.get("polymer_entities") or []) if e])
    return out


# --------------------------------------------------------------------------------------
# Parsing helpers
# --------------------------------------------------------------------------------------
PDB_ID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
ENTITY_ID_RE = re.compile(r"^([0-9][A-Za-z0-9]{3})_([0-9]+)$")


def _r_work(refine) -> Optional[float]:
    """The working R factor, from whichever mmCIF field this depositor happened to use.

    `_refine.ls_R_factor_R_work` is the field everyone quotes, and a good fraction of the
    archive leaves it empty and puts the number in `ls_R_factor_obs` or `ls_R_factor_all`
    instead: 3K34 uses `all`, 132L uses `obs`, 1AKI uses `R_work`. Reading only the obvious
    field showed a blank R-work beside a populated R-free on 82 of carbonic anhydrase II's
    1,473 X-ray entities, which reads as "this entry did not report one" rather than as
    "we looked in one place".
    """
    if not refine:
        return None
    row = refine[0] if isinstance(refine, list) else refine
    if not row:
        return None
    for field in ("ls_R_factor_R_work", "ls_R_factor_obs", "ls_R_factor_all"):
        value = row.get(field)
        if value is not None:
            return value
    return None


def _one(seq):
    """First row of a repeatable mmCIF category. These come back as lists even when they only
    ever have one row, and the validation summaries are lists too."""
    if isinstance(seq, list):
        return (seq[0] or {}) if seq else {}
    return seq or {}


def _crystal(grow) -> Optional[dict]:
    """The crystallisation condition as deposited: structured fields plus the free text."""
    g = _one(grow)
    if not g:
        return None
    ph, temp = g.get("pH"), g.get("temp")
    return {
        "method": g.get("method"),
        "ph": float(ph) if ph is not None else None,
        # Kelvin as deposited. Converted for display, not here: a stored value that has
        # silently been through a unit conversion is how two of them end up disagreeing.
        "temp_k": float(temp) if temp is not None else None,
        "details": g.get("pdbx_details"),
    }


def _validation(entry: dict) -> Optional[dict]:
    """The wwPDB validation report's headline numbers.

    Split across three sub-objects in the schema (`pdbx_vrpt_summary` itself holds almost
    nothing useful: the geometry and diffraction summaries are where the figures live).
    """
    geo = _one(entry.get("pdbx_vrpt_summary_geometry"))
    dif = _one(entry.get("pdbx_vrpt_summary_diffraction"))
    if not geo and not dif:
        return None
    out = {
        "clashscore": geo.get("clashscore"),
        "rama_outliers": geo.get("percent_ramachandran_outliers"),
        "rota_outliers": geo.get("percent_rotamer_outliers"),
        "rsrz_outliers": dif.get("percent_RSRZ_outliers"),
        "eds_r": dif.get("EDS_R"),
        "completeness": dif.get("data_completeness"),
    }
    return out if any(v is not None for v in out.values()) else None


def parse_entity(raw: dict) -> dict:
    """Flatten one GraphQL polymer_entity into the entity + entry rows we store.

    Returns `{"entity": {...}, "entry": {...}}`. Every field here is a list-or-null in the
    source, because mmCIF categories are repeatable even when they only ever have one row.
    """
    def first(seq, key=None, default=None):
        if not seq:
            return default
        head = seq[0] if isinstance(seq, list) else seq
        if head is None:
            return default
        return head.get(key, default) if key else head

    entity_id = raw["rcsb_id"]
    entry = raw.get("entry") or {}
    pdb_id = (entry.get("rcsb_id") or entity_id.split("_")[0]).upper()

    poly = raw.get("entity_poly") or {}
    ids = raw.get("rcsb_polymer_entity_container_identifiers") or {}

    # UniProt, and only UniProt: the same block also carries GenBank/EMBL/NORINE accessions
    # for some entities, and treating one of those as the reference would quietly break the
    # cross-reference on exactly the entries where it matters.
    #
    # ALL of them, not just the first. A chimera is cross-referenced to every protein it is
    # made of, and which one comes first is not meaningful: 2RH1 (the beta-2 adrenergic
    # receptor-T4 lysozyme fusion) lists T4 lysozyme first, so taking the head of the list
    # made the construct diff run backwards, treating the 164-residue fusion partner as the
    # canonical reference and the entire receptor as a pair of unexplained overhangs.
    uniprot_ids = [
        ref.get("database_accession")
        for ref in (ids.get("reference_sequence_identifiers") or [])
        if (ref or {}).get("database_name") == "UniProt" and ref.get("database_accession")
    ]
    uniprot = uniprot_ids[0] if uniprot_ids else None

    annotations = raw.get("rcsb_polymer_entity_annotation") or []
    pfam = [{"id": a["annotation_id"], "name": a.get("name")}
            for a in annotations if a and a.get("type") == "Pfam"]
    interpro = [{"id": a["annotation_id"], "name": a.get("name")}
                for a in annotations if a and a.get("type") == "InterPro"]

    # resolution_combined is a list because a structure can be solved by more than one
    # method; the best (lowest) figure is the one anyone quotes.
    info = entry.get("rcsb_entry_info") or {}
    res_list = [r for r in (info.get("resolution_combined") or []) if r is not None]
    resolution = min(res_list) if res_list else None

    cell = entry.get("cell") or {}
    has_cell = cell.get("length_a") is not None

    ligands = []
    seen_ligands = set()
    for ne in entry.get("nonpolymer_entities") or []:
        comp = ((ne or {}).get("nonpolymer_comp") or {}).get("chem_comp") or {}
        cid = comp.get("id")
        if cid and cid not in seen_ligands:
            seen_ligands.add(cid)
            ligands.append({"id": cid, "name": comp.get("name"), "formula": comp.get("formula")})

    citation = None
    for c in entry.get("citation") or []:
        if c and (c.get("rcsb_is_primary") == "Y" or c.get("id") == "primary"):
            citation = {
                "title": c.get("title"), "journal": c.get("journal_abbrev"),
                "year": c.get("year"), "doi": c.get("pdbx_database_id_DOI"),
                "authors": c.get("rcsb_authors") or [], "volume": c.get("journal_volume"),
                "pages": c.get("page_first"), "page_last": c.get("page_last"),
            }
            break

    return {
        "entity": {
            "entity_id": entity_id,
            "pdb_id": pdb_id,
            "description": (raw.get("rcsb_polymer_entity") or {}).get("pdbx_description"),
            "chains": ids.get("auth_asym_ids") or [],
            "seq_length": poly.get("rcsb_sample_sequence_length"),
            "sequence": poly.get("pdbx_seq_one_letter_code_can"),
            "organism": first(raw.get("rcsb_entity_source_organism"), "ncbi_scientific_name"),
            "taxonomy_id": first(raw.get("rcsb_entity_source_organism"), "ncbi_taxonomy_id"),
            "host_organism": first(raw.get("rcsb_entity_host_organism"), "ncbi_scientific_name"),
            "uniprot": uniprot,
            "uniprot_ids": uniprot_ids,
            "pfam": pfam,
            "interpro": interpro,
        },
        "entry": {
            "pdb_id": pdb_id,
            "title": (entry.get("struct") or {}).get("title"),
            "method": info.get("experimental_method") or first(entry.get("exptl"), "method"),
            "resolution": resolution,
            "r_work": _r_work(entry.get("refine")),
            "r_free": first(entry.get("refine"), "ls_R_factor_R_free"),
            "space_group": (entry.get("symmetry") or {}).get("space_group_name_H_M"),
            "cell": {
                "a": cell.get("length_a"), "b": cell.get("length_b"), "c": cell.get("length_c"),
                "alpha": cell.get("angle_alpha"), "beta": cell.get("angle_beta"),
                "gamma": cell.get("angle_gamma"),
            } if has_cell else None,
            "deposit_date": (entry.get("rcsb_accession_info") or {}).get("deposit_date", "")[:10] or None,
            "release_date": (entry.get("rcsb_accession_info") or {}).get("initial_release_date", "")[:10] or None,
            "chain_count": info.get("deposited_polymer_entity_instance_count"),
            "entity_count": info.get("polymer_entity_count"),
            "assembly_count": len((entry.get("rcsb_entry_container_identifiers") or {}).get("assembly_ids") or []),
            "ligands": ligands,
            "citation": citation,
            "crystal": _crystal(entry.get("exptl_crystal_grow")),
            "validation": _validation(entry),
            "has_sf": (entry.get("rcsb_accession_info") or {}).get(
                "has_released_experimental_data") == "Y",
        },
    }
