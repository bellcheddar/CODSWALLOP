"""UniProt REST client: canonical sequence, name and organism for a seed.

Phase 1 needs only enough to seed a family and label it. The features, natural variants and
isoform machinery this API also exposes belong to Phase 2's construct diff engine.
"""

from __future__ import annotations

import re
from typing import Optional

from . import config, db, http

_FIELDS = "accession,id,protein_name,gene_names,organism_name,length,sequence,reviewed"


def _protein_name(rec: dict) -> Optional[str]:
    """Recommended name, or the first submitted name for the unreviewed entries that have
    no recommended one (most of TrEMBL)."""
    desc = rec.get("proteinDescription") or {}
    rec_name = (desc.get("recommendedName") or {}).get("fullName", {}).get("value")
    if rec_name:
        return rec_name
    for sub in desc.get("submissionNames") or []:
        val = (sub.get("fullName") or {}).get("value")
        if val:
            return val
    return None


def _alt_names(rec: dict) -> tuple[list, list]:
    """Every other name UniProt records for the protein, and the short forms.

    The `protein_name` field already asked for carries these; they were simply not parsed.
    They are what the nomenclature panel reconciles the deposited descriptions against:
    lysozyme is `Lysozyme C` officially and is also, officially,
    `1,4-beta-N-acetylmuramidase C` and `Allergen Gal d IV`, and a depositor who used either
    was not making a name up.
    """
    desc = rec.get("proteinDescription") or {}
    full, short = [], []
    groups = list(desc.get("alternativeNames") or [])
    rec_name = desc.get("recommendedName")
    if rec_name:
        # The recommended name's own short forms belong with the short list, not as
        # alternative full names.
        short += [s.get("value") for s in (rec_name.get("shortNames") or [])]
    for alt in groups:
        value = (alt.get("fullName") or {}).get("value")
        if value:
            full.append(value)
        short += [s.get("value") for s in (alt.get("shortNames") or [])]
    for sub in desc.get("submissionNames") or []:
        value = (sub.get("fullName") or {}).get("value")
        if value:
            full.append(value)
    return [f for f in full if f], [s for s in short if s]


def _condense(rec: dict) -> dict:
    genes = [g.get("geneName", {}).get("value") for g in (rec.get("genes") or [])]
    # Gene synonyms count as names a depositor could reasonably have used.
    for g in rec.get("genes") or []:
        genes += [s.get("value") for s in (g.get("synonyms") or [])]
    alt, short = _alt_names(rec)
    return {
        "accession": rec.get("primaryAccession"),
        "id": rec.get("uniProtkbId"),
        "name": _protein_name(rec),
        "alt_names": alt,
        "short_names": short,
        "genes": [g for g in genes if g],
        "organism": (rec.get("organism") or {}).get("scientificName"),
        "taxonomy_id": (rec.get("organism") or {}).get("taxonId"),
        "sequence": (rec.get("sequence") or {}).get("value"),
        "length": (rec.get("sequence") or {}).get("length"),
        "reviewed": rec.get("entryType", "").startswith("UniProtKB reviewed"),
    }


def entry(accession: str) -> Optional[dict]:
    """One accession. None if it does not exist."""
    url = config.UNIPROT_ENTRY_URL.format(accession=accession.upper())

    def fetch():
        rec = http.get_json(url, params={"fields": _FIELDS})
        return _condense(rec) if rec else None

    # PARSE_VERSION in the key, for the reason spelled out where it is defined below:
    # records are cached in *parsed* form, so adding a field to `_condense` leaves it absent
    # from every already-cached accession and the feature reads as unavailable everywhere
    # rather than as new. This key was missing it while `uniprot_both` had it.
    return db.cached(("uniprot_entry", PARSE_VERSION, accession.upper()), fetch)


def search(query: str, size: int = 12) -> list[dict]:
    """Free-form UniProt search, reviewed entries first.

    Reviewed-first ordering is deliberate: a gene name like `LYZ` matches thousands of
    TrEMBL fragments, and the handful of Swiss-Prot entries are what a structural biologist
    means by it. The disambiguation card shows those.
    """
    def fetch():
        out = []
        seen = set()
        for q in (f"({query}) AND (reviewed:true)", f"({query}) AND (reviewed:false)"):
            body = http.get_json(
                config.UNIPROT_SEARCH_URL,
                params={"query": q, "fields": _FIELDS, "format": "json", "size": size},
            ) or {}
            for rec in body.get("results") or []:
                c = _condense(rec)
                if c["accession"] and c["accession"] not in seen:
                    seen.add(c["accession"])
                    out.append(c)
            if len(out) >= size:
                break
        return out[:size]

    # Ranked on the way OUT of the cache, not on the way in, so changing the ordering rule
    # never needs a cache version bump. Every caller gets it: `by_gene` ranks again after
    # merging its two queries, which is harmless because the sort is idempotent, and the
    # free-text path gets it without having to remember to ask. That path is why this is
    # here rather than in `by_gene` alone: searching "vascular endothelial growth factor
    # receptor" put Mus musculus first and a snake third, with the human protein below both.
    return rank_candidates(db.cached(("uniprot_search", PARSE_VERSION, query, size), fetch)
                           or [], query)


def by_gene(gene: str, size: int = 12) -> list[dict]:
    """Candidates for a bare gene name, across organisms.

    A gene search alone is not enough, because the thing a person types is often not a gene
    symbol in the database's sense. `VEGFR` returned twelve unreviewed *Drosophila* isoforms
    of "PDGF- and VEGF-receptor related" and no human protein at all: no reviewed entry
    carries VEGFR as a gene name, since the human genes are KDR, FLT1 and FLT4, while
    Drosophila's Pvr lists it as a synonym. The gene query matched nothing reviewed and fell
    straight through to the unreviewed tail.

    So the protein NAME is searched too when the gene finds no reviewed entry. Searching the
    name for VEGFR finds "Vascular endothelial growth factor receptor 1/2/3", which is what
    somebody typing it meant.
    """
    hits = search(f"gene:{gene}", size)
    if not any(h.get("reviewed") for h in hits):
        seen = {h["accession"] for h in hits}
        for extra in search(f'protein_name:"{gene}"', size):
            if extra["accession"] not in seen:
                seen.add(extra["accession"])
                hits.append(extra)
    return rank_candidates(hits, gene)[:size]


def search_terms(query: str) -> str:
    """The words a person actually typed, with any field prefix stripped.

    `search` is called with `gene:ABL1` and `protein_name:"VEGFR"` as well as bare text, and
    the ranking needs the terms rather than the query syntax around them.
    """
    text = (query or "").strip()
    if ":" in text:
        head, _, tail = text.partition(":")
        if head.strip().lower() in {"gene", "protein_name", "gene_exact", "family", "organism"}:
            text = tail
    return text.strip().strip('"').strip()


def name_match(name: str, terms: str) -> int:
    """How well a protein's name answers what was typed. Lower is better.

    0 exact, 1 the name begins with it, 2 it appears as a whole word, 3 not in the name at
    all. The last case is not a non-match overall: UniProt also matches on gene names, on
    keywords and on the *domains* a protein contains, which is exactly how the accident
    below happens.
    """
    n, t = (name or "").strip().lower(), (terms or "").strip().lower()
    if not t or not n:
        return 3
    if n == t:
        return 0
    if n.startswith(t):
        return 1
    return 2 if re.search(r"\b" + re.escape(t) + r"\b", n) else 3


def rank_candidates(hits: list[dict], query: str = "") -> list[dict]:
    """Reviewed first, then how well the name matches, then human, then alphabetical.

    Reviewed-before-unreviewed is the original rule and still the important one. Human-first
    came next: a disambiguation card for a receptor that lists zebrafish above *Homo sapiens*
    is technically a list and practically the wrong one, since a structural biologist typing
    a receptor's name almost always wants the human protein and can see the rest below it.

    Name-match sits ABOVE human because human-first alone is only right when a human protein
    of that name exists. **Beta-lactamase** is the case that proves it: humans do not have
    one, so the eight human hits are proteins that merely CONTAIN a beta-lactamase-like fold
    (Apollo exonuclease, MBLAC2, COA7, dipeptidase 1, LACTB2 and so on), and human-first put
    every one of them above the actual enzymes, which are bacterial and are among the most
    deposited proteins in the archive. `Beta-lactamase` itself ranked ninth.

    The remaining term is why the top hit was "5' exonuclease Apollo": with reviewed and
    human tied across all eight, the tie-break was the NAME, alphabetically, so first place
    went to whichever protein sorted first, decided by an apostrophe. Alphabetical is a fine
    last resort for a stable order and was never meant to choose the answer.

    Human-first is preserved intact wherever it was doing its job: it still decides between
    candidates whose names match the query equally well, which is the VEGFR case, since no
    protein there is named "VEGFR" at all and all of them tie at 3.
    """
    terms = search_terms(query)

    def key(h):
        organism = (h.get("organism") or "").lower()
        return (
            0 if h.get("reviewed") else 1,
            name_match(h.get("name") or "", terms),
            0 if organism.startswith("homo sapiens") else 1,
            (h.get("name") or "").lower(),
        )
    return sorted(hits, key=key)


# Phase 4 adds the post-translational ones for the motifs panel. `ft_mod_res` is the
# phosphorylation/acetylation/methylation column, `ft_carbohyd` the glycosylation and
# `ft_lipid` the myristoylation and palmitoylation: these are curated observations on the
# real protein, which is what makes them worth more than a pattern that matches by chance.
_FEATURE_FIELDS = ("ft_act_site,ft_binding,ft_signal,ft_transmem,ft_domain,ft_disulfid,"
                   "ft_mod_res,ft_carbohyd,ft_lipid,ft_site,ft_motif")

# One request per accession, not two. The entry and its features come from the same document
# and the same endpoint, and asking twice doubled the UniProt half of a cold build: 72 calls
# and 24 seconds on a family with 36 references.
_ALL_FIELDS = _FIELDS + "," + _FEATURE_FIELDS


# Part of the cache key, for the same reason rcsb.PARSE_VERSION is. Records are cached in
# *parsed* form, so adding a requested field leaves it absent from every cached row and the
# consumer reads an empty list rather than failing: adding the PTM columns without bumping
# this returned a family with no post-translational modifications at all, which for EGFR is
# a confident and completely wrong answer.
# 4: `_condense` gained alt_names, short_names and gene synonyms for the nomenclature
#    panel, and `uniprot_entry` gained this key.
PARSE_VERSION = 4


def entry_with_features(accession: str) -> tuple:
    """Both halves of one accession in a single round trip."""
    acc = accession.upper()

    def fetch():
        rec = http.get_json(config.UNIPROT_ENTRY_URL.format(accession=acc),
                            params={"fields": _ALL_FIELDS})
        if not rec:
            return [None, {}]
        return [_condense(rec), _features_from(rec)]

    got = db.cached(("uniprot_both", PARSE_VERSION, acc), fetch) or [None, {}]
    return got[0], got[1]


def features(accession: str) -> dict:
    """Positional annotations used to *describe* a construct difference.

    Only ever used to say what a mutation sits on ("at an annotated active site"), never to
    claim what it did. Whether a substitution actually inactivated the enzyme is a result
    from an experiment this tool has not seen.
    """
    return entry_with_features(accession)[1]


def _features_from(rec: dict) -> dict:
    """Positional annotations, from an already-fetched UniProt record."""
    out: dict[str, list] = {"active_site": [], "binding_site": [], "signal_peptide": [],
                            "transmembrane": [], "domain": [], "disulphide": [],
                            "modified": [], "glycosylation": [], "lipidation": [],
                            "site": [], "motif": []}
    key = {"Active site": "active_site", "Binding site": "binding_site",
           "Signal": "signal_peptide", "Transmembrane": "transmembrane",
           "Domain": "domain", "Disulfide bond": "disulphide",
           "Modified residue": "modified", "Glycosylation": "glycosylation",
           "Lipidation": "lipidation", "Site": "site", "Short sequence motif": "motif"}
    for f in rec.get("features") or []:
        k = key.get(f.get("type"))
        if not k:
            continue
        loc = f.get("location") or {}
        beg = (loc.get("start") or {}).get("value")
        end = (loc.get("end") or {}).get("value")
        if beg is None or end is None:
            continue
        if k == "disulphide":
            # A disulphide's two numbers are the two cysteines, not the ends of a span. The
            # range between them is ordinary sequence: expanding it turned lysozyme's
            # 24-145 bond into 122 consecutive "disulphide bonds", one per residue.
            out[k].append({"start": int(beg), "end": int(end),
                           "description": f.get("description") or ""})
        elif k in ("active_site", "binding_site"):
            # Single positions: the classifier tests membership, not overlap.
            out[k].extend(range(int(beg), int(end) + 1))
        else:
            out[k].append({"start": int(beg), "end": int(end),
                           "description": f.get("description") or ""})
    return out
