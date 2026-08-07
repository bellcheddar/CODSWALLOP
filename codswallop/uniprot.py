"""UniProt REST client: canonical sequence, name and organism for a seed.

Phase 1 needs only enough to seed a family and label it. The features, natural variants and
isoform machinery this API also exposes belong to Phase 2's construct diff engine.
"""

from __future__ import annotations

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


def _condense(rec: dict) -> dict:
    genes = [g.get("geneName", {}).get("value") for g in (rec.get("genes") or [])]
    return {
        "accession": rec.get("primaryAccession"),
        "id": rec.get("uniProtkbId"),
        "name": _protein_name(rec),
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

    return db.cached(("uniprot_entry", accession.upper()), fetch)


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

    return db.cached(("uniprot_search", query, size), fetch) or []


def by_gene(gene: str, size: int = 12) -> list[dict]:
    """Candidates for a bare gene name, across organisms."""
    return search(f"gene:{gene}", size)


_FEATURE_FIELDS = "ft_act_site,ft_binding,ft_signal,ft_transmem,ft_domain,ft_disulfid"


def features(accession: str) -> dict:
    """Positional annotations used to *describe* a construct difference.

    Only ever used to say what a mutation sits on ("at an annotated active site"), never to
    claim what it did. Whether a substitution actually inactivated the enzyme is a result
    from an experiment this tool has not seen.
    """
    def fetch():
        rec = http.get_json(config.UNIPROT_ENTRY_URL.format(accession=accession.upper()),
                            params={"fields": _FEATURE_FIELDS})
        if not rec:
            return {}
        out: dict[str, list] = {"active_site": [], "binding_site": [], "signal_peptide": [],
                                "transmembrane": [], "domain": [], "disulphide": []}
        key = {"Active site": "active_site", "Binding site": "binding_site",
               "Signal": "signal_peptide", "Transmembrane": "transmembrane",
               "Domain": "domain", "Disulfide bond": "disulphide"}
        for f in rec.get("features") or []:
            k = key.get(f.get("type"))
            if not k:
                continue
            loc = f.get("location") or {}
            beg = (loc.get("start") or {}).get("value")
            end = (loc.get("end") or {}).get("value")
            if beg is None or end is None:
                continue
            if k in ("active_site", "binding_site", "disulphide"):
                # Single positions: the classifier tests membership, not overlap.
                out[k].extend(range(int(beg), int(end) + 1))
            else:
                out[k].append({"start": int(beg), "end": int(end),
                               "description": f.get("description") or ""})
        return out

    return db.cached(("uniprot_features", accession.upper()), fetch) or {}
