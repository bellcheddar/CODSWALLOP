"""The family resolver: turn whatever the user typed into a seed.

Accepts a PDB ID, a PDB entity id, a UniProt accession, a gene name, a Pfam or InterPro
accession, a raw sequence (FASTA or bare), or a free-text protein name.

Everything resolves to the same thing: a **seed sequence**. That is the one definition of
family membership in CODSWALLOP, and it is what makes the identity slider mean something
exact. The input kind only decides how the seed is found, never what a family is.

Where the input genuinely admits several readings (a gene name that exists in fourteen
organisms, a PDB entry with three different polymers in it, a free-text name), the resolver
returns candidates for a disambiguation card. It never guesses. Where the *seed* rather
than the input is a choice (which member represents a Pfam family), it picks by a stated
rule and says so, which is a different thing: the answer is not in doubt, the starting
point is.
"""

from __future__ import annotations

import re
from typing import Optional

from . import interpro, rcsb, uniprot

# A UniProt accession, per their own published pattern.
UNIPROT_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$", re.I
)
# The 20 amino acids plus the ambiguity, selenocysteine, pyrrolysine and gap codes.
AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBZXUO*-")
NT_ALPHABET = set("ACGTUN-")
MIN_SEQUENCE = 20        # below this, mmseqs2 has nothing to work with


class ResolveError(Exception):
    """The input could not be turned into a family."""


# --------------------------------------------------------------------------------------
# Input classification
# --------------------------------------------------------------------------------------
def clean_sequence(text: str) -> str:
    """Strip FASTA headers, whitespace, digits and punctuation from a pasted sequence.

    Copy-pasted sequences arrive with residue numbering, line breaks and the occasional
    stray asterisk far more often than they arrive clean.
    """
    lines = [ln for ln in text.splitlines() if not ln.strip().startswith(">")]
    return re.sub(r"[^A-Za-z*\-]", "", "".join(lines)).upper()


def looks_like_sequence(text: str) -> bool:
    if text.strip().startswith(">"):
        return True
    seq = clean_sequence(text)
    if len(seq) < MIN_SEQUENCE:
        return False
    # Nucleotide strings pass this test deliberately. They are sequences, just the wrong
    # kind, and classifying them as free text instead would send them down the full-text
    # path where the only possible answer is "nothing matches ACGTACGT...". Letting them
    # through means `resolve` can say what is actually wrong.
    return set(seq) <= AA_ALPHABET


def classify(text: str) -> str:
    """Name the kind of thing the user typed. Order matters: the narrow patterns first."""
    q = text.strip()
    if not q:
        return "empty"
    if rcsb.ENTITY_ID_RE.match(q):
        return "entity_id"
    if rcsb.PDB_ID_RE.match(q):
        return "pdb_id"
    if interpro.PFAM_RE.match(q) or interpro.INTERPRO_RE.match(q):
        return "domain"
    if looks_like_sequence(q):
        return "sequence"
    if UNIPROT_RE.match(q):
        return "uniprot"
    # A bare short token with no spaces is far more likely a gene symbol than a phrase.
    if re.match(r"^[A-Za-z0-9][A-Za-z0-9\-_.]{1,14}$", q):
        return "gene"
    return "text"


# --------------------------------------------------------------------------------------
# Seeds
# --------------------------------------------------------------------------------------
def _seed_from_entity(entity_id: str, kind: str, note: Optional[str] = None) -> Optional[dict]:
    """Build a seed from one PDB polymer entity."""
    raw = rcsb.fetch_entities([entity_id])
    if not raw:
        return None
    parsed = rcsb.parse_entity(raw[0])
    ent = parsed["entity"]
    if not ent.get("sequence"):
        return None
    return {
        "kind": kind,
        "seed": entity_id,
        "name": ent.get("description"),
        "organism": ent.get("organism"),
        "uniprot": ent.get("uniprot"),
        "sequence": ent["sequence"],
        "length": ent.get("seq_length") or len(ent["sequence"]),
        "pfam": ent.get("pfam") or [],
        "interpro": ent.get("interpro") or [],
        "note": note or f"Seeded from PDB entity {entity_id}.",
    }


def _seed_from_uniprot(rec: dict, kind: str) -> dict:
    label = rec.get("name") or rec.get("id") or rec["accession"]
    return {
        "kind": kind,
        "seed": rec["accession"],
        "name": label,
        "organism": rec.get("organism"),
        "uniprot": rec["accession"],
        "sequence": rec["sequence"],
        "length": rec.get("length") or len(rec["sequence"]),
        "pfam": [],
        "interpro": [],
        "note": f"Seeded from the UniProt canonical sequence of {rec['accession']} "
                f"({rec.get('length')} residues).",
    }


def _representative_of_annotation(accession: str, meta: dict) -> Optional[dict]:
    """Pick the member that best represents a domain family, and say which and why.

    The rule: among entities carrying the annotation, the one from the highest-resolution
    experimental structure. That is a stated choice, not a guess about what the user meant,
    so it does not need a disambiguation card: it needs a sentence in the header saying
    which entity the identities are measured against.
    """
    candidates = rcsb.entities_by_annotation(accession, limit=25)
    if not candidates:
        return None
    parsed = [rcsb.parse_entity(r) for r in rcsb.fetch_entities(candidates)]
    usable = [p for p in parsed if p["entity"].get("sequence")]
    if not usable:
        return None
    # Sort by resolution, missing resolution last (cryo-EM and NMR entries without one).
    usable.sort(key=lambda p: (p["entry"].get("resolution") is None,
                               p["entry"].get("resolution") or 999))
    best = usable[0]
    ent, entry = best["entity"], best["entry"]
    res = entry.get("resolution")
    res_txt = f"{res} Å" if res else entry.get("method", "no resolution")
    name = (meta.get("name") or accession)
    return {
        "kind": "domain",
        "seed": ent["entity_id"],
        "name": name,
        "organism": ent.get("organism"),
        "uniprot": ent.get("uniprot"),
        "sequence": ent["sequence"],
        "length": ent.get("seq_length") or len(ent["sequence"]),
        "pfam": ent.get("pfam") or [],
        "interpro": ent.get("interpro") or [],
        "domain": {"accession": meta.get("accession"), "name": meta.get("name"),
                   "type": meta.get("type"), "source": meta.get("source")},
        "note": f"{accession} ({name}). Seeded from {ent['entity_id']}, the "
                f"highest-resolution deposited entity carrying this annotation ({res_txt}); "
                f"identities below are measured against it.",
    }


# --------------------------------------------------------------------------------------
# The resolver proper
# --------------------------------------------------------------------------------------
def resolve(query: str) -> dict:
    """Resolve user input.

    Returns one of:
      {"status": "resolved",  "seed": {...}}
      {"status": "ambiguous", "prompt": str, "candidates": [{...}]}
      {"status": "not_found", "message": str}

    Every candidate carries a `pick` string that resolves unambiguously on its own, so the
    disambiguation card is just a set of links back into this same function.
    """
    q = (query or "").strip()
    kind = classify(q)

    if kind == "empty":
        return {"status": "not_found", "message": "Nothing to look up."}

    # ---- a single PDB polymer entity: unambiguous -------------------------------------
    if kind == "entity_id":
        eid = q.upper()
        seed = _seed_from_entity(eid, "entity_id")
        if not seed:
            return {"status": "not_found", "message": f"No polymer entity {eid} in the PDB."}
        return {"status": "resolved", "seed": seed}

    # ---- a PDB entry: ambiguous when it holds more than one polymer --------------------
    if kind == "pdb_id":
        pdb = q.upper()
        entity_ids = rcsb.entities_of_entry(pdb)
        if not entity_ids:
            return {"status": "not_found",
                    "message": f"{pdb} is not an experimental PDB entry with a polymer in it."}
        if len(entity_ids) == 1:
            seed = _seed_from_entity(entity_ids[0], "pdb_id",
                                     note=f"Seeded from the single polymer of PDB entry {pdb}.")
            if not seed:
                return {"status": "not_found", "message": f"{pdb} has no usable protein sequence."}
            return {"status": "resolved", "seed": seed}
        # A complex. Which chain is the family? Only the user knows.
        parsed = [rcsb.parse_entity(r) for r in rcsb.fetch_entities(entity_ids)]
        candidates = [
            {
                "pick": p["entity"]["entity_id"],
                "label": p["entity"].get("description") or p["entity"]["entity_id"],
                "sublabel": p["entity"].get("organism"),
                "meta": f"{p['entity']['entity_id']} · chains "
                        f"{', '.join(p['entity'].get('chains') or []) or '?'} · "
                        f"{p['entity'].get('seq_length') or '?'} aa",
            }
            for p in parsed if p["entity"].get("sequence")
        ]
        if not candidates:
            return {"status": "not_found", "message": f"{pdb} has no protein polymer to seed from."}
        if len(candidates) == 1:
            return resolve(candidates[0]["pick"])
        return {
            "status": "ambiguous",
            "prompt": f"PDB {pdb} contains {len(candidates)} different polymers. "
                      f"Which one is the family?",
            "candidates": candidates,
        }

    # ---- Pfam / InterPro ---------------------------------------------------------------
    if kind == "domain":
        meta = interpro.entry(q)
        if not meta:
            return {"status": "not_found", "message": f"No InterPro or Pfam entry {q.upper()}."}
        seed = _representative_of_annotation(q.upper(), meta)
        if not seed:
            return {"status": "not_found",
                    "message": f"{q.upper()} ({meta.get('name')}) has no deposited structures."}
        return {"status": "resolved", "seed": seed}

    # ---- a pasted sequence: as unambiguous as input gets --------------------------------
    if kind == "sequence":
        seq = clean_sequence(q)
        if len(seq) < MIN_SEQUENCE:
            return {"status": "not_found",
                    "message": f"That sequence is {len(seq)} residues. "
                               f"{MIN_SEQUENCE} is the shortest that can be searched."}
        if set(seq) <= NT_ALPHABET:
            return {"status": "not_found",
                    "message": "That looks like a nucleotide sequence. "
                               "CODSWALLOP searches protein sequences."}
        return {"status": "resolved", "seed": {
            "kind": "sequence", "seed": None, "name": None, "organism": None, "uniprot": None,
            "sequence": seq, "length": len(seq), "pfam": [], "interpro": [],
            "note": f"Seeded from the {len(seq)}-residue sequence you pasted.",
        }}

    # ---- UniProt accession --------------------------------------------------------------
    if kind == "uniprot":
        rec = uniprot.entry(q)
        if rec and rec.get("sequence"):
            return {"status": "resolved", "seed": _seed_from_uniprot(rec, "uniprot")}
        # Not an accession after all: the pattern also matches some ordinary words.
        kind = "gene"

    # ---- gene symbol or free text -------------------------------------------------------
    hits = uniprot.by_gene(q) if kind == "gene" else uniprot.search(q)
    hits = [h for h in hits if h.get("sequence")]

    if len(hits) == 1:
        return {"status": "resolved", "seed": _seed_from_uniprot(hits[0], kind)}

    if len(hits) > 1:
        what = "gene" if kind == "gene" else "name"
        return {
            "status": "ambiguous",
            "prompt": f"That {what} matches {len(hits)} UniProt entries. Which protein?",
            "candidates": [
                {
                    "pick": h["accession"],
                    "label": h.get("name") or h["accession"],
                    "sublabel": h.get("organism"),
                    "meta": f"{h['accession']} · {h.get('length')} aa"
                            f"{' · reviewed' if h.get('reviewed') else ''}"
                            + (f" · {', '.join(h['genes'][:3])}" if h.get("genes") else ""),
                }
                for h in hits
            ],
        }

    # ---- nothing in UniProt: fall back to the PDB's own full-text index ------------------
    entity_ids = rcsb.entities_by_text(q, limit=12)
    if not entity_ids:
        return {"status": "not_found",
                "message": f"Nothing in UniProt or the PDB matches “{q}”."}
    parsed = [rcsb.parse_entity(r) for r in rcsb.fetch_entities(entity_ids)]
    candidates = [
        {
            "pick": p["entity"]["entity_id"],
            "label": p["entity"].get("description") or p["entity"]["entity_id"],
            "sublabel": p["entity"].get("organism"),
            "meta": f"{p['entity']['entity_id']} · {p['entity'].get('seq_length') or '?'} aa",
        }
        for p in parsed if p["entity"].get("sequence")
    ]
    if not candidates:
        return {"status": "not_found", "message": f"Nothing usable matches “{q}”."}
    if len(candidates) == 1:
        return resolve(candidates[0]["pick"])
    return {
        "status": "ambiguous",
        "prompt": f"No UniProt entry matches “{q}”. These PDB entities do. Which one?",
        "candidates": candidates,
    }
