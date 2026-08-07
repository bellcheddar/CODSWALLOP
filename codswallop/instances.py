"""Per-chain data: which residues were actually seen, and which domains sit where.

Both come from the same batched GraphQL call against `polymer_entity_instances`, because
both are instance-level facts: disorder is a property of a chain in a crystal, not of a
sequence, and a domain assignment is made against a deposited chain.

This is what separates Phase 2's coverage census from Phase 1's. Phase 1 answers "was this
residue in anybody's construct", which for a well-studied protein is yes everywhere. This
answers the harder and more useful one: of the residues that *were* in a construct, which has
nobody ever seen density for.
"""

from __future__ import annotations

from typing import Iterable

from . import config, db, http

# Bumped when the parsed shape changes, like rcsb.PARSE_VERSION and for the same reason.
PARSE_VERSION = 1

# One chain per entity, not every chain. A family of 1,500 entries has 3,000 to 6,000 chains,
# and the aggregate census is dominated by variation between entries rather than between the
# copies of one entity in one asymmetric unit. Stated in the UI rather than left implicit.
_QUERY = """
query($ids: [String!]!) {
  polymer_entity_instances(instance_ids: $ids) {
    rcsb_id
    rcsb_polymer_entity_instance_container_identifiers { entity_id auth_asym_id }
    rcsb_polymer_instance_feature {
      type name feature_id
      feature_positions { beg_seq_id end_seq_id }
    }
  }
}
"""

# The feature types worth carrying. The API returns two dozen per chain (clashes, ASA,
# torsion outliers, secondary structure), and pulling the lot through the cache for a
# 1,500-chain family is megabytes of data nothing reads.
_WANTED = {"UNOBSERVED_RESIDUE_XYZ", "CATH", "SCOP2B_SUPERFAMILY", "ECOD"}


def _chunks(items: list, n: int) -> Iterable[list]:
    for i in range(0, len(items), n):
        yield items[i:i + n]


def _parse(payload: dict) -> dict:
    """Reduce one GraphQL batch to `{instance_id: {unobserved, domains}}`."""
    out: dict[str, dict] = {}
    for inst in payload.get("polymer_entity_instances") or []:
        if not inst:
            continue
        ident = inst.get("rcsb_polymer_entity_instance_container_identifiers") or {}
        rec: dict = {"entity_id": ident.get("entity_id"),
                     "chain": ident.get("auth_asym_id"),
                     "unobserved": [], "domains": []}
        for f in inst.get("rcsb_polymer_instance_feature") or []:
            ftype = (f or {}).get("type")
            if ftype not in _WANTED:
                continue
            spans = [(p["beg_seq_id"], p["end_seq_id"])
                     for p in (f.get("feature_positions") or [])
                     if p and p.get("beg_seq_id") and p.get("end_seq_id")]
            if not spans:
                continue
            if ftype == "UNOBSERVED_RESIDUE_XYZ":
                rec["unobserved"].extend(spans)
            else:
                rec["domains"].append({
                    "source": {"SCOP2B_SUPERFAMILY": "SCOP2B"}.get(ftype, ftype),
                    "id": f.get("feature_id"), "name": f.get("name"),
                    "spans": spans,
                })
        out[inst["rcsb_id"]] = rec
    return out


def fetch(instance_ids: list[str]) -> dict:
    """Batch-fetch per-chain features. Each batch is cached independently."""
    out: dict[str, dict] = {}
    for batch in _chunks(instance_ids, config.GRAPHQL_BATCH):
        key_ids = sorted(batch)
        part = db.cached(
            ("instances", PARSE_VERSION, key_ids),
            lambda ids=key_ids: _parse(
                http.graphql(config.RCSB_GRAPHQL_URL, _QUERY, {"ids": ids})),
        )
        out.update(part)
    return out


def representative_ids(members: list[dict]) -> list[str]:
    """One instance id per entity: `<PDBID>.<first auth chain>`.

    The first chain rather than the best-ordered one. Picking the best-ordered chain would
    bias the census optimistic (it is the definition of the chain with least disorder) and
    picking the worst would bias it the other way; the first is arbitrary with respect to
    disorder, which is the property being measured.
    """
    ids = []
    for m in members:
        chains = m.get("chains") or []
        if chains:
            ids.append(f"{m['pdb_id']}.{chains[0]}")
    return ids
