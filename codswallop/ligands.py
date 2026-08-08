"""Chemical components, and an honest verdict on which of them anyone meant to be there.

Phase 1 shipped a 56-entry exclusion list so the amber "ligand-bound" halo was roughly right
on day one. This replaces it with a per-component classification: **ligand**, **cofactor**,
**ion**, **buffer**, **cryoprotectant**, **precipitant** or **solvent**.

The distinction matters more than it sounds. "How many structures of this family are
ligand-bound" is a headline number on the overview, and if glycerol and sulfate count as
ligands then it is the number of structures that were frozen, which is all of them.
"""

from __future__ import annotations

from typing import Iterable, Optional

from . import config, db, http

PARSE_VERSION = 1

# --------------------------------------------------------------------------------------
# The classification
# --------------------------------------------------------------------------------------
# Explicit CCD codes first: these are unambiguous and a name match could not improve on them.
WATER = {"HOH", "DOD", "D8U"}

IONS = {
    "NA", "K", "LI", "CS", "RB", "MG", "CA", "SR", "BA", "MN", "FE", "FE2", "CO", "NI",
    "CU", "CU1", "ZN", "CD", "HG", "AG", "AU", "PT", "PB", "AL", "GA", "IN", "TL",
    "CL", "BR", "IOD", "F", "3CO", "4MO", "6MO", "W", "V", "CR", "YB", "SM", "EU", "GD",
    "LU", "TB", "HO", "ER", "PR", "ND", "CE", "LA", "Y", "SC", "ZR", "NB", "TA", "RE",
    "OS", "IR", "RU", "RH", "PD", "SE", "TE", "AS", "SB", "BI", "SN",
    "NH4", "NO3", "CO3", "SO4", "PO4", "SO3", "PO3", "AZI", "CN", "SCN", "OH", "O",
    "BEF", "ALF", "MOO", "WO4", "VO4", "CAC", "PER", "CLO",
}

CRYOPROTECTANTS = {
    "GOL",   # glycerol
    "EDO",   # ethylene glycol
    "MPD", "MRD",
    "PEG", "PG4", "PGE", "1PE", "2PE", "P6G", "P33", "PE3", "PE4", "PE5", "PE8", "XPE",
    "7PE", "12P", "15P", "M2M", "MXE", "DEG", "TEG", "PGO", "PDO", "PGR", "SPD",
    "SUC", "TRE", "GLC", "XYL", "MAN", "GAL", "FRU",
    "DMS",   # DMSO
}

BUFFERS = {
    "TRS",   # Tris
    "EPE",   # HEPES
    "MES", "MPO", "PIN", "TAM", "BTB", "BIS", "CIT", "FLC", "CAC", "IMD", "HEZ",
    "ACT", "ACY", "FMT", "MLA", "MLI", "TAR", "TLA", "MAE", "SIN", "GLY", "BES",
    "CXS", "CHE", "ADE", "POP", "PPI", "BCT", "DTT", "DTU", "BME", "TCE", "MRC",
    "BCN", "BIC",   # bicine: read as a ligand until it was spotted in carbonic anhydrase II
    "TAPS", "TAPSO", "MOPS", "TES", "ADA", "AMP0", "EPPS", "HEPPS", "CHAPS", "NHE",
}

# Common cofactors and covalently-relevant chemistry: these ARE meant to be there, and
# calling them buffer would be as wrong as calling glycerol a ligand.
COFACTORS = {
    "NAD", "NAI", "NAP", "NDP", "FAD", "FMN", "FMNH", "ADP", "ATP", "AMP", "GDP", "GTP",
    "GMP", "CDP", "CTP", "UDP", "UTP", "COA", "ACO", "SAM", "SAH", "TPP", "PLP", "B12",
    "BTN", "HEM", "HEC", "HEA", "HEB", "SRM", "CLA", "BCL", "PQQ", "MGD", "MOS", "F43",
    "H4B", "THG", "FOL", "5GP", "ANP", "AGS", "ACP", "GNP", "GSP", "ADN", "NAG", "BMA",
    "MAN", "FUC", "SIA", "GLA", "NDG", "A2G", "XYP",
}

# Detergents and lipids: not a ligand in the "somebody designed this to bind" sense, but not
# noise either, because for a membrane protein they are the reason the crystal exists.
LIPIDS_DETERGENTS = {
    "OLC", "OLA", "OLB", "PLM", "MYR", "STE", "PEE", "PEF", "PC1", "PCW", "LMT", "LMU",
    "DDQ", "BOG", "BNG", "OCT", "C8E", "NG6", "HTG", "F09", "UND", "TWT", "CHS", "CLR",
    "Y01", "D10", "D12", "LDA", "LPP", "PGV", "PGT", "CDL", "9PE", "PX4", "TRD",
    "2CV", "HEG", "MEG", "P4C", "SDS", "TTX", "NOM", "L2C", "L3P", "LMN", "GLE",
}

CLASS_ORDER = ["ligand", "cofactor", "lipid/detergent", "ion", "buffer",
               "cryoprotectant", "solvent", "water"]

# What "ligand-bound" should mean on the overview: something somebody put there on purpose.
COUNTS_AS_BOUND = {"ligand", "cofactor"}

# A caveat this module cannot resolve on its own, and states rather than papers over: whether
# a metal is a structural ion or the catalytic cofactor is a property of the PROTEIN, not of
# the component. Zinc is filed here as an ion, which is right for the many proteins that
# merely have one bound and wrong for carbonic anhydrase, where it is the catalytic centre.
# Deciding per family needs the UniProt binding-site annotations (already fetched for the
# construct diff) cross-referenced against the component's contacts, which is Phase 3's
# interaction work. Until then the class is shown on every component so a reader can see the
# call that was made rather than only its consequence.
METAL_CAVEAT = ("Metals are classified as ions. Whether a metal is structural or catalytic "
                "depends on the protein, not the component: the zinc in carbonic anhydrase "
                "is a catalytic cofactor filed here as an ion.")


def classify(comp_id: str, name: Optional[str] = None,
             formula: Optional[str] = None, atom_count: Optional[int] = None) -> str:
    """Verdict for one chemical component."""
    cid = (comp_id or "").upper()
    if cid in WATER:
        return "water"
    if cid in COFACTORS:
        return "cofactor"
    if cid in LIPIDS_DETERGENTS:
        return "lipid/detergent"
    if cid in CRYOPROTECTANTS:
        return "cryoprotectant"
    if cid in BUFFERS:
        return "buffer"
    if cid in IONS:
        return "ion"

    text = (name or "").lower()
    # A name-based fallback for the long tail. Deliberately narrow: these phrases appear in
    # the component's own CCD name, not inferred from context.
    if any(k in text for k in ("polyethylene glycol", "peg ", "glycol", "glycerol")):
        return "cryoprotectant"
    if any(k in text for k in ("buffer", "tris(", "acetate ion", "citrate", "cacodylate")):
        return "buffer"
    if text.endswith(" ion") or " ion" in text and (atom_count or 99) <= 4:
        return "ion"
    # A very small non-ion component is almost never a designed ligand.
    if atom_count is not None and atom_count <= 3:
        return "solvent"
    return "ligand"


# --------------------------------------------------------------------------------------
# Fetch
# --------------------------------------------------------------------------------------
_QUERY = """
query($ids: [String!]!) {
  chem_comps(comp_ids: $ids) {
    rcsb_id
    chem_comp { id name formula formula_weight type }
    rcsb_chem_comp_descriptor { SMILES InChIKey }
    rcsb_chem_comp_info { atom_count }
  }
}
"""


def _chunks(items: list, n: int) -> Iterable[list]:
    for i in range(0, len(items), n):
        yield items[i:i + n]


def details(comp_ids: list[str]) -> dict:
    """Batch-fetch component detail and classify each one."""
    out: dict[str, dict] = {}
    for batch in _chunks(sorted(set(c.upper() for c in comp_ids if c)), config.GRAPHQL_BATCH):
        data = db.cached(
            ("chemcomps", PARSE_VERSION, batch),
            lambda ids=batch: http.graphql(config.RCSB_GRAPHQL_URL, _QUERY, {"ids": ids}),
        )
        for c in (data.get("chem_comps") or []):
            if not c:
                continue
            cc = c.get("chem_comp") or {}
            desc = c.get("rcsb_chem_comp_descriptor") or {}
            info = c.get("rcsb_chem_comp_info") or {}
            cid = (cc.get("id") or c.get("rcsb_id") or "").upper()
            out[cid] = {
                "id": cid,
                "name": cc.get("name"),
                "formula": cc.get("formula"),
                "weight": cc.get("formula_weight"),
                "type": cc.get("type"),
                "smiles": desc.get("SMILES"),
                "inchikey": desc.get("InChIKey"),
                "atom_count": info.get("atom_count"),
                "klass": classify(cid, cc.get("name"), cc.get("formula"),
                                  info.get("atom_count")),
            }
    return out


def summarise(entries: list[dict]) -> dict:
    """Every component in the family, classified, most frequent first."""
    from collections import Counter
    counts: Counter = Counter()
    best: dict[str, float] = {}
    for e in entries:
        for lig in e.get("ligands") or []:
            cid = (lig.get("id") or "").upper()
            if not cid:
                continue
            counts[cid] += 1
            r = e.get("resolution")
            if r is not None and (cid not in best or r < best[cid]):
                best[cid] = r
    if not counts:
        return {"components": [], "by_class": {}, "n": 0}

    detail = details(list(counts))
    comps = []
    for cid, n in counts.most_common():
        d = detail.get(cid) or {"id": cid, "klass": classify(cid)}
        comps.append({**d, "count": n, "best_resolution": best.get(cid)})

    by_class: dict[str, int] = {}
    for c in comps:
        by_class[c["klass"]] = by_class.get(c["klass"], 0) + 1
    return {"components": comps, "by_class": by_class, "n": len(comps)}
