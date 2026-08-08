"""Interaction fingerprints: PLIP, run family-wide.

**Workstation only**, like `embed.py`. This shells out to PLIP and OpenBabel, and the droplet
has neither. The web app reads the JSON artefact and nothing else.

The mmCIF-to-PDB conversion is **not** reimplemented here. `pipeline/cif2plip.py` already does
it correctly, including the one non-obvious part: `pdb_tidy` leaves a one-number serial gap at
each chain/TER boundary, and OpenBabel maps CONECT serials positionally, so every CONECT
record past the gap is silently discarded. The result is a hybrid distance/CONECT bond
perception that garbles the ligand: wrong SMILES, fictitious bonds, and interactions that were
never there. That script renumbers contiguously and regenerates the ligand CONECT records
against the new serials.

What is different here is what PLIP is asked for. `cif2plip` requests a PyMOL session and PNG
images, which need PyMOL importable and are useless to a web page; this asks only for the XML
and parses that.
"""

from __future__ import annotations

import importlib.util
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional

from . import config

logger = logging.getLogger(__name__)

CONTACT_DIR = config.DATA_DIR / "contacts"
STRUCT_DIR = config.DATA_DIR / "structures"
CIF2PLIP = config.ROOT_DIR / "pipeline" / "cif2plip.py"

# PLIP's XML groups its interactions under plural tags; these are the singulars used
# throughout the app and the UI.
_SINGULAR = {
    "hydrophobic_interactions": "hydrophobic_interaction",
    "hydrogen_bonds": "hydrogen_bond",
    "water_bridges": "water_bridge",
    "salt_bridges": "salt_bridge",
    "pi_stacks": "pi_stack",
    "pi_cation_interactions": "pi_cation_interaction",
    "halogen_bonds": "halogen_bond",
    "metal_complexes": "metal_complex",
}

# The interaction types PLIP reports, in the order a reader cares about them.
INTERACTION_TYPES = [
    "hydrogen_bond", "hydrophobic_interaction", "salt_bridge", "pi_stack",
    "pi_cation_interaction", "water_bridge", "halogen_bond", "metal_complex",
]


# Same single-definition rule as embed/embed_io.
from .contacts_io import VERSION, artefact_path, load  # noqa: F401


# --------------------------------------------------------------------------------------
# Running PLIP
# --------------------------------------------------------------------------------------
def _load_cif2plip():
    """Import the vendored converter as a module.

    By spec, not by faking an entry in sys.modules: a module object without a real __spec__
    breaks anything that later inspects it, and the failure surfaces far from the cause.
    """
    if not CIF2PLIP.exists():
        raise FileNotFoundError(f"{CIF2PLIP} is missing (vendored from Marc's cif2plip)")
    spec = importlib.util.spec_from_file_location("cif2plip", CIF2PLIP)
    module = importlib.util.module_from_spec(spec)
    sys.modules["cif2plip"] = module
    spec.loader.exec_module(module)
    return module


def _venv_path_env() -> dict:
    """PLIP and pdb-tools install their console scripts into the venv's bin, which is not on
    PATH when Python is invoked by absolute path. cif2plip shells out to `pdb_tidy` by name,
    so the child needs a PATH that can find it."""
    env = dict(os.environ)
    bin_dir = str(Path(sys.executable).parent)
    env["PATH"] = bin_dir + os.pathsep + env.get("PATH", "")
    return env


def contacts_for(pdb_id: str, workdir: Path) -> Optional[list[dict]]:
    """Convert one entry and run PLIP over it, returning a flat list of contacts."""
    c2p = _load_cif2plip()
    cif = STRUCT_DIR / f"{pdb_id.upper()}.cif"
    if not cif.exists():
        from . import http
        if http.download(f"https://files.rcsb.org/download/{pdb_id.upper()}.cif", cif) is None:
            return None

    out = workdir / pdb_id.upper()
    out.mkdir(parents=True, exist_ok=True)
    pdb_path = out / f"{pdb_id.upper()}.pdb"

    try:
        # (structure, chain_map, resname_map, ligand_resnames)
        st, _chain_map, _resname_map, ligands = c2p.remap_structure(str(cif))
        c2p.write_pdb(st, str(pdb_path))
        subprocess.run(["pdb_tidy", str(pdb_path)], check=True, env=_venv_path_env(),
                       stdout=open(str(pdb_path) + ".tidy", "w"),
                       stderr=subprocess.DEVNULL)
        shutil.move(str(pdb_path) + ".tidy", pdb_path)
        # The serial-gap fix, which is the whole reason for reusing this script.
        c2p.renumber_contiguous(str(pdb_path))
        c2p.add_ligand_conect(str(pdb_path), ligands)
    except Exception:
        logger.warning("conversion failed for %s", pdb_id, exc_info=True)
        return None

    # XML only. No -y (PyMOL session) and no -p (images): both need PyMOL importable and
    # neither is any use to a web page.
    try:
        subprocess.run(["plip", "-f", str(pdb_path), "-x", "-o", str(out)],
                       check=True, env=_venv_path_env(), timeout=300,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
        logger.warning("PLIP failed for %s", pdb_id, exc_info=True)
        return None

    # PLIP names the report after the input file (`1DKK_report.xml`), not `report.xml`.
    reports = sorted(out.glob("*report.xml"))
    return _parse_report(reports[0], pdb_id) if reports else None


def _parse_report(xml_path: Path, pdb_id: str) -> list[dict]:
    """Flatten a PLIP XML report into one row per contact."""
    rows: list[dict] = []
    try:
        root = ET.parse(xml_path).getroot()
    except ET.ParseError:
        return rows

    for site in root.iter("bindingsite"):
        ident = site.find("identifiers")
        lig = (ident.findtext("hetid") or "").strip() if ident is not None else ""
        chain = (ident.findtext("chain") or "").strip() if ident is not None else ""
        if not lig:
            continue
        interactions = site.find("interactions")
        if interactions is None:
            continue
        for group in interactions:
            # An explicit map, not str.rstrip("s"): rstrip removes CHARACTERS, not a suffix,
            # so <metal_complexes> came out as "metal_complexe". An unrecognised tag keeps
            # its own name rather than being mangled into a near-miss.
            kind = _SINGULAR.get(group.tag, group.tag)
            for node in group:
                resnr = node.findtext("resnr")
                restype = (node.findtext("restype") or "").strip()
                if not resnr:
                    continue
                rows.append({
                    "pdb_id": pdb_id.upper(), "ligand": lig, "lig_chain": chain,
                    "type": kind, "resnr": int(resnr), "restype": restype,
                    "reschain": (node.findtext("reschain") or "").strip(),
                })
    return rows


# --------------------------------------------------------------------------------------
# Family-wide
# --------------------------------------------------------------------------------------
def build(fam: dict, max_entries: int = 60, progress=None) -> Optional[dict]:
    """Run PLIP across a family's representative ligand-bound entries.

    Only entries that actually carry a designed ligand: running PLIP over a structure whose
    only heteroatoms are glycerol and sulfate costs the same and tells you about the
    cryoprotectant. The ligand classification the Ligands panel already computed is what
    decides, so the two panels agree on what a ligand is.
    """
    CONTACT_DIR.mkdir(parents=True, exist_ok=True)
    STRUCT_DIR.mkdir(parents=True, exist_ok=True)

    holo = [m for m in fam["members"] if m.get("has_ligand")]
    # Best resolution first: a contact seen at 1.2 A is worth more than one at 3.5 A.
    holo.sort(key=lambda m: (m.get("resolution") is None, m.get("resolution") or 99))
    chosen, seen_constructs = [], set()
    for m in holo:
        # One entry per construct, so a construct solved 200 times does not dominate the
        # fingerprint with 200 copies of the same contact.
        key = m.get("seq_id")
        if key in seen_constructs:
            continue
        seen_constructs.add(key)
        chosen.append(m)
        if len(chosen) >= max_entries:
            break
    if not chosen:
        return None

    offsets = {m["pdb_id"]: (m.get("query_beg") or 1) - 1 for m in chosen}
    all_rows: list[dict] = []
    ok = failed = 0
    t0 = time.time()

    with tempfile.TemporaryDirectory(prefix="codswallop-plip-") as tmp:
        for i, m in enumerate(chosen):
            if progress:
                progress(i + 1, len(chosen), m["pdb_id"])
            rows = contacts_for(m["pdb_id"], Path(tmp))
            if rows is None:
                failed += 1
                continue
            ok += 1
            off = offsets.get(m["pdb_id"], 0)
            for r in rows:
                # PLIP reports author residue numbering; the rest of the app speaks seed
                # coordinates. Mapped through the same alignment offset the density census
                # uses, so a hot residue here is the same residue as in the conservation
                # track and the domain ribbon.
                r["seed_pos"] = r["resnr"] + off
                all_rows.append(r)

    if not all_rows:
        return None

    # ---- the family fingerprint --------------------------------------------------------
    per_residue: Counter = Counter()
    per_type: Counter = Counter()
    lig_res: dict = defaultdict(Counter)
    restype_of: dict[int, str] = {}
    for r in all_rows:
        per_residue[r["seed_pos"]] += 1
        per_type[r["type"]] += 1
        lig_res[r["ligand"]][r["seed_pos"]] += 1
        restype_of.setdefault(r["seed_pos"], r["restype"])

    # Metal coordination, kept separate from the rest of the fingerprint. Whether a metal is
    # a structural ion or the catalytic centre is a property of the PROTEIN, not of the
    # component, and this is the evidence that decides it per family: which residues
    # coordinate it, and in how many entries.
    metal_coord: dict = defaultdict(Counter)
    metal_entries: dict = defaultdict(set)
    for r in all_rows:
        if r["type"] != "metal_complex":
            continue
        metal_coord[r["ligand"]][r["seed_pos"]] += 1
        metal_entries[r["ligand"]].add(r["pdb_id"])

    # Per-residue breakdowns, not just a total. "Residue 199 makes 47 contacts" is a
    # ranking; "37 of them are hydrogen bonds and 6 are metal coordination" is what tells a
    # reader whether it is a catalytic residue or a wall of the pocket. Built here because
    # `all_rows` exists only during the build.
    res_types: dict = defaultdict(Counter)
    res_ligands: dict = defaultdict(Counter)
    res_entries: dict = defaultdict(set)
    for r in all_rows:
        res_types[r["seed_pos"]][r["type"]] += 1
        res_ligands[r["seed_pos"]][r["ligand"]] += 1
        res_entries[r["seed_pos"]].add(r["pdb_id"])

    hot = [{"pos": pos, "restype": restype_of.get(pos, ""), "contacts": n,
            "entries": len(res_entries[pos]),
            "types": res_types[pos].most_common(),
            "ligands": res_ligands[pos].most_common(12)}
           for pos, n in per_residue.most_common(60)]

    top_ligands = [lig for lig, _ in Counter(r["ligand"] for r in all_rows).most_common(24)]
    fingerprint = {
        lig: {str(pos): lig_res[lig][pos] for pos, _ in per_residue.most_common(40)}
        for lig in top_ligands
    }

    artefact = {
        "version": VERSION,
        "slug": fam["slug"],
        "built_at": int(time.time()),
        "entries_analysed": ok,
        "entries_failed": failed,
        "seconds": round(time.time() - t0, 1),
        "n_contacts": len(all_rows),
        "by_type": per_type.most_common(),
        "hot_residues": hot,
        "ligands": top_ligands,
        "fingerprint": fingerprint,
        "positions": [pos for pos, _ in per_residue.most_common(40)],
        "metal_coordination": {
            lig: {
                "residues": [{"pos": pos, "restype": restype_of.get(pos, ""), "n": n}
                             for pos, n in coords.most_common(12)],
                "entries": len(metal_entries[lig]),
            }
            for lig, coords in metal_coord.items()
        },
    }
    artefact_path(fam["slug"]).write_text(json.dumps(artefact, separators=(",", ":")))
    return artefact
