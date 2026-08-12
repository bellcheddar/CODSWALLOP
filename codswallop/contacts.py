"""Interaction fingerprints: PLIP, run family-wide.

Runs wherever the artefacts are built, which is now the droplet. It shells out to PLIP and
OpenBabel; the web app imports none of this and reads the JSON artefact through
`contacts_io`, which imports nothing beyond the standard library.

This is the one artefact that still needs the WHOLE deposited structure. Interaction
detection is every atom, every ligand and every chain, so the single-chain fetch that made
the embedding runnable on a small box does not apply here, and an entry above
MAX_STRUCTURE_MB is skipped and counted rather than parsed.

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


def _observed_chains(st) -> dict:
    """Chain -> [(author residue number, one-letter code)], in order.

    Read off the SAME gemmi structure PLIP is about to be given, after `remap_structure` has
    renamed the chains, so the chain letters here are the ones PLIP will report. Taking them
    from the original CIF instead would look right and mis-key every entry whose chains got
    renamed.

    `one_letter_code` is lower case for a modified residue (MSE -> 'm'), which is exactly the
    information the aligner wants folded away, so it is upper-cased.
    """
    import gemmi
    out: dict = {}
    for ch in st[0]:
        seq = []
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if info is None or not info.is_amino_acid():
                continue
            code = (info.one_letter_code or "X").upper()
            seq.append((res.seqid.num, code if code.isalpha() else "X"))
        if seq:
            out[ch.name] = seq
    return out


# A backstop only. Which chains belong to the family member is READ from the entity record
# rather than inferred, because inferring it does not work: an unrelated 78-residue partner
# aligned to carbonic anhydrase at 37.8 % identity over 74 columns, which is inside the
# twilight zone and above any threshold that still admits a legitimate 35 %-identity
# orthologue. Identity cannot separate those two cases, so it is not asked to.
MIN_CHAIN_IDENTITY = 0.2

# PLIP is the one artefact that genuinely needs the WHOLE structure: interaction detection is
# every atom, every ligand, every chain, so the per-chain fetch the embedding uses does not
# apply and the memory problem comes straight back. 8GLV is 453 MB and peaks at 6.3 GB to
# parse, on a droplet with about 2.2 GB free and no swap.
#
# So an entry too large to convert safely is skipped, and the artefact records how many were
# and why. That loses interaction data for the very largest complexes, which is the honest
# trade: the alternative is not slower contacts, it is an OOM that takes the other eight apps
# on the box down with it. The cap is generous next to what a family actually holds, where
# the median entry is about 1 MB.
MAX_STRUCTURE_MB = 30


def _chain_mappings(observed: dict, seed_sequence: str,
                    only: Optional[set] = None) -> dict:
    """(chain, author residue number) -> seed position, for the member's own chains.

    Aligned, not offset, for the reason `topology.map_to_seed` gives: a deposited entry may
    be numbered on the mature protein, on the construct from 1, or on the canonical, and an
    assumed offset is right often enough to look like it works.

    `only` is the set of chains the family member actually occupies, taken from the entity
    record and translated through cif2plip's renaming. Everything else in the file is a
    partner, a tag or another copy of something, and its contacts are not this protein's.
    """
    if not seed_sequence:
        return {}
    from Bio.Align import PairwiseAligner, substitution_matrices
    from .constructs import _alignable

    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score, aligner.extend_gap_score = -11, -1
    aligner.mode = "global"

    seed = _alignable(seed_sequence)
    out: dict = {}
    for chain, residues in (observed or {}).items():
        if only is not None and chain not in only:
            continue
        if not residues:
            continue
        obs = _alignable("".join(code for _, code in residues))
        try:
            aln = aligner.align(seed, obs)[0]
        except Exception:                       # noqa: BLE001 - an unalignable chain
            continue
        pairs, matches = 0, 0
        chain_map = {}
        for (s0, s1), (o0, o1) in zip(aln.aligned[0], aln.aligned[1]):
            for k in range(s1 - s0):
                pairs += 1
                if seed[s0 + k] == obs[o0 + k]:
                    matches += 1
                chain_map[(chain, residues[o0 + k][0])] = int(s0 + k) + 1
        if not pairs:
            continue
        if matches / pairs < MIN_CHAIN_IDENTITY:
            continue
        out.update(chain_map)
    return out


def contacts_for(pdb_id: str, workdir: Path) -> Optional[tuple]:
    """Convert one entry and run PLIP over it.

    Returns `(rows, observed_chains, chain_map)`: the flat contact list, what each chain
    actually contains, and the original-to-PLIP chain renaming. The last two are what put the
    contacts into seed coordinates without an assumed offset.
    """
    c2p = _load_cif2plip()
    cif = STRUCT_DIR / f"{pdb_id.upper()}.cif"
    if not cif.exists():
        from . import http
        # Aborted at the cap rather than downloaded and then declined: a HEAD on
        # files.rcsb.org times out, so the size cannot be known in advance, and streaming
        # until the limit is passed costs at most the limit.
        if http.download(f"https://files.rcsb.org/download/{pdb_id.upper()}.cif", cif,
                         max_bytes=int(MAX_STRUCTURE_MB * 1e6)) is None:
            logger.info("skipping %s for contacts: over %d MB", pdb_id, MAX_STRUCTURE_MB)
            return "too_big"
    if cif.stat().st_size / 1e6 > MAX_STRUCTURE_MB:
        return "too_big"

    out = workdir / pdb_id.upper()
    out.mkdir(parents=True, exist_ok=True)
    pdb_path = out / f"{pdb_id.upper()}.pdb"

    try:
        # (structure, chain_map, resname_map, ligand_resnames)
        st, chain_map, _resname_map, ligands = c2p.remap_structure(str(cif))
        observed = _observed_chains(st)
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
    if not reports:
        return None
    return _parse_report(reports[0], pdb_id), observed, chain_map


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
def _empty_artefact(fam: dict, holo_entries: int = 0, analysed: int = 0, failed: int = 0,
                    too_big: int = 0, seconds: float = 0.0) -> dict:
    """A complete, current artefact recording that there is no fingerprint to draw.

    Same shape as a full one, with the collections empty, so every reader of this artefact
    keeps working without learning a second format. `n_contacts` of 0 is the thing to test.
    """
    artefact = {
        "version": VERSION,
        "slug": fam["slug"],
        "built_at": int(time.time()),
        "entries_analysed": analysed,
        "entries_failed": failed,
        "entries_too_big": too_big,
        "max_structure_mb": MAX_STRUCTURE_MB,
        "seconds": seconds,
        "n_contacts": 0,
        # How many entries carry a designed ligand at all: 0 means the family is apo and no
        # amount of compute changes that, which is the distinction the panel needs to make.
        "holo_entries": holo_entries,
        "by_type": [],
        "hot_residues": [],
        "ligands": [],
        "fingerprint": {},
        "positions": [],
        "metal_coordination": {},
    }
    CONTACT_DIR.mkdir(parents=True, exist_ok=True)
    artefact_path(fam["slug"]).write_text(json.dumps(artefact, separators=(",", ":")))
    return artefact


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
        # Not a failure: this family genuinely has nothing ligand-bound, so there is no
        # fingerprint to compute and there never will be one. Recorded as an artefact rather
        # than returned as None, because None was indistinguishable from "PLIP is missing"
        # and the panel said the latter. Beta-lactamase inhibitory protein has 0 ligand-bound
        # entities of 8, and its Contacts tab advised running the job on a workstation.
        #
        # Writing it also settles the queue. The worker marks a family served whether or not
        # contacts succeeded, so a family that reports failure here is never retried and sat
        # on that message permanently; an artefact makes the answer durable and honest.
        return _empty_artefact(fam, holo_entries=len(holo))

    all_rows: list[dict] = []
    ok = failed = too_big = 0
    t0 = time.time()

    with tempfile.TemporaryDirectory(prefix="codswallop-plip-") as tmp:
        for i, m in enumerate(chosen):
            if progress:
                progress(i + 1, len(chosen), m["pdb_id"])
            got = contacts_for(m["pdb_id"], Path(tmp))
            if got == "too_big":
                # Counted apart from failures. A conversion that broke is a bug to chase; a
                # structure too large to hold in memory is a decision this made, and the two
                # should not be added together into one number nobody can act on.
                too_big += 1
                continue
            if got is None:
                failed += 1
                continue
            rows, observed, chain_map = got
            ok += 1
            # PLIP reports author residue numbering; the rest of the app speaks seed
            # coordinates. Mapped by ALIGNING each chain to the seed, per entry.
            #
            # This used to be `seed_pos = resnr + (query_beg - 1)`, and that is wrong twice
            # over. `query_beg` is where the entity's aligned region starts in the seed, so
            # it converts an entity sequence index; `resnr` is not an entity index, it is the
            # author number, which in a well-annotated entry already follows the canonical
            # numbering. The offset was therefore counted twice: on JAK1 (query_beg 879,
            # seed 1,154 residues) the hot residues came out at 1,340 and 2,110. It was
            # invisible on the family it was validated against, carbonic anhydrase, because
            # its query_beg is 1 and the offset is zero, so His94/His96/Thr199/Thr200 were
            # right for a reason that had nothing to do with the mapping being correct.
            # 52 of 71 built families were affected.
            # The member's own author chains, renamed the way cif2plip renamed them.
            member_chains = {chain_map.get(c, c) for c in (m.get("chains") or [])}
            mapping = _chain_mappings(observed, fam.get("seed_sequence") or "",
                                      only=member_chains or None)
            for r in rows:
                pos = mapping.get((r.get("reschain"), r["resnr"]))
                if pos is None:
                    # A chain that is not this protein, or a residue outside the aligned
                    # region. Dropped rather than guessed: a contact placed on the wrong
                    # residue is worse than a contact nobody counted.
                    continue
                r["seed_pos"] = pos
                all_rows.append(r)

    if not all_rows:
        # Entries were analysed and produced no contact we could place on the seed. Kept
        # apart from the case above by its counts, which are what say which happened: 0
        # analysed and 0 too big means there was nothing to do, while 40 analysed and 0
        # contacts means something is wrong and should be chased.
        return _empty_artefact(fam, holo_entries=len(holo), analysed=ok, failed=failed,
                               too_big=too_big, seconds=round(time.time() - t0, 1))

    # ---- the family fingerprint --------------------------------------------------------
    per_residue: Counter = Counter()
    per_type: Counter = Counter()
    lig_res: dict = defaultdict(Counter)
    # The residue named for a seed position is the SEED's residue there. Anything else is
    # inconsistent by construction: the position is a coordinate on the seed, so the letter
    # at it is not a matter of opinion. This was `setdefault`, which named the position after
    # whichever entry was processed first, and a single engineered variant could therefore
    # label it; a majority vote across the entries was tried and is worse, because a family
    # holds orthologues and mutants and the majority of them frequently do not carry the
    # seed's residue at all. What varies between structures is exactly what the conservation
    # track beside this panel already reports.
    restype_votes: dict = defaultdict(Counter)
    for r in all_rows:
        per_residue[r["seed_pos"]] += 1
        per_type[r["type"]] += 1
        lig_res[r["ligand"]][r["seed_pos"]] += 1
        restype_votes[r["seed_pos"]][r["restype"]] += 1

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

    _THREE = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
              "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
              "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
              "Y": "TYR", "V": "VAL"}
    seed_seq = fam.get("seed_sequence") or ""
    restype_of = {}
    for pos in restype_votes:
        letter = seed_seq[pos - 1] if 1 <= pos <= len(seed_seq) else ""
        # Falling back to what the structures reported only where the seed cannot say, which
        # after the mapping fix should be nowhere.
        restype_of[pos] = _THREE.get(letter) or restype_votes[pos].most_common(1)[0][0]
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
        # Entries skipped because the deposited file is too large to convert on this
        # machine. Reported so the panel can say the fingerprint is missing the biggest
        # complexes rather than implying they had no interactions.
        "entries_too_big": too_big,
        "max_structure_mb": MAX_STRUCTURE_MB,
        "seconds": round(time.time() - t0, 1),
        "n_contacts": len(all_rows),
        # Present on every artefact, empty or not, so a reader never has to treat its
        # absence as a value.
        "holo_entries": len(holo),
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
