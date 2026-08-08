"""Which of a family's bound components are actually drugs, and how far each one got.

The family already knows every chemical component bound in every entry. This asks the next
question, which is the one a person starting on a target actually has: **has anyone put a
drug in this pocket, and did it reach patients?**

Three sources, each doing the one thing it is best at:

* **DrugBank, via the RCSB's own chem_comp annotation.** The link from a PDB chemical
  component to a drug already exists in the archive, batched and cached alongside everything
  else this app fetches, so no new service is needed to answer "is this a drug". It also
  carries the generic name, the brand names and the ATC codes.
* **ATC, the WHO's own classification**, for the therapeutic grouping. Inventing a set of
  categories would mean defending it; ATC is the standard the field already uses, and its
  top level is close to what anyone would have invented anyway. J and L are split at the
  second level, because "anti-infective" hides the difference between an antibacterial and
  an antiviral, and "antineoplastic and immunomodulating" hides the one between an oncology
  drug and an immunosuppressant.
* **ClinicalTrials.gov**, for the trial phase of anything not yet approved. DrugBank says
  whether a drug is approved, investigational or experimental, and those are three useful
  words, but none of them is a phase number.

The stage ladder is deliberately ordered so that the strongest evidence wins: an approved
drug is approved whatever its trials say, and a compound nobody has ever taken into a trial
is preclinical rather than unknown.
"""

from __future__ import annotations

import logging
import urllib.parse
from collections import Counter
from typing import Optional

from . import config, db, http

logger = logging.getLogger(__name__)

DRUG_VERSION = 1

# How many drugs get an external lookup per family. The DrugBank annotation is one batched
# call for all of them; these are one request each, so they are the cost.
TRIAL_LOOKUP_CAP = 40

CTGOV_URL = "https://clinicaltrials.gov/api/v2/studies"
NCATS_URL = "https://drugs.ncats.io/api/v1/substances/search"

_DRUG_QUERY = """
query($ids: [String!]!) {
  chem_comps(comp_ids: $ids) {
    rcsb_id
    chem_comp { id name formula }
    rcsb_chem_comp_target { name interaction_type provenance_source target_actions }
    drugbank { drugbank_info {
      drugbank_id name description drug_groups drug_categories synonyms atc_codes
      brand_names indication mechanism_of_action } }
  }
}
"""

# ATC level 1, with J and L split at level 2. The labels are the ones a structural biologist
# would use rather than the WHO's own wording, which is written for pharmacy stock control.
_ATC_LEVEL1 = {
    "A": "Metabolic & gastrointestinal",
    "B": "Blood & coagulation",
    "C": "Cardiovascular",
    "D": "Dermatological",
    "G": "Genito-urinary & sex hormones",
    "H": "Hormonal (systemic)",
    "M": "Musculoskeletal",
    "N": "Neurological & psychiatric",
    "P": "Antiparasitic",
    "R": "Respiratory",
    "S": "Sensory organs",
    "V": "Various",
}
_ATC_LEVEL2 = {
    "J01": "Antibacterial", "J02": "Antifungal", "J04": "Antimycobacterial",
    "J05": "Antiviral", "J06": "Immune sera", "J07": "Vaccines",
    "L01": "Oncology", "L02": "Oncology (endocrine)",
    "L03": "Immunology (stimulant)", "L04": "Immunology (suppressant)",
}
UNCLASSED = "Unclassified"


def therapeutic_classes(atc_codes: Optional[list]) -> list[str]:
    """Every therapeutic class a drug belongs to, from its ATC codes.

    All of them, not the first. A drug carries several ATC codes because it is genuinely
    used for several things: acetazolamide is both a glaucoma drug (S01EC01) and a urinary
    agent (G01AE10), and taking whichever came first filed the best-known carbonic anhydrase
    inhibitor in the world under "genito-urinary". A drug that belongs in two classes appears
    in two, which is the true answer and the one a reader can check against the codes shown
    beside it.
    """
    out: list[str] = []
    for code in atc_codes or []:
        code = (code or "").upper().strip()
        name = None
        if len(code) >= 3 and code[:3] in _ATC_LEVEL2:
            name = _ATC_LEVEL2[code[:3]]
        elif code[:1] == "J":
            name = "Anti-infective"
        elif code[:1] == "L":
            name = "Oncology & immunology"
        elif code[:1] in _ATC_LEVEL1:
            name = _ATC_LEVEL1[code[:1]]
        if name and name not in out:
            out.append(name)
    return out[:3] or [UNCLASSED]


def _base_protein(name: str) -> str:
    """A protein name with its isoform number stripped, for comparing across a family.

    "Carbonic anhydrase 2" and "Carbonic anhydrase 12" are different proteins and must not
    match each other, so the number is removed rather than ignored: comparing the full
    strings by substring makes 2 a match for 12.
    """
    import re as _re
    n = _re.sub(r"[^a-z0-9 ]+", " ", (name or "").lower())
    n = _re.sub(r"\b(isoform|type|subunit|chain|precursor)\b", " ", n)
    n = _re.sub(r"\b[ivx]+\b|\b\d+\b", " ", n)          # isoform numbers and numerals
    return " ".join(n.split())


def on_target(targets: list, fam: dict, text: str = "") -> Optional[str]:
    """How this drug is tied to the family's protein, or None.

    Returns "annotated" when DrugBank lists the protein as a target, "described" when the
    drug's own mechanism or indication text names it, and None otherwise.

    The second signal exists because the first has a coverage gap that matters. DrugBank's
    target lists lag for recently approved drugs, and the RCSB snapshot lags again: across
    all 465 candidate components in the KRAS family, not one is annotated against KRas, on a
    protein with two approved G12C inhibitors. A panel reporting "0 on target" for KRAS
    would be worse than no panel. The text match is weaker and is labelled as such rather
    than being folded into the same word.
    """
    want = _base_protein(fam.get("name") or "")
    if not want:
        return None
    for t in targets or []:
        got = _base_protein(t.get("name") or "")
        if got and (got == want or got in want or want in got):
            return "annotated"
    if want and len(want) > 4 and want in _base_protein(text or ""):
        return "described"
    return None


def _on_target_legacy(targets: list, fam: dict) -> bool:
    """Is this family's own protein among the drug's annotated targets?

    The difference between a drug that was designed for this pocket and a molecule that
    happens to be sitting in it. Acetazolamide's targets are carbonic anhydrases 1, 2, 4 and
    12, so it is on target for a carbonic anhydrase family; glucose's are glycogen
    phosphorylase and the glucosidases, so in the same family it is a cryoprotectant that
    DrugBank happens to know as an intravenous solution.

    Compared at the base protein name, because a family assembled at 30 % identity spans
    isoforms on purpose: a drug targeting CA12 is on target for a page about CA2.
    """
    want = _base_protein(fam.get("name") or "")
    if not want:
        return False
    for t in targets or []:
        got = _base_protein(t.get("name") or "")
        if got and (got == want or got in want or want in got):
            return True
    return False


# The stages, strongest evidence first. `withdrawn` sits above the trial phases deliberately:
# a drug that reached the market and was pulled is a more useful thing to know about a pocket
# than the phase of the trials it ran on the way there.
STAGE_ORDER = ["Approved", "Withdrawn", "Phase 4", "Phase 3", "Phase 2", "Phase 1",
               "Preclinical", "Unknown"]

_PHASE_LABEL = {
    "PHASE4": "Phase 4", "PHASE3": "Phase 3", "PHASE2": "Phase 2",
    "PHASE1": "Phase 1", "EARLY_PHASE1": "Phase 1",
}


def highest_trial_phase(name: str) -> Optional[str]:
    """The highest phase any registered trial has taken this drug to.

    Free-text on the intervention name, which is what the registry indexes. It is a ceiling
    rather than a status: a drug appearing in one phase 3 trial has reached phase 3, and
    that is all the claim being made.
    """
    if not name:
        return None

    def fetch():
        url = (f"{CTGOV_URL}?query.intr={urllib.parse.quote(name)}"
               "&fields=protocolSection.designModule.phases&pageSize=200")
        try:
            body = http.get_json(url)
        except Exception:                       # noqa: BLE001
            logger.info("no trial data for %r", name)
            return None
        seen = Counter()
        for st in (body or {}).get("studies") or []:
            design = ((st.get("protocolSection") or {}).get("designModule") or {})
            for ph in design.get("phases") or []:
                if ph in _PHASE_LABEL:
                    seen[_PHASE_LABEL[ph]] += 1
        if not seen:
            return None
        best = min(seen, key=lambda s: STAGE_ORDER.index(s))
        return {"phase": best, "n_studies": sum(seen.values()), "by_phase": seen.most_common()}

    return db.cached(("ctgov_phase", DRUG_VERSION, name.lower()), fetch)


def fda_record(name: str) -> Optional[dict]:
    """The NIH/NCATS substance record, which carries the FDA approval identifier.

    Used only to *confirm* an approval that DrugBank already asserts, and to give the reader
    a UNII to look up. It is not used to decide the stage: the search is by name, and a name
    search is not evidence about a molecule.
    """
    def fetch():
        try:
            body = http.get_json(f"{NCATS_URL}?q={urllib.parse.quote(name)}&top=1")
        except Exception:                       # noqa: BLE001
            return None
        rows = (body or {}).get("content") or []
        if not rows:
            return None
        r = rows[0]
        return {"unii": r.get("approvalID"), "approved_by": r.get("approvedBy"),
                "status": r.get("status"), "name": (r.get("_name") or "")}

    return db.cached(("ncats", DRUG_VERSION, name.lower()), fetch)


def annotate(comp_ids: list[str]) -> dict:
    """DrugBank and target annotation for a set of chemical components, batched."""
    out: dict = {}
    ids = [c for c in comp_ids if c]
    for i in range(0, len(ids), config.GRAPHQL_BATCH):
        batch = sorted(ids[i:i + config.GRAPHQL_BATCH])

        def fetch(batch=batch):
            body = http.graphql(config.RCSB_GRAPHQL_URL, _DRUG_QUERY, {"ids": batch})
            rows = {}
            for c in (body or {}).get("chem_comps") or []:
                info = ((c.get("drugbank") or {}).get("drugbank_info")) or {}
                rows[c["rcsb_id"]] = {
                    "drugbank_id": info.get("drugbank_id"),
                    "generic": info.get("name"),
                    "groups": info.get("drug_groups") or [],
                    "atc": info.get("atc_codes") or [],
                    "brands": info.get("brand_names") or [],
                    "categories": (info.get("drug_categories") or [])[:8],
                    "indication": (info.get("indication") or "")[:400] or None,
                    "mechanism": (info.get("mechanism_of_action") or "")[:400] or None,
                    "targets": [
                        {"name": t.get("name"), "action": (t.get("target_actions") or [None])[0],
                         "source": t.get("provenance_source")}
                        for t in (c.get("rcsb_chem_comp_target") or [])[:6] if t.get("name")
                    ],
                }
            return rows

        out.update(db.cached(("drugbank", DRUG_VERSION, batch), fetch) or {})
    return out


def _stage(groups: list, trial: Optional[dict]) -> tuple:
    """(stage, why). The ladder in one place, so the panel never has to reason about it."""
    g = {str(x).lower() for x in groups or []}
    if "approved" in g:
        return "Approved", "DrugBank lists it as approved"
    if "withdrawn" in g:
        return "Withdrawn", "approved once and withdrawn"
    if trial and trial.get("phase"):
        return trial["phase"], f"{trial['n_studies']} registered trials, highest {trial['phase']}"
    if "investigational" in g:
        return "Unknown", "investigational, with no registered trial phase found"
    if "experimental" in g:
        return "Preclinical", "DrugBank lists it as experimental"
    if "vet_approved" in g:
        return "Approved", "approved for veterinary use only"
    return "Unknown", "no development stage could be established"


def build(fam: dict, max_drugs: int = 120) -> dict:
    """The drugs bound anywhere in this family, grouped by therapeutic class."""
    comps = ((fam.get("ligands") or {}).get("components")) or []
    # Only the components that could plausibly be a drug. A buffer or a cryoprotectant is
    # never one, and asking DrugBank about 800 sulfate ions wastes everybody's time.
    candidates = [c for c in comps
                  if (c.get("klass") or "").split("/")[0] in ("ligand", "cofactor")]
    candidates.sort(key=lambda c: -(c.get("count") or 0))
    candidates = candidates[:max_drugs]
    if not candidates:
        return {"n": 0, "classes": [], "drugs": []}

    ann = annotate([c["id"] for c in candidates])

    drugs = []
    for c in candidates:
        a = ann.get(c["id"]) or {}
        if not a.get("drugbank_id"):
            continue                            # a bound molecule, but not a known drug
        trial = None
        # Only for the ones whose stage is still open: an approved drug's phase history is
        # not what the panel is for, and this is a request per drug.
        groups = {str(x).lower() for x in a.get("groups") or []}
        if "approved" not in groups and "withdrawn" not in groups:
            trial = highest_trial_phase(a.get("generic") or "")
        stage, why = _stage(a.get("groups"), trial)
        classes = therapeutic_classes(a.get("atc"))
        tie = on_target(a.get("targets") or [], fam,
                        (a.get("mechanism") or "") + " " + (a.get("indication") or ""))
        fda = fda_record(a["generic"]) if stage == "Approved" and a.get("generic") else None

        # Brand names are a long tail of pack sizes and national repackagings ("Diamox Tablets
        # 250mg", "Novo-zolamide Tab 250mg"). Deduplicated on the leading word, which is the
        # brand, and capped: the panel wants Gleevec, not forty ways of writing it.
        brands, seen_brand = [], set()
        for b in a.get("brands") or []:
            key = (b or "").split()[0].lower() if b else ""
            if key and key not in seen_brand and key != (a.get("generic") or "").lower():
                seen_brand.add(key)
                brands.append(b.split(" ")[0] if len(b.split()) > 1 else b)

        drugs.append({
            "id": c["id"],
            "generic": a.get("generic") or c.get("name"),
            "brands": brands[:8],
            "drugbank_id": a["drugbank_id"],
            "stage": stage,
            "stage_why": why,
            "klass": classes[0],
            "classes": classes,
            "on_target": bool(tie),
            "tie": tie,
            "atc": (a.get("atc") or [])[:4],
            "entries": c.get("count") or 0,
            "best_resolution": c.get("best_resolution"),
            "targets": a.get("targets") or [],
            "indication": a.get("indication"),
            "mechanism": a.get("mechanism"),
            "categories": a.get("categories") or [],
            "unii": (fda or {}).get("unii"),
            "trials": trial,
        })

    drugs.sort(key=lambda d: (STAGE_ORDER.index(d["stage"]), -d["entries"]))

    # A drug in two classes appears in two: the grouping reflects what the drug is used for,
    # not a choice made on the reader's behalf.
    by_class: dict = {}
    for d in drugs:
        for k in d["classes"]:
            by_class.setdefault(k, []).append(d)
    classes = [{"name": k, "n": len(v),
                "approved": sum(1 for x in v if x["stage"] == "Approved"),
                "drugs": v}
               for k, v in sorted(by_class.items(),
                                  key=lambda kv: (kv[0] == UNCLASSED, -len(kv[1])))]

    return {
        "n": len(drugs),
        "n_approved": sum(1 for d in drugs if d["stage"] == "Approved"),
        "n_on_target": sum(1 for d in drugs if d["on_target"]),
        "n_annotated": sum(1 for d in drugs if d.get("tie") == "annotated"),
        "n_candidates": len(candidates),
        "by_stage": Counter(d["stage"] for d in drugs).most_common(),
        "classes": classes,
        "drugs": drugs,
    }
