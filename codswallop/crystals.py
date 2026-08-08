"""Crystallisation intelligence: what actually worked, across the whole family.

`_exptl_crystal_grow.pdbx_details` is free text a depositor typed, so this is parsing, not
reading a field. The text is genuinely varied ("100 mM MES pH 6.2-6.7, 40-100 mM ammonium
phosphate dibasic, 18-24% PEG400", "hanging drop against 2.4M ammonium sulfate, 0.1M Tris
pH 8.5"), and the parser's job is to recognise the chemistry rather than to normalise
everything into one schema.

Conservative by design. A component is reported only when a known name matches; anything
unrecognised stays in the verbatim text, which is always shown. A confidently wrong
precipitant would send somebody to set up the wrong screen.
"""

from __future__ import annotations

import re
from collections import Counter
from typing import Optional

PARSE_VERSION = 2

# --------------------------------------------------------------------------------------
# Chemistry
# --------------------------------------------------------------------------------------
# Precipitants, grouped so a family's answer reads "PEG" rather than fragmenting across
# eleven molecular weights. The PEG pattern is a regex because the archive spells it
# PEG400, PEG 400, PEG-400, peg4000, "polyethylene glycol 3350" and "PEG MME 550".
PRECIPITANTS: list[tuple[str, str, str]] = [
    ("PEG",              r"\bpeg[\s\-]?(?:mme[\s\-]?)?\d{3,5}\b|polyethylene\s*glycol", "polymer"),
    ("PEG (unspecified)", r"\bpeg\b(?![\s\-]?\d)",                                       "polymer"),
    ("Ammonium sulfate", r"ammonium\s*sulfate|ammonium\s*sulphate|\(nh4\)2so4",           "salt"),
    ("Sodium chloride",  r"sodium\s*chloride|\bnacl\b",                                   "salt"),
    ("Lithium sulfate",  r"lithium\s*sulfate|lithium\s*sulphate|li2so4",                   "salt"),
    ("Sodium citrate",   r"sodium\s*citrate",                                              "salt"),
    ("Ammonium phosphate", r"ammonium\s*phosphate",                                        "salt"),
    ("Magnesium chloride", r"magnesium\s*chloride|\bmgcl2\b",                              "salt"),
    ("Calcium chloride", r"calcium\s*chloride|\bcacl2\b",                                  "salt"),
    ("Sodium formate",   r"sodium\s*formate",                                              "salt"),
    ("Sodium malonate",  r"sodium\s*malonate",                                             "salt"),
    ("Tacsimate",        r"tacsimate",                                                     "salt"),
    ("MPD",              r"\bmpd\b|methylpentanediol|2-methyl-2,4-pentanediol",           "organic"),
    ("Isopropanol",      r"isopropanol|2-propanol|\bipa\b",                                "organic"),
    ("Ethanol",          r"\bethanol\b",                                                   "organic"),
    ("Jeffamine",        r"jeffamine",                                                     "polymer"),
    ("Dioxane",          r"dioxane",                                                       "organic"),
]

# Buffers, with the pH range each is actually useful over. The range is not used to correct
# a deposited pH: it is shown so a reader can see when a deposit used a buffer outside its
# range, which is common and worth noticing.
BUFFERS: list[tuple[str, str, tuple[float, float]]] = [
    ("Tris",       r"\btris\b(?!\w)",                       (7.0, 9.0)),
    ("HEPES",      r"\bhepes\b",                            (6.8, 8.2)),
    ("MES",        r"\bmes\b(?!\w)",                        (5.5, 6.7)),
    ("Sodium acetate", r"sodium\s*acetate|\bnaoac\b",       (3.6, 5.6)),
    ("Sodium cacodylate", r"cacodylate",                    (5.0, 7.4)),
    ("Bis-Tris",   r"bis[\s\-]?tris(?!\s*propane)",         (5.8, 7.2)),
    ("Bis-Tris propane", r"bis[\s\-]?tris\s*propane",       (6.3, 9.5)),
    ("Citrate",    r"\bcitrate\b|citric\s*acid",            (3.0, 6.2)),
    ("Phosphate",  r"phosphate\s*buffer|\bkh2po4\b|sodium\s*phosphate|potassium\s*phosphate", (5.8, 8.0)),
    ("Imidazole",  r"\bimidazole\b",                        (6.2, 7.8)),
    ("CHES",       r"\bches\b",                             (8.6, 10.0)),
    ("CAPS",       r"\bcaps\b",                             (9.7, 11.1)),
    ("Glycine",    r"\bglycine\b",                          (8.6, 10.6)),
    ("ADA",        r"\bada\b",                              (6.0, 7.2)),
    ("Succinate",  r"\bsuccinate\b",                        (4.8, 6.5)),
]

ADDITIVES: list[tuple[str, str]] = [
    ("Glycerol", r"\bglycerol\b"), ("Ethylene glycol", r"ethylene\s*glycol"),
    ("DTT", r"\bdtt\b|dithiothreitol"), ("TCEP", r"\btcep\b"),
    ("beta-mercaptoethanol", r"mercaptoethanol|\bbme\b"),
    ("Zinc", r"\bzinc\b|\bzncl2\b"), ("Magnesium", r"\bmagnesium\b|\bmgso4\b"),
    ("Calcium", r"\bcalcium\b"), ("Detergent", r"\bdetergent\b|\bldao\b|\bddm\b|\bc8e4\b|\bnonyl\b|octyl"),
    ("Monoolein", r"monoolein|\bmag\b|lipidic"), ("Cholesterol", r"cholesterol"),
    ("Spermine", r"spermine"), ("Sucrose", r"\bsucrose\b"), ("Trehalose", r"trehalose"),
]

# Setup methods, normalised. Depositors write these a dozen ways.
METHODS: list[tuple[str, str]] = [
    ("Vapour diffusion, hanging drop", r"hanging\s*drop"),
    ("Vapour diffusion, sitting drop", r"sitting\s*drop"),
    ("Vapour diffusion", r"vapor\s*diffusion|vapour\s*diffusion"),
    ("Lipidic cubic phase", r"lipidic\s*cubic|\blcp\b|in\s*meso"),
    ("Microbatch", r"microbatch|micro\s*batch"),
    ("Batch", r"\bbatch\b"),
    ("Dialysis", r"\bdialysis\b"),
    ("Free interface diffusion", r"free\s*interface|counter[\s\-]?diffusion"),
    ("Seeding", r"\bseeding\b|micro\s*seed"),
]

# A concentration immediately before a component name: "2.4 M ammonium sulfate", "18% PEG400".
_CONC = re.compile(
    r"(\d+(?:\.\d+)?)\s*(?:-\s*\d+(?:\.\d+)?\s*)?\s*(m\b|mm\b|%\s*(?:w/v|v/v)?|molar)",
    re.I)


def _match_all(text: str, table) -> list[str]:
    found = []
    for entry in table:
        name, pattern = entry[0], entry[1]
        if re.search(pattern, text, re.I):
            found.append(name)
    return found


def parse(crystal: Optional[dict]) -> Optional[dict]:
    """Turn one deposited crystallisation condition into structured chemistry."""
    if not crystal:
        return None
    text = (crystal.get("details") or "")
    method = crystal.get("method") or ""
    blob = f"{method} {text}"

    precipitants = _match_all(blob, PRECIPITANTS)
    # "PEG (unspecified)" only counts when no sized PEG was found: otherwise every
    # "20% PEG 3350" reports both.
    if "PEG" in precipitants and "PEG (unspecified)" in precipitants:
        precipitants.remove("PEG (unspecified)")

    buffers = _match_all(blob, BUFFERS)
    # Three buffers whose names nest inside each other, resolved most-specific-first.
    #
    # "Bis-Tris propane" contains "Bis-Tris" contains "Tris", and `\btris\b` matches inside
    # "Bis-Tris" because a word boundary falls between the hyphen and the T. Left alone, every
    # Bis-Tris condition was also tallied as Tris: on the beta-2 adrenergic receptor that put
    # about twelve conditions into the wrong row of the "what worked" table.
    if "Bis-Tris propane" in buffers and "Bis-Tris" in buffers:
        buffers.remove("Bis-Tris")
    if ("Bis-Tris propane" in buffers or "Bis-Tris" in buffers) and "Tris" in buffers:
        # Only if plain Tris is not ALSO named separately somewhere in the text.
        if not re.search(r"(?<!bis[\s\-])\btris\b(?!\w)", blob, re.I):
            buffers.remove("Tris")

    methods = _match_all(blob, METHODS)
    if "Vapour diffusion" in methods and any(
            m.startswith("Vapour diffusion,") for m in methods):
        methods.remove("Vapour diffusion")

    ph = crystal.get("ph")
    temp_k = crystal.get("temp_k")
    return {
        "ph": ph,
        "temp_k": temp_k,
        "temp_c": round(temp_k - 273.15, 1) if temp_k else None,
        "method": methods[0] if methods else (method.title() if method else None),
        "methods": methods,
        "precipitants": precipitants,
        "precipitant_classes": sorted({
            cls for name, _, cls in PRECIPITANTS if name in precipitants}),
        "buffers": buffers,
        "additives": _match_all(blob, ADDITIVES),
        "details": text,
        # Nothing recognised at all: the condition is shown verbatim and counted separately,
        # rather than being silently absent from every tally.
        "parsed": bool(precipitants or buffers or methods),
    }


def summarise(entries: list[dict]) -> dict:
    """What worked, across the family."""
    parsed = [p for p in (parse(e.get("crystal")) for e in entries) if p]
    if not parsed:
        return {"n": 0, "n_parsed": 0}

    def tally(key, res_of):
        """Count occurrences and carry the best resolution achieved with each."""
        counts: Counter = Counter()
        best: dict[str, float] = {}
        for p, e in zip(parsed, res_of):
            for v in p.get(key) or []:
                counts[v] += 1
                r = e.get("resolution")
                if r is not None and (v not in best or r < best[v]):
                    best[v] = r
        return [{"name": k, "count": n, "best_resolution": best.get(k)}
                for k, n in counts.most_common(24)]

    with_crystal = [e for e in entries if e.get("crystal")]
    phs = sorted(p["ph"] for p in parsed if p["ph"] is not None)
    temps = sorted(p["temp_c"] for p in parsed if p["temp_c"] is not None)

    return {
        "n": len(with_crystal),
        "n_parsed": sum(1 for p in parsed if p["parsed"]),
        "ph_min": phs[0] if phs else None,
        "ph_max": phs[-1] if phs else None,
        "ph_median": phs[len(phs) // 2] if phs else None,
        "temp_median": temps[len(temps) // 2] if temps else None,
        "precipitants": tally("precipitants", with_crystal),
        "buffers": tally("buffers", with_crystal),
        "additives": tally("additives", with_crystal),
        "methods": tally("methods", with_crystal),
        # One point per entry for the pH-vs-precipitant scatter, coloured by resolution.
        "points": [
            {"pdb_id": e["pdb_id"], "ph": p["ph"], "temp_c": p["temp_c"],
             "resolution": e.get("resolution"),
             "precipitant": (p["precipitants"] or ["unclassified"])[0],
             "buffer": (p["buffers"] or [None])[0]}
            for p, e in zip(parsed, with_crystal) if p["ph"] is not None
        ],
    }
