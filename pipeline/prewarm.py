#!/usr/bin/env python3
"""Pre-warm a list of targets: family, embedding and PLIP contacts, one at a time.

Resumable by design. A target whose artefacts are already current is skipped in
milliseconds, so restarting after an interruption costs nothing and no work is repeated.
The first version of this was a heredoc that died silently at target 7 of 22 with no error,
no traceback and no crash report, and restarting it would have redone everything.

Sequential on purpose: each stage already parallelises its own HTTP, and TM-align is
CPU-bound, so running families concurrently only contends for the same two resources.

    python pipeline/prewarm.py [--reps 50] [--contacts 20]
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from codswallop import artefacts, contacts, db, embed, family, resolve  # noqa: E402

# Every family worth having warm before anyone asks for it.
#
# Two kinds, and the first kind is not optional. The landing page offers five example
# buttons, and a reader who presses one and lands on a placeholder map has been invited down
# a path the app did not prepare: 4HHB is on the front page, and going in that way built
# `hemoglobin-subunit-alpha-4hhb-1`, a different family from the pre-warmed
# `hemoglobin-subunit-alpha-p69905`, because a slug derives from its seed. Same protein, two
# addresses, artefacts on only one.
#
# The ambiguous examples are listed by the pick a reader would choose, not by the string on
# the button: "4HHB" and "TEM-1 beta-lactamase" both answer with a disambiguation card, and
# pre-warming the question warms nothing.
# Keyed by the exact text on the button, so a test can check the template against it and
# fail the moment an example is added to the front page without being warmed here.
EXAMPLE_SEEDS = {
    "1AKI": ["1AKI_1"],
    "P00918": ["P00918"],
    "PF00062": ["135L_1"],
    # Both of these answer with a disambiguation card rather than a family, so the seed is
    # the pick a reader would make. Warming the question warms nothing.
    "TEM-1 beta-lactamase": ["P62593"],
    "4HHB": ["4HHB_1", "4HHB_2"],
}

LANDING_EXAMPLES = [(seed, f"{label} (example)")
                    for label, seeds in EXAMPLE_SEEDS.items() for seed in seeds]

TARGETS = LANDING_EXAMPLES + [
    # -- Marc's own work, and the PET-degrading enzymes --------------------------------
    ("A0A0K8P6T7", "IsPETase"), ("A0A0K8P8E7", "MHETase"),

    # -- classic teaching structures: the ones people arrive already knowing -----------
    ("P00374", "DHFR"), ("P61823", "RNase A"), ("P02185", "Myoglobin"),
    ("P69905", "Haemoglobin alpha"), ("P02794", "Ferritin H"), ("P68431", "Histone H3.1"),
    ("P42212", "GFP"), ("P00760", "Trypsin"), ("P01308", "Insulin"),
    ("P0CG48", "Ubiquitin"), ("P0DP23", "Calmodulin"), ("P22629", "Streptavidin"),
    ("P02768", "Serum albumin"), ("P68133", "Actin"), ("Q71U36", "Tubulin alpha"),

    # -- kinases: the most-prosecuted target class there is -----------------------------
    ("P24941", "CDK2"), ("P23458", "JAK1"), ("P15056", "BRAF"), ("P00519", "ABL1"),
    ("P00533", "EGFR"), ("Q16539", "p38 MAPK"), ("P31749", "AKT1"), ("P11362", "FGFR1"),
    ("P04629", "TrkA"), ("Q06124", "SHP2"),

    # -- membrane proteins and receptors: the hard ones, and the pretty ones ------------
    ("P02945", "Bacteriorhodopsin"), ("P29274", "A2A receptor"), ("P08100", "Rhodopsin"),
    ("P07550", "beta-2 adrenergic"), ("B4ZY91", "GLP-1R"), ("P13569", "CFTR"),
    ("P42345", "mTOR"), ("Q9Y5Y9", "Nav1.8"),

    # -- oncology and epigenetics -------------------------------------------------------
    ("P07900", "HSP90-alpha"), ("P09874", "PARP1"), ("P04637", "p53"), ("P01116", "KRAS"),
    ("O60885", "BRD4"), ("O00255", "Menin"), ("P10275", "Androgen receptor"),
    ("P03372", "ER-alpha"),

    # -- infectious disease -------------------------------------------------------------
    ("P0DTD1", "SARS-CoV-2 nsp5"), ("P0DTC2", "SARS-CoV-2 spike"),
    ("P62593", "TEM-1 beta-lactamase"), ("P00722", "beta-galactosidase"),
    ("P9WGR1", "M. tuberculosis InhA"),

    # -- metabolic, cardiovascular and neurodegeneration ---------------------------------
    ("P56817", "BACE1"), ("P00734", "Thrombin"), ("P22303", "AChE"), ("P00742", "Factor Xa"),
    ("P27487", "DPP-4"), ("P08684", "CYP3A4"), ("P00918", "Carbonic anhydrase II"),
    ("P37840", "alpha-synuclein"), ("P10636", "Tau"), ("P00441", "SOD1"),
    ("P02766", "Transthyretin"), ("P05067", "APP"),
]

# Deduplicated, keeping the first mention: several proteins are legitimately in two groups
# above (TEM-1 is both a landing example and an infectious-disease target) and building a
# family twice is only slower, never different.
_seen = set()
TARGETS = [t for t in TARGETS if not (t[0] in _seen or _seen.add(t[0]))]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=50)
    ap.add_argument("--contacts", type=int, default=20)
    ap.add_argument("--force", action="store_true", help="rebuild even when current")
    args = ap.parse_args()

    db.init()
    done = skipped = failed = 0
    for i, (acc, label) in enumerate(TARGETS, 1):
        t0 = time.time()
        try:
            r = resolve.resolve(acc)
            if r["status"] != "resolved":
                print(f"[{i:>2}/{len(TARGETS)}] {label:<20} SKIP "
                      f"({r.get('message') or 'ambiguous'})", flush=True)
                skipped += 1
                continue
            slug = family.slug_for(r["seed"])

            # The cheap part: if both artefacts are current, this target is finished.
            st = artefacts.status(slug)
            if not args.force and st["embedding"]["current"] and st["contacts"]["current"]:
                print(f"[{i:>2}/{len(TARGETS)}] {label:<20} already current", flush=True)
                skipped += 1
                continue

            fam = family.get_or_build(r["seed"], acc)
            built = time.time() - t0

            a = embed.load(slug) if st["embedding"]["current"] and not args.force else None
            if a is None:
                a = embed.build(fam, max_representatives=args.reps)
            c = contacts.load(slug) if st["contacts"]["current"] and not args.force else None
            if c is None:
                c = contacts.build(fam, max_entries=args.contacts)

            af = (a or {}).get("alphafold")
            print(f"[{i:>2}/{len(TARGETS)}] {label:<20} {fam['stats']['entries']:>5} entries  "
                  f"build {built:>4.0f}s  emb {(a or {}).get('n_representatives', 0):>3} reps"
                  + (f" AF {af['tm']:<5}" if af else " AF -    ")
                  + f"  contacts {(c or {}).get('n_contacts', 0):>5}"
                  f"  total {time.time() - t0:>4.0f}s", flush=True)
            done += 1
        except Exception as exc:                      # one bad target must not stop the rest
            print(f"[{i:>2}/{len(TARGETS)}] {label:<20} FAILED {exc}", flush=True)
            failed += 1

    print(f"\n{done} built, {skipped} already current, {failed} failed.", flush=True)
    return 1 if failed and not done else 0


if __name__ == "__main__":
    raise SystemExit(main())
