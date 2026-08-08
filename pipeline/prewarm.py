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

TARGETS = [
    ("A0A0K8P6T7", "IsPETase"), ("A0A0K8P8E7", "MHETase"), ("P23458", "JAK1"),
    ("P00374", "DHFR"), ("P61823", "RNase A"), ("P02185", "Myoglobin"),
    ("P69905", "Haemoglobin alpha"), ("P02794", "Ferritin H"), ("P68431", "Histone H3.1"),
    ("P42212", "GFP"), ("P24941", "CDK2"), ("P02945", "Bacteriorhodopsin"),
    ("P29274", "A2A receptor"), ("P08100", "Rhodopsin"), ("P07900", "HSP90-alpha"),
    ("P56817", "BACE1"), ("P03372", "ER-alpha"), ("P00734", "Thrombin"),
    ("P15056", "BRAF"), ("P22303", "AChE"), ("P09874", "PARP1"), ("P00519", "ABL1"),
]


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
