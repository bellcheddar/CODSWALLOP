#!/usr/bin/env python3
"""CODSWALLOP command line.

    python CODSWALLOP.py init                 create the cache database
    python CODSWALLOP.py serve [--port 8006]  run the development server
    python CODSWALLOP.py build <query>        assemble a family and print a summary
    python CODSWALLOP.py stats                cache statistics
    python CODSWALLOP.py purge                drop expired cache entries

`build` exists so a family can be warmed from a shell (or a cron job) without anyone having
to sit and watch a browser do it.
"""

from __future__ import annotations

import argparse
import sys

from codswallop import config, db


def cmd_init(_args) -> int:
    config.ensure_dirs()
    db.init()
    print(f"Initialised {config.DB_PATH}")
    return 0


def cmd_serve(args) -> int:
    from codswallop.webapp import create_app

    host, _, port = args.bind.rpartition(":")
    create_app().run(host=host or "127.0.0.1", port=int(port), debug=args.debug)
    return 0


def cmd_build(args) -> int:
    from codswallop import family, resolve

    db.init()
    result = resolve.resolve(args.query)

    if result["status"] == "ambiguous":
        print(f"{result['prompt']}\n")
        for c in result["candidates"]:
            print(f"  {c['pick']:<12} {c['label']}"
                  f"{'  [' + c['sublabel'] + ']' if c.get('sublabel') else ''}")
        print("\nRe-run with one of the identifiers in the first column.")
        return 2
    if result["status"] != "resolved":
        print(result.get("message", "Could not resolve that."), file=sys.stderr)
        return 1

    seed = result["seed"]
    print(f"Seed: {seed.get('seed') or 'pasted sequence'} "
          f"({seed.get('name') or '?'}, {seed['length']} aa)")
    print(f"      {seed['note']}\n")

    fam = family.get_or_build(seed, args.query, force=args.refresh)
    s = fam["stats"]
    cov = s["coverage"]

    print(f"/f/{fam['slug']}")
    print(f"  {s['entries']:,} entries, {s['entities']:,} polymer entities, "
          f"{s['constructs']:,} distinct constructs")
    print(f"  {s['organisms']} organisms, {s['ligands']:,} distinct ligands, "
          f"{s['holo_entries']:,} ligand-bound")
    print(f"  identity {s['identity_min']}-{s['identity_max']} %, "
          f"resolution {s['best_resolution']}-{s['median_resolution']} A (best/median)")
    print(f"  released {s['first_release']} to {s['latest_release']}")
    print(f"  methods: {', '.join(f'{m} {n:,}' for m, n in s['methods'])}")
    print(f"  median construct covers {cov['median_coverage']} % of the seed; "
          f"{cov['thin_pct']} % of it is thinly covered")
    if cov["gaps"]:
        g = cov["gaps"][0]
        print(f"  poorest region: residues {g['start']}-{g['end']} "
              f"(in {g['depth']} of {s['entities']:,} constructs)")
    if fam["truncated"]:
        print(f"  NOTE: truncated -- the PDB holds {fam['total_hits']:,} matching entities, "
              f"the closest {s['entities']:,} were filed")
    return 0


def cmd_stats(_args) -> int:
    db.init()
    for k, v in db.stats().items():
        print(f"  {k:<18} {v:,}" if isinstance(v, int) else f"  {k:<18} {v}")
    return 0


def cmd_purge(_args) -> int:
    db.init()
    print(f"Dropped {db.purge_expired():,} expired cache entries.")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="CODSWALLOP.py", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("init", help="create the cache database").set_defaults(fn=cmd_init)

    sp = sub.add_parser("serve", help="run the development server")
    sp.add_argument("--bind", default=config.BIND_ADDR)
    sp.add_argument("--debug", action="store_true")
    sp.set_defaults(fn=cmd_serve)

    sp = sub.add_parser("build", help="assemble a family and print a summary")
    sp.add_argument("query", help="PDB ID, UniProt accession, gene, Pfam/InterPro, sequence or name")
    sp.add_argument("--refresh", action="store_true", help="rebuild even if cached")
    sp.set_defaults(fn=cmd_build)

    sub.add_parser("stats", help="cache statistics").set_defaults(fn=cmd_stats)
    sub.add_parser("purge", help="drop expired cache entries").set_defaults(fn=cmd_purge)

    args = p.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
