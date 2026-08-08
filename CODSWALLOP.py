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
import time

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


def cmd_embed(args) -> int:
    """Compute the pairwise TM-score matrix and the map positions it implies.

    Workstation only: this downloads mmCIF files and does real numerical work. The droplet
    reads the JSON artefact it writes; push it with deploy/push_embeddings.sh.
    """
    from codswallop import embed, family, resolve

    db.init()
    result = resolve.resolve(args.query)
    if result["status"] != "resolved":
        print(result.get("message") or result.get("prompt"), file=sys.stderr)
        return 1
    fam = family.get_or_build(result["seed"], args.query)
    print(f"{fam['name']}: {fam['stats']['entries']} entries, "
          f"{len(fam['constructs'])} distinct constructs")

    last = [0.0]

    def progress(stage, i, n, label):
        if stage == "fetch":
            if time.time() - last[0] > 1.0 or i == n:
                last[0] = time.time()
                print(f"\r  fetching structures {i}/{n} {label:<8}", end="", flush=True)
        else:
            print(f"\r  aligning {n:,} pairs...{' ' * 20}", end="", flush=True)

    t0 = time.time()
    art = embed.build(fam, max_representatives=args.max or embed.MAX_REPRESENTATIVES,
                      progress=progress)
    print()
    if not art:
        print("  Not enough usable structures to embed.", file=sys.stderr)
        return 1
    print(f"  {art['n_representatives']} representatives, {art['n_pairs']:,} pairs, "
          f"median TM {art['median_tm']}, {time.time() - t0:.0f}s total")
    print(f"  wrote {embed.artefact_path(fam['slug'])}")
    return 0


def cmd_contacts(args) -> int:
    """Run PLIP over a family's ligand-bound entries and store the fingerprint.

    Workstation only: this shells out to PLIP and OpenBabel. Push the artefact with
    deploy/push_embeddings.sh, which ships data/contacts/ alongside data/embeddings/.
    """
    from codswallop import contacts, family, resolve

    db.init()
    result = resolve.resolve(args.query)
    if result["status"] != "resolved":
        print(result.get("message") or result.get("prompt"), file=sys.stderr)
        return 1
    fam = family.get_or_build(result["seed"], args.query)
    print(f"{fam['name']}: {fam['stats']['holo_entries']} ligand-bound entities")

    last = [0.0]

    def progress(i, n, pdb_id):
        if time.time() - last[0] > 1.0 or i == n:
            last[0] = time.time()
            print(f"\r  PLIP {i}/{n} {pdb_id:<8}", end="", flush=True)

    t0 = time.time()
    art = contacts.build(fam, max_entries=args.max, progress=progress)
    print()
    if not art:
        print("  No contacts found (nothing ligand-bound, or every conversion failed).",
              file=sys.stderr)
        return 1
    print(f"  {art['entries_analysed']} entries analysed ({art['entries_failed']} failed), "
          f"{art['n_contacts']:,} contacts, {time.time() - t0:.0f}s")
    print(f"  types: {', '.join(f'{k} {v}' for k, v in art['by_type'][:5])}")
    hot = art["hot_residues"][:6]
    print(f"  hottest residues: {', '.join(h['restype'] + str(h['pos']) for h in hot)}")
    print(f"  wrote {contacts.artefact_path(fam['slug'])}")
    return 0


def cmd_warm(args) -> int:
    """Pre-build families, so a cold cache is paid for by cron rather than by a reader.

    With no arguments it refreshes every family already in the cache, which is what the
    weekly timer wants: the PDB releases on Wednesdays, and a family that was assembled last
    week is stale in exactly the way a rebuild fixes.
    """
    from codswallop import family, resolve

    db.init()
    queries = args.queries or [f["query"] for f in db.recent_families(limit=200) if f.get("query")]
    if not queries:
        # recent_families does not carry the query, so fall back to the stored one per slug.
        conn = db.connect()
        queries = [r["query"] for r in conn.execute("SELECT query FROM family") if r["query"]]
    if not queries:
        print("Nothing filed yet, and no queries given.")
        return 0

    ok = failed = 0
    for q in queries:
        t0 = time.time()
        try:
            result = resolve.resolve(q)
            if result["status"] != "resolved":
                print(f"  skip  {q!r}: {result.get('message') or 'ambiguous'}")
                continue
            fam = family.get_or_build(result["seed"], q, force=args.force)
            print(f"  warm  {fam['slug']:<44} {fam['stats']['entries']:>5} entries  "
                  f"{time.time() - t0:>5.1f}s")
            ok += 1
        except Exception as exc:            # one bad family must not stop the rest
            print(f"  FAIL  {q!r}: {exc}")
            failed += 1
    print(f"\n{ok} warmed, {failed} failed.")

    rebuilt = _warm_artefacts(args)
    _report_artefacts()
    return 1 if (failed and not ok) or rebuilt is False else 0


def _pipeline_available() -> bool:
    """Whether this machine can build artefacts at all.

    The droplet cannot: biotite, tmtools, PLIP and OpenBabel are deliberately absent, which
    is the point of the split. So `warm` reports staleness everywhere and only rebuilds
    where the pipeline exists.
    """
    import importlib.util
    return all(importlib.util.find_spec(m) for m in ("biotite", "tmtools", "numpy"))


def _warm_artefacts(args) -> bool:
    """Rebuild any embedding (and, if asked, contacts) that is missing or out of date.

    This is the step that stops an artefact version bump from silently leaving families on
    the placeholder map. Three bumps in one afternoon each did exactly that, and nothing
    anywhere said so: the page renders either way.
    """
    from codswallop import artefacts

    if args.no_artefacts:
        return True
    if not _pipeline_available():
        return True

    from codswallop import embed, family, resolve

    todo = artefacts.stale("embedding")
    if not todo:
        return True
    print(f"\n{len(todo)} embedding(s) missing or out of date "
          f"(pipeline v{artefacts.embed_io.VERSION}); rebuilding:")
    every_ok = True
    for s in todo:
        if not s["query"]:
            print(f"  skip  {s['slug']}: no stored query to rebuild from")
            continue
        t0 = time.time()
        try:
            r = resolve.resolve(s["query"])
            if r["status"] != "resolved":
                print(f"  skip  {s['slug']}: {r.get('message') or 'ambiguous'}")
                continue
            fam = family.get_or_build(r["seed"], s["query"])
            art = embed.build(fam, max_representatives=args.max_reps)
            if not art:
                print(f"  skip  {s['slug']}: not enough usable structures")
                continue
            af = art.get("alphafold")
            print(f"  embed {s['slug']:<44} {art['n_representatives']:>3} reps"
                  + (f", AF TM {af['tm']}" if af else ", no AF model")
                  + f"  {time.time() - t0:>5.0f}s")
        except Exception as exc:
            print(f"  FAIL  {s['slug']}: {exc}")
            every_ok = False
    return every_ok


def _report_artefacts() -> None:
    """Say what is still on the placeholder. Runs on the droplet too, where nothing can be
    rebuilt, because a silent fallback is the failure this exists to prevent."""
    from codswallop import artefacts

    summary = artefacts.summary()
    print(f"\nArtefacts: {summary['embeddings_current']}/{summary['families']} embeddings "
          f"current (v{summary['embedding_version']}), "
          f"{summary['contacts_current']}/{summary['families']} contacts "
          f"(v{summary['contacts_version']}).")
    if summary["on_placeholder"]:
        where = "" if _pipeline_available() else " (this machine cannot rebuild them: " \
                                                 "run `CODSWALLOP.py warm` on a workstation)"
        print(f"WARNING: {summary['on_placeholder']} famil"
              f"{'y is' if summary['on_placeholder'] == 1 else 'ies are'} showing the "
              f"placeholder map{where}.")


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

    sp = sub.add_parser("embed", help="compute the structural embedding (workstation only)")
    sp.add_argument("query", help="the family to embed")
    sp.add_argument("--max", type=int, default=None,
                    help="cap on representative structures (default 260)")
    sp.set_defaults(fn=cmd_embed)

    sp = sub.add_parser("contacts", help="run PLIP family-wide (workstation only)")
    sp.add_argument("query", help="the family to profile")
    sp.add_argument("--max", type=int, default=60, help="cap on entries analysed")
    sp.set_defaults(fn=cmd_contacts)

    sp = sub.add_parser("warm", help="pre-build families so the first visitor waits for nothing")
    sp.add_argument("queries", nargs="*", help="queries to warm (default: everything already filed)")
    sp.add_argument("--force", action="store_true", help="rebuild even if still fresh")
    sp.add_argument("--no-artefacts", action="store_true",
                    help="skip rebuilding stale embeddings (report only)")
    sp.add_argument("--max-reps", type=int, default=80,
                    help="cap on representatives when rebuilding an embedding")
    sp.set_defaults(fn=cmd_warm)

    sub.add_parser("stats", help="cache statistics").set_defaults(fn=cmd_stats)
    sub.add_parser("purge", help="drop expired cache entries").set_defaults(fn=cmd_purge)

    args = p.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
