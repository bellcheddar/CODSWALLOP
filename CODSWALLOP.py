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
import json
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
                print(f"  skip  {q!r}: {result.get('message') or 'ambiguous'}", flush=True)
                continue
            fam = family.get_or_build(result["seed"], q, force=args.force)
            # The drug enrichment is fetched here and only here. A page request reads it
            # from the cache and never fetches: these are one third-party call per drug, and
            # ABL1's 65 drug-like components put a page load past two minutes on the droplet.
            try:
                from codswallop import drugs as drug_engine
                drug_engine.build(fam, fetch_missing=True)
            except Exception as exc:            # noqa: BLE001
                print(f"        (drug enrichment skipped: {exc})", flush=True)
            print(f"  warm  {fam['slug']:<44} {fam['stats']['entries']:>5} entries  "
                  f"{time.time() - t0:>5.1f}s", flush=True)
            ok += 1
        except Exception as exc:            # one bad family must not stop the rest
            print(f"  FAIL  {q!r}: {exc}", flush=True)
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

    every_ok = _rebuild_embeddings(artefacts, args)
    # After the embeddings, deliberately: topology is drawn on the reference structure the
    # embedding chose, so a family embedded in this same pass gets its diagram on the
    # structure it was actually superposed onto.
    every_ok = _rebuild_topology(artefacts) and every_ok
    if args.contacts:
        every_ok = _rebuild_contacts(artefacts, args) and every_ok
    return every_ok


def _rebuild_topology(artefacts) -> bool:
    """Assign secondary structure for any family missing it.

    Cheap next to the embedding: one structure, one DSSP run, a few seconds. It goes after
    the embedding because it uses the same reference structure, so a family embedded in this
    same pass gets its topology drawn on the structure it was superposed onto.
    """
    from codswallop import family, resolve, topology

    todo = artefacts.stale("topology")
    if not todo:
        return True
    print(f"\n{len(todo)} topology artefact(s) missing or out of date "
          f"(v{artefacts.topology_io.VERSION}); rebuilding:")
    every_ok = True
    for s in todo:
        if not s["query"]:
            print(f"  skip  {s['slug']}: no stored query to rebuild from", flush=True)
            continue
        t0 = time.time()
        try:
            r = resolve.resolve(s["query"])
            if r["status"] != "resolved":
                print(f"  skip  {s['slug']}: {r.get('message') or 'ambiguous'}", flush=True)
                continue
            fam = family.get_or_build(r["seed"], s["query"])
            art = topology.write(fam)
            if not art:
                print(f"  skip  {s['slug']}: no secondary structure could be assigned",
                      flush=True)
                continue
            print(f"  dssp  {s['slug']:<44} {art['n_strands']:>3} strands "
                  f"{art['n_helices']:>3} helices  {art['method']:<6} "
                  f"{time.time() - t0:>5.0f}s", flush=True)
        except Exception as exc:                    # noqa: BLE001
            every_ok = False
            print(f"  FAIL  {s['slug']}: {exc}", flush=True)
    return every_ok


def _rebuild_embeddings(artefacts, args) -> bool:
    from codswallop import embed, family, resolve

    todo = artefacts.stale("embedding")
    if not todo:
        return True
    print(f"\n{len(todo)} embedding(s) missing or out of date "
          f"(pipeline v{artefacts.embed_io.VERSION}); rebuilding:")
    every_ok = True
    for s in todo:
        if not s["query"]:
            print(f"  skip  {s['slug']}: no stored query to rebuild from", flush=True)
            continue
        t0 = time.time()
        try:
            r = resolve.resolve(s["query"])
            if r["status"] != "resolved":
                print(f"  skip  {s['slug']}: {r.get('message') or 'ambiguous'}", flush=True)
                continue
            fam = family.get_or_build(r["seed"], s["query"])
            # A heartbeat, on its own line every half minute. Without it a single large
            # family is completely silent for as long as it takes: spike's embedding runs
            # for 32 minutes, and an unattended supervisor watching for silence cannot tell
            # that from a wedged socket. It killed and restarted the pass every 25 minutes
            # and would have spent the night doing so.
            # Seeded in the past, so the first notification is not swallowed by the very
            # guard meant to throttle the later ones: that left the align stage silent from
            # start to finish, which is exactly the stretch this exists to cover.
            beat = [0.0]

            def heartbeat(stage, i, n, label):
                if time.time() - beat[0] < 30:
                    return
                beat[0] = time.time()
                if stage == "fetch":
                    print(f"        {s['slug']}: fetching {i}/{n} {label}", flush=True)
                else:
                    pct = (100 * i / n) if n else 0
                    print(f"        {s['slug']}: aligning {i:,}/{n:,} pairs ({pct:.0f}%)",
                          flush=True)

            art = embed.build(fam, max_representatives=args.max_reps, progress=heartbeat)
            if not art:
                print(f"  skip  {s['slug']}: not enough usable structures", flush=True)
                continue
            af = art.get("alphafold")
            print(f"  embed {s['slug']:<44} {art['n_representatives']:>3} reps"
                  + (f", AF TM {af['tm']}" if af else ", no AF model")
                  + f"  {time.time() - t0:>5.0f}s", flush=True)
        except Exception as exc:
            print(f"  FAIL  {s['slug']}: {exc}", flush=True)
            every_ok = False
    return every_ok


def _rebuild_contacts(artefacts, args) -> bool:
    """PLIP for any family whose fingerprint is missing or stale.

    Opt-in (`--contacts`) rather than automatic, unlike the embeddings. PLIP is minutes per
    family where an embedding is seconds-to-minutes, and a weekly run that quietly grew to
    an hour is a weekly run somebody turns off. The report still names what is missing
    either way.
    """
    from codswallop import contacts as contact_engine, family, resolve

    todo = artefacts.stale("contacts")
    if not todo:
        return True
    print(f"\n{len(todo)} contact profile(s) missing or out of date "
          f"(pipeline v{artefacts.contacts_io.VERSION}); rebuilding:")
    every_ok = True
    for s in todo:
        if not s["query"]:
            print(f"  skip  {s['slug']}: no stored query to rebuild from", flush=True)
            continue
        t0 = time.time()
        try:
            r = resolve.resolve(s["query"])
            if r["status"] != "resolved":
                print(f"  skip  {s['slug']}: {r.get('message') or 'ambiguous'}", flush=True)
                continue
            fam = family.get_or_build(r["seed"], s["query"])
            art = contact_engine.build(fam, max_entries=args.max_contacts)
            if not art:
                print(f"  skip  {s['slug']}: nothing ligand-bound to profile", flush=True)
                continue
            print(f"  plip  {s['slug']:<44} {art['n_contacts']:>6,} contacts from "
                  f"{art['entries_analysed']:>3} entries  {time.time() - t0:>5.0f}s",
                  flush=True)
        except Exception as exc:
            print(f"  FAIL  {s['slug']}: {exc}", flush=True)
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


def cmd_queue(args) -> int:
    """List, or resolve, the families a reader asked for that no workstation has built.

    `--json` is what `deploy/drain_queue.sh` reads over SSH, so its shape is an interface:
    one object per line rather than one array, so the reading end can stream it and a
    truncated transfer loses a family instead of the whole queue.
    """
    db.init()
    if args.served:
        db.mark_request_served(args.served)
        print(f"  marked {args.served} as served")
        return 0

    rows = db.open_requests(limit=args.limit)
    if args.json:
        for r in rows:
            print(json.dumps({k: r[k] for k in ("slug", "query", "kind", "hits", "n_entries")}))
        return 0
    if not rows:
        print("  the queue is empty: every family a reader has opened has its artefacts")
        return 0
    print(f"  {len(rows)} famil{'y' if len(rows) == 1 else 'ies'} waiting on a workstation:\n")
    for r in rows:
        print("  %-42s %-12s %-9s %4d hit%s %6s entries" % (
            r["slug"], (r["query"] or "-")[:12], r["kind"], r["hits"],
            " " if r["hits"] == 1 else "s", r["n_entries"] or "?"))
    print("\n  Drain them with:  bash deploy/drain_queue.sh")
    return 0


def cmd_artefacts(args) -> int:
    """Report artefact state. Runs anywhere, including the droplet, which cannot build them.

    This is what closes the loop on families a READER created: they exist only in the
    droplet's database, so a workstation never knows to build artefacts for them and they
    sit on the placeholder map indefinitely. `deploy/push_embeddings.sh` calls this over SSH
    after every push so the gap is named rather than discovered months later.
    """
    from codswallop import artefacts

    db.init()
    rows = artefacts.survey()
    for s in rows:
        e, c = s["embedding"], s["contacts"]
        if args.missing and e["current"] and c["current"]:
            continue
        print("  %-42s %-12s emb %-5s %-8s contacts %-5s %s" % (
            s["slug"], (s["query"] or "-")[:12],
            ("v%s" % e["version"]) if e["version"] else "none",
            "ok" if e["current"] else "MISSING",
            ("v%s" % c["version"]) if c["version"] else "none",
            "ok" if c["current"] else "missing"))
    if args.missing and not any(
            not (s["embedding"]["current"] and s["contacts"]["current"]) for s in rows):
        print("  none: every family has current artefacts")
    summary = artefacts.summary()
    print(f"\n  {summary['embeddings_current']}/{summary['families']} embeddings current, "
          f"{summary['contacts_current']}/{summary['families']} contacts. "
          f"{summary['on_placeholder']} on the placeholder map.")
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
    sp.add_argument("--contacts", action="store_true",
                    help="also run PLIP for families whose fingerprint is missing "
                         "(minutes per family; off by default)")
    sp.add_argument("--max-contacts", type=int, default=40,
                    help="cap on entries per PLIP run")
    sp.set_defaults(fn=cmd_warm)

    sp = sub.add_parser("artefacts", help="which families have a current embedding/contacts")
    sp.add_argument("--missing", action="store_true", help="list only the stale ones")
    sp.set_defaults(fn=cmd_artefacts)

    sp = sub.add_parser("queue", help="families a reader opened that need a workstation")
    sp.add_argument("--json", action="store_true", help="one JSON object per line, for scripts")
    sp.add_argument("--limit", type=int, default=50, help="how many to list")
    sp.add_argument("--served", metavar="SLUG", help="mark one as built and pushed")
    sp.set_defaults(fn=cmd_queue)

    sub.add_parser("stats", help="cache statistics").set_defaults(fn=cmd_stats)
    sub.add_parser("purge", help="drop expired cache entries").set_defaults(fn=cmd_purge)

    args = p.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
