# NEXT: where CODSWALLOP is, and what to do next

Handover note, 2026-08-10. Everything described here is committed, pushed and deployed.
Live at [codswallop.mdeller.com](https://codswallop.mdeller.com).

---

## The decision that changed the architecture

**All compute now runs on the droplet. The Mac is not used for anything.**

It could not before, for one measured reason: parsing a whole deposited assembly to extract
one chain's alpha carbons peaked at **6.3 GB RSS** (8GLV is a 453 MB file that yields 426
alpha carbons). The droplet has 2 cores, ~2.2 GB free and **no swap**, so that was not
"slow", it was an OOM that takes the other eight apps down with it.

Fetching the single chain from the RCSB **Model Server** instead:

| | whole file | one chain |
|---|---|---|
| 8GLV transfer | 453 MB | **533 kB** |
| Peak RSS | 6.3 GB | **102 MB** |
| Per structure | ~0.8 s | **~0.4-0.6 s** |
| 232-chain family | n/a | ~90 s of fetching |

gemmi was tried and is 12x faster on the same file while still needing ~6 GB. The parser was
never the problem; the file really is that big.

`embed.structure_path(pdb_id, chain)` is the **one** place anything asks for a structure
file. It prefers the chain, falls back to the whole entry, and refuses anything over
`MAX_WHOLE_FILE_MB` (12) rather than parsing it.

### Proven on the droplet

VEGFR-1 (P17948, 1,328 entries) was queued by a reader, picked up by the timer, and embedded
in **6,142 s (~1h42m)** on one `nice`d core. Peak RSS 102 MB; free memory never dipped below
2.2 GB even with two families building at once. Structure cache after the run: **20 kB**.

---

## What is left

**Nothing from the migration.** All three steps are done and verified on the droplet:

1. **Embeddings** through the Model Server, one chain at a time.
2. **Topology** through the same route, `CODSWALLOP.py topology <query>`, biotite P-SEA
   where `mkdssp` is absent (it ships with CCP4 and is not pip-installable) with the artefact
   recording which method drew it.
3. **Contacts** with a 30 MB structure cap, skips counted apart from failures and reported on
   the panel.

`deploy/sync.sh`, `drain_queue.sh`, `push_embeddings.sh`, `scripts/overnight.sh` and the
launchd timer are deleted, and the launchd job is unloaded. The Mac builds nothing.

Proven end to end on VEGFR-1 (P17948, 1,328 entries), queued by a reader and built entirely
on the droplet: embedding 1h42m at 102 MB peak, topology by P-SEA in 1 s, PLIP 6 entries in
21 s. Free memory never dipped below 2.2 GB.

### Ideas, not commitments

- **`mkdssp` on the droplet** would upgrade every topology from P-SEA to real DSSP, which is
  the better assignment and the only one that gives bridge partners. It is a CCP4 component,
  so it needs a system package rather than a wheel.
- **Contacts is the slow one**: 1-3 min per entry over up to 60 entries, so hours per family
  against minutes for the rest. It runs last, at nice 19, and `CONTACTS=0` disables it.
- The dossier still has no **lego plots**, and the structure still is the RCSB's render
  rather than one drawn from the family's own superposition.

## How the droplet worker works

- `deploy/worker.sh`, driven by `codswallop-worker.timer` **every 15 minutes**.
- Drains `artefact_request`, the queue the web app already writes when a reader opens a
  family with no current artefact.
- **One family at a time**, `nice 15`, `ionice -c3`, `MemoryMax=900M`, `flock` so a timer
  firing mid-job does nothing. Two cores shared with eight apps, and TM-align is
  single-threaded, so a second job buys nothing and costs the web its core.
- A family is marked served **only after its artefact exists**. Marking on command
  completion would drop a failed family from the queue and leave it on the placeholder for
  ever with nothing recording that it still needs one.
- Chain files are cleared per family: they are a cache on an 18 GB disk shared with
  everything else, and re-fetching one costs under half a second.
- Log: `/var/log/codswallop-worker.log`.
- Compute deps are in `requirements-compute.txt` (tmtools, biotite, gemmi). All ship prebuilt
  x86_64 wheels, so the droplet needs no compiler. **The web process still imports none of
  them**, so a web worker cannot be killed by parsing a structure.

---

## Traps worth not rediscovering

**A side effect another module relies on is not an interface.** Moving the fetch to the Model
Server silently broke topology: it read `embed._cif_path()` and depended on `ca_trace` having
downloaded the whole entry as a side effect. The failure wore the costume of a normal answer,
because "no structure file" was already a legitimate "no 2D layout for this one".

**The same offset applied twice.** `contacts.py` had `seed_pos = resnr + (query_beg - 1)`,
adding the seed offset to PLIP's *author* residue number which already carries it. 52 of 71
families were wrong. It was invisible on carbonic anhydrase, the family it was validated on,
because that family's `query_beg` is 1 so the offset is zero. Now mapped by **aligning**, as
`topology.map_to_seed` already did. **If you touch a seed-coordinate conversion, check the
mapped residue identity against the seed sequence**: that check found it and would have
prevented it.

**Verify in a real browser, with real input events.** Three shipped features had never worked
and looked fine: `focusResidue` (built its selection with `lib.structure.Script`, absent from
this Mol\* build, so it threw and returned a bare `false` for its whole existence); a cluster
labeller whose preferred input was empty on every real request; and the map's node clicks.
Synthetic clicks dispatched at an element bypass the browser's own hit-testing and cannot
fail. The CDP harness also needs `--use-angle=swiftshader --enable-unsafe-swiftshader` or
every 3D screenshot is a blank panel reading "WebGL does not seem to be available".

**Cache the parse version.** `uniprot.py` has now needed `PARSE_VERSION` in a cache key three
separate times. Records are cached in *parsed* form, so adding a field leaves it absent from
every already-cached row and the feature reads as unavailable everywhere rather than as new.

**Per-family logging made a healthy run look dead.** PLIP takes 1-3 minutes per entry, so a
40-entry family printed nothing for 35 minutes and I called it a stall. `warm` now emits a
heartbeat per entry, throttled to one line a minute.

---

## Open, not urgent

- The 2D map holds only **33-100%** of the variance (median 65%) and 1&nbsp;&minus;&nbsp;TM is
  not Euclidean (negative eigenvalue mass a median 19%, up to 49%). Documented in the README
  under "What the map is not", not on the panel.
- `family.py` fetches UniProt features for PDB-seeded families by aligning the canonical onto
  the seed. Works (lysozyme gives Glu35/Asp52), but the 0.6 agreement floor is shared with the
  conservation colouring and has not been tuned separately.
