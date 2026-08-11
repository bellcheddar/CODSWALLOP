#!/usr/bin/env bash
# The droplet's own artefact worker.
#
# This replaces `drain_queue.sh`, which ran on Marc's Mac, read the queue over SSH and pushed
# the results back. That made "this family has been queued for it" mean "when a laptop next
# wakes up", which for a laptop that is off for a week means never. The droplet is now the
# only machine that computes anything, so the queue drains where it is written.
#
# Deliberately modest about what it takes on:
#
#   * ONE family at a time. Two cores, shared with eight other apps, and TM-align is
#     single-threaded: a second job buys nothing and costs the web workers their core.
#   * `nice`d and `ionice`d, because a reader waiting on a page matters more than an
#     artefact that nobody is watching.
#   * Embedding, then topology, then contacts. The first two need a single chain, which the
#     Model Server returns in a few hundred kB. Contacts needs whole structures, so it caps
#     the file size and counts what it skipped; it runs last and at nice 19 because it is
#     hours where the others are minutes, and CONTACTS=0 disables it.
#   * Chain files are deleted after each family. They are a cache, the disk is 18 GB shared
#     with everything else on the box, and re-fetching one is under half a second.
#
# Runs from a systemd timer. Safe to run twice: `warm` and `embed` both skip what is already
# current, and the queue is marked served only after the artefact exists.
set -uo pipefail

APP=/opt/codswallop
cd "$APP" || exit 1
PY="$APP/.venv/bin/python"
LOG=/var/log/codswallop-worker.log

say() { printf '[%s] %s\n' "$(date '+%F %H:%M:%S')" "$*" >> "$LOG"; }

# One worker, always. A timer that fires while the last run is still going would otherwise
# put two TM-aligns on two cores and leave nothing for the web.
exec 9>/run/codswallop-worker.lock
flock -n 9 || { say "another worker holds the lock; leaving it to finish"; exit 0; }

# How many families to take in one pass. The queue is a priority order (most-requested
# first), so stopping early costs the least-wanted family a wait rather than losing it.
MAX_FAMILIES=${MAX_FAMILIES:-3}

mapfile_compat() { while IFS= read -r line; do printf '%s\n' "$line"; done; }

queue=$("$PY" CODSWALLOP.py queue --json --limit "$MAX_FAMILIES" 2>/dev/null)
if [[ -z "$queue" ]]; then
  say "queue empty"
  exit 0
fi

count=0
while IFS= read -r row; do
  [[ -z "$row" ]] && continue
  slug=$(printf '%s' "$row" | "$PY" -c 'import json,sys; print(json.load(sys.stdin)["slug"])')
  query=$(printf '%s' "$row" | "$PY" -c 'import json,sys; print(json.load(sys.stdin)["query"])')
  [[ -z "$slug" || -z "$query" ]] && continue

  count=$((count + 1))
  say "start ${slug} (${query})"
  t0=$(date +%s)

  # nice 15 and idle I/O: this must never be the reason a page is slow.
  if nice -n 15 ionice -c3 "$PY" -u CODSWALLOP.py embed "$query" >> "$LOG" 2>&1; then
    say "embedded ${slug} in $(( $(date +%s) - t0 ))s"

    # Topology after the embedding, deliberately: the diagram is drawn on the reference
    # structure the EMBEDDING chose, so running it first would draw on a structure the
    # family is not superposed onto. Seconds rather than minutes, and it needs one chain.
    #
    # Its failure does NOT hold the family in the queue. The PDBe publishes a 2D layout for
    # most of the archive and not all of it, so "no topology" is a legitimate answer for
    # this structure rather than work still owed; leaving the family queued for it would
    # mean retrying a thing that will never succeed, every fifteen minutes, for ever.
    t1=$(date +%s)
    if nice -n 15 ionice -c3 "$PY" -u CODSWALLOP.py topology "$query" >> "$LOG" 2>&1; then
      say "topology for ${slug} in $(( $(date +%s) - t1 ))s"
    else
      say "topology failed for ${slug}; the embedding stands"
    fi

    # Contacts last, and optional. PLIP is one to three minutes per entry over up to sixty
    # entries, so this is hours where the embedding is minutes: it must not stand between a
    # reader and their map. CONTACTS=0 turns it off entirely if the box gets busy.
    #
    # It needs the WHOLE structure, unlike the other two, so entries above the size cap are
    # skipped and counted rather than parsed. Its failure does not hold the family either;
    # the artefact it would write is an addition to a page that already works.
    if [[ "${CONTACTS:-1}" == "1" ]]; then
      t2=$(date +%s)
      if nice -n 19 ionice -c3 "$PY" -u CODSWALLOP.py contacts "$query" >> "$LOG" 2>&1; then
        say "contacts for ${slug} in $(( $(date +%s) - t2 ))s"
      else
        say "contacts failed for ${slug}; the map and diagram stand"
      fi
    fi

    # Only now: marking it served on a failed build would drop it from the queue and the
    # family would sit on the placeholder for ever with nothing recording that it needs one.
    "$PY" CODSWALLOP.py queue --served "$slug" >> "$LOG" 2>&1
  else
    say "FAILED ${slug}; left in the queue for the next pass"
  fi

  # The chain files are a cache, not an artefact. Clearing them per family keeps a runaway
  # queue from filling a disk that eight other apps share.
  rm -f "$APP"/data/structures/*.cif 2>/dev/null
done <<< "$queue"

say "pass done: ${count} famil$([[ $count -eq 1 ]] && echo y || echo ies)"
