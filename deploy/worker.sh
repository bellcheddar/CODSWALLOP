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
#   * Embeddings and topology only. Both need a single chain, which the Model Server
#     returns in a few hundred kB; PLIP needs whole structures and is not run here yet.
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
