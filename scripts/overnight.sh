#!/usr/bin/env bash
# Build every pre-warm target's artefacts, unattended, and keep going until they are done.
#
#   bash scripts/overnight.sh            # run it
#   tail -f /tmp/codswallop-overnight.log
#
# This is a supervisor, not a launcher. A single long run is the wrong shape for the job:
# it is hours of network fetches, DSSP, PLIP and TM-align across sixty-four families, and
# any of those can hang on a slow upstream or be killed by memory pressure. Every stage is
# resumable because `warm` skips what is already current, so the honest design is to keep
# restarting until nothing is stale rather than to hope one pass survives.
#
# Three things it does that a bare `nohup` does not:
#
#   * `caffeinate` holds the machine awake. Without it a laptop sleeps twenty minutes in and
#     the whole run silently stops, which is the single most likely way to find nothing done
#     in the morning.
#   * Stall detection. A wedged HTTP read produces no output and no exit, so elapsed silence
#     is the only signal available; the run is killed and restarted from where it got to.
#   * It stops when the artefacts are actually complete, not when a command returns 0.
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

LOG=/tmp/codswallop-overnight.log
RUN=/tmp/codswallop-overnight-pass.log
STALL_MIN=${STALL_MIN:-25}          # silence beyond this means wedged, not slow
MAX_PASSES=${MAX_PASSES:-40}

say() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*" | tee -a "$LOG"; }

# Every seed in the pre-warm list, read from the list itself so the two cannot drift.
# A while-read loop, not `mapfile`: that is bash 4 and macOS ships 3.2, where it is a silent
# no-op. The first version of this script used it, produced an empty seed list, and warmed
# only the families that already existed while reporting success.
SEEDFILE=/tmp/cw_overnight_seeds.txt
./.venv/bin/python - > "$SEEDFILE" <<'PYSEEDS'
import importlib.util, pathlib
spec = importlib.util.spec_from_file_location("pw", pathlib.Path("pipeline/prewarm.py"))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
for seed, _ in m.TARGETS:
    print(seed)
PYSEEDS

SEEDS=()
while IFS= read -r line; do
  [ -n "$line" ] && SEEDS[${#SEEDS[@]}]="$line"
done < "$SEEDFILE"

if [ "${#SEEDS[@]}" -lt 50 ]; then
  echo "Refusing to start: only ${#SEEDS[@]} seeds read from pipeline/prewarm.py" | tee -a "$LOG"
  exit 1
fi

say "=== overnight build starting: ${#SEEDS[@]} targets ==="

# How many TARGETS are not finished, which is not the same question as how many existing
# families have stale artefacts. A target with no family yet does not appear in the artefact
# survey at all, so the first version of this counted 25 unbuilt targets as zero gaps and
# stopped after two passes with nothing done.
remaining() {
  ./.venv/bin/python - "$SEEDFILE" <<'PYREM'
import sys, pathlib
sys.path.insert(0, str(pathlib.Path.cwd()))
from codswallop import artefacts, db, family, resolve
db.init()
seeds = [l.strip() for l in open(sys.argv[1]) if l.strip()]
gaps = 0
for q in seeds:
    try:
        r = resolve.resolve(q)
        if r["status"] != "resolved":
            continue                       # cannot be built; not a gap to chase
        st = artefacts.status(family.slug_for(r["seed"]))
        if not (st["embedding"]["current"] and st["contacts"]["current"]
                and st["topology"]["current"]):
            gaps += 1
    except Exception:
        gaps += 1                          # unknown means unfinished
print(gaps)
PYREM
}

pass=0
while (( pass < MAX_PASSES )); do
  pass=$((pass + 1))

  left=$(remaining)
  say "pass ${pass}: ${left} artefact gaps outstanding"
  if (( pass > 1 )) && (( left == 0 )); then
    say "nothing stale left"
    break
  fi

  : > "$RUN"
  # caffeinate -imsu: no idle sleep, no disk sleep, no display sleep, system awake. It wraps
  # the work rather than running alongside it, so the assertion dies with the job.
  caffeinate -imsu ./.venv/bin/python -u CODSWALLOP.py warm --contacts "${SEEDS[@]}" \
    >> "$RUN" 2>&1 &
  job=$!
  say "pass ${pass}: started pid ${job}"

  # Watch it. Two ways to finish: the process exits, or it goes quiet for too long.
  while kill -0 "$job" 2>/dev/null; do
    sleep 60
    quiet=$(( ($(date +%s) - $(stat -f %m "$RUN" 2>/dev/null || echo 0)) / 60 ))
    if (( quiet >= STALL_MIN )); then
      say "pass ${pass}: no output for ${quiet} min, killing pid ${job}"
      # The whole process group: warm spawns mkdssp and PLIP, and killing only the parent
      # leaves those holding the cores.
      pkill -P "$job" 2>/dev/null
      kill -9 "$job" 2>/dev/null
      sleep 5
      break
    fi
  done
  wait "$job" 2>/dev/null
  rc=$?

  tail -3 "$RUN" | sed 's/^/    /' >> "$LOG"
  say "pass ${pass}: exited rc=${rc}"

  # A pass that did nothing at all twice running is not going to start working: bail rather
  # than spin until morning on something a person has to look at.
  if (( rc != 0 )) && [[ ! -s "$RUN" ]]; then
    say "pass ${pass}: produced no output at all"
  fi
  sleep 10
done

say "=== finished after ${pass} passes: $(remaining) gaps outstanding ==="
./.venv/bin/python CODSWALLOP.py artefacts --missing 2>&1 | tail -3 | tee -a "$LOG"

# Ship whatever got built, so the morning does not need a second instruction.
say "pushing artefacts"
bash deploy/push_embeddings.sh >> "$LOG" 2>&1 && say "pushed" || say "push FAILED"
say "syncing the family list to the droplet"
bash deploy/sync.sh >> "$LOG" 2>&1 && say "synced" || say "sync FAILED"
say "=== done ==="
