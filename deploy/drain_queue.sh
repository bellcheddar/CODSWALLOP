#!/usr/bin/env bash
# Build, on this workstation, the artefacts that readers of the live site have asked for.
#
#   bash deploy/drain_queue.sh              # build every open request
#   bash deploy/drain_queue.sh --max 3      # build the three most-wanted and stop
#   bash deploy/drain_queue.sh --dry-run    # show the queue without building anything
#
# The gap this closes: the embedding and the interaction fingerprint need biotite, tmtools,
# PLIP and OpenBabel, none of which are installed on the droplet and none of which belong on
# a box with two cores shared between eight apps. So a family that a reader assembles for
# the first time on the live site has an artefact on neither machine, and its map falls back
# to the sequence-identity placeholder. The droplet is the only machine that knows this
# happened and the one machine that cannot do anything about it, so it writes the request
# down and this script comes and collects.
#
# Idempotent, and safe to interrupt: a family is marked served on the droplet only after its
# artefact has actually been pushed, so a run that dies half way leaves the rest queued.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

if [[ -f .env ]]; then set -a; source .env; set +a; fi
DROPLET_SSH="${DROPLET_SSH:-}"
DROPLET_PATH="${DROPLET_PATH:-/opt/codswallop}"
SSH_KEY="${SSH_KEY:-}"
MAX=0
DRY=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --max) MAX="$2"; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    *) echo "unknown option: $1"; exit 2 ;;
  esac
done

[[ -z "$DROPLET_SSH" ]] && { echo "DROPLET_SSH is not set (see .env.example)."; exit 1; }
[[ -x .venv/bin/python ]] || { echo "No venv here. Run: python -m venv .venv"; exit 1; }

SSH_CMD=(ssh); [[ -n "$SSH_KEY" ]] && SSH_CMD=(ssh -i "${SSH_KEY/#\~/$HOME}")

echo "==> Reading the queue from ${DROPLET_SSH}"
# One JSON object per line rather than one array, so a truncated transfer costs a family
# rather than the whole queue.
queue=$("${SSH_CMD[@]}" "$DROPLET_SSH" \
  "cd ${DROPLET_PATH} && ./.venv/bin/python CODSWALLOP.py queue --json --limit 200" || true)

if [[ -z "${queue//[[:space:]]/}" ]]; then
  echo "    The queue is empty: every family a reader has opened has its artefacts."
  exit 0
fi

n=$(printf '%s\n' "$queue" | grep -c . || true)
echo "    ${n} famil$([[ $n -eq 1 ]] && echo y || echo ies) waiting, most-wanted first"
printf '%s\n' "$queue" | .venv/bin/python -c '
import json, sys
for line in sys.stdin:
    if line.strip():
        r = json.loads(line)
        print("      %-42s %-10s %4s hits  %s entries"
              % (r["slug"], r["kind"], r["hits"], r.get("n_entries") or "?"))
'

if [[ "$DRY" -eq 1 ]]; then
  echo "==> --dry-run: stopping before building anything"
  exit 0
fi

built=0
while IFS= read -r line; do
  [[ -z "${line//[[:space:]]/}" ]] && continue
  slug=$(printf '%s' "$line" | .venv/bin/python -c 'import json,sys; print(json.load(sys.stdin)["slug"])')
  query=$(printf '%s' "$line" | .venv/bin/python -c 'import json,sys; print(json.load(sys.stdin)["query"])')
  kind=$(printf '%s' "$line" | .venv/bin/python -c 'import json,sys; print(json.load(sys.stdin)["kind"])')

  if [[ "$MAX" -gt 0 && "$built" -ge "$MAX" ]]; then
    echo "==> Reached --max ${MAX}; the rest stay queued"
    break
  fi

  echo
  echo "==> ${slug}  (${query}, needs ${kind})"
  # Not `set -e`-fatal: one family that cannot be built must not strand the others. It stays
  # in the queue, which is the correct record of the situation.
  if ! ./.venv/bin/python -u CODSWALLOP.py embed "$query"; then
    echo "    embedding failed; leaving it queued"
    continue
  fi
  if [[ "$kind" == "both" || "$kind" == "contacts" ]]; then
    ./.venv/bin/python -u CODSWALLOP.py contacts "$query" || \
      echo "    contacts failed; the embedding still ships"
  fi
  built=$((built + 1))
done <<< "$queue"

if [[ "$built" -eq 0 ]]; then
  echo; echo "==> Nothing built, so nothing to push."
  exit 0
fi

echo
echo "==> Pushing ${built} newly built famil$([[ $built -eq 1 ]] && echo y || echo ies)"
bash deploy/push_embeddings.sh

# Marked served only now, after the artefacts are actually on the droplet. Marking them when
# the build finished would have quietly emptied the queue on a failed push, which is the one
# state where the record matters most.
echo
echo "==> Clearing the served entries"
while IFS= read -r line; do
  [[ -z "${line//[[:space:]]/}" ]] && continue
  slug=$(printf '%s' "$line" | .venv/bin/python -c 'import json,sys; print(json.load(sys.stdin)["slug"])')
  if [[ -f "data/embeddings/${slug}.json" ]]; then
    "${SSH_CMD[@]}" "$DROPLET_SSH" \
      "cd ${DROPLET_PATH} && ./.venv/bin/python CODSWALLOP.py queue --served '${slug}'" \
      >/dev/null && echo "    cleared ${slug}"
  fi
done <<< "$queue"

echo
echo "==> Done. Remaining queue:"
"${SSH_CMD[@]}" "$DROPLET_SSH" "cd ${DROPLET_PATH} && ./.venv/bin/python CODSWALLOP.py queue" || true
