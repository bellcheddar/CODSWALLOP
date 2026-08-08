#!/usr/bin/env bash
# Bring the workstation and the droplet into step, in both directions.
#
#   bash deploy/sync.sh              # the whole cycle
#   bash deploy/sync.sh --dry-run    # say what it would do
#
# Why this exists. The two machines hold different halves of the same app and neither can do
# the other's job:
#
#   * The DROPLET is the only one that knows what readers have asked for. Families they
#     assemble live there and nowhere else, because the cache database is deliberately
#     excluded from every rsync: it is the server's own state and a code deploy must never
#     clobber it.
#   * The WORKSTATION is the only one that can build an embedding, an interaction fingerprint
#     or a topology. Those need biotite, tmtools, PLIP, OpenBabel and DSSP, none of which
#     belong on a two-core box shared with eight other apps.
#
# So the family *set* has to travel one way and the artefacts the other, and until this
# script existed neither happened on its own: the droplet had 13 families, the workstation
# had 34, and the landing page showed a third of the archive.
#
# Idempotent and safe to run at any time. Everything it does is a rebuild of something
# already derivable, so an interrupted run costs time and nothing else.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

if [[ -f .env ]]; then set -a; source .env; set +a; fi
DROPLET_SSH="${DROPLET_SSH:-}"
DROPLET_PATH="${DROPLET_PATH:-/opt/codswallop}"
SSH_KEY="${SSH_KEY:-}"
DRY=0
[[ "${1:-}" == "--dry-run" ]] && DRY=1

[[ -z "$DROPLET_SSH" ]] && { echo "DROPLET_SSH is not set (see .env.example)."; exit 1; }
[[ -x .venv/bin/python ]] || { echo "No venv here."; exit 1; }

SSH_CMD=(ssh); [[ -n "$SSH_KEY" ]] && SSH_CMD=(ssh -i "${SSH_KEY/#\~/$HOME}")
SCP_CMD=(scp); [[ -n "$SSH_KEY" ]] && SCP_CMD=(scp -i "${SSH_KEY/#\~/$HOME}")

say() { printf '==> %s\n' "$*"; }

# ---------------------------------------------------------------------------------------
# 1. What each machine has filed, by the query it was built from.
#
# Keyed on the query rather than the slug: the slug is derived from the seed and can change
# if a protein is renamed upstream, while the query is what a person actually typed and is
# the only thing either machine can rebuild from.
# ---------------------------------------------------------------------------------------
say "Reading both family lists"
./.venv/bin/python - <<'PY' > /tmp/cw_local_queries.txt
from codswallop import db
db.init()
for r in db.connect().execute("SELECT query FROM family ORDER BY built_at"):
    if r["query"]:
        print(r["query"])
PY

"${SSH_CMD[@]}" "$DROPLET_SSH" \
  "cd ${DROPLET_PATH} && ./.venv/bin/python -c \"
from codswallop import db
db.init()
for r in db.connect().execute('SELECT query FROM family ORDER BY built_at'):
    print(r['query']) if r['query'] else None
\"" > /tmp/cw_droplet_queries.txt || true

local_n=$(grep -c . /tmp/cw_local_queries.txt || true)
drop_n=$(grep -c . /tmp/cw_droplet_queries.txt || true)
sort -u /tmp/cw_local_queries.txt /tmp/cw_droplet_queries.txt > /tmp/cw_all_queries.txt
all_n=$(grep -c . /tmp/cw_all_queries.txt || true)
say "workstation ${local_n}, droplet ${drop_n}, union ${all_n}"

# Only the ones each side is missing, so a run with nothing to do says so and stops.
comm -13 <(sort -u /tmp/cw_local_queries.txt) <(sort -u /tmp/cw_droplet_queries.txt) \
  > /tmp/cw_missing_here.txt || true
comm -23 <(sort -u /tmp/cw_local_queries.txt) <(sort -u /tmp/cw_droplet_queries.txt) \
  > /tmp/cw_missing_there.txt || true
here_n=$(grep -c . /tmp/cw_missing_here.txt || true)
there_n=$(grep -c . /tmp/cw_missing_there.txt || true)
say "${here_n} families to pull in here, ${there_n} to seed over there"

if [[ "$DRY" -eq 1 ]]; then
  say "--dry-run: stopping"
  [[ "$here_n" -gt 0 ]] && { echo "  would build here:"; sed 's/^/    /' /tmp/cw_missing_here.txt; }
  [[ "$there_n" -gt 0 ]] && { echo "  would seed there:"; sed 's/^/    /' /tmp/cw_missing_there.txt; }
  exit 0
fi

# ---------------------------------------------------------------------------------------
# 2. Warm everything here, which assembles anything new and builds every missing artefact.
#
# The whole union, not just the difference: `warm` skips what is already current, and a
# family that went stale for its own reasons should be caught by the same pass.
# ---------------------------------------------------------------------------------------
say "Warming ${all_n} families on this workstation (artefacts included)"
# shellcheck disable=SC2046
./.venv/bin/python -u CODSWALLOP.py warm $(tr '\n' ' ' < /tmp/cw_all_queries.txt) || \
  say "warm reported a problem; continuing so the artefacts that did build still ship"

# ---------------------------------------------------------------------------------------
# 3. Artefacts to the droplet, and the family list back the other way.
# ---------------------------------------------------------------------------------------
say "Pushing artefacts"
bash deploy/push_embeddings.sh

say "Seeding the droplet with the full family list"
"${SCP_CMD[@]}" /tmp/cw_all_queries.txt "${DROPLET_SSH}:/tmp/cw_all_queries.txt" >/dev/null
# --no-artefacts: the droplet cannot build one and must not spend an hour discovering that.
"${SSH_CMD[@]}" "$DROPLET_SSH" \
  "cd ${DROPLET_PATH} && ./.venv/bin/python -u CODSWALLOP.py warm --no-artefacts \
     \$(tr '\n' ' ' < /tmp/cw_all_queries.txt)" || \
  say "the droplet's warm reported a problem"

# ---------------------------------------------------------------------------------------
# 4. Say where both machines ended up, because a sync that silently half-worked is the one
#    failure mode worth catching.
# ---------------------------------------------------------------------------------------
echo
say "Workstation"
./.venv/bin/python CODSWALLOP.py artefacts --missing 2>&1 | tail -2
say "Droplet"
"${SSH_CMD[@]}" "$DROPLET_SSH" \
  "cd ${DROPLET_PATH} && ./.venv/bin/python CODSWALLOP.py artefacts --missing 2>&1 | tail -2; \
   echo; ./.venv/bin/python CODSWALLOP.py queue 2>&1 | tail -3"
