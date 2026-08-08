#!/usr/bin/env bash
# Push locally computed structural embeddings to the droplet.
#
#   python CODSWALLOP.py embed P00918      # on your Mac: downloads structures, runs TM-align
#   bash deploy/push_embeddings.sh         # ship the artefacts
#
# Separate from deploy.sh on purpose. deploy.sh excludes data/ so a code deploy can never
# clobber the server's cache, and these artefacts live in data/embeddings/. They are also
# the only thing in data/ that a workstation produces and the droplet cannot: everything
# else there is a cache the server rebuilds for itself.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

if [[ -f .env ]]; then set -a; source .env; set +a; fi
DROPLET_SSH="${DROPLET_SSH:-}"
DROPLET_PATH="${DROPLET_PATH:-/opt/codswallop}"
SSH_KEY="${SSH_KEY:-}"

[[ -z "$DROPLET_SSH" ]] && { echo "DROPLET_SSH is not set (see .env.example)."; exit 1; }
[[ -d data/embeddings ]] || { echo "No embeddings yet. Run: python CODSWALLOP.py embed <query>"; exit 1; }

n=$(find data/embeddings -name '*.json' | wc -l | tr -d ' ')
echo "==> Pushing ${n} embedding(s) to ${DROPLET_SSH}:${DROPLET_PATH}/data/embeddings/"

SSH_OPTS=()
[[ -n "$SSH_KEY" ]] && SSH_OPTS=(-e "ssh -i ${SSH_KEY/#\~/$HOME}")

# No --delete: an embedding on the server that this workstation has not computed is not
# stale, it was computed somewhere else. Removing artefacts is a deliberate act, not a
# side effect of pushing from a different machine.
rsync -az ${SSH_OPTS[@]+"${SSH_OPTS[@]}"} \
  --include '*/' --include '*.json' --exclude '*' \
  data/embeddings/ "${DROPLET_SSH}:${DROPLET_PATH}/data/embeddings/"

SSH_CMD=(ssh); [[ -n "$SSH_KEY" ]] && SSH_CMD=(ssh -i "${SSH_KEY/#\~/$HOME}")
"${SSH_CMD[@]}" "$DROPLET_SSH" \
  "chown -R codswallop:codswallop ${DROPLET_PATH}/data/embeddings && \
   ls -1 ${DROPLET_PATH}/data/embeddings | wc -l | xargs echo '   embeddings on the droplet:'"
echo "==> Done. The map switches to the embedding on the next page load."
