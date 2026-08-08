#!/usr/bin/env bash
# Push CODSWALLOP from your Mac to the droplet and restart the web service.
# Run from the repo root:  bash deploy/deploy.sh
#
# Reads DROPLET_SSH / DROPLET_PATH from .env (see .env.example). Idempotent; excludes the
# venv, the cache database and secrets so the server's state is never clobbered.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

if [[ -f .env ]]; then set -a; source .env; set +a; fi
DROPLET_SSH="${DROPLET_SSH:-}"
DROPLET_PATH="${DROPLET_PATH:-/opt/codswallop}"
SSH_KEY="${SSH_KEY:-}"

if [[ -z "$DROPLET_SSH" ]]; then
  echo "DROPLET_SSH is not set. Copy .env.example to .env and fill it in."; exit 1
fi

SSH_OPTS=()
[[ -n "$SSH_KEY" ]] && SSH_OPTS=(-e "ssh -i ${SSH_KEY/#\~/$HOME}")

# A pipeline version bump invalidates every artefact built before it, and the reader sees
# that as the map silently reverting to the placeholder and superposition refusing to run.
# Deploying the bump BEFORE rebuilding is how that reaches the live site, so say so here
# rather than letting it be discovered by looking at the page.
if [[ -x .venv/bin/python ]]; then
  stale=$(.venv/bin/python CODSWALLOP.py artefacts --missing 2>/dev/null \
          | grep -c "MISSING" || true)
  if [[ "${stale:-0}" -gt 0 ]]; then
    echo "WARNING: ${stale} local famil$([[ $stale -eq 1 ]] && echo y || echo ies) have an"
    echo "         out-of-date embedding. Deploying now puts them on the placeholder map"
    echo "         until you run:  python CODSWALLOP.py warm  &&  bash deploy/push_embeddings.sh"
    echo
  fi
fi

echo "==> Syncing code to ${DROPLET_SSH}:${DROPLET_PATH}"
# ${arr[@]+"${arr[@]}"} expands to nothing when empty without tripping `set -u` (needed for
# macOS's bash 3.2, where "${arr[@]}" on an empty array is an error).
rsync -az --delete ${SSH_OPTS[@]+"${SSH_OPTS[@]}"} \
  --exclude '.venv/' --exclude 'data/' --exclude '__pycache__/' \
  --exclude '*.pyc' --exclude '.git/' --exclude '.env' \
  --exclude 'codswallop.db' --exclude 'codswallop.db-*' \
  ./ "${DROPLET_SSH}:${DROPLET_PATH}/"

echo "==> Installing dependencies + restarting the service on the droplet"
SSH_CMD=(ssh)
[[ -n "$SSH_KEY" ]] && SSH_CMD=(ssh -i "${SSH_KEY/#\~/$HOME}")
"${SSH_CMD[@]}" "$DROPLET_SSH" bash -s <<REMOTE
set -euo pipefail
cd "${DROPLET_PATH}"
if [[ ! -x .venv/bin/python ]]; then
  echo "No venv yet -- run deploy/provision.sh as root first."; exit 1
fi

# PIP_NO_CACHE_DIR silences pip's warning about the root-owned ~/.cache/pip under the
# codswallop user; the cache buys nothing here, since deps are already installed after the
# first provision and reinstalls are near-instant regardless.
sudo -u codswallop env PIP_NO_CACHE_DIR=1 ./.venv/bin/pip install --quiet -r requirements.txt

# rsync (run as root) leaves new code files root-owned, so chown them back. data/ may not
# exist on a first deploy, and \`find -path ... -prune\` on a missing directory is fine, but
# the whole block is kept off the failure path deliberately: a sibling app's deploy script
# aborted here on a missing data/ directory, which left the new code on disk while the
# service carried on running the old build, with the deploy reporting success.
mkdir -p "${DROPLET_PATH}/data"
sudo find "${DROPLET_PATH}" \\
  -path "${DROPLET_PATH}/.venv" -prune -o \\
  -name 'codswallop.db*' -prune -o \\
  -exec chown codswallop:codswallop {} + || echo "  (chown reported a problem; continuing to restart)"

sudo systemctl restart codswallop-web.service
sleep 1
sudo systemctl --no-pager --lines=3 status codswallop-web.service || true
REMOTE

echo "==> Deployed. Checking health..."
curl -fsS "https://${SERVER_NAME:-codswallop.mdeller.com}/healthz" && echo || \
  echo "    (health check failed -- if this is the first deploy, run provision.sh)"
