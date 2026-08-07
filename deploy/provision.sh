#!/usr/bin/env bash
# One-time droplet provisioning for CODSWALLOP. Run as root ON the droplet, AFTER the code
# is present at /opt/codswallop (push it first with deploy/deploy.sh from your Mac).
#
#   sudo SERVER_NAME=codswallop.mdeller.com bash /opt/codswallop/deploy/provision.sh
#
# Idempotent: safe to re-run. Installs system packages, a service user, the Python venv, the
# systemd unit, the nginx site, and a Let's Encrypt certificate.
set -euo pipefail

APP_DIR=/opt/codswallop
APP_USER=codswallop
BIND_ADDR="${BIND_ADDR:-127.0.0.1:8006}"

# Pull SERVER_NAME/BIND_ADDR from .env if not passed in the environment.
if [[ -f "$APP_DIR/.env" ]]; then
  set -a; # shellcheck disable=SC1091
  source "$APP_DIR/.env"; set +a
fi
SERVER_NAME="${SERVER_NAME:-codswallop.mdeller.com}"

echo "==> CODSWALLOP provisioning for ${SERVER_NAME} on ${BIND_ADDR}"

if [[ $EUID -ne 0 ]]; then echo "Run as root (sudo)."; exit 1; fi
if [[ ! -f "$APP_DIR/wsgi.py" ]]; then
  echo "No code at $APP_DIR. Push it first: bash deploy/deploy.sh (from your Mac)."; exit 1
fi

echo "==> Installing system packages"
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq
apt-get install -y -qq python3-venv python3-pip nginx certbot python3-certbot-nginx rsync

echo "==> Creating service user '${APP_USER}'"
id -u "$APP_USER" &>/dev/null || useradd --system --home "$APP_DIR" --shell /usr/sbin/nologin "$APP_USER"
mkdir -p "$APP_DIR/data"
chown -R "$APP_USER:$APP_USER" "$APP_DIR"

echo "==> Building Python venv"
if [[ ! -x "$APP_DIR/.venv/bin/python" ]]; then
  sudo -u "$APP_USER" python3 -m venv "$APP_DIR/.venv"
fi
sudo -u "$APP_USER" "$APP_DIR/.venv/bin/pip" install --quiet --upgrade pip
sudo -u "$APP_USER" "$APP_DIR/.venv/bin/pip" install --quiet -r "$APP_DIR/requirements.txt"

echo "==> Initialising the cache database"
sudo -u "$APP_USER" "$APP_DIR/.venv/bin/python" "$APP_DIR/CODSWALLOP.py" init

echo "==> Installing the shared long-cache snippet"
# Shared across every vhost on this droplet. Written here too so CODSWALLOP can be
# provisioned onto a fresh droplet without mdeller-landing having gone first.
mkdir -p /etc/nginx/snippets
if [[ ! -f /etc/nginx/snippets/long-cache.conf ]]; then
  cat > /etc/nginx/snippets/long-cache.conf <<'SNIP'
# Long-lived caching for rarely-changing static assets (images/audio/fonts/js/css).
# Include inside any location block that serves such files directly from disk.
# Shared across vhosts on the droplet.
add_header Cache-Control "public, max-age=31536000, immutable" always;

# Logging stays ON here, deliberately. These asset requests are the cheapest honest
# evidence that a real browser rendered a page -- a scanner fetches the HTML and stops --
# so mdeller.com's per-app hit counter reads exactly this traffic (see apps.json's
# "beacon" fields). `access_log off;` here silently zeroes the count for every app whose
# beacon is a static file, with nothing to show that it happened.
SNIP
fi

echo "==> Installing the systemd unit"
cp "$APP_DIR/deploy/codswallop-web.service" /etc/systemd/system/
systemctl daemon-reload
systemctl enable --now codswallop-web.service

echo "==> Installing the nginx site"
sed -e "s|__SERVER_NAME__|${SERVER_NAME}|g" -e "s|__BIND_ADDR__|${BIND_ADDR}|g" \
  "$APP_DIR/deploy/nginx-codswallop.conf" > /etc/nginx/sites-available/codswallop
ln -sf /etc/nginx/sites-available/codswallop /etc/nginx/sites-enabled/codswallop
nginx -t && systemctl reload nginx

echo "==> Requesting/deploying the TLS certificate (certbot)"
# Always run certbot --nginx, never skip on "a certificate already exists". The site file
# above is unconditionally re-templated from source on every run (plain HTTP only, no SSL
# block), so a skip-if-cert-exists check would leave a freshly re-templated vhost with no
# SSL directives at all. nginx then answers port 443 for this hostname with *some other*
# vhost's certificate: a real incident on a sibling app's re-provision, where
# boltzmaker.mdeller.com served AlphaFraud's certificate. certbot's own --nginx plugin
# already reuses a still-valid certificate rather than re-requesting one, so running it
# unconditionally is both safe and correctly idempotent.
#
# Retried, because certbot's droplet-wide renewal timer can hold certbot's lock at the
# moment this runs ("Another instance of Certbot is already running"), hit during that
# same incident.
certbot_ok=0
for attempt in 1 2 3; do
  if certbot --nginx -d "$SERVER_NAME" --non-interactive --agree-tos \
       -m "${CERTBOT_EMAIL:-marc@marcdeller.com}" --redirect; then
    certbot_ok=1
    break
  fi
  echo "    certbot attempt ${attempt}/3 failed, retrying in 10s..."
  sleep 10
done
if [[ "$certbot_ok" -ne 1 ]]; then
  echo "    certbot failed after 3 attempts (DNS not pointed yet? renewal-timer lock"
  echo "    contention?). ${SERVER_NAME} has NO SSL config right now -- re-run:"
  echo "    certbot --nginx -d ${SERVER_NAME}"
  exit 1
fi

# Certbot's `listen 443 ssl;` lines do not enable HTTP/2 on nginx 1.24, so add it here,
# idempotently, which also fixes an already-provisioned site on a re-run.
if grep -q "listen.*443 ssl" /etc/nginx/sites-available/codswallop && \
   ! grep -q "listen.*443 ssl http2" /etc/nginx/sites-available/codswallop; then
  echo "==> Enabling HTTP/2"
  python3 - <<'PYEOF'
import re
p = "/etc/nginx/sites-available/codswallop"
text = open(p).read()
text = re.sub(
    r'listen ((?:\[::\]:)?443) ssl( ipv6only=on)?;',
    lambda m: f'listen {m.group(1)} ssl http2{m.group(2) or ""};',
    text,
)
open(p, "w").write(text)
PYEOF
  nginx -t && systemctl reload nginx
fi

echo "==> Done. Status:"
systemctl --no-pager --lines=3 status codswallop-web || true
echo "    Site:   https://${SERVER_NAME}/"
echo "    Health: https://${SERVER_NAME}/healthz"
