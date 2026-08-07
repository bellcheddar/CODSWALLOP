"""Gunicorn config for the CODSWALLOP web service. Referenced by codswallop-web.service."""
import os

bind = os.environ.get("BIND_ADDR", "127.0.0.1:8006")
workers = int(os.environ.get("WEB_WORKERS", "3"))
worker_class = "sync"
# Assembling a cold family means a sequence search plus tens of batched GraphQL calls
# against someone else's API. Twelve minutes is generous, but the alternative is a worker
# killed halfway through a build that would have been cached for a week.
timeout = 720
graceful_timeout = 30
keepalive = 5
# Log to stdout/stderr so journald captures everything (journalctl -u codswallop-web).
accesslog = "-"
errorlog = "-"
loglevel = os.environ.get("LOG_LEVEL", "info")
proc_name = "codswallop-web"
