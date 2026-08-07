"""Gunicorn entrypoint:  gunicorn wsgi:app  (see deploy/codswallop-web.service)."""

from codswallop.webapp import create_app

app = create_app()
