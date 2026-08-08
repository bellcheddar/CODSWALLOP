"""Reading the PLIP contacts artefact.

Split from `contacts.py` for the same reason `embed_io` is split from `embed`: that module
shells out to PLIP and OpenBabel and imports gemmi, none of which exist on the droplet. This
one imports json and nothing else, and is what the web app calls.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Optional

from . import config

logger = logging.getLogger(__name__)

# 3: hot residues carry their own contact-type and ligand breakdowns.
VERSION = 3
CONTACT_DIR = config.DATA_DIR / "contacts"


def artefact_path(slug: str) -> Path:
    return CONTACT_DIR / f"{slug}.json"


def load(slug: str) -> Optional[dict]:
    p = artefact_path(slug)
    if not p.exists():
        return None
    try:
        data = json.loads(p.read_text())
    except (OSError, ValueError):
        logger.warning("contacts for %s are unreadable", slug)
        return None
    return data if data.get("version") == VERSION else None
