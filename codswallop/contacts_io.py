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
# 4: seed positions are mapped by ALIGNING each of the member's own chains to the seed,
#    replacing `resnr + (query_beg - 1)`. That added the seed offset to PLIP's author
#    residue number, which already carries it, so it was counted twice: JAK1's hot residues
#    came out at 1,340 and 2,110 on a 1,154-residue seed, and 52 of 71 families were wrong.
#    Every version-3 artefact is therefore untrustworthy and must be rebuilt, not migrated:
#    the raw residue numbers are not kept, so the double offset cannot be undone after the
#    fact. Contacts from chains that are not the family member are no longer counted at all.
VERSION = 4
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
