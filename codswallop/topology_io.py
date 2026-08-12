"""Reading the topology artefact.

Split from `topology.py` for the same reason `embed_io` is split from `embed`: that module
shells out to DSSP and imports biotite, and neither is installed on the droplet. This one
imports nothing beyond the standard library and is what the web app calls.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Optional

from . import config

logger = logging.getLogger(__name__)

# 2: carries the PDBe 2D fold diagram alongside the DSSP elements.
# 3: the PDBe fold-diagram elements carry seed coordinates alongside the entry's own
# numbering, so the diagram's labels agree with every other panel's.
VERSION = 3
TOPOLOGY_DIR = config.DATA_DIR / "topology"


def artefact_path(slug: str) -> Path:
    return TOPOLOGY_DIR / f"{slug}.json"


def load(slug: str) -> Optional[dict]:
    """The topology artefact for a family, or None when it is absent or from an older run."""
    path = artefact_path(slug)
    if not path.exists():
        return None
    try:
        data = json.loads(path.read_text())
    except (OSError, ValueError):
        logger.warning("unreadable topology artefact for %s", slug, exc_info=True)
        return None
    if data.get("version") != VERSION:
        logger.info("topology artefact for %s is v%s, wanted v%s",
                    slug, data.get("version"), VERSION)
        return None
    return data
