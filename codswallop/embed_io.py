"""Reading the structural embedding artefact.

Split from `embed.py` deliberately: that module imports numpy, biotite and tmtools and runs
only on a workstation, while this one imports nothing beyond the standard library and is
what the web app calls. Without the split, a droplet without the pipeline installed would
fail to import `family` at all.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Optional

from . import config

logger = logging.getLogger(__name__)

# 6: the superposition reference and the AlphaFold accession come from the family's seed
# rather than from a majority vote among its members, which the subject was losing.
VERSION = 6
EMBED_DIR = config.DATA_DIR / "embeddings"


def artefact_path(slug: str) -> Path:
    return EMBED_DIR / f"{slug}.json"


def load(slug: str) -> Optional[dict]:
    """The embedding for one family, or None when there is not one (or it is stale)."""
    p = artefact_path(slug)
    if not p.exists():
        return None
    try:
        data = json.loads(p.read_text())
    except (OSError, ValueError):
        logger.warning("embedding for %s is unreadable; using the placeholder", slug)
        return None
    if data.get("version") != VERSION:
        logger.info("embedding for %s is version %s, expected %s; using the placeholder",
                    slug, data.get("version"), VERSION)
        return None
    return data
