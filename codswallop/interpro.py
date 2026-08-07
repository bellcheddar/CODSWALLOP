"""InterPro REST client: resolve a Pfam or InterPro accession to a name and type.

Used to label a family seeded from a domain identifier, and to give the Domains section of
the divider rail something real to show before Phase 2 fills it in properly.
"""

from __future__ import annotations

import re
from typing import Optional

from . import config, db, http

PFAM_RE = re.compile(r"^PF\d{5}$", re.I)
INTERPRO_RE = re.compile(r"^IPR\d{6}$", re.I)


def entry(accession: str) -> Optional[dict]:
    """Look up a PFxxxxx or IPRxxxxxx accession. None if it does not exist."""
    acc = accession.upper()
    if PFAM_RE.match(acc):
        source = "pfam"
    elif INTERPRO_RE.match(acc):
        source = "interpro"
    else:
        return None

    def fetch():
        body = http.get_json(config.INTERPRO_ENTRY_URL.format(db=source, accession=acc))
        if not body:
            return None
        meta = body.get("metadata") or {}
        return {
            "accession": meta.get("accession"),
            "source": source,
            "name": (meta.get("name") or {}).get("name"),
            "short_name": (meta.get("name") or {}).get("short"),
            "type": meta.get("type"),
            # Pfam entries name the InterPro entry they are integrated into; carrying it
            # through lets the header cite both identifiers a reader might recognise.
            "integrated": meta.get("integrated"),
        }

    return db.cached(("interpro_entry", acc), fetch)
