"""Which families are missing a current embedding or contacts artefact.

Imports nothing heavy, because both sides need it and they are not the same machine:

* the **workstation** uses it to decide what to rebuild, then runs the pipeline;
* the **droplet** uses it to *report* what is stale, because it cannot rebuild anything:
  biotite, tmtools, PLIP and OpenBabel are deliberately not installed there.

That asymmetry is the whole reason this exists as a separate module. Three artefact version
bumps in one afternoon each silently invalidated every artefact, and a family whose embedding
is rejected falls back to the placeholder map with no error anywhere: the page still renders,
it just quietly stops being a measurement. Detection has to live where the serving happens
even though the fix cannot.
"""

from __future__ import annotations

import json
from typing import Optional

from . import config, contacts_io, db, embed_io, topology_io


def _version_of(path) -> Optional[int]:
    try:
        return json.loads(path.read_text()).get("version")
    except (OSError, ValueError):
        return None


def status(slug: str) -> dict:
    """The artefact state of one family."""
    emb, con = embed_io.artefact_path(slug), contacts_io.artefact_path(slug)
    top = topology_io.artefact_path(slug)
    emb_v = _version_of(emb) if emb.exists() else None
    con_v = _version_of(con) if con.exists() else None
    top_v = _version_of(top) if top.exists() else None
    return {
        "slug": slug,
        "embedding": {
            "present": emb.exists(),
            "version": emb_v,
            "current": emb_v == embed_io.VERSION,
            "wanted": embed_io.VERSION,
        },
        "contacts": {
            "present": con.exists(),
            "version": con_v,
            "current": con_v == contacts_io.VERSION,
            "wanted": contacts_io.VERSION,
        },
        "topology": {
            "present": top.exists(),
            "version": top_v,
            "current": top_v == topology_io.VERSION,
            "wanted": topology_io.VERSION,
        },
    }


def survey() -> list[dict]:
    """Every filed family, with its artefact state and the query that rebuilds it."""
    rows = db.connect().execute(
        "SELECT slug, query, name FROM family ORDER BY built_at DESC").fetchall()
    out = []
    for r in rows:
        st = status(r["slug"])
        st["query"] = r["query"]
        st["name"] = r["name"]
        out.append(st)
    return out


def stale(kind: str = "embedding") -> list[dict]:
    """Families whose `kind` artefact is missing or built by an older pipeline."""
    return [s for s in survey() if not s[kind]["current"]]


def summary() -> dict:
    """Counts, for /api/stats and for the weekly warm to log.

    `on_placeholder` is the number that matters: those families are rendering a
    sequence-identity placeholder where the map claims to show structural distance.
    """
    rows = survey()
    return {
        "families": len(rows),
        "embeddings_current": sum(1 for s in rows if s["embedding"]["current"]),
        "contacts_current": sum(1 for s in rows if s["contacts"]["current"]),
        "topology_current": sum(1 for s in rows if s["topology"]["current"]),
        "on_placeholder": sum(1 for s in rows if not s["embedding"]["current"]),
        "embedding_version": embed_io.VERSION,
        "contacts_version": contacts_io.VERSION,
        "topology_version": topology_io.VERSION,
    }
