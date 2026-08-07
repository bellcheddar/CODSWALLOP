"""SQLite cache: assembled families, their entries and their entities, plus a raw-response
cache in front of every external API.

Two layers, deliberately:

* `http_cache` sits directly in front of `http.py`. It is keyed by the request itself, so a
  second family that happens to share a GraphQL batch, a UniProt lookup or a Pfam name with
  the first pays nothing for it.
* `family` / `entry` / `entity` hold an assembled family. The second query for a family is
  a handful of indexed SELECTs, which is what "must be instant" means in practice.

Both carry a TTL. The PDB releases weekly, so a week is the natural staleness bound for
each: a family cannot gain a member mid-week.
"""

from __future__ import annotations

import hashlib
import json
import sqlite3
import threading
import time
from typing import Any, Optional

from . import config

_local = threading.local()

SCHEMA = """
CREATE TABLE IF NOT EXISTS http_cache (
    key        TEXT PRIMARY KEY,
    body       TEXT NOT NULL,
    fetched_at INTEGER NOT NULL
);
CREATE INDEX IF NOT EXISTS ix_http_cache_fetched ON http_cache(fetched_at);

CREATE TABLE IF NOT EXISTS family (
    slug               TEXT PRIMARY KEY,
    query              TEXT NOT NULL,      -- the raw user input, verbatim
    kind               TEXT NOT NULL,      -- pdb_id|uniprot|gene|pfam|interpro|sequence|text
    seed               TEXT,               -- resolved seed identifier
    name               TEXT,
    organism           TEXT,
    seed_sequence      TEXT,
    seed_length        INTEGER,
    pfam               TEXT,               -- JSON [{id, name}]
    interpro           TEXT,               -- JSON [{id, name}]
    identity_threshold INTEGER NOT NULL,
    total_hits         INTEGER,            -- hits the search reported, before the cap
    truncated          INTEGER NOT NULL DEFAULT 0,
    built_at           INTEGER NOT NULL
);

CREATE TABLE IF NOT EXISTS entry (
    slug           TEXT NOT NULL,
    pdb_id         TEXT NOT NULL,
    title          TEXT,
    method         TEXT,
    resolution     REAL,
    r_work         REAL,
    r_free         REAL,
    space_group    TEXT,
    cell           TEXT,                   -- JSON {a,b,c,alpha,beta,gamma}
    deposit_date   TEXT,
    release_date   TEXT,
    chain_count    INTEGER,
    entity_count   INTEGER,
    assembly_count INTEGER,
    ligands        TEXT,                   -- JSON [{id, name, formula}]
    citation       TEXT,                   -- JSON {title, journal, year, doi, authors, ...}
    PRIMARY KEY (slug, pdb_id)
);
CREATE INDEX IF NOT EXISTS ix_entry_slug ON entry(slug);

CREATE TABLE IF NOT EXISTS entity (
    slug           TEXT NOT NULL,
    entity_id      TEXT NOT NULL,          -- RCSB polymer entity id, e.g. 132L_1
    pdb_id         TEXT NOT NULL,
    description    TEXT,
    chains         TEXT,                   -- JSON ["A", "B"]
    seq_length     INTEGER,
    sequence       TEXT,
    organism       TEXT,
    taxonomy_id    INTEGER,
    host_organism  TEXT,
    uniprot        TEXT,                   -- primary (first-listed) UniProt accession
    uniprot_ids    TEXT,                   -- JSON list: a chimera has more than one
    identity       REAL,                   -- % identity to the family seed
    aligned_length INTEGER,
    query_beg      INTEGER,                -- span of the SEED this member aligns to; the
    query_end      INTEGER,                -- union across members is the coverage census
    is_fusion      INTEGER NOT NULL DEFAULT 0,
    is_orthologue  INTEGER NOT NULL DEFAULT 0,
    pfam           TEXT,                   -- JSON [{id, name}]
    interpro       TEXT,                   -- JSON [{id, name}]
    PRIMARY KEY (slug, entity_id)
);
CREATE INDEX IF NOT EXISTS ix_entity_slug ON entity(slug);
CREATE INDEX IF NOT EXISTS ix_entity_pdb ON entity(slug, pdb_id);
"""


def connect() -> sqlite3.Connection:
    """Per-thread connection. gunicorn sync workers are one thread apiece, but the CLI and
    Flask's dev server are not, and a SQLite connection is not safe to share across threads.
    """
    conn = getattr(_local, "conn", None)
    if conn is None:
        config.DB_PATH.parent.mkdir(parents=True, exist_ok=True)
        conn = sqlite3.connect(config.DB_PATH, timeout=30)
        conn.row_factory = sqlite3.Row
        # WAL lets the readers serving page requests carry on while a family is being
        # assembled and written, which is the whole concurrency story of this app.
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA synchronous=NORMAL")
        conn.execute("PRAGMA foreign_keys=ON")
        _local.conn = conn
    return conn


# Columns added after the first release. `CREATE TABLE IF NOT EXISTS` silently does nothing
# to an existing table, so a new column in SCHEMA above never reaches a database that
# already exists: the app then reads the column back as NULL and computes a plausible wrong
# answer instead of failing. Every column added from here on gets a line here too.
_MIGRATIONS: list[tuple[str, str, str]] = [
    ("entity", "query_beg", "INTEGER"),
    ("entity", "query_end", "INTEGER"),
    ("entity", "uniprot_ids", "TEXT"),
]


def init() -> None:
    """Create the schema and apply any column additions. Idempotent."""
    conn = connect()
    conn.executescript(SCHEMA)
    for table, column, decl in _MIGRATIONS:
        existing = {r["name"] for r in conn.execute(f"PRAGMA table_info({table})")}
        if column not in existing:
            conn.execute(f"ALTER TABLE {table} ADD COLUMN {column} {decl}")
    conn.commit()


# --------------------------------------------------------------------------------------
# Raw response cache
# --------------------------------------------------------------------------------------
def cache_key(*parts: Any) -> str:
    """Stable key for a request. JSON payloads are dumped with sorted keys so that two
    logically identical queries built in a different dict order share a cache entry."""
    blob = "\x1f".join(
        json.dumps(p, sort_keys=True, separators=(",", ":")) if isinstance(p, (dict, list)) else str(p)
        for p in parts
    )
    return hashlib.sha256(blob.encode()).hexdigest()


def cache_get(key: str, ttl_hours: Optional[int] = None) -> Optional[Any]:
    ttl = config.CACHE_TTL_HOURS if ttl_hours is None else ttl_hours
    row = connect().execute(
        "SELECT body, fetched_at FROM http_cache WHERE key = ?", (key,)
    ).fetchone()
    if row is None:
        return None
    if ttl and time.time() - row["fetched_at"] > ttl * 3600:
        return None
    # A cached `null` is a real answer (a 404 from UniProt, an entry with no ligands), so it
    # is stored and returned as such. Callers distinguish "no cache entry" from "cached
    # nothing" via the sentinel wrapper below.
    return json.loads(row["body"])


def cache_put(key: str, value: Any) -> None:
    conn = connect()
    conn.execute(
        "INSERT INTO http_cache(key, body, fetched_at) VALUES (?, ?, ?) "
        "ON CONFLICT(key) DO UPDATE SET body = excluded.body, fetched_at = excluded.fetched_at",
        (key, json.dumps(value, separators=(",", ":")), int(time.time())),
    )
    conn.commit()


def cached(key_parts: tuple, fetch, ttl_hours: Optional[int] = None) -> Any:
    """Return a cached value or call `fetch()` and store the result.

    The stored value is wrapped in a one-element list so that a legitimately cached `None`
    (UniProt 404, RCSB 204 no-hits) is distinguishable from a cache miss. Without the
    wrapper every empty answer would be re-fetched on every request, which is exactly the
    query most worth caching: the one that found nothing.
    """
    key = cache_key(*key_parts)
    hit = cache_get(key, ttl_hours)
    if hit is not None:
        return hit[0]
    value = fetch()
    cache_put(key, [value])
    return value


# --------------------------------------------------------------------------------------
# Families
# --------------------------------------------------------------------------------------
def family_fresh(slug: str, ttl_hours: Optional[int] = None) -> bool:
    ttl = config.FAMILY_TTL_HOURS if ttl_hours is None else ttl_hours
    row = connect().execute("SELECT built_at FROM family WHERE slug = ?", (slug,)).fetchone()
    if row is None:
        return False
    return not ttl or (time.time() - row["built_at"]) <= ttl * 3600


def save_family(fam: dict, entries: list[dict], entities: list[dict]) -> None:
    """Write an assembled family, replacing any previous build of the same slug.

    One transaction: a half-written family that the next request reads as complete is worse
    than no cached family at all.
    """
    conn = connect()
    slug = fam["slug"]
    with conn:
        conn.execute("DELETE FROM entry  WHERE slug = ?", (slug,))
        conn.execute("DELETE FROM entity WHERE slug = ?", (slug,))
        conn.execute(
            "INSERT INTO family(slug, query, kind, seed, name, organism, seed_sequence, "
            "  seed_length, pfam, interpro, identity_threshold, total_hits, truncated, built_at) "
            "VALUES (:slug, :query, :kind, :seed, :name, :organism, :seed_sequence, "
            "  :seed_length, :pfam, :interpro, :identity_threshold, :total_hits, :truncated, :built_at) "
            "ON CONFLICT(slug) DO UPDATE SET "
            "  query=excluded.query, kind=excluded.kind, seed=excluded.seed, name=excluded.name, "
            "  organism=excluded.organism, seed_sequence=excluded.seed_sequence, "
            "  seed_length=excluded.seed_length, pfam=excluded.pfam, interpro=excluded.interpro, "
            "  identity_threshold=excluded.identity_threshold, total_hits=excluded.total_hits, "
            "  truncated=excluded.truncated, built_at=excluded.built_at",
            {
                "slug": slug,
                "query": fam["query"],
                "kind": fam["kind"],
                "seed": fam.get("seed"),
                "name": fam.get("name"),
                "organism": fam.get("organism"),
                "seed_sequence": fam.get("seed_sequence"),
                "seed_length": fam.get("seed_length"),
                "pfam": json.dumps(fam.get("pfam") or []),
                "interpro": json.dumps(fam.get("interpro") or []),
                "identity_threshold": fam.get("identity_threshold", config.IDENTITY_DEFAULT),
                "total_hits": fam.get("total_hits"),
                "truncated": 1 if fam.get("truncated") else 0,
                "built_at": int(time.time()),
            },
        )
        conn.executemany(
            "INSERT INTO entry(slug, pdb_id, title, method, resolution, r_work, r_free, "
            "  space_group, cell, deposit_date, release_date, chain_count, entity_count, "
            "  assembly_count, ligands, citation) "
            "VALUES (:slug, :pdb_id, :title, :method, :resolution, :r_work, :r_free, "
            "  :space_group, :cell, :deposit_date, :release_date, :chain_count, :entity_count, "
            "  :assembly_count, :ligands, :citation)",
            [
                {
                    "slug": slug, "pdb_id": e["pdb_id"], "title": e.get("title"),
                    "method": e.get("method"), "resolution": e.get("resolution"),
                    "r_work": e.get("r_work"), "r_free": e.get("r_free"),
                    "space_group": e.get("space_group"),
                    "cell": json.dumps(e.get("cell")) if e.get("cell") else None,
                    "deposit_date": e.get("deposit_date"), "release_date": e.get("release_date"),
                    "chain_count": e.get("chain_count"), "entity_count": e.get("entity_count"),
                    "assembly_count": e.get("assembly_count"),
                    "ligands": json.dumps(e.get("ligands") or []),
                    "citation": json.dumps(e.get("citation")) if e.get("citation") else None,
                }
                for e in entries
            ],
        )
        conn.executemany(
            "INSERT INTO entity(slug, entity_id, pdb_id, description, chains, seq_length, "
            "  sequence, organism, taxonomy_id, host_organism, uniprot, identity, "
            "  aligned_length, query_beg, query_end, is_fusion, is_orthologue, pfam, interpro, uniprot_ids) "
            "VALUES (:slug, :entity_id, :pdb_id, :description, :chains, :seq_length, "
            "  :sequence, :organism, :taxonomy_id, :host_organism, :uniprot, :identity, "
            "  :aligned_length, :query_beg, :query_end, :is_fusion, :is_orthologue, :pfam, :interpro, "
            "  :uniprot_ids)",
            [
                {
                    "slug": slug, "entity_id": t["entity_id"], "pdb_id": t["pdb_id"],
                    "description": t.get("description"),
                    "chains": json.dumps(t.get("chains") or []),
                    "seq_length": t.get("seq_length"), "sequence": t.get("sequence"),
                    "organism": t.get("organism"), "taxonomy_id": t.get("taxonomy_id"),
                    "host_organism": t.get("host_organism"), "uniprot": t.get("uniprot"),
                    "uniprot_ids": json.dumps(t.get("uniprot_ids") or []),
                    "identity": t.get("identity"), "aligned_length": t.get("aligned_length"),
                    "query_beg": t.get("query_beg"), "query_end": t.get("query_end"),
                    "is_fusion": 1 if t.get("is_fusion") else 0,
                    "is_orthologue": 1 if t.get("is_orthologue") else 0,
                    "pfam": json.dumps(t.get("pfam") or []),
                    "interpro": json.dumps(t.get("interpro") or []),
                }
                for t in entities
            ],
        )


def load_family(slug: str) -> Optional[dict]:
    """Read an assembled family back as the shape the API and templates expect."""
    conn = connect()
    row = conn.execute("SELECT * FROM family WHERE slug = ?", (slug,)).fetchone()
    if row is None:
        return None

    entries = []
    for r in conn.execute("SELECT * FROM entry WHERE slug = ? ORDER BY pdb_id", (slug,)):
        e = dict(r)
        e.pop("slug", None)
        e["cell"] = json.loads(e["cell"]) if e["cell"] else None
        e["ligands"] = json.loads(e["ligands"] or "[]")
        e["citation"] = json.loads(e["citation"]) if e["citation"] else None
        entries.append(e)

    entities = []
    for r in conn.execute("SELECT * FROM entity WHERE slug = ? ORDER BY entity_id", (slug,)):
        t = dict(r)
        t.pop("slug", None)
        t["chains"] = json.loads(t["chains"] or "[]")
        t["uniprot_ids"] = json.loads(t["uniprot_ids"] or "[]")
        t["pfam"] = json.loads(t["pfam"] or "[]")
        t["interpro"] = json.loads(t["interpro"] or "[]")
        t["is_fusion"] = bool(t["is_fusion"])
        t["is_orthologue"] = bool(t["is_orthologue"])
        entities.append(t)

    fam = dict(row)
    fam["pfam"] = json.loads(fam["pfam"] or "[]")
    fam["interpro"] = json.loads(fam["interpro"] or "[]")
    fam["truncated"] = bool(fam["truncated"])
    fam["entries"] = entries
    fam["entities"] = entities
    return fam


def recent_families(limit: int = 12) -> list[dict]:
    """Most recently assembled families, for the landing page's drawer of previous work."""
    rows = connect().execute(
        "SELECT f.slug, f.name, f.organism, f.kind, f.built_at, "
        "       (SELECT COUNT(*) FROM entry e WHERE e.slug = f.slug) AS n_entries "
        "FROM family f ORDER BY f.built_at DESC LIMIT ?",
        (limit,),
    ).fetchall()
    return [dict(r) for r in rows]


# --------------------------------------------------------------------------------------
# Housekeeping
# --------------------------------------------------------------------------------------
def purge_expired() -> int:
    """Drop http_cache rows past their TTL. Returns the number removed."""
    conn = connect()
    cutoff = int(time.time() - config.CACHE_TTL_HOURS * 3600)
    with conn:
        cur = conn.execute("DELETE FROM http_cache WHERE fetched_at < ?", (cutoff,))
    return cur.rowcount


def db_size_bytes() -> int:
    try:
        return config.DB_PATH.stat().st_size
    except OSError:
        return 0


def stats() -> dict:
    conn = connect()
    one = lambda q: conn.execute(q).fetchone()[0]  # noqa: E731
    return {
        "families": one("SELECT COUNT(*) FROM family"),
        "entries": one("SELECT COUNT(*) FROM entry"),
        "entities": one("SELECT COUNT(*) FROM entity"),
        "cached_responses": one("SELECT COUNT(*) FROM http_cache"),
        "db_mb": round(db_size_bytes() / 1048576, 1),
    }
