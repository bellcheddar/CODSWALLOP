"""Input classification, the cache's sentinel behaviour, and schema migration."""
import pytest
from codswallop import config, db, resolve


@pytest.mark.parametrize("text,kind", [
    ("4HHB", "pdb_id"),
    ("1abc_1", "entity_id"),
    ("P00918", "uniprot"),
    ("PF00062", "domain"),
    ("IPR001916", "domain"),
    ("LYZ", "gene"),
    ("carbonic anhydrase II", "text"),
    ("MKVFGRCELAAAMKRHGLDNYRGYSLGNW", "sequence"),
    ("", "empty"),
])
def test_input_classification(text, kind):
    assert resolve.classify(text) == kind


def test_nucleotide_sequence_is_classified_as_a_sequence():
    """So that `resolve` can say what is actually wrong.

    Classifying it as free text sent it down the full-text path, where the only possible
    answer was "nothing matches ACGTACGT...", and the specific message could never fire.
    """
    dna = "ACGT" * 8
    assert resolve.classify(dna) == "sequence"


def test_pasted_sequence_is_cleaned_of_fasta_furniture():
    raw = ">sp|P00698|LYSC_CHICK Lysozyme C\nKVFGRCELAA 10\nAMKRHGLDNY 20\n"
    assert resolve.clean_sequence(raw) == "KVFGRCELAAAMKRHGLDNY"


# ---- cache -----------------------------------------------------------------------------
@pytest.fixture
def fresh_db(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "t.db")
    import threading
    db._local = threading.local()
    db.init()
    return db


def test_a_cached_none_is_not_refetched(fresh_db):
    """The query most worth caching is the one that found nothing.

    Without the sentinel wrapper a legitimately empty answer (a UniProt 404, an RCSB 204)
    is indistinguishable from a cache miss and is re-fetched on every single request.
    """
    calls = []

    def fetch():
        calls.append(1)
        return None

    assert fresh_db.cached(("k",), fetch) is None
    assert fresh_db.cached(("k",), fetch) is None
    assert len(calls) == 1, "a cached None must count as an answer"


def test_cache_key_is_order_independent_for_dicts(fresh_db):
    a = fresh_db.cache_key({"x": 1, "y": 2})
    b = fresh_db.cache_key({"y": 2, "x": 1})
    assert a == b


def test_a_parse_version_change_orphans_the_old_entry(fresh_db):
    """Searches are cached in PARSED form, so a new field is absent from every cached row.

    Adding the alignment spans without bumping the version left the coverage figure reading
    100 % and the fusion count reading 0: both look like answers, not like missing inputs.
    """
    fresh_db.cached(("search", 1, "abc"), lambda: {"old": True})
    got = fresh_db.cached(("search", 2, "abc"), lambda: {"new": True})
    assert got == {"new": True}


def test_migrations_add_columns_to_an_existing_table(fresh_db):
    """CREATE TABLE IF NOT EXISTS silently does nothing to a table that already exists."""
    conn = fresh_db.connect()
    conn.execute("DROP TABLE entity")
    conn.execute("CREATE TABLE entity (slug TEXT, entity_id TEXT, pdb_id TEXT)")
    conn.commit()
    fresh_db.init()
    cols = {r["name"] for r in conn.execute("PRAGMA table_info(entity)")}
    assert {"query_beg", "query_end", "uniprot_ids"} <= cols


# The artefact queue. A family a reader assembles on the live site has an embedding on
# neither machine, and the droplet is the only one that knows while being the only one that
# cannot fix it.
def test_a_request_is_recorded_once_and_counted_thereafter(tmp_path, monkeypatch):
    from codswallop import config, db
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "q.db")
    monkeypatch.setattr(db._local, "conn", None, raising=False)
    db.init()
    db.request_artefact("gfp-p42212", "P42212", "both", "GFP", 1138)
    db.request_artefact("gfp-p42212", "P42212", "both", "GFP", 1138)
    rows = db.open_requests()
    assert len(rows) == 1 and rows[0]["hits"] == 2, "the queue is a priority order, not a log"

    db.mark_request_served("gfp-p42212")
    assert db.open_requests() == []


def test_a_family_that_needs_more_than_it_did_is_requeued(tmp_path, monkeypatch):
    """Served means served for what was asked. A family whose contacts were built and whose
    embedding was later invalidated by a pipeline bump must come back, or a version bump
    silently strands every family that had already been through the queue once."""
    from codswallop import config, db
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "q2.db")
    monkeypatch.setattr(db._local, "conn", None, raising=False)
    db.init()
    db.request_artefact("abl1", "P00519", "contacts", "ABL1", 1988)
    db.mark_request_served("abl1")
    assert db.open_requests() == []
    db.request_artefact("abl1", "P00519", "both", "ABL1", 1988)
    assert [r["slug"] for r in db.open_requests()] == ["abl1"]


def test_recording_a_request_never_breaks_the_page(monkeypatch):
    """It runs in the request path. A queue that cannot be written is a missing rebuild; an
    exception there would be a 500 on a page that is otherwise perfectly serviceable."""
    from codswallop import db
    monkeypatch.setattr(db, "connect", lambda: (_ for _ in ()).throw(db.sqlite3.Error("nope")))
    db.request_artefact("x", "y", "both")      # must not raise


def test_a_family_built_by_an_older_parser_is_not_fresh(tmp_path, monkeypatch):
    """PARSE_VERSION is part of every HTTP cache key, so bumping it re-fetches and re-parses
    from the API. But family/entry/entity are a second cache downstream of that one, and the
    bump never reached them: adding the assembly fields left every already-filed family
    reporting no biological assembly at all. The panel rendered, said "no assembly
    annotation for this family", and looked like a property of the data rather than a stale
    row. A warm then "refreshed" all 32 families without rebuilding one of them."""
    from codswallop import config, db, rcsb
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "pv.db")
    monkeypatch.setattr(db._local, "conn", None, raising=False)
    db.init()
    fam = {"slug": "x", "query": "P1", "kind": "uniprot", "seed": "P1", "name": "X",
           "organism": "o", "seed_sequence": "AAA", "seed_length": 3, "pfam": [],
           "interpro": [], "identity_threshold": 30, "total_hits": 1, "truncated": False}
    db.save_family(fam, [], [])
    assert db.family_fresh("x") is True

    monkeypatch.setattr(rcsb, "PARSE_VERSION", rcsb.PARSE_VERSION + 1)
    assert db.family_fresh("x") is False, \
        "a family built by an older parser is stale however recently it was written"


def test_a_request_whose_artefact_arrived_is_no_longer_open(tmp_path, monkeypatch):
    """The queue is cleared by drain_queue.sh when it builds something, but an artefact
    pushed from a workstation that was not draining the queue left the row behind: four
    families sat in the queue asking for contacts they already had."""
    from codswallop import artefacts, config, db
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "q3.db")
    monkeypatch.setattr(db._local, "conn", None, raising=False)
    db.init()
    db.request_artefact("abl1", "P00519", "contacts", "ABL1", 2000)
    assert [r["slug"] for r in db.open_requests()] == ["abl1"]

    # The artefact turns up by some other route.
    monkeypatch.setattr(artefacts, "status", lambda slug: {
        "embedding": {"current": True}, "contacts": {"current": True},
        "topology": {"current": True}})
    assert db.open_requests() == []
    # And it is marked served, so the next call does not have to work it out again.
    row = db.connect().execute(
        "SELECT served_at FROM artefact_request WHERE slug='abl1'").fetchone()
    assert row["served_at"] is not None
