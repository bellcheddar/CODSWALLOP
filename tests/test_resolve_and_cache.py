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
