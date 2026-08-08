"""Route behaviour, without touching the network."""
import pytest
from codswallop import config, db


@pytest.fixture
def client(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "t.db")
    import threading
    db._local = threading.local()
    from codswallop.webapp import create_app
    app = create_app()
    app.config["TESTING"] = True
    return app.test_client()


def test_landing_page_renders(client):
    r = client.get("/")
    assert r.status_code == 200
    assert b"CODSWALLOP" in r.data


def test_healthz(client):
    assert client.get("/healthz").get_json()["ok"] is True


def test_stats_endpoint_exists_because_it_is_the_hit_beacon(client):
    """mdeller.com counts a page view by watching for this request in the access log.

    It was declared as the beacon in apps.json before any page requested it, which would
    have left the launcher's hit count reading zero forever with nothing to show why.
    """
    assert client.get("/api/stats").status_code == 200


def test_every_page_actually_requests_the_beacon():
    js = (config.PACKAGE_DIR / "static" / "app.js").read_text()
    assert "/api/stats" in js, "app.js loads on every page and must fire the beacon"


def test_empty_lookup_redirects_home(client):
    assert client.get("/lookup?q=").status_code == 302


def test_cold_permalink_with_no_query_is_a_404_not_a_guess(client):
    r = client.get("/f/nothing-is-filed-here")
    assert r.status_code == 404


def test_api_refuses_to_build_without_a_query(client):
    assert client.get("/api/family/unknown-slug").status_code == 404


def test_html_is_not_heuristically_cacheable(client):
    """Flask sends no Cache-Control on a template, and heuristic caching then pins the old
    ?v= asset URLs, making a CSS or JS deploy invisible to anyone who has visited before."""
    r = client.get("/")
    assert "no-cache" in r.headers.get("Cache-Control", "")


def test_stats_reports_families_left_on_the_placeholder(client):
    """A stale artefact version leaves the map silently non-structural.

    Three artefact version bumps in one afternoon each invalidated every artefact, and the
    page renders either way: the map just quietly stops being a measurement. This is the
    only place that says so without a shell on the server.
    """
    art = client.get("/api/stats").get_json()["artefacts"]
    assert art is not None
    assert "on_placeholder" in art
    assert art["embedding_version"] >= 1


# ---- the dossier ---------------------------------------------------------------------
# Patched at the family boundary rather than by assembling a real one: decorate() reaches
# UniProt for the construct references, and no test here is allowed near the network. What
# these tests are about is the document, not the assembly that feeds it.
FAKE = {
    "slug": "lysozyme-c-p00698", "name": "Lysozyme C", "organism": "Gallus gallus",
    "seed": "P00698", "seed_length": 147, "seed_sequence": "KVFGRCELAA",
    "identity_threshold": 30, "truncated": False, "total_hits": 1686,
    "stats": {"entries": 1686, "entities": 1687, "constructs": 280, "organisms": 26,
              "holo_entries": 330, "best_resolution": 0.65, "median_resolution": 1.73,
              "tagged": 1, "engineered": 1606, "fusions": 0},
    "constructs": [{"n_entities": 1239, "length": 129, "best_resolution": 0.65,
                    "best_pdb_id": "2VB1", "summary": "residues 19-147"}],
    "msa": {"engineered": []}, "assemblies": {"n": 0}, "domains": {"domains": []},
    "ligands": {"components": []}, "crystals": {"n": 0, "n_parsed": 0},
    "quality": {"n": 0}, "orthologues": [], "citations": {},
}


@pytest.fixture
def dossier_client(client, monkeypatch):
    from codswallop import db as db_mod, family as family_mod, webapp
    monkeypatch.setattr(webapp.db, "family_fresh", lambda slug, *a, **k: slug == FAKE["slug"])
    monkeypatch.setattr(webapp.db, "load_family", lambda slug: dict(FAKE))
    monkeypatch.setattr(family_mod, "decorate", lambda fam: dict(FAKE))
    return client


def test_dossier_is_self_contained(dossier_client):
    """The point of it is to outlive the session: something to attach to a grant appendix or
    hand over. So it must load nothing. No script, stylesheet, font, image, or url() in the
    CSS, or it degrades to a broken page the moment it is opened offline."""
    r = dossier_client.get("/f/lysozyme-c-p00698/dossier")
    assert r.status_code == 200
    html = r.get_data(as_text=True)
    for tag in ("<script", "<link", "<img", "<iframe", "url("):
        assert tag not in html, f"the dossier must not pull in {tag!r}"
    assert "@media print" in html, "a dossier that cannot become a PDF is not a dossier"


def test_dossier_refuses_to_assemble_a_family_on_demand(dossier_client):
    """Assembly takes up to ninety seconds and this URL is exactly what a link checker
    fetches. A family that is not filed gets a 404 with an explanation, never a build."""
    assert dossier_client.get("/f/not-a-real-family-x9/dossier").status_code == 404


def test_dossier_survives_a_citation_with_no_year(dossier_client, monkeypatch):
    """Undated papers are common in the archive, and Jinja's sort filter compares None to an
    int directly, so one of them took the whole document down with a TypeError."""
    from codswallop import family as family_mod
    fam = dict(FAKE)
    fam["citations"] = {"a": {"title": "Undated paper", "year": None, "journal": "J"},
                        "b": {"title": "Dated paper", "year": 2007, "journal": "J"}}
    monkeypatch.setattr(family_mod, "decorate", lambda f: fam)
    r = dossier_client.get("/f/lysozyme-c-p00698/dossier")
    assert r.status_code == 200
    body = r.get_data(as_text=True)
    assert "Undated paper" in body and "Dated paper" in body
    assert body.index("Dated paper") < body.index("Undated paper"), \
        "dated papers lead; undated ones go last rather than crashing the sort"


def test_a_stale_permalink_rebuilds_instead_of_404ing(client, monkeypatch):
    """Stale is not absent. Families go stale weekly by design and a pipeline version bump
    stales every one at once, so without this a shared link rots on a timer: every filed
    family on the live site returned "that drawer is empty" the moment the parser moved."""
    from codswallop import db as db_mod, webapp
    monkeypatch.setattr(webapp.db, "family_fresh", lambda slug, *a, **k: False)
    monkeypatch.setattr(webapp.db, "family_row",
                        lambda slug: {"slug": slug, "query": "P00698", "name": "Lysozyme C"})
    r = client.get("/f/lysozyme-c-p00698")
    assert r.status_code == 200, "a family that has been filed before must rebuild itself"
    assert b"Lysozyme C" in r.data


def test_a_slug_never_filed_still_says_so(client, monkeypatch):
    from codswallop import webapp
    monkeypatch.setattr(webapp.db, "family_fresh", lambda slug, *a, **k: False)
    monkeypatch.setattr(webapp.db, "family_row", lambda slug: None)
    assert client.get("/f/never-heard-of-it-x9").status_code == 404
