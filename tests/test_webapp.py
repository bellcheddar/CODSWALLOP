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
