"""The pre-warm list against the front page.

The landing page offers example buttons. A reader who presses one and lands on a placeholder
map has been invited down a path the app did not prepare, and the map then says the family
has been queued, which reads as the feature being broken rather than as it being honest.

That happened: 4HHB is on the front page, and going in that way builds
`hemoglobin-subunit-alpha-4hhb-1`, a different family from the pre-warmed
`hemoglobin-subunit-alpha-p69905`, because a slug derives from its seed. Same protein, two
addresses, artefacts on only one of them.
"""
import importlib.util
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def _prewarm():
    spec = importlib.util.spec_from_file_location("prewarm", ROOT / "pipeline" / "prewarm.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _example_buttons():
    html = (ROOT / "codswallop" / "templates" / "landing.html").read_text()
    return [m.strip() for m in
            re.findall(r'<button[^>]*class="eg"[^>]*>(.*?)</button>', html, re.S)]


def test_every_landing_example_is_pre_warmed():
    pw = _prewarm()
    buttons = _example_buttons()
    assert buttons, "the landing page should still offer examples"
    missing = [b for b in buttons if b not in pw.EXAMPLE_SEEDS]
    assert not missing, f"front-page examples with no pre-warm entry: {missing}"


def test_no_stale_example_entries():
    """The other direction: an entry here for a button that no longer exists is a family
    being rebuilt every night for a route nobody can take."""
    pw = _prewarm()
    buttons = set(_example_buttons())
    stale = [k for k in pw.EXAMPLE_SEEDS if k not in buttons]
    assert not stale, f"pre-warm entries for examples not on the page: {stale}"


def test_example_seeds_reach_the_target_list():
    pw = _prewarm()
    seeds = {s for v in pw.EXAMPLE_SEEDS.values() for s in v}
    listed = {t[0] for t in pw.TARGETS}
    assert seeds <= listed, f"example seeds missing from TARGETS: {seeds - listed}"


def test_targets_are_deduplicated():
    """A protein can legitimately belong to two groups in the list; building it twice is
    only slower, never different."""
    pw = _prewarm()
    seeds = [t[0] for t in pw.TARGETS]
    assert len(seeds) == len(set(seeds)), "duplicate seeds in TARGETS"
