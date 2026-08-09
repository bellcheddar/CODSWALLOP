"""The browser-side superposition geometry, exercised in node.

It lives in JavaScript because it runs in the reader's browser, and it is checked here
because a wrong rotation matrix does not look wrong: the model lands near the structure and
reads as roughly aligned. Two separate bugs in it each gave 21 A RMSD on two clouds related
by an exact rotation, and neither was visible on screen.
"""
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


@pytest.mark.skipif(shutil.which("node") is None, reason="node is not installed")
def test_kabsch_recovers_a_known_rotation():
    r = subprocess.run(
        ["node", str(ROOT / "tests" / "js" / "kabsch_check.js"),
         str(ROOT / "codswallop" / "static" / "viewer.js")],
        capture_output=True, text=True, timeout=60,
    )
    assert r.returncode == 0, r.stderr or r.stdout


@pytest.mark.skipif(shutil.which("node") is None, reason="node is not installed")
def test_one_outlier_does_not_set_the_scale_for_the_whole_family():
    """Hen lysozyme's amyloid fibrils are the same sequence as the seed and TM 0.15-0.20 to
    every other representative, so they sit at x = -1.0 while the rest of the family spans
    -0.036 to +0.046. Fitting the panel to them left 1,686 of 1,688 nodes inside 5.4 % of its
    width, which is a map of one structure and a smudge."""
    r = subprocess.run(
        ["node", str(ROOT / "tests" / "js" / "fit_check.js"),
         str(ROOT / "codswallop" / "static" / "map.js")],
        capture_output=True, text=True, timeout=60,
    )
    assert r.returncode == 0, r.stderr or r.stdout


@pytest.mark.skipif(shutil.which("node") is None, reason="node is not installed")
def test_depth_does_not_override_the_dim_class():
    """Written as an inline opacity, the depth cue beat `.node.dim` and silently disabled
    every filter on the map: the entity count changed and nothing on the map did, which
    reads as the toggles being broken rather than as a styling fault."""
    r = subprocess.run(
        ["node", str(ROOT / "tests" / "js" / "depth_check.js"),
         str(ROOT / "codswallop" / "static" / "map.js"),
         str(ROOT / "codswallop" / "static" / "theme.css")],
        capture_output=True, text=True, timeout=60,
    )
    assert r.returncode == 0, r.stderr or r.stdout
