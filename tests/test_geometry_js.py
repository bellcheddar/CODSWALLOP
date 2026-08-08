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
