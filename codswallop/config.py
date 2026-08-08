"""Central configuration for CODSWALLOP.

All paths are anchored to the repository root (the parent of this package) so the code runs
identically from a laptop checkout and from /opt/codswallop on the droplet. The deploy
target comes from a gitignored .env (see .env.example); everything else has a sensible
default here.
"""

from __future__ import annotations

import os
from pathlib import Path

try:
    from dotenv import load_dotenv
except ImportError:  # dotenv is optional at import time (e.g. bare python3 running `init`)
    def load_dotenv(*_args, **_kwargs):  # type: ignore
        return False

# --------------------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------------------
PACKAGE_DIR = Path(__file__).resolve().parent
ROOT_DIR = PACKAGE_DIR.parent
DATA_DIR = ROOT_DIR / "data"
DB_PATH = ROOT_DIR / "codswallop.db"
STATIC_DIR = PACKAGE_DIR / "static"

load_dotenv(ROOT_DIR / ".env")

VERSION = "0.1.0"

# --------------------------------------------------------------------------------------
# External API endpoints (all public, none require a key)
# --------------------------------------------------------------------------------------
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"
INTERPRO_ENTRY_URL = "https://www.ebi.ac.uk/interpro/api/entry/{db}/{accession}"

# A descriptive User-Agent is the price of admission for hammering someone else's free API.
USER_AGENT = f"CODSWALLOP/{VERSION} (+https://codswallop.mdeller.com; marc@marcdeller.com)"

HTTP_TIMEOUT = 60          # seconds per request
# Concurrent requests to any one upstream. Low on purpose: these are free public APIs
# and the aim is to stop wasting round-trip latency, not to hammer them.
HTTP_WORKERS = int(os.environ.get('HTTP_WORKERS', '6'))
HTTP_MAX_RETRIES = 5

# --------------------------------------------------------------------------------------
# Family assembly
# --------------------------------------------------------------------------------------
# The identity-threshold slider's range and where it starts. 30 % is the conventional
# "same fold, plausibly same family" floor; the default sits higher so a first look at a
# family is not swamped by distant relatives.
IDENTITY_MIN = 30
IDENTITY_MAX = 100
IDENTITY_DEFAULT = 30

# Hard cap on entities pulled into one family. A free-text query like "kinase" would
# otherwise try to assemble a third of the PDB; past this we truncate and say so in the UI
# rather than melt the droplet.
#
# Truncation keeps the *closest* members, because the search returns them in score order.
# That is the right end to keep (they are the family proper) but it does narrow the identity
# range on a very large family, so the UI reports both the cap and the identity range it
# actually ended up with rather than implying the slider spans 30-100 %.
MAX_FAMILY_ENTITIES = int(os.environ.get("MAX_FAMILY_ENTITIES", "2500"))

# GraphQL batch size. RCSB tolerates larger, but 50 keeps each response small enough that a
# single flaky request costs little to retry.
GRAPHQL_BATCH = 50

# An entity whose deposited sequence runs this many residues longer than its aligned region
# is treated as carrying a fusion partner (MBP, GST, SUMO, BRIL and friends are all well
# over 100 residues). A blunt Phase 1 heuristic: Phase 2's construct diff engine replaces it
# with a real alignment against the UniProt canonical.
FUSION_EXCESS_RESIDUES = 80

# --------------------------------------------------------------------------------------
# Cache
# --------------------------------------------------------------------------------------
# Raw API responses. The PDB releases weekly, so a week is the natural staleness bound.
CACHE_TTL_HOURS = int(os.environ.get("CACHE_TTL_HOURS", "168"))
# Assembled families. Same reasoning: a family cannot gain a member mid-week.
FAMILY_TTL_HOURS = int(os.environ.get("FAMILY_TTL_HOURS", "168"))

# --------------------------------------------------------------------------------------
# Deploy / serving (from .env)
# --------------------------------------------------------------------------------------
DROPLET_SSH = os.environ.get("DROPLET_SSH", "")
DROPLET_PATH = os.environ.get("DROPLET_PATH", "/opt/codswallop")
SERVER_NAME = os.environ.get("SERVER_NAME", "codswallop.mdeller.com")
BIND_ADDR = os.environ.get("BIND_ADDR", "127.0.0.1:8006")


def ensure_dirs() -> None:
    """Create the runtime directory tree if missing. Safe to call repeatedly."""
    for d in (DATA_DIR, STATIC_DIR):
        d.mkdir(parents=True, exist_ok=True)
