"""Shared retrying HTTP session for the four external services (RCSB Search, RCSB GraphQL,
UniProt REST, InterPro REST).

One session gives connection reuse; tenacity gives polite backoff on the transient 5xx and
rate-limit responses these public APIs occasionally return. Every one of them is free and
maintained by someone else, so the rules here are: identify ourselves, batch what can be
batched, and back off when told to.
"""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Callable, Iterable, Optional

import requests
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

from . import config

logger = logging.getLogger(__name__)

_SESSION = requests.Session()
_SESSION.headers.update({"User-Agent": config.USER_AGENT, "Accept": "application/json"})

# Retry only on connection errors and the handful of server-side statuses worth retrying.
_RETRYABLE_STATUS = {429, 500, 502, 503, 504}


class RetryableStatus(Exception):
    """Raised for a retryable HTTP status so tenacity re-attempts the request."""

    def __init__(self, message: str, retry_after: Optional[float] = None):
        super().__init__(message)
        self.retry_after = retry_after


def _wait(retry_state) -> float:
    """Exponential backoff, but honour a server-supplied Retry-After when there is one.

    tenacity's wait_exponential has no idea the server just told us exactly how long to
    wait. On a 429 the header is the authoritative answer and guessing shorter is how an
    IP gets blocked.
    """
    exc = retry_state.outcome.exception() if retry_state.outcome else None
    if isinstance(exc, RetryableStatus) and exc.retry_after:
        return min(exc.retry_after, 120)
    return wait_exponential(multiplier=1, min=1, max=30)(retry_state)


_retry = retry(
    reraise=True,
    stop=stop_after_attempt(config.HTTP_MAX_RETRIES),
    wait=_wait,
    retry=retry_if_exception_type((requests.ConnectionError, requests.Timeout, RetryableStatus)),
)


def _check(resp: requests.Response, what: str) -> requests.Response:
    if resp.status_code in _RETRYABLE_STATUS:
        after = resp.headers.get("Retry-After")
        try:
            after_s = float(after) if after else None
        except ValueError:                  # Retry-After may be an HTTP-date; fall back to backoff
            after_s = None
        raise RetryableStatus(f"{resp.status_code} for {what}", after_s)
    return resp


@_retry
def get(url: str, **kwargs) -> requests.Response:
    return _check(_SESSION.get(url, timeout=config.HTTP_TIMEOUT, **kwargs), f"GET {url}")


@_retry
def post_json(url: str, payload: dict[str, Any], **kwargs) -> requests.Response:
    resp = _SESSION.post(url, json=payload, timeout=config.HTTP_TIMEOUT, **kwargs)
    return _check(resp, f"POST {url}")


def get_json(url: str, **kwargs) -> Optional[Any]:
    """GET returning parsed JSON, or None on 404 (a miss, not an error, for these APIs)."""
    resp = get(url, **kwargs)
    if resp.status_code == 404:
        return None
    resp.raise_for_status()
    return resp.json()


def get_bytes(url: str, **kwargs) -> Optional[bytes]:
    """GET returning the raw body, or None on 404.

    For the images the dossier embeds. Same retry and rate-limit behaviour as everything
    else here, because the RCSB's image CDN is the same host being a good citizen towards.
    """
    resp = get(url, **kwargs)
    if resp.status_code == 404:
        return None
    resp.raise_for_status()
    return resp.content


def post_search(url: str, payload: dict[str, Any]) -> Optional[dict]:
    """POST a search query, returning parsed JSON.

    RCSB Search answers a query with no hits with **204 No Content**, not an empty result
    set, so `resp.json()` on it raises. That is a legitimate outcome for a search (nothing
    matched your sequence), so it maps to None rather than an exception.
    """
    resp = post_json(url, payload)
    if resp.status_code == 204:
        return None
    resp.raise_for_status()
    return resp.json()


def post_form(url: str, fields: dict[str, Any]) -> Optional[Any]:
    """POST form-encoded, returning parsed JSON.

    ScanProsite is a CGI script from an older web and takes a form POST rather than JSON. It
    also 308-redirects its documented URL to a new path, and urllib re-issues a redirected
    POST as a GET, which loses the body and returns an HTML page that json refuses. requests
    keeps the method across a 308, so the redirect is followed correctly here, but the
    endpoint is pinned to the post-redirect URL anyway rather than relying on that.
    """
    # Accept is set to application/json on the shared session; ScanProsite ignores it and
    # keys off the `output` field instead, so it is left alone rather than special-cased.
    resp = _check(_SESSION.post(url, data=fields, timeout=config.HTTP_TIMEOUT),
                  f"POST {url}")
    if resp.status_code == 204 or not resp.content:
        return None
    try:
        return resp.json()
    except ValueError:
        # A CGI that answers with an HTML error page rather than a status code. Not fatal:
        # the panel it feeds is optional and a family without motifs still renders.
        logger.warning("%s answered with non-JSON (%d bytes)", url, len(resp.content))
        return None


def graphql(url: str, query: str, variables: Optional[dict] = None) -> dict:
    """POST a GraphQL query; return the `data` object (raises on transport/GraphQL error)."""
    resp = post_json(url, {"query": query, "variables": variables or {}})
    resp.raise_for_status()
    body = resp.json()
    if body.get("errors"):
        raise RuntimeError(f"GraphQL errors: {body['errors']}")
    return body.get("data") or {}


def parallel_map(fn: Callable, items: list, workers: Optional[int] = None) -> list:
    """Run `fn` over `items` concurrently, preserving order.

    A cold family is ~160 HTTP calls to two APIs, and run one at a time they were 57 of the
    75 seconds a large family took to assemble. The batches are independent by construction
    (each is a disjoint set of ids), so the only thing serialising them was the loop.

    Bounded deliberately, and low. These are free public APIs maintained by other people: the
    point is to stop wasting round-trip latency, not to open thirty sockets to the RCSB. Six
    concurrent requests is roughly what a browser opens to one host.

    Exceptions propagate, so a failed batch still fails the build rather than silently
    yielding a family with a hole in it.
    """
    if not items:
        return []
    n = workers or config.HTTP_WORKERS
    if n <= 1 or len(items) == 1:
        return [fn(x) for x in items]
    with ThreadPoolExecutor(max_workers=min(n, len(items))) as pool:
        return list(pool.map(fn, items))


def download(url: str, dest, skip_if_exists: bool = True, params: Optional[dict] = None,
             max_bytes: Optional[int] = None):
    """Stream a URL to `dest` (a Path). Returns dest, or None on 404.

    Written to a unique temporary name and moved into place, so a half-written file can
    never be mistaken for a complete one: the cache check is `exists()`, and an interrupted
    download that left a truncated mmCIF behind would be served as though it were whole for
    as long as the file sat there.
    """
    import os
    import threading
    from pathlib import Path

    dest = Path(dest)
    if skip_if_exists and dest.exists() and dest.stat().st_size > 0:
        return dest
    resp = get(url, stream=True, params=params) if params else get(url, stream=True)
    if resp.status_code == 404:
        return None
    resp.raise_for_status()
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(f"{dest.suffix}.{os.getpid()}.{threading.get_ident()}.part")
    try:
        written = 0
        with open(tmp, "wb") as fh:
            for chunk in resp.iter_content(chunk_size=1 << 16):
                written += len(chunk)
                # Abort rather than finish and refuse. A HEAD on files.rcsb.org times out,
                # so the size cannot be known in advance; streaming until the cap is passed
                # costs at most the cap, where downloading 453 MB and then declining to
                # parse it costs 453 MB of transfer and disk on a box that has neither to
                # spare.
                if max_bytes is not None and written > max_bytes:
                    resp.close()
                    return None
                fh.write(chunk)
        tmp.replace(dest)
    finally:
        if tmp.exists():
            tmp.unlink(missing_ok=True)
    return dest
