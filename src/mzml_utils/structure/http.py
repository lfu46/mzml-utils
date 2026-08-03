"""One HTTP session for every structure-database call.

Before this module the lab's structure fetchers made bare ``requests.get`` calls
with no retry, no backoff, and (in the CIF fetcher) no session at all.  Both EBI
and RCSB rate-limit by IP and answer with 429; without ``Retry-After`` handling a
batch over a few thousand accessions degrades into silent failures.

Policy here, in one place:

* **Retry** 4 times on 429/500/502/503/504 with exponential backoff, honouring
  ``Retry-After``.  ``raise_on_status=False`` so callers still see 400 and 404 and
  can tell "malformed identifier" from "no such model" -- that distinction is the
  whole point of ``FetchStatus``.
* **Throttle** to a few requests per second per host.  RCSB publishes no numeric
  limit and advises "a handful of requests per second"; EBI limits by IP.
* **Redirects followed.**  UniProt answers a merged accession with 303, not 404,
  so a client with redirects off silently loses those entries.
* **Descriptive User-Agent.**  Both EBI and RCSB ask for one, and UniProt sends
  ``Vary: User-Agent``.  Set ``WU_LAB_CONTACT`` to append a contact address; no
  personal address is baked into this source.
"""
from __future__ import annotations

import os
import threading
import time
from typing import Any, Dict, Optional, Tuple
from urllib.parse import urlsplit

import requests
from requests.adapters import HTTPAdapter

try:  # urllib3 v2 and v1 disagree on where Retry lives / what it is called
    from urllib3.util.retry import Retry
except ImportError:  # pragma: no cover - very old urllib3
    from requests.packages.urllib3.util.retry import Retry  # type: ignore

from .records import StructureFetchError

#: Requests per second, per host. Override with WU_LAB_STRUCTURE_RPS.
DEFAULT_RATE_LIMIT = float(os.environ.get("WU_LAB_STRUCTURE_RPS", "3"))
#: (connect, read) seconds -- matches the tuple already used by the InterPro fetchers.
DEFAULT_TIMEOUT: Tuple[float, float] = (10.0, 45.0)

RETRY_STATUSES = (429, 500, 502, 503, 504)

_session: Optional[requests.Session] = None
_session_lock = threading.Lock()
_last_request: Dict[str, float] = {}
_throttle_lock = threading.Lock()


def user_agent() -> str:
    """A descriptive UA, with an optional contact from ``WU_LAB_CONTACT``."""
    from mzml_utils import __version__

    ua = f"mzml-utils-structure/{__version__} (+https://github.com/longpingfu/mzml-utils)"
    contact = os.environ.get("WU_LAB_CONTACT", "").strip()
    return f"{ua} ({contact})" if contact else ua


def get_session() -> requests.Session:
    """The process-wide session. Connection pooling plus the retry policy."""
    global _session
    if _session is None:
        with _session_lock:
            if _session is None:
                s = requests.Session()
                retry = Retry(
                    total=4,
                    connect=4,
                    read=4,
                    status=4,
                    backoff_factor=0.5,
                    status_forcelist=RETRY_STATUSES,
                    allowed_methods=frozenset(["GET", "POST", "HEAD"]),
                    respect_retry_after_header=True,
                    raise_on_status=False,  # let callers see 400/404 themselves
                )
                adapter = HTTPAdapter(
                    max_retries=retry, pool_connections=16, pool_maxsize=16
                )
                s.mount("https://", adapter)
                s.mount("http://", adapter)
                s.headers.update(
                    {"User-Agent": user_agent(), "Accept-Encoding": "gzip, deflate"}
                )
                _session = s
    return _session


def reset_session() -> None:
    """Drop the cached session. For tests, and after changing WU_LAB_CONTACT."""
    global _session
    with _session_lock:
        if _session is not None:
            _session.close()
        _session = None


def _throttle(url: str, rate: float) -> None:
    """Space requests to the same host at most ``rate`` per second."""
    if rate <= 0:
        return
    host = urlsplit(url).netloc
    min_gap = 1.0 / rate
    with _throttle_lock:
        now = time.monotonic()
        wait = _last_request.get(host, 0.0) + min_gap - now
        if wait > 0:
            time.sleep(wait)
            now = time.monotonic()
        _last_request[host] = now


def request(
    method: str,
    url: str,
    *,
    params: Optional[Dict[str, Any]] = None,
    data: Optional[Any] = None,
    json_body: Optional[Any] = None,
    headers: Optional[Dict[str, str]] = None,
    timeout: Tuple[float, float] = DEFAULT_TIMEOUT,
    rate: float = DEFAULT_RATE_LIMIT,
    stream: bool = False,
) -> requests.Response:
    """Throttled, retrying request. Raises ``StructureFetchError`` on transport failure.

    A non-2xx response is *returned*, not raised -- status codes carry meaning
    here (400 malformed vs 404 absent vs 200-with-empty-body).
    """
    _throttle(url, rate)
    try:
        return get_session().request(
            method,
            url,
            params=params,
            data=data,
            json=json_body,
            headers=headers,
            timeout=timeout,
            allow_redirects=True,  # UniProt merges answer 303
            stream=stream,
        )
    except requests.RequestException as exc:
        raise StructureFetchError(f"{method} {url} failed: {exc}") from exc


def get_json(
    url: str,
    *,
    params: Optional[Dict[str, Any]] = None,
    timeout: Tuple[float, float] = DEFAULT_TIMEOUT,
    rate: float = DEFAULT_RATE_LIMIT,
) -> Tuple[int, Any]:
    """GET returning ``(status_code, parsed_json_or_None)``.

    Returns ``None`` for the body when the response is not 2xx or is not valid
    JSON, so callers branch on the status code rather than on an exception.
    """
    resp = request("GET", url, params=params, timeout=timeout, rate=rate)
    if not resp.ok:
        return resp.status_code, None
    try:
        return resp.status_code, resp.json()
    except ValueError:
        return resp.status_code, None


def post_json(
    url: str,
    payload: Any,
    *,
    timeout: Tuple[float, float] = DEFAULT_TIMEOUT,
    rate: float = DEFAULT_RATE_LIMIT,
) -> Tuple[int, Any]:
    """POST a JSON body, returning ``(status_code, parsed_json_or_None)``."""
    resp = request("POST", url, json_body=payload, timeout=timeout, rate=rate)
    if not resp.ok:
        return resp.status_code, None
    try:
        return resp.status_code, resp.json()
    except ValueError:
        return resp.status_code, None


def get_bytes(
    url: str,
    *,
    timeout: Tuple[float, float] = DEFAULT_TIMEOUT,
    rate: float = DEFAULT_RATE_LIMIT,
) -> Tuple[int, Optional[bytes]]:
    """GET a binary/text payload, returning ``(status_code, content_or_None)``."""
    resp = request("GET", url, timeout=timeout, rate=rate, stream=True)
    if not resp.ok:
        resp.close()
        return resp.status_code, None
    try:
        return resp.status_code, resp.content
    finally:
        resp.close()
