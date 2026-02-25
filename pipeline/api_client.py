"""IPD REST API client with pagination, caching, and rate limiting."""

import json
import logging
import time
import hashlib
from pathlib import Path
from typing import Iterator, Optional
from urllib.parse import urlparse, parse_qs

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .config import API_BASE, PAGE_LIMIT, REQUEST_DELAY, MAX_RETRIES, RETRY_BACKOFF

logger = logging.getLogger(__name__)


class IPDClient:
    """Client for the IPD REST API."""

    def __init__(
        self,
        cache_dir: Optional[Path] = None,
        use_cache: bool = True,
        delay: float = REQUEST_DELAY,
    ):
        self.session = requests.Session()
        self.session.headers["Accept"] = "application/json"

        # Configure retry with exponential backoff
        retry = Retry(
            total=MAX_RETRIES,
            backoff_factor=RETRY_BACKOFF,
            status_forcelist=[500, 502, 503, 504],
            allowed_methods=["GET"],
        )
        self.session.mount("https://", HTTPAdapter(max_retries=retry))

        self.cache_dir = Path(cache_dir) if cache_dir else None
        self.use_cache = use_cache and self.cache_dir is not None
        self.delay = delay
        self._last_request_time = 0.0

    def _rate_limit(self) -> None:
        """Enforce minimum delay between requests."""
        elapsed = time.monotonic() - self._last_request_time
        if elapsed < self.delay:
            time.sleep(self.delay - elapsed)
        self._last_request_time = time.monotonic()

    def _cache_key(self, url: str, params: dict) -> Path:
        """Generate a filesystem-safe cache key."""
        raw = f"{url}|{json.dumps(params, sort_keys=True)}"
        digest = hashlib.sha256(raw.encode()).hexdigest()[:16]
        return self.cache_dir / f"{digest}.json"

    def _get(self, url: str, params: dict) -> dict:
        """GET with caching and rate limiting."""
        if self.use_cache:
            cache_path = self._cache_key(url, params)
            if cache_path.exists():
                return json.loads(cache_path.read_text())

        self._rate_limit()
        resp = self.session.get(url, params=params, timeout=60)
        resp.raise_for_status()
        data = resp.json()

        if self.use_cache:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            cache_path = self._cache_key(url, params)
            cache_path.write_text(json.dumps(data))

        return data

    def list_alleles(
        self,
        project: str,
        query: Optional[str] = None,
        fields: Optional[str] = None,
        limit: int = PAGE_LIMIT,
    ) -> Iterator[dict]:
        """
        Paginate through all alleles matching a query.

        Yields individual allele records (dicts).
        """
        url = f"{API_BASE}/allele"
        params = {"project": project, "limit": limit}
        if query:
            params["query"] = query
        if fields:
            params["fields"] = fields

        while True:
            data = self._get(url, params)
            yield from data.get("data", [])

            next_cursor = data.get("meta", {}).get("next")
            if not next_cursor:
                break

            # Parse next cursor URL to extract params
            parsed = parse_qs(urlparse(next_cursor).query)
            params = {k: v[0] for k, v in parsed.items()}
            params["project"] = project

    def get_allele(self, accession: str, project: str) -> dict:
        """Fetch a single full allele record by accession."""
        url = f"{API_BASE}/allele/{accession}"
        return self._get(url, {"project": project})

    def get_total_count(
        self, project: str, query: Optional[str] = None
    ) -> int:
        """Get total allele count without fetching all records."""
        url = f"{API_BASE}/allele"
        params = {"project": project, "limit": 1, "fields": "accession"}
        if query:
            params["query"] = query
        data = self._get(url, params)
        return data.get("meta", {}).get("total", 0)
