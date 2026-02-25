"""Tests for the IPD REST API client."""

import json
import time

import pytest
import responses

from pipeline.api_client import IPDClient
from pipeline.config import API_BASE


@pytest.fixture
def client(tmp_path):
    """API client with cache in tmp dir and no rate limiting."""
    return IPDClient(cache_dir=tmp_path / "cache", use_cache=True, delay=0)


@pytest.fixture
def client_no_cache(tmp_path):
    """API client with caching disabled."""
    return IPDClient(cache_dir=tmp_path / "cache", use_cache=False, delay=0)


class TestListAlleles:
    @responses.activate
    def test_single_page(self, client):
        """List alleles returns all records from a single page."""
        responses.add(
            responses.GET,
            f"{API_BASE}/allele",
            json={
                "data": [
                    {"accession": "NHP00001", "name": "Test-A*001"},
                    {"accession": "NHP00002", "name": "Test-A*002"},
                ],
                "meta": {"next": None, "total": 2},
            },
        )
        result = list(client.list_alleles("MHC"))
        assert len(result) == 2
        assert result[0]["accession"] == "NHP00001"
        assert result[1]["accession"] == "NHP00002"

    @responses.activate
    def test_pagination(self, client_no_cache):
        """List alleles follows pagination cursors."""
        responses.add(
            responses.GET,
            f"{API_BASE}/allele",
            json={
                "data": [{"accession": "NHP00001"}],
                "meta": {"next": "?limit=1&next=cursor123", "total": 2},
            },
        )
        responses.add(
            responses.GET,
            f"{API_BASE}/allele",
            json={
                "data": [{"accession": "NHP00002"}],
                "meta": {"next": None, "total": 2},
            },
        )
        result = list(client_no_cache.list_alleles("MHC", limit=1))
        assert len(result) == 2
        assert result[0]["accession"] == "NHP00001"
        assert result[1]["accession"] == "NHP00002"

    @responses.activate
    def test_with_query_and_fields(self, client):
        """Query and fields parameters are passed through."""
        responses.add(
            responses.GET,
            f"{API_BASE}/allele",
            json={"data": [{"accession": "NHP00001"}], "meta": {"next": None}},
        )
        result = list(client.list_alleles(
            "MHC",
            query="eq(organism.name,Mamu)",
            fields="accession,name",
        ))
        assert len(result) == 1
        assert "project" in responses.calls[0].request.params
        assert responses.calls[0].request.params["project"] == "MHC"

    @responses.activate
    def test_empty_response(self, client):
        """Empty data returns no results."""
        responses.add(
            responses.GET,
            f"{API_BASE}/allele",
            json={"data": [], "meta": {"next": None, "total": 0}},
        )
        result = list(client.list_alleles("NHKIR"))
        assert result == []


class TestGetAllele:
    @responses.activate
    def test_fetch_single(self, client):
        """Fetch a single allele by accession."""
        record = {"accession": "NHP01224", "name": "Mamu-A1*026:01:01:01"}
        responses.add(
            responses.GET,
            f"{API_BASE}/allele/NHP01224",
            json=record,
        )
        result = client.get_allele("NHP01224", "MHC")
        assert result["accession"] == "NHP01224"

    @responses.activate
    def test_http_error(self, client_no_cache):
        """HTTP errors are raised."""
        responses.add(
            responses.GET,
            f"{API_BASE}/allele/NHP99999",
            json={"error": "not found"},
            status=404,
        )
        with pytest.raises(Exception):
            client_no_cache.get_allele("NHP99999", "MHC")


class TestGetTotalCount:
    @responses.activate
    def test_returns_total(self, client):
        """Total count is extracted from meta."""
        responses.add(
            responses.GET,
            f"{API_BASE}/allele",
            json={"data": [{"accession": "NHP00001"}], "meta": {"total": 9558}},
        )
        assert client.get_total_count("MHC") == 9558


class TestCaching:
    @responses.activate
    def test_cache_hit(self, client):
        """Second call uses cache, not network."""
        responses.add(
            responses.GET,
            f"{API_BASE}/allele/NHP00001",
            json={"accession": "NHP00001"},
        )
        client.get_allele("NHP00001", "MHC")
        client.get_allele("NHP00001", "MHC")
        assert len(responses.calls) == 1  # Only one HTTP call

    @responses.activate
    def test_cache_disabled(self, client_no_cache):
        """With cache disabled, every call hits network."""
        responses.add(
            responses.GET,
            f"{API_BASE}/allele/NHP00001",
            json={"accession": "NHP00001"},
        )
        client_no_cache.get_allele("NHP00001", "MHC")
        client_no_cache.get_allele("NHP00001", "MHC")
        assert len(responses.calls) == 2


class TestRateLimiting:
    @responses.activate
    def test_delay_between_requests(self, tmp_path):
        """Rate limiting enforces delay between requests."""
        client = IPDClient(cache_dir=tmp_path / "cache", use_cache=False, delay=0.1)
        responses.add(
            responses.GET,
            f"{API_BASE}/allele/NHP00001",
            json={"accession": "NHP00001"},
        )
        responses.add(
            responses.GET,
            f"{API_BASE}/allele/NHP00002",
            json={"accession": "NHP00002"},
        )
        start = time.monotonic()
        client.get_allele("NHP00001", "MHC")
        client.get_allele("NHP00002", "MHC")
        elapsed = time.monotonic() - start
        assert elapsed >= 0.1
