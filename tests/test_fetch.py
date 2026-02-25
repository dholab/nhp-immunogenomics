"""Tests for fetch orchestration."""

import responses

from pipeline.api_client import IPDClient
from pipeline.config import API_BASE
from pipeline.fetch import (
    detect_changed_alleles,
    fetch_allele_listing,
    fetch_full_records,
    group_by_species_locus,
)


class TestFetchAlleleListing:
    @responses.activate
    def test_returns_alleles(self, tmp_path):
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
        client = IPDClient(cache_dir=tmp_path / "cache", delay=0)
        result = fetch_allele_listing(client, "MHC")
        assert len(result) == 2


class TestFetchFullRecords:
    @responses.activate
    def test_fetches_individual_records(self, tmp_path):
        for acc in ["NHP00001", "NHP00002"]:
            responses.add(
                responses.GET,
                f"{API_BASE}/allele/{acc}",
                json={"accession": acc, "name": f"Test-{acc}"},
            )
        client = IPDClient(cache_dir=tmp_path / "cache", delay=0)
        result = fetch_full_records(client, "MHC", ["NHP00001", "NHP00002"])
        assert len(result) == 2

    @responses.activate
    def test_skips_failed_records(self, tmp_path):
        responses.add(
            responses.GET,
            f"{API_BASE}/allele/NHP00001",
            json={"accession": "NHP00001"},
        )
        responses.add(
            responses.GET,
            f"{API_BASE}/allele/NHP00002",
            status=500,
            json={"error": "server error"},
        )
        client = IPDClient(cache_dir=tmp_path / "cache", use_cache=False, delay=0)
        result = fetch_full_records(client, "MHC", ["NHP00001", "NHP00002"])
        assert len(result) == 1

    @responses.activate
    def test_progress_callback(self, tmp_path):
        responses.add(
            responses.GET,
            f"{API_BASE}/allele/NHP00001",
            json={"accession": "NHP00001"},
        )
        client = IPDClient(cache_dir=tmp_path / "cache", delay=0)
        calls = []
        fetch_full_records(client, "MHC", ["NHP00001"],
                          progress_callback=lambda i, t: calls.append((i, t)))
        assert calls == [(1, 1)]


class TestDetectChangedAlleles:
    def test_detects_new(self):
        listing = [
            {"accession": "NHP00001", "date_modified": "2025-01-01"},
            {"accession": "NHP00002", "date_modified": "2025-01-01"},
        ]
        known = {"NHP00001": "2025-01-01"}
        new, modified = detect_changed_alleles(listing, known)
        assert new == ["NHP00002"]
        assert modified == []

    def test_detects_modified(self):
        listing = [
            {"accession": "NHP00001", "date_modified": "2025-06-01"},
        ]
        known = {"NHP00001": "2025-01-01"}
        new, modified = detect_changed_alleles(listing, known)
        assert new == []
        assert modified == ["NHP00001"]

    def test_no_changes(self):
        listing = [
            {"accession": "NHP00001", "date_modified": "2025-01-01"},
        ]
        known = {"NHP00001": "2025-01-01"}
        new, modified = detect_changed_alleles(listing, known)
        assert new == []
        assert modified == []


class TestGroupBySpeciesLocus:
    def test_groups_correctly(self):
        listing = [
            {"accession": "NHP00001", "organism.name": "Mamu", "locus": "A1"},
            {"accession": "NHP00002", "organism.name": "Mamu", "locus": "A1"},
            {"accession": "NHP00003", "organism.name": "Mamu", "locus": "B"},
            {"accession": "NHP00004", "organism.name": "Mafa", "locus": "A1"},
        ]
        groups = group_by_species_locus(listing)
        assert ("Mamu", "A1") in groups
        assert len(groups[("Mamu", "A1")]) == 2
        assert ("Mamu", "B") in groups
        assert len(groups[("Mamu", "B")]) == 1
        assert ("Mafa", "A1") in groups
