"""Tests for release version tracking."""

import json

import responses

from pipeline.release_tracker import (
    has_new_release,
    load_version_file,
    resolve_version,
    save_version_file,
    scrape_current_release,
)


class TestScrapeCurrentRelease:
    @responses.activate
    def test_parses_version(self):
        responses.add(
            responses.GET,
            "https://www.ebi.ac.uk/ipd/mhc/",
            body='<html>Release 3.16.0.0 (2026-01)</html>',
        )
        result = scrape_current_release("MHC")
        assert result["version"] == "3.16.0.0"

    @responses.activate
    def test_version_not_found(self):
        responses.add(
            responses.GET,
            "https://www.ebi.ac.uk/ipd/mhc/",
            body="<html>No version here</html>",
        )
        result = scrape_current_release("MHC")
        assert result["version"] == "unknown"


class TestVersionFile:
    def test_save_and_load(self, tmp_path):
        data = {"MHC": {"version": "3.16.0.0", "allele_count": 10125}}
        save_version_file(tmp_path, data)
        loaded = load_version_file(tmp_path)
        assert loaded["MHC"]["version"] == "3.16.0.0"

    def test_load_missing(self, tmp_path):
        assert load_version_file(tmp_path) == {}


class TestHasNewRelease:
    @responses.activate
    def test_detects_new_release(self, tmp_path):
        save_version_file(tmp_path, {"MHC": {"version": "3.15.0.0"}})
        responses.add(
            responses.GET,
            "https://www.ebi.ac.uk/ipd/mhc/",
            body='<html>Release 3.16.0.0</html>',
        )
        responses.add(
            responses.GET,
            "https://www.ebi.ac.uk/ipd/nhkir/",
            body='<html>Release 1.7.0.0</html>',
        )
        save_version_file(tmp_path, {
            "MHC": {"version": "3.15.0.0"},
            "NHKIR": {"version": "1.7.0.0"},
        })
        result = has_new_release(tmp_path)
        assert result["MHC"] is True
        assert result["NHKIR"] is False

    @responses.activate
    def test_no_new_release(self, tmp_path):
        save_version_file(tmp_path, {
            "MHC": {"version": "3.16.0.0"},
            "NHKIR": {"version": "1.7.0.0"},
        })
        responses.add(
            responses.GET,
            "https://www.ebi.ac.uk/ipd/mhc/",
            body='<html>Release 3.16.0.0</html>',
        )
        responses.add(
            responses.GET,
            "https://www.ebi.ac.uk/ipd/nhkir/",
            body='<html>Release 1.7.0.0</html>',
        )
        result = has_new_release(tmp_path)
        assert result["MHC"] is False
        assert result["NHKIR"] is False

    @responses.activate
    def test_detects_new_release_with_ipd_version_field(self, tmp_path):
        """has_new_release works with the newer ipd_version field."""
        save_version_file(tmp_path, {
            "MHC": {"ipd_version": "3.15.0.0", "version": "3.15.0.0.2", "revision": 2},
            "NHKIR": {"ipd_version": "1.7.0.0", "version": "1.7.0.0", "revision": 0},
        })
        responses.add(
            responses.GET,
            "https://www.ebi.ac.uk/ipd/mhc/",
            body='<html>Release 3.16.0.0</html>',
        )
        responses.add(
            responses.GET,
            "https://www.ebi.ac.uk/ipd/nhkir/",
            body='<html>Release 1.7.0.0</html>',
        )
        result = has_new_release(tmp_path)
        assert result["MHC"] is True
        assert result["NHKIR"] is False


class TestResolveVersion:
    def test_new_upstream_release(self):
        stored = {"ipd_version": "3.15.0.0", "version": "3.15.0.0.3", "revision": 3}
        assert resolve_version("3.16.0.0", stored, False) == "3.16.0.0"

    def test_new_upstream_release_ignores_allele_changes(self):
        stored = {"ipd_version": "3.15.0.0", "version": "3.15.0.0", "revision": 0}
        assert resolve_version("3.16.0.0", stored, True) == "3.16.0.0"

    def test_silent_update_first_revision(self):
        stored = {"ipd_version": "3.16.0.0", "version": "3.16.0.0", "revision": 0}
        assert resolve_version("3.16.0.0", stored, True) == "3.16.0.0.1"

    def test_silent_update_increments_revision(self):
        stored = {"ipd_version": "3.16.0.0", "version": "3.16.0.0.2", "revision": 2}
        assert resolve_version("3.16.0.0", stored, True) == "3.16.0.0.3"

    def test_no_changes_returns_stored_version(self):
        stored = {"ipd_version": "3.16.0.0", "version": "3.16.0.0.1", "revision": 1}
        assert resolve_version("3.16.0.0", stored, False) == "3.16.0.0.1"

    def test_backwards_compat_no_ipd_version(self):
        """Falls back to 'version' field when ipd_version is missing."""
        stored = {"version": "3.16.0.0", "revision": 0}
        assert resolve_version("3.16.0.0", stored, True) == "3.16.0.0.1"

    def test_empty_stored(self):
        """First run with no stored data treats as new release."""
        assert resolve_version("3.16.0.0", {}, False) == "3.16.0.0"
