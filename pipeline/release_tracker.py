"""Track IPD release versions and detect updates."""

import json
import logging
import re
from pathlib import Path

import requests

from .config import VERSION_FILE

logger = logging.getLogger(__name__)

RELEASE_PAGES = {
    "MHC": "https://www.ebi.ac.uk/ipd/mhc/",
    "NHKIR": "https://www.ebi.ac.uk/ipd/nhkir/",
}


def scrape_current_release(project: str) -> dict:
    """
    Scrape the current IPD release version from the project homepage.

    Returns {"version": "3.16.0.0", "date": "2026-01"} or similar.
    """
    url = RELEASE_PAGES[project]
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()

    match = re.search(r"Release\s+([\d.]+)", resp.text)
    if match:
        return {"version": match.group(1)}
    return {"version": "unknown"}


def load_version_file(repo_root: Path) -> dict:
    """Load version.json from the repo root."""
    path = repo_root / VERSION_FILE
    if path.exists():
        return json.loads(path.read_text())
    return {}


def save_version_file(repo_root: Path, data: dict) -> None:
    """Write version.json."""
    path = repo_root / VERSION_FILE
    path.write_text(json.dumps(data, indent=2) + "\n")


def has_new_release(repo_root: Path) -> dict[str, bool]:
    """
    Check if either IPD database has a newer release.

    Returns {"MHC": True/False, "NHKIR": True/False}.
    """
    current = load_version_file(repo_root)
    result = {}
    for project in ("MHC", "NHKIR"):
        try:
            remote = scrape_current_release(project)
            local_base = current.get(project, {}).get("ipd_version", "")
            # Fall back to "version" for backwards compatibility with
            # version.json files written before the ipd_version field existed.
            if not local_base:
                local_base = current.get(project, {}).get("version", "")
            result[project] = remote["version"] != local_base
        except Exception as e:
            logger.error("Failed to check %s release: %s", project, e)
            result[project] = False
    return result


def resolve_version(
    ipd_version: str,
    stored: dict,
    has_allele_changes: bool,
) -> str:
    """
    Determine the display version, appending a revision suffix for silent updates.

    If the upstream IPD version changed, return it directly (e.g., "3.17.0.0").
    If the IPD version is unchanged but alleles changed, increment a local
    revision suffix (e.g., "3.16.0.0.1", "3.16.0.0.2").
    If nothing changed, return the current stored version as-is.

    Args:
        ipd_version: Version string scraped from IPD homepage.
        stored: The existing version.json entry for this project.
        has_allele_changes: Whether detect_changed_alleles found differences.

    Returns:
        The version string to store and display.
    """
    stored_ipd = stored.get("ipd_version", "")
    # Backwards compat: older version.json files lack ipd_version
    if not stored_ipd:
        stored_ipd = stored.get("version", "")

    if ipd_version != stored_ipd:
        # New upstream release — reset revision
        return ipd_version

    if has_allele_changes:
        # Silent update — increment revision
        current_revision = stored.get("revision", 0)
        return f"{ipd_version}.{current_revision + 1}"

    # No changes
    return stored.get("version", ipd_version)
