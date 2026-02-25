"""Orchestrates fetching alleles from the IPD API."""

import logging
from typing import Optional

from .api_client import IPDClient
from .config import LISTING_FIELDS

logger = logging.getLogger(__name__)


def fetch_allele_listing(
    client: IPDClient,
    project: str,
    query: Optional[str] = None,
) -> list[dict]:
    """
    Fetch metadata listing for all NHP alleles in a project.

    Returns list of dicts with flat dot-notation keys
    (e.g., 'organism.name', 'organism.scientificName').
    """
    if query is None:
        if project == "MHC":
            # MHC contains all taxonomic groups; filter to NHP only
            query = "and(eq(status,public),eq(organism.group,NHP))"
        else:
            query = "eq(status,public)"
    alleles = list(client.list_alleles(project, query=query, fields=LISTING_FIELDS))
    logger.info("Fetched listing: %d alleles from %s", len(alleles), project)
    return alleles


def fetch_full_records(
    client: IPDClient,
    project: str,
    accessions: list[str],
    progress_callback: Optional[callable] = None,
) -> list[dict]:
    """
    Fetch full records (with sequences) for a list of accessions.

    Returns list of complete allele record dicts.
    """
    records = []
    total = len(accessions)
    for i, acc in enumerate(accessions):
        try:
            record = client.get_allele(acc, project)
            records.append(record)
        except Exception as e:
            logger.error("Failed to fetch %s: %s", acc, e)
        if progress_callback:
            progress_callback(i + 1, total)
    logger.info("Fetched %d/%d full records from %s", len(records), total, project)
    return records


def detect_changed_alleles(
    current_listing: list[dict],
    known_alleles: dict[str, str],
) -> tuple[list[str], list[str]]:
    """
    Compare current API listing against known alleles.

    Args:
        current_listing: Alleles from fetch_allele_listing()
        known_alleles: Mapping of accession -> date_modified from previous build

    Returns:
        (new_accessions, modified_accessions)
    """
    new = []
    modified = []
    for allele in current_listing:
        acc = allele["accession"]
        date_mod = allele.get("date_modified", "")
        if acc not in known_alleles:
            new.append(acc)
        elif date_mod != known_alleles[acc]:
            modified.append(acc)
    return new, modified


def group_by_species_locus(
    listing: list[dict],
) -> dict[tuple[str, str], list[str]]:
    """
    Group allele accessions by (species_prefix, locus).

    Returns {("Mamu", "A1"): ["NHP01224", "NHP01225", ...], ...}
    """
    groups: dict[tuple[str, str], list[str]] = {}
    for allele in listing:
        species = allele.get("organism.name", allele.get("organism", {}).get("name", "unknown"))
        locus = allele.get("locus", "unknown")
        acc = allele["accession"]
        groups.setdefault((species, locus), []).append(acc)
    return groups
