"""CLI entry points for the NHP immunogenomics pipeline."""

import argparse
import json
import logging
import sys
from pathlib import Path

from .api_client import IPDClient
from .config import CACHE_DIR, DATA_DIR, DOCS_DIR, LISTING_FIELDS
from .fetch import (
    detect_changed_alleles,
    fetch_allele_listing,
    fetch_full_records,
    group_by_species_locus,
)
from .genbank_builder import write_genbank_file
from .metadata_index import build_metadata_index
from .release_tracker import (
    has_new_release,
    load_version_file,
    resolve_version,
    save_version_file,
    scrape_current_release,
)

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="NHP Immunogenomics data pipeline"
    )
    parser.add_argument(
        "--repo-root", type=Path, default=Path(__file__).resolve().parent.parent,
        help="Repository root directory",
    )
    sub = parser.add_subparsers(dest="command")

    build_cmd = sub.add_parser("build", help="Full rebuild from API")
    build_cmd.add_argument("--no-cache", action="store_true")
    build_cmd.add_argument(
        "--project", choices=["MHC", "NHKIR", "all"], default="all"
    )

    sub.add_parser("check", help="Check if IPD has new releases")
    sub.add_parser("update", help="Incremental update")
    sub.add_parser("index", help="Rebuild metadata index only")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    repo_root = args.repo_root

    if args.command == "build":
        run_full_build(repo_root, args)
    elif args.command == "check":
        run_check(repo_root)
    elif args.command == "update":
        run_incremental_update(repo_root)
    elif args.command == "index":
        run_build_index(repo_root)
    else:
        parser.print_help()
        sys.exit(1)


def run_full_build(repo_root: Path, args) -> None:
    """Full rebuild: fetch all alleles, generate all files."""
    client = IPDClient(
        cache_dir=repo_root / CACHE_DIR,
        use_cache=not args.no_cache,
    )

    projects = ["MHC", "NHKIR"] if args.project == "all" else [args.project]
    all_listings: dict[str, list] = {}

    for project in projects:
        logger.info("=== Processing %s ===", project)

        listing = fetch_allele_listing(client, project)
        all_listings[project] = listing

        groups = group_by_species_locus(listing)
        logger.info("Grouped into %d species/locus combinations", len(groups))

        for (species, locus), accessions in groups.items():
            logger.info("  %s/%s: %d alleles", species, locus, len(accessions))

            records = fetch_full_records(client, project, accessions)

            data_dir = repo_root / DATA_DIR / project / species
            gb_path = data_dir / f"{locus}.gb"
            n = write_genbank_file(records, gb_path, project)
            logger.info("    Wrote %d GenBank records to %s", n, gb_path)

        # Update version
        try:
            release = scrape_current_release(project)
            ipd_version = release.get("version", "unknown")
            version_data = load_version_file(repo_root)
            version_data[project] = {
                "ipd_version": ipd_version,
                "version": ipd_version,
                "revision": 0,
                "allele_count": len(listing),
            }
            save_version_file(repo_root, version_data)
        except Exception as e:
            logger.warning("Could not update version info for %s: %s", project, e)

    # Build metadata index
    run_build_index(repo_root, all_listings)
    logger.info("=== Full build complete ===")


def run_check(repo_root: Path) -> None:
    """Check for new IPD releases. Exit 0 if update needed, 1 if current."""
    updates = has_new_release(repo_root)
    if any(updates.values()):
        for project, needs_update in updates.items():
            if needs_update:
                logger.info("%s has a new release available", project)
        sys.exit(0)
    else:
        logger.info("No new releases detected.")
        sys.exit(1)


def run_incremental_update(repo_root: Path) -> None:
    """Fetch only new/modified alleles and regenerate affected files."""
    client = IPDClient(cache_dir=repo_root / CACHE_DIR)
    all_listings: dict[str, list] = {}
    project_changed: dict[str, bool] = {}
    total_new = 0
    total_mod = 0
    any_changes = False

    # Load existing index to get known alleles
    index_path = repo_root / DOCS_DIR / "alleles.json"
    if not index_path.exists():
        logger.error("No existing index found. Run 'build' first.")
        sys.exit(1)

    existing_index = json.loads(index_path.read_text())
    known = {a["a"]: a.get("dm", "") for a in existing_index.get("alleles", [])}

    for project in ("MHC", "NHKIR"):
        logger.info("=== Checking %s for updates ===", project)
        listing = fetch_allele_listing(client, project)
        all_listings[project] = listing

        new_accs, mod_accs = detect_changed_alleles(listing, known)
        if not new_accs and not mod_accs:
            logger.info("  No changes in %s", project)
            project_changed[project] = False
            continue

        any_changes = True
        project_changed[project] = True
        total_new += len(new_accs)
        total_mod += len(mod_accs)
        logger.info("  %d new, %d modified alleles in %s", len(new_accs), len(mod_accs), project)

        # Determine affected species/locus groups
        changed_set = set(new_accs + mod_accs)
        groups = group_by_species_locus(listing)
        affected = {
            k: v for k, v in groups.items()
            if any(acc in changed_set for acc in v)
        }

        for (species, locus), accessions in affected.items():
            logger.info("  Rebuilding %s/%s (%d alleles)", species, locus, len(accessions))
            records = fetch_full_records(client, project, accessions)

            data_dir = repo_root / DATA_DIR / project / species
            write_genbank_file(records, data_dir / f"{locus}.gb", project)

    # Update version info (with revision suffix for silent updates)
    version_data = load_version_file(repo_root)
    for project in ("MHC", "NHKIR"):
        try:
            release = scrape_current_release(project)
            ipd_version = release.get("version", "unknown")
            stored = version_data.get(project, {})
            changed = project_changed.get(project, False)

            new_version = resolve_version(ipd_version, stored, changed)

            # Compute revision number from the resolved version
            stored_ipd = stored.get("ipd_version", stored.get("version", ""))
            if ipd_version != stored_ipd:
                revision = 0
            elif changed:
                revision = stored.get("revision", 0) + 1
            else:
                revision = stored.get("revision", 0)

            version_data[project] = {
                "ipd_version": ipd_version,
                "version": new_version,
                "revision": revision,
                "allele_count": len(all_listings.get(project, [])),
            }
            logger.info("  %s version: %s (IPD %s, revision %d)",
                        project, new_version, ipd_version, revision)
        except Exception as e:
            logger.warning("Could not update version info for %s: %s", project, e)

    save_version_file(repo_root, version_data)
    run_build_index(repo_root, all_listings)

    if any_changes:
        logger.info("=== Incremental update complete: %d new, %d modified ===", total_new, total_mod)
    else:
        logger.info("=== No changes detected ===")


def run_build_index(
    repo_root: Path,
    listings: dict[str, list] | None = None,
) -> None:
    """Build docs/alleles.json from allele listings."""
    if listings is None:
        client = IPDClient(cache_dir=repo_root / CACHE_DIR)
        listings = {
            "MHC": fetch_allele_listing(client, "MHC"),
            "NHKIR": fetch_allele_listing(client, "NHKIR"),
        }

    version_data = load_version_file(repo_root)
    index_path = repo_root / DOCS_DIR / "alleles.json"
    index = build_metadata_index(
        listings.get("MHC", []),
        listings.get("NHKIR", []),
        index_path,
        mhc_version=version_data.get("MHC", {}).get("version", ""),
        nhkir_version=version_data.get("NHKIR", {}).get("version", ""),
    )
    logger.info("Built metadata index: %d alleles at %s", len(index["alleles"]), index_path)
