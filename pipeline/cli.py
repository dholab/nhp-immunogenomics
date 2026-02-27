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
from .metadata_index import build_metadata_index, scan_genbank_seq_info
from .provisional import (
    add_provisional_alleles,
    backfill_rationale,
    build_genbank,
    build_metadata,
    build_sequence_index,
    check_retirements,
    load_manifest,
    load_sequences,
    retire_alleles,
    validate,
    validate_no_ipd_collisions,
)
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

    sub.add_parser("provisional-build", help="Build provisional GenBank files and update metadata index")
    pv_cmd = sub.add_parser("provisional-validate", help="Validate provisional alleles")

    pa_cmd = sub.add_parser("provisional-add", help="Add provisional allele(s)")
    pa_cmd.add_argument("--species", required=True, help="Species prefix (e.g., Mamu)")
    pa_cmd.add_argument("--locus", required=True, help="Locus name (e.g., A1)")
    pa_cmd.add_argument("--class", dest="allele_class", default="", help="MHC class (I, II)")
    pa_cmd.add_argument("--submitter", required=True, help="Name of the submitter")
    pa_cmd.add_argument("--fasta", type=Path, help="Path to FASTA file (one or more sequences)")
    pa_cmd.add_argument("--fasta-text", default="", help="FASTA content as string (alternative to --fasta)")
    pa_cmd.add_argument("--name", default="", help="Override auto-assigned name (single sequence only)")
    pa_cmd.add_argument("--seq-type", default="coding", choices=["coding", "genomic"],
                        help="Sequence type (default: coding)")
    pa_cmd.add_argument("--notes", default="", help="Optional notes")

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
    elif args.command == "provisional-build":
        _process_provisionals(repo_root)
        run_build_index(repo_root)
    elif args.command == "provisional-validate":
        run_provisional_validate(repo_root)
    elif args.command == "provisional-add":
        run_provisional_add(repo_root, args)
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

    # Process provisional alleles
    _process_provisionals(repo_root)

    # Build metadata index (includes provisional alleles)
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

    # Process provisional alleles + auto-retirement
    _process_provisionals(repo_root, run_retirement=True)

    # Build metadata index (includes provisional alleles)
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

    # Scan GenBank files for sequence length and type
    data_dir = repo_root / DATA_DIR
    seq_info = scan_genbank_seq_info(data_dir)
    logger.info("Scanned sequence info for %d alleles from GenBank files", len(seq_info))

    # Load provisional alleles for inclusion in the index
    prov_alleles = _get_provisional_metadata(repo_root)

    index_path = repo_root / DOCS_DIR / "alleles.json"
    index = build_metadata_index(
        listings.get("MHC", []),
        listings.get("NHKIR", []),
        index_path,
        mhc_version=version_data.get("MHC", {}).get("version", ""),
        nhkir_version=version_data.get("NHKIR", {}).get("version", ""),
        provisional_alleles=prov_alleles if prov_alleles else None,
        seq_info=seq_info,
    )
    logger.info("Built metadata index: %d alleles at %s", len(index["alleles"]), index_path)


def _get_provisional_metadata(repo_root: Path) -> list[dict]:
    """Load provisional manifest and build metadata records for the index."""
    manifest = load_manifest(repo_root)
    if not manifest:
        return []
    sequences = load_sequences(repo_root, manifest)
    return build_metadata(manifest, sequences)


def _get_species_map(repo_root: Path) -> dict[str, dict]:
    """Build a species map from the existing alleles.json index."""
    index_path = repo_root / DOCS_DIR / "alleles.json"
    if index_path.exists():
        index = json.loads(index_path.read_text())
        return index.get("species", {})
    return {}


def _process_provisionals(repo_root: Path, run_retirement: bool = False) -> None:
    """Build provisional GenBank files and optionally run auto-retirement."""
    manifest = load_manifest(repo_root)
    if not manifest:
        logger.info("No provisional alleles found.")
        return

    sequences = load_sequences(repo_root, manifest)

    # Auto-retirement (during updates only)
    if run_retirement:
        logger.info("=== Checking provisional allele retirements ===")
        ipd_index = build_sequence_index(repo_root)
        retirements = check_retirements(manifest, sequences, ipd_index)
        if retirements:
            retire_alleles(repo_root, retirements)
            # Reload after retirement
            manifest = load_manifest(repo_root)
            sequences = load_sequences(repo_root, manifest)
        else:
            logger.info("  No provisional alleles to retire.")

    if not manifest:
        return

    # Backfill rationale for entries that lack it
    backfill_rationale(repo_root)
    # Reload in case backfill updated the manifest
    manifest = load_manifest(repo_root)
    sequences = load_sequences(repo_root, manifest)

    # Build provisional GenBank files
    species_map = _get_species_map(repo_root)
    n = build_genbank(repo_root, manifest, sequences, species_map)
    logger.info("Built %d provisional GenBank records", n)


def run_provisional_validate(repo_root: Path) -> None:
    """Validate provisional manifest and sequences."""
    errors = validate(repo_root)

    # Also check for IPD name collisions if we have an existing index
    manifest = load_manifest(repo_root)
    if manifest:
        index_path = repo_root / DOCS_DIR / "alleles.json"
        if index_path.exists():
            index = json.loads(index_path.read_text())
            ipd_names = {
                a["n"] for a in index.get("alleles", [])
                if not a.get("prov")
            }
            errors.extend(validate_no_ipd_collisions(manifest, ipd_names))

    if errors:
        for err in errors:
            logger.error("  %s", err)
        logger.error("Validation failed with %d error(s).", len(errors))
        sys.exit(1)
    else:
        logger.info("Provisional alleles are valid.")


def run_provisional_add(repo_root: Path, args) -> None:
    """Add provisional allele(s) from a FASTA file or string."""
    from io import StringIO

    from Bio import SeqIO

    # Parse records from --fasta file or --fasta-text string
    if args.fasta_text:
        records = list(SeqIO.parse(StringIO(args.fasta_text), "fasta"))
    elif args.fasta:
        if not args.fasta.exists():
            logger.error("FASTA file not found: %s", args.fasta)
            sys.exit(1)
        records = list(SeqIO.parse(str(args.fasta), "fasta"))
    else:
        logger.error("Must provide either --fasta or --fasta-text")
        sys.exit(1)

    if not records:
        logger.error("No sequences found in FASTA input")
        sys.exit(1)

    try:
        names = add_provisional_alleles(
            repo_root=repo_root,
            records=records,
            species=args.species,
            locus=args.locus,
            allele_class=args.allele_class,
            seq_type=args.seq_type,
            submitter=args.submitter,
            notes=args.notes,
            name_override=args.name,
        )
    except ValueError as exc:
        logger.error("%s", exc)
        sys.exit(1)

    rel_fasta = f"{args.species}/{args.locus}.fasta"
    for name in names:
        logger.info("  Manifest: provisional/manifest.tsv")
        logger.info("  Sequence: provisional/sequences/%s", rel_fasta)
