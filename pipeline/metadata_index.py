"""Build the JSON metadata index consumed by the GitHub Pages site."""

import json
from datetime import datetime, timezone
from pathlib import Path


def build_metadata_index(
    mhc_listing: list[dict],
    nhkir_listing: list[dict],
    output_path: Path,
    mhc_version: str = "",
    nhkir_version: str = "",
) -> dict:
    """
    Build docs/alleles.json from allele listings.

    Listings use flat dot-notation keys from the API (e.g., 'organism.name').
    """
    species_map = {}
    all_alleles = []
    loci: dict[str, set] = {"MHC": set(), "NHKIR": set()}

    for allele in mhc_listing:
        _process_allele(allele, "MHC", species_map, all_alleles, loci)

    for allele in nhkir_listing:
        _process_allele(allele, "NHKIR", species_map, all_alleles, loci)

    index = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "mhc_version": mhc_version,
        "nhkir_version": nhkir_version,
        "species": species_map,
        "alleles": all_alleles,
        "loci": {k: sorted(v) for k, v in loci.items()},
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(index, separators=(",", ":")))

    return index


def _process_allele(
    allele: dict,
    project: str,
    species_map: dict,
    all_alleles: list,
    loci: dict[str, set],
) -> None:
    """Process a single allele listing record into the index."""
    sp = allele.get("organism.name", "")
    sci = allele.get("organism.scientificName", "")
    common = allele.get("organism.commonName", "")
    taxon = allele.get("organism.taxon")

    if sp and sp not in species_map:
        species_map[sp] = {
            "scientificName": sci,
            "commonName": common,
            "taxon": taxon,
        }

    locus = allele.get("locus", "")
    if locus:
        loci[project].add(locus)

    entry = {
        "a": allele.get("accession", ""),
        "n": allele.get("name", ""),
        "l": locus,
        "c": allele.get("class", ""),
        "s": sp,
        "p": project,
        "da": allele.get("date_assigned", ""),
        "dm": allele.get("date_modified", ""),
    }
    prev = allele.get("previous", [])
    if prev:
        entry["prev"] = prev
    all_alleles.append(entry)
