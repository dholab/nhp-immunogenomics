"""Build the JSON metadata index consumed by the GitHub Pages site."""

import json
import logging
from datetime import datetime, timezone
from pathlib import Path

from Bio import SeqIO

logger = logging.getLogger(__name__)


def scan_genbank_seq_info(data_dir: Path) -> dict[str, dict]:
    """Scan all GenBank files to extract sequence length and type per accession.

    Returns {accession: {"len": int, "st": str}} where st is
    "genomic", "coding", or "partial".
    """
    info: dict[str, dict] = {}

    for project_dir in ("MHC", "NHKIR"):
        project_path = data_dir / project_dir
        if not project_path.exists():
            continue

        for gb_file in project_path.rglob("*.gb"):
            try:
                for record in SeqIO.parse(str(gb_file), "genbank"):
                    acc = record.id
                    seq_len = len(record.seq)

                    # Determine sequence type from features
                    mol_type = ""
                    has_introns = False
                    has_cds_join = False

                    for feat in record.features:
                        if feat.type == "source":
                            mol_type = feat.qualifiers.get("mol_type", [""])[0]
                        elif feat.type == "intron":
                            has_introns = True
                        elif feat.type == "CDS":
                            if hasattr(feat.location, "parts") and len(feat.location.parts) > 1:
                                has_cds_join = True

                    if has_introns or (has_cds_join and "genomic" in mol_type.lower()):
                        seq_type = "genomic"
                    elif "mRNA" in mol_type or "cDNA" in mol_type:
                        seq_type = "coding"
                    elif "genomic" in mol_type.lower():
                        seq_type = "genomic"
                    else:
                        seq_type = "coding"

                    info[acc] = {"len": seq_len, "st": seq_type}
            except Exception as e:
                logger.warning("Error scanning %s: %s", gb_file, e)

    return info


def build_metadata_index(
    mhc_listing: list[dict],
    nhkir_listing: list[dict],
    output_path: Path,
    mhc_version: str = "",
    nhkir_version: str = "",
    provisional_alleles: list[dict] | None = None,
    seq_info: dict[str, dict] | None = None,
) -> dict:
    """
    Build docs/alleles.json from allele listings.

    Listings use flat dot-notation keys from the API (e.g., 'organism.name').
    Provisional alleles (if provided) are appended with their pre-built metadata.
    seq_info maps accession -> {"len": int, "st": str} from GenBank scanning.
    """
    species_map = {}
    all_alleles = []
    loci: dict[str, set] = {"MHC": set(), "NHKIR": set()}
    si = seq_info or {}

    for allele in mhc_listing:
        _process_allele(allele, "MHC", species_map, all_alleles, loci, si)

    for allele in nhkir_listing:
        _process_allele(allele, "NHKIR", species_map, all_alleles, loci, si)

    # Append provisional alleles (already in index-ready format)
    if provisional_alleles:
        for prov in provisional_alleles:
            all_alleles.append(prov)
            locus = prov.get("l", "")
            if locus:
                # Infer database from locus for loci tracking
                db = "NHKIR" if locus.startswith("KIR") else "MHC"
                loci[db].add(locus)

    index = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "mhc_version": mhc_version,
        "nhkir_version": nhkir_version,
        "species": species_map,
        "alleles": all_alleles,
        "loci": {k: sorted(v) for k, v in loci.items()},
    }

    if provisional_alleles:
        index["provisional_count"] = len(provisional_alleles)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(index, separators=(",", ":")))

    return index


def _process_allele(
    allele: dict,
    project: str,
    species_map: dict,
    all_alleles: list,
    loci: dict[str, set],
    seq_info: dict[str, dict] | None = None,
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

    acc = allele.get("accession", "")
    entry = {
        "a": acc,
        "n": allele.get("name", ""),
        "l": locus,
        "c": allele.get("class", ""),
        "s": sp,
        "p": project,
        "da": allele.get("date_assigned", ""),
        "dm": allele.get("date_modified", ""),
    }

    # Add sequence length and type from GenBank scan
    si = seq_info or {}
    if acc in si:
        entry["len"] = si[acc]["len"]
        entry["st"] = si[acc]["st"]

    prev = allele.get("previous", [])
    if prev:
        entry["prev"] = prev
    all_alleles.append(entry)
