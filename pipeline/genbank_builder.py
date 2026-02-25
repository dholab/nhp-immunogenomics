"""Convert IPD API JSON allele records to GenBank format."""

import logging
from datetime import datetime
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, Reference

logger = logging.getLogger(__name__)


def allele_to_seqrecord(record: dict, project: str) -> SeqRecord:
    """
    Convert a single IPD API allele JSON to a BioPython SeqRecord.

    Uses the genomic sequence as the primary sequence when it differs from
    coding (indicating introns are present). Falls back to coding sequence.
    """
    seq_data = record.get("sequence", {})
    genomic_seq = seq_data.get("genomic", "")
    coding_seq = seq_data.get("coding", "")
    protein_seq = seq_data.get("protein", "")
    codon_start = seq_data.get("codon_start", 1)

    use_genomic = bool(genomic_seq and coding_seq and genomic_seq != coding_seq)
    primary_seq = genomic_seq if use_genomic else coding_seq

    if not primary_seq:
        raise ValueError(f"No sequence data for {record.get('accession', '?')}")

    accession = record["accession"]
    allele_name = record.get("name", accession)
    locus_name = record.get("locus", "")
    allele_class = record.get("class", "")
    date_mod = record.get("date_modified", "")
    previous = record.get("previous", [])

    organism = record.get("organism", {})
    sci_name = organism.get("scientificName", "unknown organism")
    common_name = organism.get("commonName", "")
    lineage = organism.get("lineage", "")
    taxon = organism.get("taxon")
    species_prefix = organism.get("name", "")

    sr = SeqRecord(
        seq=Seq(primary_seq),
        id=accession,
        name=accession,
        description=f"{allele_name}, {locus_name} locus allele",
    )

    # Annotations
    sr.annotations["molecule_type"] = "DNA"
    sr.annotations["topology"] = "linear"
    sr.annotations["data_file_division"] = "PRI"

    if date_mod:
        try:
            dt = datetime.strptime(date_mod, "%Y-%m-%d")
            sr.annotations["date"] = dt.strftime("%d-%b-%Y").upper()
        except ValueError:
            pass

    sr.annotations["organism"] = sci_name
    sr.annotations["source"] = (
        f"{sci_name} ({common_name})" if common_name else sci_name
    )
    sr.annotations["taxonomy"] = [
        t.strip() for t in lineage.rstrip(";").split(";") if t.strip()
    ]
    sr.annotations["accessions"] = [accession]
    sr.annotations["sequence_version"] = 1

    # Keywords
    keywords = []
    if project == "MHC":
        keywords.append("IPD-MHC")
    else:
        keywords.append("IPD-NHKIR")
    if allele_class and allele_class != "unknown":
        keywords.append(f"MHC class {allele_class}")
    keywords.append(locus_name)
    sr.annotations["keywords"] = [k for k in keywords if k]

    # Cross-references
    insdc = record.get("insdc", [])
    sr.dbxrefs = [f"INSDC:{acc}" for acc in insdc]

    # Comment with previous designations
    comment_lines = [
        f"IPD-{'MHC' if project == 'MHC' else 'NHKIR'} allele: {allele_name}",
        f"IPD accession: {accession}",
        f"Species prefix: {species_prefix}",
    ]
    if allele_class and allele_class != "unknown":
        comment_lines.append(f"MHC class: {allele_class}")
    if previous:
        comment_lines.append(f"Previous designations: {'; '.join(previous)}")
    sr.annotations["comment"] = "\n".join(comment_lines)

    # References
    bio_refs = []
    for ref_data in record.get("reference", []):
        bio_ref = Reference()
        bio_ref.authors = ref_data.get("authors", "")
        bio_ref.title = ref_data.get("title", "")
        journal = ref_data.get("journal", "")
        if ref_data.get("volume"):
            journal += f" {ref_data['volume']}"
        if ref_data.get("year"):
            journal += f" ({ref_data['year']})"
        bio_ref.journal = journal
        if ref_data.get("pmid"):
            bio_ref.pubmed_id = str(ref_data["pmid"])
        bio_refs.append(bio_ref)
    sr.annotations["references"] = bio_refs

    # Features
    features = record.get("feature", [])
    exons = sorted(
        [f for f in features if f.get("name") == "exon"],
        key=lambda x: x["start"],
    )
    introns = sorted(
        [f for f in features if f.get("name") == "intron"],
        key=lambda x: x["start"],
    )

    seq_len = len(primary_seq)

    # Source feature
    source_qualifiers = {
        "organism": [sci_name],
        "mol_type": ["genomic DNA" if use_genomic else "mRNA"],
    }
    if taxon:
        source_qualifiers["db_xref"] = [f"taxon:{taxon}"]
    sr.features.append(SeqFeature(
        location=FeatureLocation(0, seq_len),
        type="source",
        qualifiers=source_qualifiers,
    ))

    # Gene feature
    gene_qualifiers = {
        "gene": [locus_name],
        "allele": [allele_name],
    }
    if previous:
        gene_qualifiers["note"] = [
            f"previous designations: {'; '.join(previous)}"
        ]
    sr.features.append(SeqFeature(
        location=FeatureLocation(0, seq_len),
        type="gene",
        qualifiers=gene_qualifiers,
    ))

    # Exon features
    for i, exon in enumerate(exons, 1):
        start = exon["start"]
        size = exon["size"]
        sr.features.append(SeqFeature(
            location=FeatureLocation(start, start + size),
            type="exon",
            qualifiers={"gene": [locus_name], "number": [str(i)]},
        ))

    # Intron features (only when genomic sequence)
    if use_genomic:
        for i, intron in enumerate(introns, 1):
            start = intron["start"]
            size = intron["size"]
            sr.features.append(SeqFeature(
                location=FeatureLocation(start, start + size),
                type="intron",
                qualifiers={"gene": [locus_name], "number": [str(i)]},
            ))

    # CDS feature
    if exons and coding_seq:
        if use_genomic:
            exon_locations = [
                FeatureLocation(e["start"], e["start"] + e["size"])
                for e in exons
            ]
            cds_location = (
                exon_locations[0] if len(exon_locations) == 1
                else CompoundLocation(exon_locations)
            )
        else:
            cds_location = FeatureLocation(0, seq_len)

        cds_qualifiers = {
            "gene": [locus_name],
            "allele": [allele_name],
            "codon_start": [str(codon_start)],
            "product": [_infer_product(locus_name, allele_class, project)],
        }
        if protein_seq:
            cds_qualifiers["translation"] = [protein_seq]

        sr.features.append(SeqFeature(
            location=cds_location,
            type="CDS",
            qualifiers=cds_qualifiers,
        ))

    return sr


def _infer_product(locus: str, allele_class: str, project: str) -> str:
    """Infer the gene product description from locus and class."""
    if project == "NHKIR":
        return f"killer cell immunoglobulin-like receptor {locus}"
    if allele_class == "I":
        return f"MHC class I {locus} antigen"
    elif allele_class == "II":
        return f"MHC class II {locus} antigen"
    return f"MHC {locus} antigen"


def provisional_to_seqrecord(
    name: str,
    sequence: str,
    species: str,
    locus: str,
    allele_class: str,
    seq_type: str,
    species_info: dict,
    date_added: str,
    submitter: str,
    accession: str = "",
) -> SeqRecord:
    """Convert a provisional allele to a BioPython SeqRecord.

    Simpler than allele_to_seqrecord: no exon/intron features from API,
    no references, no INSDC cross-references.

    Args:
        name: Provisional allele name (e.g., "Mamu-A1*026:new01").
        sequence: Uppercase nucleotide sequence.
        species: Species prefix (e.g., "Mamu").
        locus: Locus name (e.g., "A1").
        allele_class: MHC class ("I", "II", or "").
        seq_type: "coding" or "genomic".
        species_info: Dict with scientificName, commonName, taxon.
        date_added: Date string (YYYY-MM-DD).
        submitter: Name of the person who submitted the allele.
        accession: Provisional accession (PROVxxxxx).
    """
    sci_name = species_info.get("scientificName", "unknown organism")
    common_name = species_info.get("commonName", "")
    taxon = species_info.get("taxon")
    is_genomic = seq_type == "genomic"

    sr = SeqRecord(
        seq=Seq(sequence),
        id=accession,
        name=accession,
        description=f"{name}, {locus} locus allele (provisional)",
    )

    sr.annotations["molecule_type"] = "DNA"
    sr.annotations["topology"] = "linear"
    sr.annotations["data_file_division"] = "PRI"
    sr.annotations["organism"] = sci_name
    sr.annotations["source"] = (
        f"{sci_name} ({common_name})" if common_name else sci_name
    )
    sr.annotations["accessions"] = [accession]
    sr.annotations["sequence_version"] = 1

    if date_added:
        try:
            dt = datetime.strptime(date_added, "%Y-%m-%d")
            sr.annotations["date"] = dt.strftime("%d-%b-%Y").upper()
        except ValueError:
            pass

    # Keywords
    keywords = ["Provisional"]
    if allele_class and allele_class != "unknown":
        keywords.append(f"MHC class {allele_class}")
    keywords.append(locus)
    sr.annotations["keywords"] = keywords

    # Comment
    sr.annotations["comment"] = (
        f"Provisional allele submitted by {submitter}. Not yet in IPD.\n"
        f"Species prefix: {species}\n"
        f"Sequence type: {seq_type}"
    )

    seq_len = len(sequence)

    # Source feature
    source_qualifiers = {
        "organism": [sci_name],
        "mol_type": ["genomic DNA" if is_genomic else "mRNA"],
    }
    if taxon:
        source_qualifiers["db_xref"] = [f"taxon:{taxon}"]
    sr.features.append(SeqFeature(
        location=FeatureLocation(0, seq_len),
        type="source",
        qualifiers=source_qualifiers,
    ))

    # Gene feature
    sr.features.append(SeqFeature(
        location=FeatureLocation(0, seq_len),
        type="gene",
        qualifiers={
            "gene": [locus],
            "allele": [name],
            "note": ["provisional allele, not yet in IPD"],
        },
    ))

    # CDS feature (for coding sequences)
    if not is_genomic:
        product = _infer_product(
            locus, allele_class,
            "NHKIR" if locus.startswith("KIR") else "MHC",
        )
        sr.features.append(SeqFeature(
            location=FeatureLocation(0, seq_len),
            type="CDS",
            qualifiers={
                "gene": [locus],
                "allele": [name],
                "codon_start": ["1"],
                "product": [product],
            },
        ))

    return sr


def write_genbank_file(
    records: list[dict],
    output_path: Path,
    project: str,
) -> int:
    """
    Convert a list of API allele records to a multi-entry GenBank file.

    Returns the number of records written.
    """
    seq_records = []
    for rec in records:
        try:
            sr = allele_to_seqrecord(rec, project)
            seq_records.append(sr)
        except Exception as e:
            logger.warning("Skipping %s: %s", rec.get("accession", "?"), e)

    if seq_records:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as fh:
            SeqIO.write(seq_records, fh, "genbank")

    return len(seq_records)
