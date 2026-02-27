"""Manage provisional (lab-submitted) alleles not yet in IPD."""

import csv
import hashlib
import json
import logging
import re
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner

from .config import (
    DATA_DIR,
    PROVISIONAL_DATA_DIR,
    PROVISIONAL_MANIFEST,
    PROVISIONAL_RETIRED,
    PROVISIONAL_SEQ_INDEX,
    PROVISIONAL_SEQUENCES_DIR,
)

logger = logging.getLogger(__name__)

MANIFEST_COLUMNS = (
    "name", "species", "locus", "class", "seq_type",
    "sequence_file", "submitter", "date_added", "notes",
    "ref_nt", "ref_nt_pct", "ref_aa", "ref_aa_pct",
    "cds_diffs", "aa_diffs", "fasta_header",
)


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load_manifest(repo_root: Path) -> list[dict]:
    """Load provisional/manifest.tsv into a list of dicts.

    Returns an empty list if the manifest is empty or missing.
    """
    manifest_path = repo_root / PROVISIONAL_MANIFEST
    if not manifest_path.exists():
        return []

    rows = []
    with open(manifest_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            # Skip completely blank rows
            if not any(row.values()):
                continue
            rows.append(row)
    return rows


def load_sequences(
    repo_root: Path,
    manifest: list[dict],
) -> dict[str, str]:
    """Load all provisional FASTA sequences referenced by the manifest.

    Returns {allele_name: uppercase_nucleotide_sequence}.
    """
    sequences: dict[str, str] = {}
    seq_dir = repo_root / PROVISIONAL_SEQUENCES_DIR

    # Collect all FASTA files referenced by manifest entries
    files_to_load: dict[str, list[str]] = {}
    for entry in manifest:
        sf = entry.get("sequence_file", "")
        name = entry.get("name", "")
        if sf and name:
            files_to_load.setdefault(sf, []).append(name)

    for rel_path, expected_names in files_to_load.items():
        fasta_path = seq_dir / rel_path
        if not fasta_path.exists():
            logger.warning("Missing FASTA file: %s", fasta_path)
            continue

        # Parse all records from the FASTA file
        records_in_file: dict[str, str] = {}
        for rec in SeqIO.parse(str(fasta_path), "fasta"):
            records_in_file[rec.id] = str(rec.seq).upper()

        for name in expected_names:
            if name in records_in_file:
                sequences[name] = records_in_file[name]
            else:
                logger.warning(
                    "Allele %s not found in %s (records: %s)",
                    name, fasta_path, list(records_in_file.keys()),
                )

    return sequences


# ---------------------------------------------------------------------------
# Accessions
# ---------------------------------------------------------------------------

def compute_accession(sequence: str) -> str:
    """Compute a deterministic provisional accession from a nucleotide sequence.

    Format: PROV + first 5 hex chars of SHA-256 of uppercase sequence.
    """
    h = hashlib.sha256(sequence.upper().encode()).hexdigest()
    return f"PROV{h[:5]}"


def sequence_hash(sequence: str, length: int = 12) -> str:
    """Compute a SHA-256 hash prefix for sequence matching."""
    return hashlib.sha256(sequence.upper().encode()).hexdigest()[:length]


# ---------------------------------------------------------------------------
# Nearest-neighbor naming
# ---------------------------------------------------------------------------

def find_nearest_neighbor(
    sequence: str,
    locus_gb_path: Path,
) -> tuple[str, float]:
    """Find the closest IPD allele by pairwise alignment.

    Args:
        sequence: Uppercase nucleotide sequence of the provisional allele.
        locus_gb_path: Path to the GenBank file for the relevant locus
            (e.g., data/MHC/Mamu/A1.gb).

    Returns:
        (best_allele_name, percent_identity). Returns ("", 0.0) if no
        alleles found at the locus.
    """
    if not locus_gb_path.exists():
        return ("", 0.0)

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    best_name = ""
    best_pct = 0.0
    query = sequence.upper()
    query_len = len(query)

    for record in SeqIO.parse(str(locus_gb_path), "genbank"):
        target = str(record.seq).upper()
        if not target:
            continue

        score = aligner.score(query, target)
        # Percent identity approximation: score / max possible score
        max_len = max(query_len, len(target))
        pct = (score / max_len) * 100 if max_len > 0 else 0.0

        allele_name = record.description.split(",")[0].strip()
        if pct > best_pct:
            best_pct = pct
            best_name = allele_name

    return (best_name, best_pct)


def extract_allele_group(allele_name: str) -> str:
    """Extract the allele group number from an IPD allele name.

    E.g., "Mamu-A1*026:01:01:01" → "026"
          "Mamu-B*017:01" → "017"
          "Mafa-KIR3DL01*001:01" → "001"
    """
    match = re.search(r"\*(\d+)", allele_name)
    return match.group(1) if match else ""


def assign_name(
    species: str,
    locus: str,
    nearest_name: str,
    manifest: list[dict],
) -> str:
    """Assign a provisional allele name based on the nearest IPD neighbor.

    Format: {Species}-{Locus}*{Group}:new{NN}

    Increments NN to avoid collisions with existing provisional names.
    """
    group = extract_allele_group(nearest_name)
    if not group:
        group = "000"

    # Find existing provisional names in the same group
    prefix = f"{species}-{locus}*{group}:new"
    existing_nums = []
    for entry in manifest:
        name = entry.get("name", "")
        if name.startswith(prefix):
            suffix = name[len(prefix):]
            if suffix.isdigit():
                existing_nums.append(int(suffix))

    next_num = max(existing_nums, default=0) + 1
    return f"{prefix}{next_num:02d}"


_NUCLEOTIDE_RE = re.compile(r"^[ACGTNRYSWKMBDHVacgtnryswkmbdhv]+$")


def add_provisional_alleles(
    repo_root: Path,
    records: list,
    species: str,
    locus: str,
    allele_class: str,
    seq_type: str,
    submitter: str,
    notes: str = "",
    name_override: str = "",
) -> list[str]:
    """Add one or more provisional alleles from SeqRecords.

    Args:
        repo_root: Repository root directory.
        records: List of Bio.SeqRecord objects (from FASTA parsing).
        species: Species prefix (e.g., "Mamu").
        locus: Locus name (e.g., "A1").
        allele_class: MHC class ("I", "II", or "").
        seq_type: "coding" or "genomic".
        submitter: Name of the submitter.
        notes: Optional notes.
        name_override: If set, use this name instead of auto-assigning.
            Only valid for single-record input.

    Returns:
        List of assigned allele names.

    Raises:
        ValueError: On validation errors (empty sequences, duplicates, etc.).
    """
    from datetime import date

    if not records:
        raise ValueError("No sequences provided")

    if name_override and len(records) > 1:
        raise ValueError(
            f"--name cannot be used with multiple sequences ({len(records)} records)"
        )

    # Pre-validate all sequences before writing anything
    # (record_id, uppercase_seq, fasta_header)
    sequences_to_add: list[tuple[str, str, str]] = []
    batch_hashes: dict[str, str] = {}  # hash -> record_id

    for i, record in enumerate(records):
        seq = str(record.seq).upper()
        if not seq:
            raise ValueError(f"Record {i + 1} ('{record.id}') has an empty sequence")
        if not _NUCLEOTIDE_RE.match(seq):
            raise ValueError(
                f"Record {i + 1} ('{record.id}') contains invalid nucleotide characters"
            )

        h = sequence_hash(seq)
        if h in batch_hashes:
            raise ValueError(
                f"Record {i + 1} ('{record.id}') has the same sequence as "
                f"record '{batch_hashes[h]}'"
            )
        batch_hashes[h] = record.id
        sequences_to_add.append((record.id, seq, record.description))

    # Check for duplicates against existing manifest
    manifest = load_manifest(repo_root)
    existing_seqs = load_sequences(repo_root, manifest)
    existing_hashes = {sequence_hash(s): n for n, s in existing_seqs.items()}

    for rec_id, seq, _hdr in sequences_to_add:
        h = sequence_hash(seq)
        if h in existing_hashes:
            raise ValueError(
                f"Record '{rec_id}' has the same sequence as existing "
                f"provisional allele '{existing_hashes[h]}'"
            )

    # All validation passed — now process each record
    from .namer import ReferenceAllele, name_provisional_allele

    assigned_names: list[str] = []
    rel_fasta = f"{species}/{locus}.fasta"
    today = date.today().isoformat()

    # Resolve GenBank path for reference alleles
    gb_path = repo_root / DATA_DIR / "MHC" / species / f"{locus}.gb"
    if not gb_path.exists():
        gb_path = repo_root / DATA_DIR / "NHKIR" / species / f"{locus}.gb"

    # Collect existing provisional names at this locus for collision avoidance
    existing_names = [
        e.get("name", "")
        for e in manifest
        if e.get("species") == species and e.get("locus") == locus
    ]

    # Track previously-named provisional alleles in this batch for
    # intra-batch synonymous detection
    provisional_refs: list[ReferenceAllele] = []

    for i, (rec_id, seq, fasta_header) in enumerate(sequences_to_add):
        if name_override:
            allele_name = name_override
            relationship = ""
        else:
            result = name_provisional_allele(
                sequence=seq,
                species=species,
                locus=locus,
                seq_type=seq_type,
                locus_gb_path=gb_path,
                existing_names=existing_names,
                provisional_refs=provisional_refs,
            )
            allele_name = result.provisional_name
            relationship = result.relationship

            if result.closest_nt_allele:
                logger.info(
                    "Record %d ('%s'): nearest IPD allele: %s (%.1f%% nt identity), "
                    "relationship: %s, closest protein: %s (%.1f%%)",
                    i + 1, rec_id, result.closest_nt_allele, result.closest_nt_identity,
                    result.relationship,
                    result.closest_protein_allele or "none",
                    result.closest_protein_identity,
                )
            else:
                logger.warning(
                    "Record %d ('%s'): no IPD alleles found at %s/%s",
                    i + 1, rec_id, species, locus,
                )

            for warn in result.warnings:
                logger.warning("Record %d ('%s'): %s", i + 1, rec_id, warn)

            # Track name for next record's collision avoidance
            existing_names.append(allele_name)

            # Track protein for intra-batch synonymous detection
            if result.extracted_protein:
                provisional_refs.append(
                    ReferenceAllele(
                        name=allele_name,
                        accession="",
                        sequence=seq,
                        mol_type="genomic DNA" if seq_type == "genomic" else "mRNA",
                        protein=result.extracted_protein,
                        coding_seq=result.extracted_cds,
                        exon_coords=[],
                        intron_coords=[],
                        codon_start=1,
                    )
                )

        # Build notes: include relationship and any user-provided notes
        entry_notes_parts = []
        if relationship:
            entry_notes_parts.append(relationship)
        if notes:
            entry_notes_parts.append(notes)
        entry_notes = "; ".join(entry_notes_parts)

        # Extract rationale fields from NamingResult
        ref_nt = ""
        ref_nt_pct = ""
        ref_aa = ""
        ref_aa_pct = ""
        cds_diffs = ""
        aa_diffs = ""
        if not name_override and result:
            ref_nt = result.closest_nt_allele
            ref_nt_pct = f"{result.closest_nt_identity:.1f}" if result.closest_nt_allele else ""
            ref_aa = result.closest_protein_allele
            ref_aa_pct = f"{result.closest_protein_identity:.1f}" if result.closest_protein_allele else ""
            if result.cds_diffs >= 0:
                cds_diffs = str(result.cds_diffs)
            if result.aa_diffs >= 0:
                aa_diffs = str(result.aa_diffs)

        entry = {
            "name": allele_name,
            "species": species,
            "locus": locus,
            "class": allele_class,
            "seq_type": seq_type,
            "sequence_file": rel_fasta,
            "submitter": submitter,
            "date_added": today,
            "notes": entry_notes,
            "ref_nt": ref_nt,
            "ref_nt_pct": ref_nt_pct,
            "ref_aa": ref_aa,
            "ref_aa_pct": ref_aa_pct,
            "cds_diffs": cds_diffs,
            "aa_diffs": aa_diffs,
            "fasta_header": fasta_header,
        }

        append_to_manifest(repo_root, entry)
        append_to_fasta(repo_root, rel_fasta, allele_name, seq)

        # Track in-memory so next iteration sees this entry
        manifest.append(entry)
        assigned_names.append(allele_name)

        logger.info("Added provisional allele %s", allele_name)

    logger.info("Added %d provisional allele(s)", len(assigned_names))
    return assigned_names


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate(repo_root: Path) -> list[str]:
    """Validate the provisional manifest and sequences.

    Returns a list of error messages. Empty list means valid.
    """
    errors: list[str] = []
    manifest = load_manifest(repo_root)

    if not manifest:
        return errors  # Empty manifest is valid

    # Check required columns
    for i, entry in enumerate(manifest, 1):
        for col in ("name", "species", "locus", "class", "seq_type", "sequence_file", "submitter", "date_added"):
            if not entry.get(col, "").strip():
                errors.append(f"Row {i}: missing required field '{col}'")

        # Validate seq_type
        st = entry.get("seq_type", "").strip()
        if st and st not in ("coding", "genomic"):
            errors.append(f"Row {i}: seq_type must be 'coding' or 'genomic', got '{st}'")

        # Validate date format
        da = entry.get("date_added", "").strip()
        if da and not re.match(r"\d{4}-\d{2}-\d{2}$", da):
            errors.append(f"Row {i}: date_added must be YYYY-MM-DD, got '{da}'")

    # Check FASTA files exist and contain the expected alleles
    sequences = load_sequences(repo_root, manifest)
    for entry in manifest:
        name = entry.get("name", "")
        if name and name not in sequences:
            errors.append(f"Allele '{name}' not found in FASTA file '{entry.get('sequence_file', '')}'")

    # Check for valid nucleotide sequences
    for name, seq in sequences.items():
        if not seq:
            errors.append(f"Allele '{name}' has an empty sequence")
        elif not re.match(r"^[ACGTNRYSWKMBDHV]+$", seq, re.IGNORECASE):
            errors.append(f"Allele '{name}' contains invalid nucleotide characters")

    # Check for duplicate names
    names = [e.get("name", "") for e in manifest]
    seen = set()
    for n in names:
        if n in seen:
            errors.append(f"Duplicate provisional allele name: '{n}'")
        seen.add(n)

    # Check for duplicate sequences
    seq_hashes: dict[str, str] = {}
    for name, seq in sequences.items():
        h = sequence_hash(seq)
        if h in seq_hashes:
            errors.append(
                f"Allele '{name}' has the same sequence as '{seq_hashes[h]}'"
            )
        else:
            seq_hashes[h] = name

    return errors


def validate_no_ipd_collisions(
    manifest: list[dict],
    ipd_names: set[str],
) -> list[str]:
    """Check that no provisional names collide with IPD allele names."""
    errors = []
    for entry in manifest:
        name = entry.get("name", "")
        if name in ipd_names:
            errors.append(f"Provisional name '{name}' collides with an existing IPD allele")
    return errors


# ---------------------------------------------------------------------------
# Metadata (for alleles.json)
# ---------------------------------------------------------------------------

def build_metadata(
    manifest: list[dict],
    sequences: dict[str, str],
) -> list[dict]:
    """Convert provisional manifest entries to alleles.json-compatible records."""
    records = []
    for entry in manifest:
        name = entry.get("name", "")
        seq = sequences.get(name, "")
        if not seq:
            continue

        rec = {
            "a": compute_accession(seq),
            "n": name,
            "l": entry.get("locus", ""),
            "c": entry.get("class", ""),
            "s": entry.get("species", ""),
            "p": "provisional",
            "da": entry.get("date_added", ""),
            "dm": entry.get("date_added", ""),
            "prov": True,
            "sub": entry.get("submitter", ""),
            "len": len(seq),
            "st": entry.get("seq_type", "coding"),
        }

        rationale = _build_rationale(entry)
        if rationale:
            rec["nr"] = rationale

        fh = entry.get("fasta_header", "")
        if fh:
            rec["fh"] = fh

        records.append(rec)
    return records


def _build_rationale(entry: dict) -> str:
    """Generate a human-readable naming rationale from manifest fields."""
    notes = entry.get("notes", "")
    ref_aa = entry.get("ref_aa", "")
    ref_aa_pct = entry.get("ref_aa_pct", "")
    ref_nt = entry.get("ref_nt", "")
    ref_nt_pct = entry.get("ref_nt_pct", "")
    name = entry.get("name", "")
    cds_diffs_str = entry.get("cds_diffs", "")
    aa_diffs_str = entry.get("aa_diffs", "")

    cds_diffs = int(cds_diffs_str) if cds_diffs_str else None
    aa_diffs = int(aa_diffs_str) if aa_diffs_str else None

    # Determine relationship from the notes field
    # Check non_synonymous before synonymous to avoid substring match
    relationship = ""
    for rel in ("extension", "non_synonymous", "synonymous", "fallback"):
        if rel in notes:
            relationship = rel
            break

    if not relationship:
        return ""

    if relationship == "extension" and ref_nt:
        parts = [f"Extension of {ref_nt}."]
        if ref_nt_pct:
            parts.append(f"Provisional sequence contains the full IPD nucleotide sequence ({ref_nt_pct}% nucleotide identity) plus additional flanking regions.")
        return " ".join(parts)

    if relationship == "synonymous" and ref_aa:
        # Check if synonymous to a provisional allele (ref_aa contains "new")
        if "new" in ref_aa:
            parts = [
                f"Amino acid sequence identical to provisional allele {ref_aa}.",
            ]
            if cds_diffs is not None and cds_diffs > 0:
                parts.append(f"{cds_diffs} synonymous CDS nucleotide difference{'s' if cds_diffs != 1 else ''}.")
            elif cds_diffs == 0:
                parts.append("Identical CDS; differs in intronic/UTR regions.")
            else:
                parts.append("Differs at synonymous nucleotide positions.")
            parts.append("Both encode a novel protein not yet in IPD.")
            return " ".join(parts)
        else:
            parts = [f"Amino acid sequence identical to {ref_aa}."]
            if cds_diffs is not None and cds_diffs > 0:
                parts.append(f"{cds_diffs} synonymous CDS nucleotide difference{'s' if cds_diffs != 1 else ''}.")
            elif cds_diffs == 0:
                parts.append("Identical CDS; differs in intronic/UTR regions.")
            else:
                parts.append("Differs at synonymous nucleotide positions.")
            return " ".join(parts)

    if relationship == "non_synonymous":
        parts = ["Novel protein sequence not yet in IPD."]
        if ref_aa and ref_aa_pct:
            diff_detail = ""
            if aa_diffs is not None and aa_diffs > 0:
                diff_detail = f", {aa_diffs} amino acid difference{'s' if aa_diffs != 1 else ''}"
            if cds_diffs is not None and cds_diffs > 0:
                diff_detail += f", {cds_diffs} CDS nucleotide difference{'s' if cds_diffs != 1 else ''}"
            parts.append(
                f"Closest match: {ref_aa} ({ref_aa_pct}% amino acid identity{diff_detail})."
            )
        # Extract group from the name to explain assignment
        group_match = re.search(r"\*(\d+):new", name)
        if group_match:
            parts.append(
                f"Assigned to the *{group_match.group(1)} group based on protein similarity."
            )
        return " ".join(parts)

    if relationship == "fallback":
        parts = ["Could not be confidently assigned to an allele group."]
        if ref_nt and ref_nt_pct:
            parts.append(
                f"Closest nucleotide match: {ref_nt} ({ref_nt_pct}% identity)."
            )
        parts.append("Assigned placeholder group *000.")
        return " ".join(parts)

    return ""


# ---------------------------------------------------------------------------
# Sequence index (for retirement matching)
# ---------------------------------------------------------------------------

def build_sequence_index(repo_root: Path) -> dict[str, dict]:
    """Build a hash→{accession, name} index from all IPD GenBank files on disk.

    Parses ORIGIN sequences from data/MHC/ and data/NHKIR/ GenBank files.
    Saves the index to provisional/ipd_sequences.json.
    """
    index: dict[str, dict] = {}
    data_dir = repo_root / DATA_DIR

    for project_dir in (data_dir / "MHC", data_dir / "NHKIR"):
        if not project_dir.exists():
            continue
        for gb_file in project_dir.rglob("*.gb"):
            try:
                for record in SeqIO.parse(str(gb_file), "genbank"):
                    seq = str(record.seq).upper()
                    if seq:
                        h = sequence_hash(seq)
                        index[h] = {
                            "accession": record.id,
                            "name": record.description.split(",")[0].strip(),
                        }
            except Exception as e:
                logger.warning("Error parsing %s: %s", gb_file, e)

    # Save the index
    index_path = repo_root / PROVISIONAL_SEQ_INDEX
    index_path.parent.mkdir(parents=True, exist_ok=True)
    index_path.write_text(json.dumps(index, separators=(",", ":")))
    logger.info("Built sequence index: %d IPD sequences", len(index))

    return index


def load_sequence_index(repo_root: Path) -> dict[str, dict]:
    """Load the sequence index from disk."""
    index_path = repo_root / PROVISIONAL_SEQ_INDEX
    if not index_path.exists():
        return {}
    return json.loads(index_path.read_text())


# ---------------------------------------------------------------------------
# Auto-retirement
# ---------------------------------------------------------------------------

def check_retirements(
    manifest: list[dict],
    sequences: dict[str, str],
    ipd_index: dict[str, dict],
) -> list[dict]:
    """Find provisional alleles whose sequences now exist in IPD.

    Returns a list of retirement records.
    """
    retirements = []
    for entry in manifest:
        name = entry.get("name", "")
        seq = sequences.get(name, "")
        if not seq:
            continue

        h = sequence_hash(seq)
        if h in ipd_index:
            retirements.append({
                "provisional_name": name,
                "ipd_accession": ipd_index[h]["accession"],
                "ipd_name": ipd_index[h]["name"],
            })

    return retirements


def retire_alleles(
    repo_root: Path,
    retirements: list[dict],
) -> None:
    """Remove retired alleles from manifest/FASTA and log to retired.tsv."""
    if not retirements:
        return

    from datetime import date

    retired_names = {r["provisional_name"] for r in retirements}

    # 1. Remove from manifest
    manifest = load_manifest(repo_root)
    remaining = [e for e in manifest if e.get("name") not in retired_names]
    _write_manifest(repo_root, remaining)

    # 2. Remove from FASTA files
    seq_dir = repo_root / PROVISIONAL_SEQUENCES_DIR
    fasta_files: set[str] = set()
    for entry in manifest:
        if entry.get("name") in retired_names:
            fasta_files.add(entry.get("sequence_file", ""))

    for rel_path in fasta_files:
        if not rel_path:
            continue
        fasta_path = seq_dir / rel_path
        if not fasta_path.exists():
            continue

        # Keep only records NOT being retired
        kept = []
        for rec in SeqIO.parse(str(fasta_path), "fasta"):
            if rec.id not in retired_names:
                kept.append(rec)

        if kept:
            with open(fasta_path, "w") as fh:
                SeqIO.write(kept, fh, "fasta")
        else:
            fasta_path.unlink()

    # 3. Clean up empty GenBank files
    prov_data_dir = repo_root / PROVISIONAL_DATA_DIR
    if prov_data_dir.exists():
        for gb_file in prov_data_dir.rglob("*.gb"):
            # Check if any remaining manifest entries reference this locus
            species = gb_file.parent.name
            locus = gb_file.stem
            has_remaining = any(
                e.get("species") == species and e.get("locus") == locus
                for e in remaining
            )
            if not has_remaining:
                gb_file.unlink()
                logger.info("Removed empty provisional GenBank: %s", gb_file)

    # 4. Append to retired.tsv
    retired_path = repo_root / PROVISIONAL_RETIRED
    retired_path.parent.mkdir(parents=True, exist_ok=True)
    with open(retired_path, "a", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        today = date.today().isoformat()
        for r in retirements:
            writer.writerow([
                r["provisional_name"],
                r["ipd_accession"],
                r["ipd_name"],
                today,
            ])

    for r in retirements:
        logger.info(
            "Retired provisional allele %s → IPD %s (%s)",
            r["provisional_name"], r["ipd_accession"], r["ipd_name"],
        )


def _write_manifest(repo_root: Path, entries: list[dict]) -> None:
    """Write the manifest TSV from a list of dicts."""
    manifest_path = repo_root / PROVISIONAL_MANIFEST
    with open(manifest_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=MANIFEST_COLUMNS, delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(entries)


def append_to_manifest(repo_root: Path, entry: dict) -> None:
    """Append a single entry to the manifest TSV."""
    manifest_path = repo_root / PROVISIONAL_MANIFEST
    file_exists = manifest_path.exists() and manifest_path.stat().st_size > 0
    with open(manifest_path, "a", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=MANIFEST_COLUMNS, delimiter="\t",
            extrasaction="ignore",
        )
        if not file_exists:
            writer.writeheader()
        writer.writerow(entry)


# ---------------------------------------------------------------------------
# GenBank generation
# ---------------------------------------------------------------------------

def build_genbank(
    repo_root: Path,
    manifest: list[dict],
    sequences: dict[str, str],
    species_map: dict[str, dict],
) -> int:
    """Generate provisional GenBank files at data/provisional/{species}/{locus}.gb.

    Args:
        repo_root: Repository root directory.
        manifest: Loaded provisional manifest entries.
        sequences: {allele_name: uppercase_sequence} dict.
        species_map: {species_prefix: {scientificName, commonName, taxon}} dict.

    Returns:
        Number of GenBank records written.
    """
    from .genbank_builder import provisional_to_seqrecord

    # Group entries by species/locus
    groups: dict[tuple[str, str], list[dict]] = {}
    for entry in manifest:
        name = entry.get("name", "")
        if name not in sequences:
            continue
        key = (entry.get("species", ""), entry.get("locus", ""))
        groups.setdefault(key, []).append(entry)

    total = 0
    prov_data_dir = repo_root / PROVISIONAL_DATA_DIR

    for (species, locus), entries in groups.items():
        seq_records = []
        for entry in entries:
            name = entry.get("name", "")
            seq = sequences.get(name, "")
            if not seq:
                continue

            species_info = species_map.get(species, {})
            accession = compute_accession(seq)

            try:
                sr = provisional_to_seqrecord(
                    name=name,
                    sequence=seq,
                    species=species,
                    locus=locus,
                    allele_class=entry.get("class", ""),
                    seq_type=entry.get("seq_type", "coding"),
                    species_info=species_info,
                    date_added=entry.get("date_added", ""),
                    submitter=entry.get("submitter", ""),
                    accession=accession,
                )
                seq_records.append(sr)
            except Exception as e:
                logger.warning("Skipping provisional %s: %s", name, e)

        if seq_records:
            from Bio import SeqIO

            gb_path = prov_data_dir / species / f"{locus}.gb"
            gb_path.parent.mkdir(parents=True, exist_ok=True)
            with open(gb_path, "w") as fh:
                SeqIO.write(seq_records, fh, "genbank")
            total += len(seq_records)
            logger.info("Wrote %d provisional records to %s", len(seq_records), gb_path)

    return total


def backfill_rationale(repo_root: Path) -> int:
    """Re-run the namer on manifest entries missing rationale fields.

    Returns the number of entries updated.
    """
    from .namer import ReferenceAllele, name_provisional_allele

    manifest = load_manifest(repo_root)
    if not manifest:
        return 0

    # Find entries missing rationale (ref_nt or cds_diffs)
    needs_backfill = [
        e for e in manifest
        if (not e.get("ref_nt") or not e.get("cds_diffs")) and e.get("notes", "") != ""
    ]
    if not needs_backfill:
        return 0

    logger.info("Backfilling rationale for %d provisional allele(s)", len(needs_backfill))
    sequences = load_sequences(repo_root, manifest)

    # Group by species/locus so we load references once per locus
    by_locus: dict[tuple[str, str], list[dict]] = {}
    for entry in needs_backfill:
        key = (entry.get("species", ""), entry.get("locus", ""))
        by_locus.setdefault(key, []).append(entry)

    updated = 0
    for (species, locus), entries in by_locus.items():
        gb_path = repo_root / DATA_DIR / "MHC" / species / f"{locus}.gb"
        if not gb_path.exists():
            gb_path = repo_root / DATA_DIR / "NHKIR" / species / f"{locus}.gb"

        # Collect all existing names at this locus
        existing_names = [
            e.get("name", "") for e in manifest
            if e.get("species") == species and e.get("locus") == locus
        ]

        # Build provisional refs for intra-batch synonymous detection
        provisional_refs: list[ReferenceAllele] = []

        for entry in entries:
            name = entry.get("name", "")
            seq = sequences.get(name, "")
            if not seq:
                continue

            seq_type = entry.get("seq_type", "coding")
            result = name_provisional_allele(
                sequence=seq,
                species=species,
                locus=locus,
                seq_type=seq_type,
                locus_gb_path=gb_path,
                existing_names=existing_names,
                provisional_refs=provisional_refs,
            )

            entry["ref_nt"] = result.closest_nt_allele
            entry["ref_nt_pct"] = f"{result.closest_nt_identity:.1f}" if result.closest_nt_allele else ""
            entry["ref_aa"] = result.closest_protein_allele
            entry["ref_aa_pct"] = f"{result.closest_protein_identity:.1f}" if result.closest_protein_allele else ""
            entry["cds_diffs"] = str(result.cds_diffs) if result.cds_diffs >= 0 else ""
            entry["aa_diffs"] = str(result.aa_diffs) if result.aa_diffs >= 0 else ""
            updated += 1

            # Track for intra-batch detection
            if result.extracted_protein:
                provisional_refs.append(
                    ReferenceAllele(
                        name=name,
                        accession="",
                        sequence=seq,
                        mol_type="genomic DNA" if seq_type == "genomic" else "mRNA",
                        protein=result.extracted_protein,
                        coding_seq=result.extracted_cds,
                        exon_coords=[],
                        intron_coords=[],
                        codon_start=1,
                    )
                )

            logger.info(
                "  %s: ref_nt=%s (%.1f%%), ref_aa=%s (%.1f%%)",
                name,
                result.closest_nt_allele, result.closest_nt_identity,
                result.closest_protein_allele or "none", result.closest_protein_identity,
            )

    if updated:
        _write_manifest(repo_root, manifest)
        logger.info("Backfilled rationale for %d entries", updated)

    return updated


def append_to_fasta(repo_root: Path, rel_path: str, name: str, sequence: str) -> None:
    """Append a sequence to a provisional FASTA file, creating it if needed."""
    fasta_path = repo_root / PROVISIONAL_SEQUENCES_DIR / rel_path
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    with open(fasta_path, "a") as fh:
        fh.write(f">{name}\n")
        # Wrap at 70 chars
        for i in range(0, len(sequence), 70):
            fh.write(sequence[i:i + 70] + "\n")
