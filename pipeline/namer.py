"""Improved provisional allele naming using protein-level comparison.

This module extracts CDS from genomic submissions by transferring exon
structure from the closest IPD reference allele, translates to protein,
and classifies the relationship (identical protein, non-synonymous, new
group) to produce names like:

    Mamu-E*02:01:new01   (synonymous variant of 02:01)
    Mamu-E*02:new01      (non-synonymous variant in the 02 group)

Algorithm overview:
    1. Auto-detect sequence type (genomic vs CDS vs partial)
    2. Find closest IPD reference at the nucleotide level
    3. Transfer exon coordinates from reference to query via alignment
    4. Extract CDS, translate to protein
    5. Compare protein to all IPD reference proteins at the locus
    6. Classify relationship and assign name accordingly

Dependencies: BioPython (PairwiseAligner, SeqIO, Seq, SeqFeature)
"""

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Cross-species fallback mapping
# ---------------------------------------------------------------------------
# Closely related species that can share genomic reference templates for
# exon coordinate transfer when one species lacks genomic IPD references.
# All macaques (genus Macaca) are closely related enough for reliable exon
# boundary transfer.  Mamu and Mafa have the richest IPD reference sets,
# so they are listed first as preferred fallback sources.
_MACAQUE_SPECIES = [
    "Mamu",  # Macaca mulatta (rhesus)
    "Mafa",  # Macaca fascicularis (cynomolgus)
    "Mane",  # Macaca nemestrina (pig-tailed)
    "Maar",  # Macaca arctoides (stump-tailed)
    "Maas",  # Macaca assamensis (Assamese)
    "Mafu",  # Macaca fuscata (Japanese)
    "Male",  # Macaca leonina (northern pig-tailed)
    "Malo",  # Macaca ?
    "Masi",  # Macaca sinica (toque)
    "Masp",  # Macaca sp.
    "Math",  # Macaca thibetana (Tibetan)
]
RELATED_SPECIES: dict[str, list[str]] = {
    sp: [r for r in ["Mamu", "Mafa"] + _MACAQUE_SPECIES if r != sp]
    for sp in _MACAQUE_SPECIES
}


def _find_cross_species_gb(
    locus_gb_path: Path,
    species: str,
    locus: str,
) -> Path | None:
    """Look for a GenBank file at the same locus in a related species.

    Walks the RELATED_SPECIES map and checks for a file at the same
    relative position (e.g., data/MHC/Mafa/E.gb when Mamu/E.gb is given).

    Returns the first existing path or None.
    """
    # Determine the project directory (MHC or NHKIR)
    # locus_gb_path looks like: .../data/MHC/Mamu/E.gb
    parts = locus_gb_path.parts
    if len(parts) < 3:
        return None

    # Find the species directory in the path
    species_dir_idx = None
    for i, part in enumerate(parts):
        if part == species:
            species_dir_idx = i
            break

    if species_dir_idx is None:
        return None

    relatives = RELATED_SPECIES.get(species, [])
    for rel_species in relatives:
        candidate = Path(*parts[:species_dir_idx], rel_species, *parts[species_dir_idx + 1:])
        if candidate.exists():
            return candidate

    return None


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass
class ReferenceAllele:
    """Parsed IPD reference allele with extracted features."""

    name: str
    accession: str
    sequence: str  # full primary sequence (genomic or CDS)
    mol_type: str  # "genomic DNA" or "mRNA"
    protein: str  # translation from /translation qualifier
    coding_seq: str  # CDS nucleotide sequence (exons joined)
    exon_coords: list[tuple[int, int]]  # list of (start, end) 0-based
    intron_coords: list[tuple[int, int]]
    codon_start: int  # 1-based codon start position

    @property
    def is_genomic(self) -> bool:
        return "genomic" in self.mol_type.lower()

    @property
    def allele_group(self) -> str:
        """Extract first field: e.g. '02' from 'Mamu-E*02:01'."""
        m = re.search(r"\*(\d+)", self.name)
        return m.group(1) if m else ""

    @property
    def allele_fields(self) -> list[str]:
        """Extract all colon-separated fields after the asterisk.

        E.g. 'Mamu-E*02:01:03:01' -> ['02', '01', '03', '01']
             'Mamu-E*02:01'       -> ['02', '01']
        """
        m = re.search(r"\*(.+)$", self.name)
        if not m:
            return []
        return m.group(1).split(":")

    @property
    def protein_field(self) -> str:
        """Return the first two fields joined, e.g. '02:01'."""
        fields = self.allele_fields
        if len(fields) >= 2:
            return f"{fields[0]}:{fields[1]}"
        elif len(fields) == 1:
            return fields[0]
        return ""


@dataclass
class NamingResult:
    """Result of the naming algorithm for a single submitted sequence."""

    # Assigned name
    provisional_name: str

    # Classification
    relationship: str  # "synonymous", "non_synonymous", "novel_group", "fallback", "extension"

    # Closest matches
    closest_nt_allele: str  # closest by nucleotide alignment
    closest_nt_identity: float  # percent identity at nt level
    closest_protein_allele: str  # closest by protein identity
    closest_protein_identity: float  # percent identity at aa level

    # Sequence type detection
    detected_seq_type: str  # "genomic", "coding", "partial"

    # Extracted CDS info
    cds_extracted: bool
    cds_length: int
    protein_length: int
    exon_count: int

    # Difference counts vs closest protein-level match
    cds_diffs: int = -1  # CDS nucleotide mismatches (-1 = not computed)
    aa_diffs: int = -1   # amino acid mismatches (-1 = not computed)

    # Extension detection
    is_extension: bool = False  # query contains an IPD allele as a substring

    # Extracted sequences (used by batch naming for intra-batch comparison)
    extracted_cds: str = ""
    extracted_protein: str = ""

    # Warnings
    warnings: list[str] = field(default_factory=list)


def _count_diffs(seq_a: str, seq_b: str) -> int:
    """Count mismatches between two sequences of equal or similar length.

    Compares character-by-character up to the shorter length, then adds
    the length difference as additional differences.
    """
    if not seq_a or not seq_b:
        return -1
    a = seq_a.upper().rstrip("*")
    b = seq_b.upper().rstrip("*")
    diffs = abs(len(a) - len(b))
    for i in range(min(len(a), len(b))):
        if a[i] != b[i]:
            diffs += 1
    return diffs


# ---------------------------------------------------------------------------
# Sequence type auto-detection
# ---------------------------------------------------------------------------

# Canonical splice donor/acceptor dinucleotides
_SPLICE_DONOR = "GT"
_SPLICE_ACCEPTOR = "AG"


def detect_sequence_type(
    sequence: str,
    expected_cds_lengths: list[int] | None = None,
    expected_genomic_lengths: list[int] | None = None,
) -> str:
    """Classify a submitted sequence as genomic, coding, or partial.

    Strategy:
        1. Length heuristics: MHC class I CDS is ~1,080 bp; genomic is
           ~2,800-3,200 bp. If the sequence is within ~15% of known CDS
           lengths and has no intron signatures, it is CDS.
        2. Intron signature scan: look for GT...AG bounded regions between
           exon-like coding segments. Genomic sequences will have multiple
           GT...AG intron candidates in conserved positions.
        3. ORF analysis: if the sequence translates cleanly from ATG to
           stop without interruption, it is more likely CDS.
        4. If the sequence is much shorter than expected CDS, flag as partial.

    Args:
        sequence: uppercase nucleotide sequence.
        expected_cds_lengths: typical CDS lengths for this locus (from
            reference alleles). Used for length-based classification.
        expected_genomic_lengths: typical genomic lengths for this locus.

    Returns:
        "genomic", "coding", or "partial".
    """
    seq = sequence.upper()
    seq_len = len(seq)

    # Compute typical CDS / genomic lengths from references if provided
    if expected_cds_lengths:
        median_cds = sorted(expected_cds_lengths)[len(expected_cds_lengths) // 2]
    else:
        median_cds = 1080  # MHC class I default

    if expected_genomic_lengths:
        median_genomic = sorted(expected_genomic_lengths)[
            len(expected_genomic_lengths) // 2
        ]
    else:
        median_genomic = 2900  # MHC class I default

    # Check if sequence is too short to be a full CDS
    min_expected_cds = int(median_cds * 0.5)
    if seq_len < min_expected_cds:
        return "partial"

    # Length-based heuristic: if within 20% of genomic length, likely genomic
    if expected_genomic_lengths and abs(seq_len - median_genomic) / median_genomic < 0.20:
        return "genomic"

    # If within 15% of CDS length and translates cleanly, likely CDS
    if abs(seq_len - median_cds) / median_cds < 0.15:
        if _has_clean_orf(seq):
            return "coding"

    # Intron signature scan: count GT...AG candidate introns
    # In genomic MHC sequences, introns are typically 80-700 bp long
    intron_candidates = _count_intron_signatures(seq)
    if intron_candidates >= 3:
        return "genomic"

    # If significantly longer than CDS but no clear intron signatures,
    # still likely genomic (could be divergent splice sites)
    if seq_len > median_cds * 1.5:
        return "genomic"

    # If close to CDS length
    if abs(seq_len - median_cds) / median_cds < 0.25:
        return "coding"

    # Default: if longer than expected CDS, call genomic; otherwise partial
    if seq_len >= median_cds:
        return "genomic"

    return "partial"


def _has_clean_orf(seq: str) -> bool:
    """Check if the sequence has a clean ORF from first ATG to stop."""
    start = seq.find("ATG")
    if start < 0:
        return False
    # Try all three frames starting from first ATG
    coding = seq[start:]
    # Trim to multiple of 3
    trim_len = (len(coding) // 3) * 3
    coding = coding[:trim_len]
    if len(coding) < 300:  # too short to be meaningful
        return False
    try:
        protein = str(Seq(coding).translate())
        # A clean ORF has at most 1 stop (at the end)
        internal_stops = protein[:-1].count("*")
        return internal_stops == 0
    except Exception:
        return False


def _count_intron_signatures(seq: str, min_intron: int = 60, max_intron: int = 1500) -> int:
    """Count GT...AG candidate intron sites in the sequence.

    Scans for GT dinucleotides followed by AG dinucleotides at a
    plausible intron distance (60-1500 bp), where the region between
    them would interrupt a reading frame.
    """
    count = 0
    i = 0
    while i < len(seq) - min_intron:
        if seq[i : i + 2] == _SPLICE_DONOR:
            # Look for AG acceptor at plausible distance
            for j in range(i + min_intron, min(i + max_intron, len(seq) - 1)):
                if seq[j : j + 2] == _SPLICE_ACCEPTOR:
                    count += 1
                    i = j + 2  # skip past this intron
                    break
            else:
                i += 1
        else:
            i += 1
    return count


# ---------------------------------------------------------------------------
# Reference allele parsing
# ---------------------------------------------------------------------------


def load_reference_alleles(gb_path: Path) -> list[ReferenceAllele]:
    """Parse all IPD reference alleles from a locus GenBank file.

    Extracts exon coordinates, CDS, protein translation, and other
    metadata needed for comparison.

    Args:
        gb_path: Path to a per-locus GenBank file
            (e.g., data/MHC/Mamu/E.gb).

    Returns:
        List of ReferenceAllele objects.
    """
    if not gb_path.exists():
        return []

    alleles = []
    for record in SeqIO.parse(str(gb_path), "genbank"):
        seq_str = str(record.seq).upper()
        if not seq_str:
            continue

        # Parse allele name from description
        allele_name = record.description.split(",")[0].strip()

        # Extract mol_type from source feature
        mol_type = ""
        for feat in record.features:
            if feat.type == "source":
                mol_type = feat.qualifiers.get("mol_type", [""])[0]
                break

        # Extract exon coordinates (0-based, half-open)
        exons = []
        for feat in record.features:
            if feat.type == "exon":
                exons.append((int(feat.location.start), int(feat.location.end)))
        exons.sort(key=lambda x: x[0])

        # Extract intron coordinates
        introns = []
        for feat in record.features:
            if feat.type == "intron":
                introns.append((int(feat.location.start), int(feat.location.end)))
        introns.sort(key=lambda x: x[0])

        # Extract CDS translation and codon_start
        protein = ""
        codon_start = 1
        for feat in record.features:
            if feat.type == "CDS":
                protein = feat.qualifiers.get("translation", [""])[0]
                cs = feat.qualifiers.get("codon_start", ["1"])[0]
                codon_start = int(cs)
                break

        # Build coding sequence from exon coordinates
        if exons:
            coding_parts = [seq_str[s:e] for s, e in exons]
            coding_seq = "".join(coding_parts)
        else:
            # No exon annotations: entire sequence is CDS
            coding_seq = seq_str

        # If no protein from /translation, derive it
        if not protein and coding_seq:
            protein = _translate_cds(coding_seq, codon_start)

        alleles.append(
            ReferenceAllele(
                name=allele_name,
                accession=record.id,
                sequence=seq_str,
                mol_type=mol_type,
                protein=protein,
                coding_seq=coding_seq,
                exon_coords=exons,
                intron_coords=introns,
                codon_start=codon_start,
            )
        )

    return alleles


def _translate_cds(coding_seq: str, codon_start: int = 1) -> str:
    """Translate a CDS nucleotide sequence to protein.

    Handles codon_start offset and strips trailing stop codon.
    Returns empty string on failure.
    """
    try:
        offset = codon_start - 1  # codon_start is 1-based
        trimmed = coding_seq[offset:]
        # Trim to multiple of 3
        trim_len = (len(trimmed) // 3) * 3
        trimmed = trimmed[:trim_len]
        if not trimmed:
            return ""
        protein = str(Seq(trimmed).translate())
        # Strip trailing stop
        if protein.endswith("*"):
            protein = protein[:-1]
        return protein
    except Exception:
        return ""


# ---------------------------------------------------------------------------
# CDS extraction from genomic sequence via coordinate transfer
# ---------------------------------------------------------------------------


def transfer_exon_coordinates(
    query_seq: str,
    ref: ReferenceAllele,
) -> list[tuple[int, int]] | None:
    """Transfer exon/intron coordinates from a reference allele to a query
    genomic sequence using pairwise alignment.

    This is the core of the genomic-to-CDS extraction. The algorithm:

        1. Align the query genomic sequence to the reference genomic
           sequence (or to reference CDS if reference is CDS-only).
        2. Walk the alignment to build a coordinate map:
           ref_pos -> query_pos.
        3. Map each reference exon boundary through the coordinate map
           to get query exon boundaries.
        4. Validate that transferred intron boundaries have GT...AG
           splice site consensus in the query.
        5. Fall back to splice-site scanning if alignment-based transfer
           fails validation.

    Args:
        query_seq: Uppercase genomic nucleotide sequence.
        ref: The closest ReferenceAllele (should be genomic if possible).

    Returns:
        List of (start, end) tuples for exons in the query sequence,
        0-based half-open. Returns None if extraction fails.
    """
    if not ref.exon_coords:
        return None

    ref_seq = ref.sequence

    # Build the alignment
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    try:
        alignments = aligner.align(ref_seq, query_seq)
        aln = alignments[0]  # take top alignment
    except Exception as e:
        logger.warning("Alignment failed for coordinate transfer: %s", e)
        return None

    # Build coordinate map: ref_pos -> query_pos
    # Walk the alignment character by character
    ref_to_query = _build_coordinate_map(aln, len(ref_seq), len(query_seq))

    # Map each reference exon boundary to query coordinates
    query_exons = []
    for ref_start, ref_end in ref.exon_coords:
        q_start = ref_to_query.get(ref_start)
        # For the end coordinate, we want the position that maps to
        # ref_end - 1, then add 1 (since end is exclusive)
        q_end_inclusive = ref_to_query.get(ref_end - 1)

        if q_start is None or q_end_inclusive is None:
            logger.warning(
                "Could not map exon %d-%d from ref %s to query",
                ref_start, ref_end, ref.name,
            )
            return None

        q_end = q_end_inclusive + 1
        query_exons.append((q_start, q_end))

    # Validate splice sites in the query
    query_exons = _validate_and_refine_splice_sites(query_seq, query_exons)

    return query_exons


def _build_coordinate_map(
    alignment,
    ref_len: int,
    query_len: int,
) -> dict[int, int]:
    """Build a reference-to-query coordinate map from a BioPython alignment.

    Walks through the aligned pairs to create a mapping from each
    reference position to its corresponding query position.
    """
    coord_map: dict[int, int] = {}

    try:
        # Use the aligned property which gives pairs of indices
        # aligned is a list of arrays, one per segment
        for ref_indices, query_indices in zip(
            alignment.aligned[0], alignment.aligned[1]
        ):
            ref_start, ref_end = ref_indices
            q_start, q_end = query_indices
            seg_len = ref_end - ref_start
            for offset in range(seg_len):
                coord_map[ref_start + offset] = q_start + offset
    except Exception:
        # Fallback: walk the alignment path
        # This handles edge cases in different BioPython versions
        try:
            path = alignment.path
            ref_pos = 0
            query_pos = 0
            for i in range(1, len(path)):
                prev_r, prev_q = path[i - 1]
                curr_r, curr_q = path[i]
                dr = curr_r - prev_r
                dq = curr_q - prev_q
                if dr > 0 and dq > 0:
                    # Match/mismatch block
                    for j in range(min(dr, dq)):
                        coord_map[prev_r + j] = prev_q + j
        except Exception as e:
            logger.warning("Could not build coordinate map: %s", e)

    return coord_map


def _validate_and_refine_splice_sites(
    query_seq: str,
    exons: list[tuple[int, int]],
    search_window: int = 5,
) -> list[tuple[int, int]]:
    """Validate and optionally refine exon boundaries by checking splice sites.

    For each intron (gap between consecutive exons), verify that the
    query sequence has GT at the donor site and AG at the acceptor site.
    If not, search within a small window to find the nearest GT...AG pair.

    Args:
        query_seq: Uppercase query genomic sequence.
        exons: List of (start, end) exon coordinates.
        search_window: Number of bases to search around each boundary.

    Returns:
        Refined list of exon coordinates.
    """
    if len(exons) <= 1:
        return exons

    refined = [exons[0]]

    for i in range(1, len(exons)):
        prev_end = refined[-1][1]  # end of previous exon = donor site
        curr_start = exons[i][0]  # start of current exon = acceptor site

        donor_ok = (
            prev_end < len(query_seq) - 1
            and query_seq[prev_end : prev_end + 2] == _SPLICE_DONOR
        )
        acceptor_ok = (
            curr_start >= 2
            and query_seq[curr_start - 2 : curr_start] == _SPLICE_ACCEPTOR
        )

        if donor_ok and acceptor_ok:
            refined.append(exons[i])
            continue

        # Try to refine: search for GT near donor and AG near acceptor
        best_donor = None
        best_acceptor = None

        for d in range(-search_window, search_window + 1):
            pos = prev_end + d
            if 0 <= pos < len(query_seq) - 1:
                if query_seq[pos : pos + 2] == _SPLICE_DONOR:
                    if best_donor is None or abs(d) < abs(best_donor - prev_end):
                        best_donor = pos

        for d in range(-search_window, search_window + 1):
            pos = curr_start + d
            if 2 <= pos <= len(query_seq):
                if query_seq[pos - 2 : pos] == _SPLICE_ACCEPTOR:
                    if best_acceptor is None or abs(d) < abs(
                        best_acceptor - curr_start
                    ):
                        best_acceptor = pos

        if best_donor is not None and best_acceptor is not None:
            # Adjust the previous exon's end and the current exon's start
            refined[-1] = (refined[-1][0], best_donor)
            refined.append((best_acceptor, exons[i][1]))
        else:
            # Could not find splice sites; keep original coordinates
            # (the CDS may still be usable even without perfect splicing)
            refined.append(exons[i])

    return refined


def extract_cds_from_genomic(
    query_seq: str,
    ref: ReferenceAllele,
) -> tuple[str, list[tuple[int, int]]]:
    """Extract CDS from a genomic query using reference exon structure.

    Args:
        query_seq: Uppercase genomic nucleotide sequence.
        ref: Closest ReferenceAllele with exon annotations.

    Returns:
        (coding_sequence, exon_coords_in_query). Returns ("", []) on
        failure.
    """
    if not ref.exon_coords:
        # Reference has no exon annotations (single-exon CDS-only record).
        # The query is genomic but we have no template for exon structure.
        # Fall back: try to find the longest ORF.
        logger.warning(
            "Reference %s has no exon annotations; attempting ORF-based "
            "CDS extraction from genomic query",
            ref.name,
        )
        return _extract_cds_by_orf(query_seq)

    # If reference is CDS-only (mRNA) but has exon annotations, those
    # exons are contiguous in the CDS. We need a genomic reference to
    # transfer intron positions. However, we can still align the query
    # to the reference CDS and attempt to identify exon boundaries by
    # splice site scanning.
    if not ref.is_genomic:
        logger.info(
            "Reference %s is CDS-only; will use exon sizes + splice site "
            "scanning for coordinate transfer",
            ref.name,
        )
        return _extract_cds_by_exon_sizes_and_splice_scan(
            query_seq, ref
        )

    # Primary path: reference is genomic with exon/intron annotations
    query_exons = transfer_exon_coordinates(query_seq, ref)
    if query_exons is None:
        return ("", [])

    # Join exon sequences to build CDS
    cds_parts = []
    for start, end in query_exons:
        if start < 0 or end > len(query_seq):
            logger.warning("Exon coordinates out of bounds: %d-%d", start, end)
            return ("", [])
        cds_parts.append(query_seq[start:end])

    cds = "".join(cds_parts)
    return (cds, query_exons)


def _extract_cds_by_orf(
    query_seq: str,
    max_cds_length: int = 1500,
) -> tuple[str, list[tuple[int, int]]]:
    """Extract the best ORF from a sequence as a fallback.

    For genomic sequences that contain introns, the longest ORF may span
    intron boundaries and be too long. This function collects all ORFs
    and picks the one closest to typical MHC CDS length (~1,080 bp) that
    doesn't exceed max_cds_length.

    Returns (cds, [(start, end)]) for the best ATG..stop ORF.
    """
    orfs: list[tuple[str, int, int]] = []

    for frame in range(3):
        i = frame
        while i < len(query_seq) - 2:
            codon = query_seq[i : i + 3]
            if codon == "ATG":
                orf_start = i
                j = i + 3
                while j < len(query_seq) - 2:
                    c = query_seq[j : j + 3]
                    if c in ("TAA", "TAG", "TGA"):
                        orf_end = j + 3
                        orf = query_seq[orf_start:orf_end]
                        if len(orf) >= 9:  # at least ATG + 1 codon + stop
                            orfs.append((orf, orf_start, orf_end))
                        break
                    j += 3
                i = j + 3 if j < len(query_seq) - 2 else i + 3
            else:
                i += 3

    if not orfs:
        return ("", [])

    # Filter to ORFs under max_cds_length, then pick the longest
    reasonable = [(o, s, e) for o, s, e in orfs if len(o) <= max_cds_length]
    if reasonable:
        best = max(reasonable, key=lambda x: len(x[0]))
    else:
        # All ORFs exceed max length — pick shortest (least wrong)
        best = min(orfs, key=lambda x: len(x[0]))
        logger.warning(
            "All ORFs exceed max CDS length %d bp; shortest is %d bp",
            max_cds_length, len(best[0]),
        )

    return (best[0], [(best[1], best[2])])


def _extract_cds_by_exon_sizes_and_splice_scan(
    query_seq: str,
    ref: ReferenceAllele,
) -> tuple[str, list[tuple[int, int]]]:
    """Extract CDS when reference is CDS-only but has exon size annotations.

    Uses the known exon sizes from the reference to guide a splice-site
    scan of the genomic query. Strategy:

        1. Align query to reference CDS to find approximate coding regions.
        2. Use reference exon sizes as expected block sizes.
        3. Scan for GT...AG splice sites at expected intron positions.
    """
    ref_exon_sizes = [end - start for start, end in ref.exon_coords]

    # Align query to reference CDS to find the coding region
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    try:
        alignments = aligner.align(ref.coding_seq, query_seq)
        aln = alignments[0]
    except Exception:
        return ("", [])

    # Get the approximate start position in the query
    # The alignment tells us where the coding region begins in the query
    try:
        query_aligned_blocks = aln.aligned[1]
        if not len(query_aligned_blocks):
            return ("", [])
        query_start = query_aligned_blocks[0][0]
    except Exception:
        return ("", [])

    # Now walk through the query starting from query_start, consuming
    # exon-sized blocks and looking for GT...AG introns between them
    exons = []
    pos = query_start

    for exon_idx, exon_size in enumerate(ref_exon_sizes):
        exon_end = pos + exon_size

        if exon_end > len(query_seq):
            # Truncated sequence; take what we have
            exons.append((pos, len(query_seq)))
            break

        exons.append((pos, exon_end))

        # If not the last exon, find the intron
        if exon_idx < len(ref_exon_sizes) - 1:
            # Look for GT at exon_end position
            found_intron = False
            for d in range(-3, 4):
                donor_pos = exon_end + d
                if (
                    0 <= donor_pos < len(query_seq) - 1
                    and query_seq[donor_pos : donor_pos + 2] == _SPLICE_DONOR
                ):
                    # Scan for AG acceptor
                    for intron_len in range(60, 1500):
                        acc_pos = donor_pos + intron_len
                        if acc_pos + 2 > len(query_seq):
                            break
                        if query_seq[acc_pos : acc_pos + 2] == _SPLICE_ACCEPTOR:
                            # Adjust exon end and set next exon start
                            exons[-1] = (exons[-1][0], donor_pos)
                            pos = acc_pos + 2
                            found_intron = True
                            break
                    if found_intron:
                        break

            if not found_intron:
                # Could not find intron; might be CDS-only or misaligned
                pos = exon_end
        else:
            pos = exon_end

    cds = "".join(query_seq[s:e] for s, e in exons)
    return (cds, exons)


# ---------------------------------------------------------------------------
# Protein comparison
# ---------------------------------------------------------------------------


def compare_proteins(
    query_protein: str,
    references: list[ReferenceAllele],
) -> tuple[ReferenceAllele | None, float]:
    """Find the closest IPD reference allele by protein identity.

    Uses global pairwise alignment of amino acid sequences.

    Args:
        query_protein: Translated protein from the submitted sequence.
        references: All reference alleles at the locus.

    Returns:
        (best_match, percent_identity). Returns (None, 0.0) if no
        references have protein sequences.
    """
    if not query_protein:
        return (None, 0.0)

    aligner = PairwiseAligner()
    aligner.mode = "global"
    # Use BLOSUM62-like scoring for more biologically meaningful alignment
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    best_ref = None
    best_pct = 0.0

    for ref in references:
        if not ref.protein:
            continue

        try:
            score = aligner.score(query_protein, ref.protein)
            max_len = max(len(query_protein), len(ref.protein))
            # Normalize: perfect match score would be max_len * match_score
            pct = (score / (max_len * aligner.match_score)) * 100 if max_len > 0 else 0.0
        except Exception:
            continue

        if pct > best_pct:
            best_pct = pct
            best_ref = ref

    return (best_ref, best_pct)


def find_identical_protein(
    query_protein: str,
    references: list[ReferenceAllele],
) -> list[ReferenceAllele]:
    """Find all reference alleles with identical protein sequences.

    Exact string match after stripping trailing stops.

    Returns:
        List of matching ReferenceAllele objects (may be empty).
    """
    query_clean = query_protein.rstrip("*").upper()
    if not query_clean:
        return []

    matches = []
    for ref in references:
        ref_clean = ref.protein.rstrip("*").upper()
        if ref_clean == query_clean:
            matches.append(ref)

    return matches


# ---------------------------------------------------------------------------
# Naming decision tree
# ---------------------------------------------------------------------------


def classify_relationship(
    query_protein: str,
    query_cds: str,
    references: list[ReferenceAllele],
) -> tuple[str, ReferenceAllele | None, float]:
    """Classify the relationship between a submitted sequence and IPD alleles.

    Decision tree:
        1. If query protein is identical to one or more reference proteins,
           the relationship is "synonymous" (third-field variant). The
           best match is the one whose CDS is most similar.
        2. If query protein differs but the closest match is >95% identical,
           the relationship is "non_synonymous" (second-field variant within
           the same group).
        3. If the closest protein match is <80% identical, the relationship
           is "novel_group" -- but we still assign to the nearest group.
        4. If no protein comparison is possible, return "fallback".

    Args:
        query_protein: Translated protein sequence.
        query_cds: CDS nucleotide sequence.
        references: All reference alleles at the locus.

    Returns:
        (relationship, best_match_ref, protein_identity_pct)
    """
    if not query_protein:
        return ("fallback", None, 0.0)

    # Step 1: Check for exact protein matches
    identical_matches = find_identical_protein(query_protein, references)
    if identical_matches:
        # Among identical-protein alleles, find the one with closest CDS
        # This determines the second field (protein-level designation)
        best_cds_match = _find_closest_cds_among(query_cds, identical_matches)
        if best_cds_match is not None:
            return ("synonymous", best_cds_match, 100.0)
        # Protein matches but CDS comparison failed; still synonymous
        return ("synonymous", identical_matches[0], 100.0)

    # Step 2: Find closest protein match
    best_ref, best_pct = compare_proteins(query_protein, references)
    if best_ref is None:
        return ("fallback", None, 0.0)

    # Any non-identical protein is a "non_synonymous" relationship.
    # We use the best match's group for naming.
    return ("non_synonymous", best_ref, best_pct)


def _find_closest_cds_among(
    query_cds: str,
    candidates: list[ReferenceAllele],
) -> ReferenceAllele | None:
    """Among alleles with identical protein, find the one with closest CDS.

    This is useful because alleles that share the same protein but differ
    at synonymous positions should be named relative to their CDS match.
    """
    if not query_cds or not candidates:
        return None

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    best = None
    best_score = float("-inf")

    for ref in candidates:
        if not ref.coding_seq:
            continue
        try:
            score = aligner.score(query_cds, ref.coding_seq)
            if score > best_score:
                best_score = score
                best = ref
        except Exception:
            continue

    return best


# ---------------------------------------------------------------------------
# Name assignment
# ---------------------------------------------------------------------------


def assign_provisional_name(
    species: str,
    locus: str,
    relationship: str,
    best_match: ReferenceAllele | None,
    existing_names: list[str],
) -> str:
    """Assign a provisional allele name based on the classification.

    Naming patterns:
        - extension:       {Sp}-{Locus}*{Group}:{Protein}:ext{NN}
          Example:         Mamu-E*02:28:ext01
          Meaning:         extends an existing IPD allele with additional sequence

        - synonymous:      {Sp}-{Locus}*{Group}:{Protein}:new{NN}
          Example:         Mamu-E*02:01:new01
          Meaning:         same protein as 02:01, differs at synonymous sites

        - non_synonymous:  {Sp}-{Locus}*{Group}:new{NN}
          Example:         Mamu-E*02:new01
          Meaning:         different protein, but in the 02 group

        - novel_group:     {Sp}-{Locus}*{Group}:new{NN}
          (same as non_synonymous -- we assign to nearest group)

        - fallback:        {Sp}-{Locus}*000:new{NN}
          (no match found; legacy behavior)

    Args:
        species: Species prefix (e.g., "Mamu").
        locus: Locus name (e.g., "E").
        relationship: One of "extension", "synonymous", "non_synonymous",
            "novel_group", "fallback".
        best_match: Closest ReferenceAllele (None for fallback).
        existing_names: All existing provisional names (to avoid collisions).

    Returns:
        Provisional allele name string.
    """
    if relationship == "extension" and best_match is not None:
        # Extension: name at third field level with ext suffix
        protein_field = best_match.protein_field
        if protein_field:
            prefix = f"{species}-{locus}*{protein_field}:ext"
        else:
            group = best_match.allele_group or "000"
            prefix = f"{species}-{locus}*{group}:ext"

    elif relationship == "synonymous" and best_match is not None:
        # Name at the third field level: *Group:Protein:newNN
        protein_field = best_match.protein_field
        if protein_field:
            prefix = f"{species}-{locus}*{protein_field}:new"
        else:
            # Shouldn't happen, but fall back to group level
            group = best_match.allele_group or "000"
            prefix = f"{species}-{locus}*{group}:new"

    elif relationship in ("non_synonymous", "novel_group") and best_match is not None:
        # Name at the second field level: *Group:newNN
        group = best_match.allele_group or "000"
        prefix = f"{species}-{locus}*{group}:new"

    else:
        # Fallback
        prefix = f"{species}-{locus}*000:new"

    # Find the next available number
    existing_nums = []
    for name in existing_names:
        if name.startswith(prefix):
            suffix = name[len(prefix) :]
            if suffix.isdigit():
                existing_nums.append(int(suffix))

    next_num = max(existing_nums, default=0) + 1
    return f"{prefix}{next_num:02d}"


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def name_provisional_allele(
    sequence: str,
    species: str,
    locus: str,
    seq_type: str,
    locus_gb_path: Path,
    existing_names: list[str],
    provisional_refs: list[ReferenceAllele] | None = None,
) -> NamingResult:
    """Run the full naming algorithm for a single submitted sequence.

    This is the primary entry point. It:
        1. Loads all reference alleles from the locus GenBank file.
        2. Auto-detects the sequence type if needed.
        3. Finds the closest nucleotide-level match.
        4. Extracts CDS if the sequence is genomic.
           - If no genomic references exist for this species, tries a
             related species (e.g., Mafa for Mamu) for exon coordinate
             transfer.
        5. Translates to protein.
        6. Classifies the relationship (against IPD refs + provisional_refs).
        7. Assigns the name.

    Args:
        sequence: Uppercase nucleotide sequence of the submitted allele.
        species: Species prefix (e.g., "Mamu").
        locus: Locus name (e.g., "E").
        seq_type: "coding", "genomic", or "auto".
        locus_gb_path: Path to the locus GenBank file.
        existing_names: All existing provisional + IPD names at this locus
            (to avoid collisions).
        provisional_refs: Previously-named provisional alleles in this batch.
            Used to detect intra-batch synonymous relationships (two
            submissions with the same protein but both novel to IPD).

    Returns:
        NamingResult with the assigned name and metadata.
    """
    warnings: list[str] = []

    # Step 1: Load reference alleles
    references = load_reference_alleles(locus_gb_path)
    if not references:
        return NamingResult(
            provisional_name=assign_provisional_name(
                species, locus, "fallback", None, existing_names
            ),
            relationship="fallback",
            closest_nt_allele="",
            closest_nt_identity=0.0,
            closest_protein_allele="",
            closest_protein_identity=0.0,
            detected_seq_type=seq_type,
            cds_extracted=False,
            cds_length=0,
            protein_length=0,
            exon_count=0,
            warnings=["No reference alleles found at locus"],
        )

    # Step 2: Auto-detect sequence type
    if seq_type == "auto" or not seq_type:
        cds_lengths = [len(r.coding_seq) for r in references if r.coding_seq]
        genomic_lengths = [len(r.sequence) for r in references if r.is_genomic]
        detected_type = detect_sequence_type(
            sequence, cds_lengths or None, genomic_lengths or None
        )
    else:
        detected_type = seq_type

    # Step 3: Find closest nucleotide-level match
    # When comparing genomic query to CDS-only references, use CDS-aware
    # comparison for meaningful identity scores.
    has_genomic_refs = any(r.is_genomic for r in references)
    if detected_type == "genomic" and not has_genomic_refs:
        # All references are CDS-only. Global alignment of 3,190 bp
        # genomic vs 1,080 bp CDS gives misleading low identity.
        # Use local alignment for the NT identity score.
        closest_nt_ref, nt_identity = _find_closest_nt(
            sequence, references, mode="local"
        )
    else:
        closest_nt_ref, nt_identity = _find_closest_nt(sequence, references)

    # Step 4: Extract CDS
    cds = ""
    query_exons: list[tuple[int, int]] = []

    if detected_type == "genomic":
        # Need to extract CDS from genomic sequence.
        # Prefer genomic references for coordinate transfer.
        template_ref = _select_template_reference(
            closest_nt_ref, references, query_length=len(sequence)
        )

        # If no genomic template available among same-species references,
        # try a related species for exon coordinate transfer.
        cross_species_refs: list[ReferenceAllele] = []
        if template_ref is None or (
            not template_ref.is_genomic and not template_ref.exon_coords
        ):
            cross_gb = _find_cross_species_gb(locus_gb_path, species, locus)
            if cross_gb is not None:
                cross_species_refs = load_reference_alleles(cross_gb)
                cross_genomic = [
                    r for r in cross_species_refs
                    if r.is_genomic and r.exon_coords
                ]
                if cross_genomic:
                    # Find the closest cross-species genomic reference
                    cross_template, _ = _find_closest_nt(
                        sequence, cross_genomic
                    )
                    if cross_template is not None:
                        template_ref = cross_template
                        logger.info(
                            "Using cross-species template %s from %s for "
                            "exon coordinate transfer",
                            cross_template.name, cross_gb.parent.name,
                        )

        if template_ref is not None:
            cds, query_exons = extract_cds_from_genomic(sequence, template_ref)
            if cds:
                logger.info(
                    "Extracted CDS (%d bp, %d exons) from genomic sequence "
                    "using template %s",
                    len(cds), len(query_exons), template_ref.name,
                )
            else:
                warnings.append(
                    f"CDS extraction failed using template {template_ref.name}"
                )
        else:
            warnings.append("No suitable template reference for CDS extraction")

    elif detected_type == "coding":
        # Sequence is already CDS
        cds = sequence
        query_exons = [(0, len(sequence))]

    elif detected_type == "partial":
        # Partial sequence: use what we have
        cds = sequence
        query_exons = [(0, len(sequence))]
        warnings.append("Sequence detected as partial; naming may be imprecise")

    # Step 5: Translate to protein
    protein = ""
    if cds:
        # Try codon_start=1 first
        protein = _translate_cds(cds, codon_start=1)
        if not protein or "*" in protein[:-1]:
            # Try codon_start from the closest reference
            if closest_nt_ref and closest_nt_ref.codon_start != 1:
                protein = _translate_cds(cds, closest_nt_ref.codon_start)
            if not protein:
                warnings.append("CDS could not be translated to a valid protein")

    # Step 6: Classify relationship
    # First check against IPD references
    relationship, best_protein_ref, protein_pct = classify_relationship(
        protein, cds, references
    )

    # If not synonymous to any IPD allele, check against previously-named
    # provisionals in the same batch — two submissions may share a novel
    # protein that doesn't match any IPD allele.
    if relationship != "synonymous" and provisional_refs and protein:
        prov_identical = find_identical_protein(protein, provisional_refs)
        if prov_identical:
            best_prov_cds = _find_closest_cds_among(cds, prov_identical)
            best_protein_ref = best_prov_cds or prov_identical[0]
            relationship = "synonymous"
            protein_pct = 100.0

    # Step 6.5: Check for extension
    # If any IPD reference's nucleotide sequence is contained within the
    # provisional sequence, this is an extension of an existing allele.
    # Two checks:
    #   1. ref full sequence is a substring of query (genomic-in-genomic)
    #   2. ref CDS is a substring of query's extracted CDS (CDS-only refs
    #      vs genomic provisionals — CDS won't be contiguous in genomic
    #      due to introns, but will match after CDS extraction)
    is_extension = False
    extension_ref = None
    query_upper = sequence.upper()
    cds_upper = cds.upper() if cds else ""
    for ref in references:
        if len(ref.sequence) >= len(sequence):
            continue  # provisional must be longer overall
        # Check 1: full sequence substring
        if ref.sequence.upper() in query_upper:
            if extension_ref is None or len(ref.sequence) > len(extension_ref.sequence):
                extension_ref = ref
            continue
        # Check 2: CDS substring (for CDS-only refs vs genomic queries)
        if cds_upper and ref.coding_seq:
            ref_cds = ref.coding_seq.upper()
            if ref_cds in cds_upper:
                if extension_ref is None or len(ref.sequence) > len(extension_ref.sequence):
                    extension_ref = ref

    if extension_ref is not None:
        is_extension = True
        relationship = "extension"
        best_protein_ref = extension_ref
        logger.info(
            "Extension detected: query contains %s (%d bp within %d bp)",
            extension_ref.name, len(extension_ref.sequence), len(sequence),
        )

    # Step 7: Compute difference counts vs closest protein-level match
    cds_diffs = -1
    aa_diffs = -1
    if best_protein_ref:
        if cds and best_protein_ref.coding_seq:
            cds_diffs = _count_diffs(cds, best_protein_ref.coding_seq)
        if protein and best_protein_ref.protein:
            aa_diffs = _count_diffs(protein, best_protein_ref.protein)

    # Step 8: Assign name
    # For the naming, use the protein-level best match
    naming_ref = best_protein_ref if best_protein_ref else closest_nt_ref

    provisional_name = assign_provisional_name(
        species, locus, relationship, naming_ref, existing_names
    )

    return NamingResult(
        provisional_name=provisional_name,
        relationship=relationship,
        closest_nt_allele=closest_nt_ref.name if closest_nt_ref else "",
        closest_nt_identity=nt_identity,
        closest_protein_allele=best_protein_ref.name if best_protein_ref else "",
        closest_protein_identity=protein_pct,
        detected_seq_type=detected_type,
        cds_extracted=bool(cds),
        cds_length=len(cds),
        protein_length=len(protein),
        exon_count=len(query_exons),
        cds_diffs=cds_diffs,
        aa_diffs=aa_diffs,
        is_extension=is_extension,
        extracted_cds=cds,
        extracted_protein=protein,
        warnings=warnings,
    )


def _find_closest_nt(
    query: str,
    references: list[ReferenceAllele],
    mode: str = "global",
) -> tuple[ReferenceAllele | None, float]:
    """Find the closest reference allele by nucleotide alignment.

    Args:
        query: Query nucleotide sequence.
        references: Reference alleles to compare against.
        mode: "global" (default) or "local". Use "local" when comparing
            a genomic query against CDS-only references, since global
            alignment of 3,190 bp vs 1,080 bp gives misleadingly low
            identity.

    Returns:
        (closest_ref, percent_identity).
    """
    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    best_ref = None
    best_pct = 0.0
    query_len = len(query)

    for ref in references:
        target = ref.sequence
        if not target:
            continue
        try:
            score = aligner.score(query, target)
            if mode == "local":
                # For local alignment, normalize against the shorter
                # sequence (the CDS reference) to get meaningful identity.
                norm_len = min(query_len, len(target))
            else:
                norm_len = max(query_len, len(target))
            pct = (score / norm_len) * 100 if norm_len > 0 else 0.0
        except Exception:
            continue

        if pct > best_pct:
            best_pct = pct
            best_ref = ref

    return (best_ref, best_pct)


def _select_template_reference(
    closest_nt: ReferenceAllele | None,
    references: list[ReferenceAllele],
    query_length: int = 0,
) -> ReferenceAllele | None:
    """Select the best reference to use as a template for CDS extraction.

    Prefers genomic references with exon annotations since those allow
    reliable coordinate transfer. The template must be length-compatible
    with the query: a partial genomic reference (e.g. 1,300 bp) cannot
    serve as a template for a full-length genomic query (e.g. 3,200 bp).

    Returns None if no suitable genomic template is available (the caller
    should then try cross-species fallback).

    Priority:
        1. The closest nucleotide match if it is genomic with exon
           annotations and length-compatible.
        2. Any genomic reference in the same allele group that is
           length-compatible.
        3. Any length-compatible genomic reference at the locus.
        4. None — signals caller to try cross-species fallback.
    """
    def _is_length_compatible(ref: ReferenceAllele) -> bool:
        """Check if a reference is length-compatible with the query.

        A template should be at least 60% of the query length to
        provide reliable coordinate transfer.
        """
        if query_length <= 0:
            return True
        return len(ref.sequence) >= query_length * 0.6

    # Priority 1: closest NT match, if genomic with exons and compatible
    if (
        closest_nt
        and closest_nt.is_genomic
        and closest_nt.exon_coords
        and _is_length_compatible(closest_nt)
    ):
        return closest_nt

    # Priority 2: genomic reference in the same group
    if closest_nt:
        group = closest_nt.allele_group
        for ref in references:
            if (
                ref.is_genomic
                and ref.exon_coords
                and ref.allele_group == group
                and _is_length_compatible(ref)
            ):
                return ref

    # Priority 3: any compatible genomic reference at the locus
    for ref in references:
        if ref.is_genomic and ref.exon_coords and _is_length_compatible(ref):
            return ref

    # No suitable genomic template — return None to trigger
    # cross-species fallback in name_provisional_allele()
    return None


# ---------------------------------------------------------------------------
# Batch naming (replaces the current assign_name in provisional.py)
# ---------------------------------------------------------------------------


def name_provisional_batch(
    sequences: list[tuple[str, str]],
    species: str,
    locus: str,
    seq_type: str,
    locus_gb_path: Path,
    existing_manifest: list[dict],
) -> list[NamingResult]:
    """Name a batch of provisional alleles.

    Handles intra-batch synonymous relationships: if two submissions in
    the same batch share an identical protein (both novel to IPD), the
    second is named as a synonymous variant of the first.

    Args:
        sequences: List of (record_id, uppercase_sequence) tuples.
        species: Species prefix.
        locus: Locus name.
        seq_type: "coding", "genomic", or "auto".
        locus_gb_path: Path to the locus GenBank file.
        existing_manifest: Current manifest entries (for collision avoidance).

    Returns:
        List of NamingResult, one per input sequence.
    """
    # Collect all existing provisional names at this locus
    existing_names = [
        e.get("name", "")
        for e in existing_manifest
        if e.get("species") == species and e.get("locus") == locus
    ]

    # Track previously-named provisional alleles for intra-batch
    # synonymous detection
    provisional_refs: list[ReferenceAllele] = []

    results = []
    for rec_id, sequence in sequences:
        result = name_provisional_allele(
            sequence=sequence,
            species=species,
            locus=locus,
            seq_type=seq_type,
            locus_gb_path=locus_gb_path,
            existing_names=existing_names,
            provisional_refs=provisional_refs,
        )
        # Track the newly assigned name so the next iteration avoids collision
        existing_names.append(result.provisional_name)

        # If this sequence has a valid protein, add it to provisional_refs
        # so later sequences can be classified as synonymous to it
        if result.extracted_protein:
            provisional_refs.append(
                ReferenceAllele(
                    name=result.provisional_name,
                    accession="",
                    sequence=sequence,
                    mol_type="genomic DNA" if seq_type == "genomic" else "mRNA",
                    protein=result.extracted_protein,
                    coding_seq=result.extracted_cds,
                    exon_coords=[],
                    intron_coords=[],
                    codon_start=1,
                )
            )

        results.append(result)

        logger.info(
            "Named '%s' as %s (relationship=%s, closest_protein=%s at %.1f%%)",
            rec_id,
            result.provisional_name,
            result.relationship,
            result.closest_protein_allele,
            result.closest_protein_identity,
        )

    return results
