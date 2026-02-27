"""Tests for the improved provisional allele naming module."""

import pytest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio import SeqIO

from pipeline.namer import (
    ReferenceAllele,
    NamingResult,
    detect_sequence_type,
    load_reference_alleles,
    transfer_exon_coordinates,
    extract_cds_from_genomic,
    compare_proteins,
    find_identical_protein,
    classify_relationship,
    assign_provisional_name,
    name_provisional_allele,
    name_provisional_batch,
    _has_clean_orf,
    _translate_cds,
    _find_closest_nt,
)


# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------


def _make_reference(
    name: str,
    sequence: str,
    protein: str = "",
    coding_seq: str = "",
    mol_type: str = "mRNA",
    exon_coords: list[tuple[int, int]] | None = None,
    intron_coords: list[tuple[int, int]] | None = None,
    codon_start: int = 1,
) -> ReferenceAllele:
    """Helper to create a ReferenceAllele for testing."""
    if exon_coords is None:
        exon_coords = [(0, len(sequence))]
    if intron_coords is None:
        intron_coords = []
    if not coding_seq:
        coding_seq = "".join(sequence[s:e] for s, e in exon_coords)
    if not protein and coding_seq:
        protein = _translate_cds(coding_seq, codon_start)
    return ReferenceAllele(
        name=name,
        accession="NHP00000",
        sequence=sequence,
        mol_type=mol_type,
        protein=protein,
        coding_seq=coding_seq,
        exon_coords=exon_coords,
        intron_coords=intron_coords,
        codon_start=codon_start,
    )


def _write_genbank_with_features(tmp_path, species, locus, alleles):
    """Write a GenBank file with proper exon/CDS features for testing.

    alleles: list of dicts with keys:
        name, sequence, exon_coords, protein, mol_type
    """
    gb_dir = tmp_path / "data" / "MHC" / species
    gb_dir.mkdir(parents=True, exist_ok=True)
    records = []

    for allele in alleles:
        seq_str = allele["sequence"]
        rec = SeqRecord(
            Seq(seq_str),
            id="NHP00000",
            name="NHP00000",
            description=f"{allele['name']}, {locus} locus allele",
        )
        rec.annotations["molecule_type"] = "DNA"

        mol_type = allele.get("mol_type", "mRNA")

        # Source feature
        rec.features.append(SeqFeature(
            location=FeatureLocation(0, len(seq_str)),
            type="source",
            qualifiers={"mol_type": [mol_type], "organism": ["Test"]},
        ))

        # Exon features
        exon_coords = allele.get("exon_coords", [(0, len(seq_str))])
        for i, (start, end) in enumerate(exon_coords, 1):
            rec.features.append(SeqFeature(
                location=FeatureLocation(start, end),
                type="exon",
                qualifiers={"gene": [locus], "number": [str(i)]},
            ))

        # Intron features (if genomic)
        if "genomic" in mol_type:
            for i in range(len(exon_coords) - 1):
                intron_start = exon_coords[i][1]
                intron_end = exon_coords[i + 1][0]
                rec.features.append(SeqFeature(
                    location=FeatureLocation(intron_start, intron_end),
                    type="intron",
                    qualifiers={"gene": [locus], "number": [str(i + 1)]},
                ))

        # CDS feature
        cds_qualifiers = {
            "gene": [locus],
            "allele": [allele["name"]],
            "codon_start": ["1"],
        }
        if allele.get("protein"):
            cds_qualifiers["translation"] = [allele["protein"]]

        if len(exon_coords) > 1 and "genomic" in mol_type:
            locs = [FeatureLocation(s, e) for s, e in exon_coords]
            cds_location = CompoundLocation(locs)
        else:
            cds_location = FeatureLocation(0, len(seq_str))

        rec.features.append(SeqFeature(
            location=cds_location,
            type="CDS",
            qualifiers=cds_qualifiers,
        ))

        records.append(rec)

    gb_path = gb_dir / f"{locus}.gb"
    with open(gb_path, "w") as fh:
        SeqIO.write(records, fh, "genbank")
    return gb_path


# ---------------------------------------------------------------------------
# ReferenceAllele property tests
# ---------------------------------------------------------------------------


class TestReferenceAllele:
    def test_allele_group(self):
        ref = _make_reference("Mamu-E*02:01:01", "ATGCCC")
        assert ref.allele_group == "02"

    def test_allele_group_three_digit(self):
        ref = _make_reference("Mamu-A1*026:01:01:01", "ATGCCC")
        assert ref.allele_group == "026"

    def test_allele_fields(self):
        ref = _make_reference("Mamu-E*02:01:03:01", "ATGCCC")
        assert ref.allele_fields == ["02", "01", "03", "01"]

    def test_allele_fields_short(self):
        ref = _make_reference("Mamu-E*02:01", "ATGCCC")
        assert ref.allele_fields == ["02", "01"]

    def test_protein_field(self):
        ref = _make_reference("Mamu-E*02:01:03:01", "ATGCCC")
        assert ref.protein_field == "02:01"

    def test_protein_field_short(self):
        ref = _make_reference("Mamu-E*02:01", "ATGCCC")
        assert ref.protein_field == "02:01"

    def test_is_genomic(self):
        ref = _make_reference("Test", "ATGCCC", mol_type="genomic DNA")
        assert ref.is_genomic is True

    def test_is_not_genomic(self):
        ref = _make_reference("Test", "ATGCCC", mol_type="mRNA")
        assert ref.is_genomic is False


# ---------------------------------------------------------------------------
# Sequence type detection
# ---------------------------------------------------------------------------


class TestDetectSequenceType:
    def test_cds_length_sequence(self):
        """A ~1080 bp sequence with clean ORF should be classified as coding."""
        # Build a fake CDS-length sequence with ATG...stop
        cds = "ATG" + "GCA" * 358 + "TAA"  # 1080 bp
        result = detect_sequence_type(cds, expected_cds_lengths=[1080])
        assert result == "coding"

    def test_genomic_length_sequence(self):
        """A ~2900 bp sequence should be classified as genomic."""
        seq = "A" * 2900
        result = detect_sequence_type(
            seq,
            expected_cds_lengths=[1080],
            expected_genomic_lengths=[2900],
        )
        assert result == "genomic"

    def test_short_partial(self):
        """A very short sequence should be classified as partial."""
        seq = "ATGCCC"
        result = detect_sequence_type(seq, expected_cds_lengths=[1080])
        assert result == "partial"

    def test_sequence_with_intron_signatures(self):
        """A sequence with multiple GT...AG patterns at intron distances."""
        # Build a fake genomic sequence with intron-like patterns
        exon1 = "ATG" + "GCA" * 20  # 63 bp
        intron1 = "GT" + "A" * 100 + "AG"  # 104 bp
        exon2 = "GCA" * 20  # 60 bp
        intron2 = "GT" + "C" * 100 + "AG"  # 104 bp
        exon3 = "GCA" * 20  # 60 bp
        intron3 = "GT" + "T" * 100 + "AG"  # 104 bp
        exon4 = "GCA" * 10 + "TAA"  # 33 bp

        seq = exon1 + intron1 + exon2 + intron2 + exon3 + intron3 + exon4
        result = detect_sequence_type(seq, expected_cds_lengths=[216])
        assert result == "genomic"


class TestHasCleanOrf:
    def test_clean_orf(self):
        cds = "ATG" + "GCA" * 100 + "TAA"
        assert _has_clean_orf(cds) is True

    def test_no_atg(self):
        assert _has_clean_orf("GCAGCAGCA") is False

    def test_internal_stop(self):
        # Insert a stop codon in the middle
        cds = "ATG" + "GCA" * 50 + "TAA" + "GCA" * 50 + "TAA"
        assert _has_clean_orf(cds) is False

    def test_too_short(self):
        assert _has_clean_orf("ATGGCATAA") is False


# ---------------------------------------------------------------------------
# Reference allele loading
# ---------------------------------------------------------------------------


class TestLoadReferenceAlleles:
    def test_loads_from_genbank(self, tmp_path):
        # Build a simple CDS-only allele
        cds = "ATG" + "GCA" * 100 + "TAA"
        protein = str(Seq(cds).translate())[:-1]  # strip stop

        gb_path = _write_genbank_with_features(tmp_path, "Mamu", "E", [{
            "name": "Mamu-E*02:01",
            "sequence": cds,
            "protein": protein,
            "mol_type": "mRNA",
        }])

        refs = load_reference_alleles(gb_path)
        assert len(refs) == 1
        assert refs[0].name == "Mamu-E*02:01"
        assert refs[0].protein == protein
        assert refs[0].coding_seq == cds

    def test_missing_file(self, tmp_path):
        refs = load_reference_alleles(tmp_path / "nonexistent.gb")
        assert refs == []

    def test_genomic_allele_exon_parsing(self, tmp_path):
        # Build a genomic sequence with exon/intron structure
        exon1 = "ATG" + "GCA" * 10  # 33 bp
        intron1 = "GT" + "A" * 100 + "AG"
        exon2 = "GCA" * 10 + "TAA"  # 33 bp

        genomic = exon1 + intron1 + exon2
        exon_coords = [(0, 33), (33 + 104, 33 + 104 + 33)]
        coding = exon1 + exon2
        protein = str(Seq(coding).translate())[:-1]

        gb_path = _write_genbank_with_features(tmp_path, "Mamu", "E", [{
            "name": "Mamu-E*02:01",
            "sequence": genomic,
            "protein": protein,
            "mol_type": "genomic DNA",
            "exon_coords": exon_coords,
        }])

        refs = load_reference_alleles(gb_path)
        assert len(refs) == 1
        ref = refs[0]
        assert ref.is_genomic
        assert len(ref.exon_coords) == 2
        assert ref.exon_coords[0] == (0, 33)
        assert ref.exon_coords[1] == exon_coords[1]
        assert ref.coding_seq == coding


# ---------------------------------------------------------------------------
# Protein comparison
# ---------------------------------------------------------------------------


class TestCompareProteins:
    def test_identical_match(self):
        refs = [
            _make_reference("Mamu-E*02:01", "ATG", protein="MAVMAPRTL"),
            _make_reference("Mamu-E*02:02", "ATG", protein="MAVMAPRTX"),
        ]
        best, pct = compare_proteins("MAVMAPRTL", refs)
        assert best is not None
        assert best.name == "Mamu-E*02:01"
        assert pct == 100.0

    def test_closest_match(self):
        refs = [
            _make_reference("Mamu-E*02:01", "ATG", protein="MAVMAPR"),
            _make_reference("Mamu-E*02:02", "ATG", protein="XXXXXXX"),
        ]
        best, pct = compare_proteins("MAVMAPX", refs)
        assert best is not None
        assert best.name == "Mamu-E*02:01"
        assert pct > 50.0

    def test_no_protein_refs(self):
        refs = [_make_reference("Test", "ATG", protein="")]
        best, pct = compare_proteins("MAVMAPR", refs)
        assert best is None

    def test_empty_query(self):
        refs = [_make_reference("Test", "ATG", protein="MAVMAPR")]
        best, pct = compare_proteins("", refs)
        assert best is None


class TestFindIdenticalProtein:
    def test_exact_match(self):
        refs = [
            _make_reference("Mamu-E*02:01", "ATG", protein="MAVMAPR"),
            _make_reference("Mamu-E*02:02", "ATG", protein="MAVMAPX"),
        ]
        matches = find_identical_protein("MAVMAPR", refs)
        assert len(matches) == 1
        assert matches[0].name == "Mamu-E*02:01"

    def test_multiple_matches(self):
        refs = [
            _make_reference("Mamu-E*02:01:01", "ATG", protein="MAVMAPR"),
            _make_reference("Mamu-E*02:01:02", "ATG", protein="MAVMAPR"),
        ]
        matches = find_identical_protein("MAVMAPR", refs)
        assert len(matches) == 2

    def test_trailing_stop_ignored(self):
        refs = [_make_reference("Test", "ATG", protein="MAVMAPR*")]
        matches = find_identical_protein("MAVMAPR", refs)
        assert len(matches) == 1

    def test_no_match(self):
        refs = [_make_reference("Test", "ATG", protein="XXXXXXX")]
        matches = find_identical_protein("MAVMAPR", refs)
        assert len(matches) == 0


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------


class TestClassifyRelationship:
    def test_synonymous(self):
        """Identical protein -> synonymous."""
        refs = [
            _make_reference(
                "Mamu-E*02:01",
                "ATGGCAGTTTGA",
                protein="MAV",
                coding_seq="ATGGCAGTTTGA",
            ),
        ]
        rel, best, pct = classify_relationship("MAV", "ATGGCAGTTTGA", refs)
        assert rel == "synonymous"
        assert best.name == "Mamu-E*02:01"
        assert pct == 100.0

    def test_non_synonymous(self):
        """Different protein -> non_synonymous."""
        refs = [
            _make_reference("Mamu-E*02:01", "ATG", protein="MAVMAPR"),
        ]
        rel, best, pct = classify_relationship("MAVMAPX", "ATGGCG", refs)
        assert rel == "non_synonymous"
        assert best.name == "Mamu-E*02:01"
        assert pct < 100.0

    def test_fallback_no_protein(self):
        """No query protein -> fallback."""
        refs = [_make_reference("Test", "ATG", protein="MAVMAPR")]
        rel, best, pct = classify_relationship("", "ATGGCG", refs)
        assert rel == "fallback"

    def test_fallback_no_refs(self):
        """No reference proteins -> fallback."""
        refs = [_make_reference("Test", "ATG", protein="")]
        rel, best, pct = classify_relationship("MAVMAPR", "ATGGCG", refs)
        assert rel == "fallback"


# ---------------------------------------------------------------------------
# Name assignment
# ---------------------------------------------------------------------------


class TestAssignProvisionalName:
    def test_synonymous_naming(self):
        ref = _make_reference("Mamu-E*02:01", "ATG")
        name = assign_provisional_name("Mamu", "E", "synonymous", ref, [])
        assert name == "Mamu-E*02:01:new01"

    def test_synonymous_naming_four_field_ref(self):
        ref = _make_reference("Mamu-E*02:01:03:01", "ATG")
        name = assign_provisional_name("Mamu", "E", "synonymous", ref, [])
        # Should use first two fields (protein designation)
        assert name == "Mamu-E*02:01:new01"

    def test_non_synonymous_naming(self):
        ref = _make_reference("Mamu-E*02:01", "ATG")
        name = assign_provisional_name("Mamu", "E", "non_synonymous", ref, [])
        assert name == "Mamu-E*02:new01"

    def test_fallback_naming(self):
        name = assign_provisional_name("Mamu", "E", "fallback", None, [])
        assert name == "Mamu-E*000:new01"

    def test_collision_avoidance(self):
        ref = _make_reference("Mamu-E*02:01", "ATG")
        existing = ["Mamu-E*02:01:new01", "Mamu-E*02:01:new02"]
        name = assign_provisional_name("Mamu", "E", "synonymous", ref, existing)
        assert name == "Mamu-E*02:01:new03"

    def test_non_synonymous_collision_avoidance(self):
        ref = _make_reference("Mamu-E*02:01", "ATG")
        existing = ["Mamu-E*02:new01"]
        name = assign_provisional_name("Mamu", "E", "non_synonymous", ref, existing)
        assert name == "Mamu-E*02:new02"

    def test_kir_naming(self):
        ref = _make_reference("Mamu-KIR3DS03*001:01", "ATG")
        name = assign_provisional_name("Mamu", "KIR3DS03", "synonymous", ref, [])
        assert name == "Mamu-KIR3DS03*001:01:new01"

    def test_kir_non_synonymous(self):
        ref = _make_reference("Mamu-KIR3DS03*001:01", "ATG")
        name = assign_provisional_name("Mamu", "KIR3DS03", "non_synonymous", ref, [])
        assert name == "Mamu-KIR3DS03*001:new01"


# ---------------------------------------------------------------------------
# CDS extraction (integration-style tests)
# ---------------------------------------------------------------------------


class TestExtractCdsFromGenomic:
    def test_with_genomic_reference(self):
        """Test CDS extraction when reference is genomic."""
        # Construct a simple 2-exon gene
        exon1 = "ATGGCAGCA"  # 9 bp
        intron = "GT" + "A" * 80 + "AG"  # 84 bp
        exon2 = "GCAGCATGA"  # 9 bp (includes stop)

        ref_genomic = exon1 + intron + exon2
        ref_exons = [(0, 9), (93, 102)]

        ref = _make_reference(
            "Mamu-E*02:01",
            ref_genomic,
            mol_type="genomic DNA",
            exon_coords=ref_exons,
        )

        # Query: same structure, slightly different sequence
        q_exon1 = "ATGGCAGCC"  # 1 synonymous change
        q_intron = "GT" + "A" * 80 + "AG"  # same intron
        q_exon2 = "GCAGCATGA"  # same exon2
        query = q_exon1 + q_intron + q_exon2

        cds, exons = extract_cds_from_genomic(query, ref)

        assert cds != ""
        assert len(exons) == 2
        # CDS should be exon1 + exon2
        expected_cds = q_exon1 + q_exon2
        assert cds == expected_cds

    def test_no_exon_coords_fallback_to_orf(self):
        """When reference has no exon coords, should attempt ORF extraction."""
        ref = _make_reference(
            "Mamu-E*02:01",
            "ATGGCAGCATGA",  # simple CDS
            exon_coords=[],
        )

        # Query is genomic but reference has no exons
        query = "NNNNATGGCAGCATGANNN"
        cds, exons = extract_cds_from_genomic(query, ref)

        # Should fall back to ORF detection
        # The longest ORF starting with ATG should be found
        assert "ATG" in cds


# ---------------------------------------------------------------------------
# Full naming pipeline
# ---------------------------------------------------------------------------


class TestNameProvisionalAllele:
    def test_coding_submission_identical_protein(self, tmp_path):
        """CDS submission that encodes same protein as reference -> synonymous."""
        # Reference: Mamu-E*02:01 with protein MAV
        cds_ref = "ATGGCAGTTTGA"  # MAV*
        protein = "MAV"

        gb_path = _write_genbank_with_features(tmp_path, "Mamu", "E", [{
            "name": "Mamu-E*02:01",
            "sequence": cds_ref,
            "protein": protein,
            "mol_type": "mRNA",
        }])

        # Submission: same protein, different synonymous codon
        cds_query = "ATGGCCGTCTGA"  # MAV* (GCA->GCC synonymous change in A)
        # Actually: ATG GCC GTC TGA -> M A V * ... wait, let me fix
        # ATG = M, GCC = A, GTC = V, TGA = stop -> MAV
        # Reference: ATG GCA GTT TGA -> M A V stop -> MAV
        # Different codons, same protein

        result = name_provisional_allele(
            sequence=cds_query,
            species="Mamu",
            locus="E",
            seq_type="coding",
            locus_gb_path=gb_path,
            existing_names=[],
        )

        assert result.relationship == "synonymous"
        assert "02:01:new" in result.provisional_name
        assert result.closest_protein_identity == 100.0

    def test_coding_submission_different_protein(self, tmp_path):
        """CDS submission with different protein -> non_synonymous."""
        cds_ref = "ATGGCAGTTTGA"  # MAV*
        protein = "MAV"

        gb_path = _write_genbank_with_features(tmp_path, "Mamu", "E", [{
            "name": "Mamu-E*02:01",
            "sequence": cds_ref,
            "protein": protein,
            "mol_type": "mRNA",
        }])

        # Submission: different protein
        cds_query = "ATGCAAGTTTGA"  # MQV* (different second AA)

        result = name_provisional_allele(
            sequence=cds_query,
            species="Mamu",
            locus="E",
            seq_type="coding",
            locus_gb_path=gb_path,
            existing_names=[],
        )

        assert result.relationship == "non_synonymous"
        assert "02:new" in result.provisional_name
        assert result.closest_protein_identity < 100.0

    def test_no_references(self, tmp_path):
        """When no GenBank file exists, should use fallback."""
        result = name_provisional_allele(
            sequence="ATGGCAGTTTGA",
            species="Mamu",
            locus="Z9",
            seq_type="coding",
            locus_gb_path=tmp_path / "nonexistent.gb",
            existing_names=[],
        )

        assert result.relationship == "fallback"
        assert "000:new" in result.provisional_name


class TestNameProvisionalBatch:
    def test_batch_increments(self, tmp_path):
        """Multiple sequences in a batch should get incrementing numbers."""
        cds_ref = "ATGGCAGTTTGA"
        protein = "MAV"

        gb_path = _write_genbank_with_features(tmp_path, "Mamu", "E", [{
            "name": "Mamu-E*02:01",
            "sequence": cds_ref,
            "protein": protein,
            "mol_type": "mRNA",
        }])

        sequences = [
            ("seq1", "ATGGCCGTCTGA"),  # MAV (synonymous)
            ("seq2", "ATGGCTGTGTGA"),  # MAV (synonymous, different codons)
        ]

        results = name_provisional_batch(
            sequences=sequences,
            species="Mamu",
            locus="E",
            seq_type="coding",
            locus_gb_path=gb_path,
            existing_manifest=[],
        )

        assert len(results) == 2
        # Both should be synonymous variants
        names = [r.provisional_name for r in results]
        # They should have different numbers
        assert names[0] != names[1]
        # Both should reference 02:01
        for n in names:
            assert "02:01:new" in n or "02:new" in n


# ---------------------------------------------------------------------------
# Translation helper
# ---------------------------------------------------------------------------


class TestTranslateCds:
    def test_standard_translation(self):
        cds = "ATGGCAGTTTGA"  # MAV*
        protein = _translate_cds(cds)
        assert protein == "MAV"

    def test_codon_start_offset(self):
        # codon_start=2 means skip first base
        cds = "XATGGCATGA"  # skip X, then ATG GCA TGA -> MA*
        protein = _translate_cds(cds, codon_start=2)
        assert protein == "MA"

    def test_empty_input(self):
        assert _translate_cds("") == ""

    def test_strips_trailing_stop(self):
        cds = "ATGTGA"  # M*
        protein = _translate_cds(cds)
        assert protein == "M"
