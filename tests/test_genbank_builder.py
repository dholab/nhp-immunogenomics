"""Tests for GenBank record conversion."""

import io

from Bio import SeqIO

from pipeline.genbank_builder import allele_to_seqrecord, write_genbank_file, _infer_product


class TestAlleleToSeqrecord:
    def test_genomic_class1(self, mhc_genomic_record):
        """Full genomic Class I allele with exons, introns, CDS join."""
        sr = allele_to_seqrecord(mhc_genomic_record, "MHC")
        assert sr.id == "NHP01224"
        assert sr.name == "NHP01224"
        assert "Mamu-A1*026:01:01:01" in sr.description
        assert len(sr.seq) == 2874  # Genomic length
        assert sr.annotations["molecule_type"] == "DNA"
        assert sr.annotations["data_file_division"] == "PRI"
        assert sr.annotations["organism"] == "Macaca mulatta"
        assert "Primates" in sr.annotations["taxonomy"]
        assert "Previous designations: Mamu-A*25" in sr.annotations["comment"]

        # Check features
        feature_types = [f.type for f in sr.features]
        assert "source" in feature_types
        assert "gene" in feature_types
        assert "exon" in feature_types
        assert "intron" in feature_types
        assert "CDS" in feature_types

        # Check exon count
        exons = [f for f in sr.features if f.type == "exon"]
        assert len(exons) == 8

        # Check intron count
        introns = [f for f in sr.features if f.type == "intron"]
        assert len(introns) == 7

        # Check CDS is a compound location (join)
        cds = [f for f in sr.features if f.type == "CDS"][0]
        assert hasattr(cds.location, "parts")  # CompoundLocation
        assert cds.qualifiers["gene"] == ["A1"]
        assert cds.qualifiers["allele"] == ["Mamu-A1*026:01:01:01"]
        assert cds.qualifiers["translation"][0].startswith("M")
        assert cds.qualifiers["codon_start"] == ["1"]
        assert cds.qualifiers["product"] == ["MHC class I A1 antigen"]

        # Check source feature
        source = [f for f in sr.features if f.type == "source"][0]
        assert source.qualifiers["mol_type"] == ["genomic DNA"]
        assert source.qualifiers["db_xref"] == ["taxon:9544"]

        # Check gene feature has previous designations
        gene = [f for f in sr.features if f.type == "gene"][0]
        assert "Mamu-A*25" in gene.qualifiers["note"][0]

    def test_coding_only_single_exon(self, mhc_coding_only_record):
        """CDS-only allele with single exon and codon_start=3."""
        sr = allele_to_seqrecord(mhc_coding_only_record, "MHC")
        assert sr.id == "NHP00001"
        assert len(sr.seq) == 224  # Coding length (same as genomic)

        # Should use mRNA since no introns
        source = [f for f in sr.features if f.type == "source"][0]
        assert source.qualifiers["mol_type"] == ["mRNA"]

        # No introns
        introns = [f for f in sr.features if f.type == "intron"]
        assert len(introns) == 0

        # Single exon
        exons = [f for f in sr.features if f.type == "exon"]
        assert len(exons) == 1

        # CDS spans full length, codon_start=3
        cds = [f for f in sr.features if f.type == "CDS"][0]
        assert cds.qualifiers["codon_start"] == ["3"]
        assert cds.qualifiers["product"] == ["MHC class II DQA1 antigen"]

        # Keywords
        assert "IPD-MHC" in sr.annotations["keywords"]
        assert "MHC class II" in sr.annotations["keywords"]

    def test_nhkir_record(self, nhkir_record):
        """NHKIR KIR record with 9 exons."""
        sr = allele_to_seqrecord(nhkir_record, "NHKIR")
        assert sr.id == "NHP00001"
        assert "Mamu-KIR3DL01*013" in sr.description
        assert "IPD-NHKIR" in sr.annotations["keywords"]

        exons = [f for f in sr.features if f.type == "exon"]
        assert len(exons) == 9

        cds = [f for f in sr.features if f.type == "CDS"][0]
        assert "killer cell immunoglobulin-like receptor" in cds.qualifiers["product"][0]

    def test_no_sequence_raises(self):
        """Allele with no sequence data raises ValueError."""
        import pytest
        record = {"accession": "NHP99999", "sequence": {}}
        with pytest.raises(ValueError, match="No sequence data"):
            allele_to_seqrecord(record, "MHC")

    def test_previous_designations_in_comment(self, mhc_genomic_record):
        """Previous designations appear in COMMENT."""
        sr = allele_to_seqrecord(mhc_genomic_record, "MHC")
        comment = sr.annotations["comment"]
        assert "Previous designations:" in comment
        assert "Mamu-A*25" in comment
        assert "Mamu-A1*026:01" in comment
        assert "Mamu-A1*02601" in comment

    def test_no_previous_designations(self, nhkir_record):
        """Record without previous designations omits that line."""
        sr = allele_to_seqrecord(nhkir_record, "NHKIR")
        assert "Previous designations:" not in sr.annotations["comment"]

    def test_references(self, mhc_genomic_record):
        """References are populated with authors, title, journal, PUBMED."""
        sr = allele_to_seqrecord(mhc_genomic_record, "MHC")
        refs = sr.annotations["references"]
        assert len(refs) >= 1
        assert refs[0].authors
        assert refs[0].title

    def test_insdc_crossrefs(self, mhc_genomic_record):
        """INSDC cross-references are in dbxrefs."""
        sr = allele_to_seqrecord(mhc_genomic_record, "MHC")
        assert any("INSDC:" in x for x in sr.dbxrefs)

    def test_genbank_roundtrip(self, mhc_genomic_record):
        """GenBank output parses back correctly."""
        sr = allele_to_seqrecord(mhc_genomic_record, "MHC")
        buf = io.StringIO()
        SeqIO.write(sr, buf, "genbank")
        buf.seek(0)
        parsed = SeqIO.read(buf, "genbank")
        assert parsed.id == "NHP01224"
        assert len(parsed.seq) == 2874
        assert parsed.annotations["organism"] == "Macaca mulatta"

    def test_coding_only_roundtrip(self, mhc_coding_only_record):
        """Coding-only record roundtrips through GenBank format."""
        sr = allele_to_seqrecord(mhc_coding_only_record, "MHC")
        buf = io.StringIO()
        SeqIO.write(sr, buf, "genbank")
        buf.seek(0)
        parsed = SeqIO.read(buf, "genbank")
        assert parsed.id == "NHP00001"
        assert len(parsed.seq) == 224


class TestWriteGenbankFile:
    def test_write_multiple(self, mhc_genomic_record, mhc_coding_only_record, tmp_path):
        """Write multi-record GenBank file."""
        output = tmp_path / "test.gb"
        n = write_genbank_file(
            [mhc_genomic_record, mhc_coding_only_record],
            output,
            "MHC",
        )
        assert n == 2
        assert output.exists()

        # Parse back all records
        records = list(SeqIO.parse(output, "genbank"))
        assert len(records) == 2
        assert records[0].id == "NHP01224"
        assert records[1].id == "NHP00001"

    def test_skip_bad_records(self, tmp_path):
        """Records with no sequence are skipped."""
        good = {"accession": "NHP00001", "name": "Test", "locus": "A",
                "sequence": {"coding": "atgccc", "genomic": "atgccc"},
                "feature": [], "organism": {"scientificName": "Test sp"}}
        bad = {"accession": "NHP99999", "sequence": {}}
        output = tmp_path / "test.gb"
        n = write_genbank_file([good, bad], output, "MHC")
        assert n == 1

    def test_creates_parent_dirs(self, tmp_path):
        """Output directories are created automatically."""
        output = tmp_path / "a" / "b" / "c" / "test.gb"
        record = {"accession": "NHP00001", "name": "Test", "locus": "A",
                  "sequence": {"coding": "atg", "genomic": "atg"},
                  "feature": [], "organism": {"scientificName": "Test sp"}}
        write_genbank_file([record], output, "MHC")
        assert output.exists()


class TestInferProduct:
    def test_class_i(self):
        assert _infer_product("A1", "I", "MHC") == "MHC class I A1 antigen"

    def test_class_ii(self):
        assert _infer_product("DRB", "II", "MHC") == "MHC class II DRB antigen"

    def test_unknown_class(self):
        assert _infer_product("X", "unknown", "MHC") == "MHC X antigen"

    def test_nhkir(self):
        result = _infer_product("KIR3DL01", "unknown", "NHKIR")
        assert "killer cell immunoglobulin-like receptor" in result
