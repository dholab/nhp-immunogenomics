"""Tests for provisional allele management."""

import csv
import json

from Bio import SeqIO

from pipeline.genbank_builder import provisional_to_seqrecord
import pytest

from pipeline.provisional import (
    add_provisional_alleles,
    append_to_fasta,
    append_to_manifest,
    assign_name,
    build_genbank,
    build_metadata,
    build_sequence_index,
    check_retirements,
    compute_accession,
    extract_allele_group,
    find_nearest_neighbor,
    load_manifest,
    load_sequences,
    retire_alleles,
    sequence_hash,
    validate,
    validate_no_ipd_collisions,
)


def _write_manifest(tmp_path, rows):
    """Helper to write a manifest TSV for testing."""
    manifest_path = tmp_path / "provisional" / "manifest.tsv"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(manifest_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "name", "species", "locus", "class", "seq_type",
            "sequence_file", "submitter", "date_added", "notes",
        ])
        for row in rows:
            writer.writerow(row)


def _write_fasta(tmp_path, rel_path, records):
    """Helper to write a FASTA file for testing."""
    fasta_path = tmp_path / "provisional" / "sequences" / rel_path
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    with open(fasta_path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")


class TestLoadManifest:
    def test_empty_manifest(self, tmp_path):
        _write_manifest(tmp_path, [])
        result = load_manifest(tmp_path)
        assert result == []

    def test_missing_manifest(self, tmp_path):
        result = load_manifest(tmp_path)
        assert result == []

    def test_loads_entries(self, tmp_path):
        _write_manifest(tmp_path, [
            ["Mamu-A1*026:new01", "Mamu", "A1", "I", "coding",
             "Mamu/A1.fasta", "J. Karl", "2026-02-15", "test allele"],
        ])
        result = load_manifest(tmp_path)
        assert len(result) == 1
        assert result[0]["name"] == "Mamu-A1*026:new01"
        assert result[0]["species"] == "Mamu"
        assert result[0]["locus"] == "A1"
        assert result[0]["class"] == "I"
        assert result[0]["seq_type"] == "coding"
        assert result[0]["submitter"] == "J. Karl"

    def test_multiple_entries(self, tmp_path):
        _write_manifest(tmp_path, [
            ["Mamu-A1*026:new01", "Mamu", "A1", "I", "coding",
             "Mamu/A1.fasta", "J. Karl", "2026-02-15", ""],
            ["Mamu-B*017:new01", "Mamu", "B", "I", "genomic",
             "Mamu/B.fasta", "J. Karl", "2026-02-15", ""],
        ])
        result = load_manifest(tmp_path)
        assert len(result) == 2


class TestLoadSequences:
    def test_loads_fasta(self, tmp_path):
        _write_manifest(tmp_path, [
            ["Mamu-A1*026:new01", "Mamu", "A1", "I", "coding",
             "Mamu/A1.fasta", "J. Karl", "2026-02-15", ""],
        ])
        _write_fasta(tmp_path, "Mamu/A1.fasta", [
            ("Mamu-A1*026:new01", "ATGGCACCCTGAAC"),
        ])
        manifest = load_manifest(tmp_path)
        seqs = load_sequences(tmp_path, manifest)
        assert "Mamu-A1*026:new01" in seqs
        assert seqs["Mamu-A1*026:new01"] == "ATGGCACCCTGAAC"

    def test_uppercase(self, tmp_path):
        _write_manifest(tmp_path, [
            ["Mamu-A1*026:new01", "Mamu", "A1", "I", "coding",
             "Mamu/A1.fasta", "J. Karl", "2026-02-15", ""],
        ])
        _write_fasta(tmp_path, "Mamu/A1.fasta", [
            ("Mamu-A1*026:new01", "atggcaccctgaac"),
        ])
        manifest = load_manifest(tmp_path)
        seqs = load_sequences(tmp_path, manifest)
        assert seqs["Mamu-A1*026:new01"] == "ATGGCACCCTGAAC"

    def test_missing_fasta(self, tmp_path):
        _write_manifest(tmp_path, [
            ["Mamu-A1*026:new01", "Mamu", "A1", "I", "coding",
             "Mamu/A1.fasta", "J. Karl", "2026-02-15", ""],
        ])
        manifest = load_manifest(tmp_path)
        seqs = load_sequences(tmp_path, manifest)
        assert seqs == {}

    def test_missing_allele_in_fasta(self, tmp_path):
        _write_manifest(tmp_path, [
            ["Mamu-A1*026:new01", "Mamu", "A1", "I", "coding",
             "Mamu/A1.fasta", "J. Karl", "2026-02-15", ""],
        ])
        _write_fasta(tmp_path, "Mamu/A1.fasta", [
            ("Mamu-A1*999:new01", "ATGGCACCCTGAAC"),  # Wrong name
        ])
        manifest = load_manifest(tmp_path)
        seqs = load_sequences(tmp_path, manifest)
        assert "Mamu-A1*026:new01" not in seqs

    def test_multi_record_fasta(self, tmp_path):
        _write_manifest(tmp_path, [
            ["Mamu-A1*026:new01", "Mamu", "A1", "I", "coding",
             "Mamu/A1.fasta", "J. Karl", "2026-02-15", ""],
            ["Mamu-A1*026:new02", "Mamu", "A1", "I", "coding",
             "Mamu/A1.fasta", "J. Karl", "2026-02-15", ""],
        ])
        _write_fasta(tmp_path, "Mamu/A1.fasta", [
            ("Mamu-A1*026:new01", "ATGGCACCCTGAAC"),
            ("Mamu-A1*026:new02", "ATGGCACCCTGAAT"),
        ])
        manifest = load_manifest(tmp_path)
        seqs = load_sequences(tmp_path, manifest)
        assert len(seqs) == 2


class TestExtractAlleleGroup:
    def test_standard_mhc(self):
        assert extract_allele_group("Mamu-A1*026:01:01:01") == "026"

    def test_short_name(self):
        assert extract_allele_group("Mamu-B*017:01") == "017"

    def test_kir(self):
        assert extract_allele_group("Mafa-KIR3DL01*001:01") == "001"

    def test_no_group(self):
        assert extract_allele_group("SomethingWeird") == ""


class TestAssignName:
    def test_basic_assignment(self):
        name = assign_name("Mamu", "A1", "Mamu-A1*026:01:01:01", [])
        assert name == "Mamu-A1*026:new01"

    def test_increments_on_existing(self):
        manifest = [{"name": "Mamu-A1*026:new01"}]
        name = assign_name("Mamu", "A1", "Mamu-A1*026:01:01:01", manifest)
        assert name == "Mamu-A1*026:new02"

    def test_multiple_existing(self):
        manifest = [
            {"name": "Mamu-A1*026:new01"},
            {"name": "Mamu-A1*026:new02"},
            {"name": "Mamu-A1*026:new03"},
        ]
        name = assign_name("Mamu", "A1", "Mamu-A1*026:01:01:01", manifest)
        assert name == "Mamu-A1*026:new04"

    def test_different_group_no_conflict(self):
        manifest = [{"name": "Mamu-A1*026:new01"}]
        name = assign_name("Mamu", "A1", "Mamu-A1*017:01", manifest)
        assert name == "Mamu-A1*017:new01"

    def test_fallback_group(self):
        name = assign_name("Mamu", "A1", "NoGroupHere", [])
        assert name == "Mamu-A1*000:new01"


class TestFindNearestNeighbor:
    def test_finds_exact_match(self, tmp_path):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        gb_dir = tmp_path / "data" / "MHC" / "Mamu"
        gb_dir.mkdir(parents=True)
        rec = SeqRecord(Seq("ATGGCACCCTGAAC"), id="NHP01224", name="NHP01224",
                        description="Mamu-A1*026:01:01:01, A1 locus allele")
        rec.annotations["molecule_type"] = "DNA"
        from Bio import SeqIO
        with open(gb_dir / "A1.gb", "w") as fh:
            SeqIO.write([rec], fh, "genbank")

        name, pct = find_nearest_neighbor("ATGGCACCCTGAAC", gb_dir / "A1.gb")
        assert name == "Mamu-A1*026:01:01:01"
        assert pct > 99.0

    def test_finds_closest_among_multiple(self, tmp_path):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        gb_dir = tmp_path / "data" / "MHC" / "Mamu"
        gb_dir.mkdir(parents=True)
        recs = [
            SeqRecord(Seq("ATGGCACCCTGAAC"), id="NHP001", name="NHP001",
                      description="Mamu-A1*026:01, A1"),
            SeqRecord(Seq("ATGGCACCCTGAAT"), id="NHP002", name="NHP002",
                      description="Mamu-A1*017:01, A1"),
        ]
        for r in recs:
            r.annotations["molecule_type"] = "DNA"
        from Bio import SeqIO
        with open(gb_dir / "A1.gb", "w") as fh:
            SeqIO.write(recs, fh, "genbank")

        # Query is closer to 026 (differs only in last base from 017)
        name, pct = find_nearest_neighbor("ATGGCACCCTGAAC", gb_dir / "A1.gb")
        assert "026" in name

    def test_missing_gb_file(self, tmp_path):
        name, pct = find_nearest_neighbor("ATGGC", tmp_path / "nonexistent.gb")
        assert name == ""
        assert pct == 0.0


class TestComputeAccession:
    def test_deterministic(self):
        seq = "ATGGCACCCTGAAC"
        assert compute_accession(seq) == compute_accession(seq)

    def test_case_insensitive(self):
        assert compute_accession("ATGGC") == compute_accession("atggc")

    def test_format(self):
        acc = compute_accession("ATGGCACCCTGAAC")
        assert acc.startswith("PROV")
        assert len(acc) == 9  # PROV + 5 hex chars

    def test_different_sequences(self):
        assert compute_accession("ATGGC") != compute_accession("ATGGA")


class TestValidate:
    def test_empty_valid(self, tmp_path):
        _write_manifest(tmp_path, [])
        assert validate(tmp_path) == []

    def test_valid_entry(self, tmp_path):
        _write_manifest(tmp_path, [
            ["Mamu-A1*026:new01", "Mamu", "A1", "I", "coding",
             "Mamu/A1.fasta", "J. Karl", "2026-02-15", ""],
        ])
        _write_fasta(tmp_path, "Mamu/A1.fasta", [
            ("Mamu-A1*026:new01", "ATGGCACCCTGAAC"),
        ])
        assert validate(tmp_path) == []

    def test_missing_required_field(self, tmp_path):
        _write_manifest(tmp_path, [
            ["", "Mamu", "A1", "I", "coding",
             "Mamu/A1.fasta", "J. Karl", "2026-02-15", ""],
        ])
        _write_fasta(tmp_path, "Mamu/A1.fasta", [])
        errors = validate(tmp_path)
        assert any("missing required field 'name'" in e for e in errors)

    def test_invalid_seq_type(self, tmp_path):
        _write_manifest(tmp_path, [
            ["Mamu-A1*026:new01", "Mamu", "A1", "I", "protein",
             "Mamu/A1.fasta", "J. Karl", "2026-02-15", ""],
        ])
        _write_fasta(tmp_path, "Mamu/A1.fasta", [
            ("Mamu-A1*026:new01", "ATGGCACCCTGAAC"),
        ])
        errors = validate(tmp_path)
        assert any("seq_type must be 'coding' or 'genomic'" in e for e in errors)

    def test_invalid_date(self, tmp_path):
        _write_manifest(tmp_path, [
            ["Mamu-A1*026:new01", "Mamu", "A1", "I", "coding",
             "Mamu/A1.fasta", "J. Karl", "02-15-2026", ""],
        ])
        _write_fasta(tmp_path, "Mamu/A1.fasta", [
            ("Mamu-A1*026:new01", "ATGGCACCCTGAAC"),
        ])
        errors = validate(tmp_path)
        assert any("date_added must be YYYY-MM-DD" in e for e in errors)

    def test_missing_sequence(self, tmp_path):
        _write_manifest(tmp_path, [
            ["Mamu-A1*026:new01", "Mamu", "A1", "I", "coding",
             "Mamu/A1.fasta", "J. Karl", "2026-02-15", ""],
        ])
        # No FASTA file created
        errors = validate(tmp_path)
        assert any("not found in FASTA" in e for e in errors)

    def test_invalid_nucleotides(self, tmp_path):
        _write_manifest(tmp_path, [
            ["Mamu-A1*026:new01", "Mamu", "A1", "I", "coding",
             "Mamu/A1.fasta", "J. Karl", "2026-02-15", ""],
        ])
        _write_fasta(tmp_path, "Mamu/A1.fasta", [
            ("Mamu-A1*026:new01", "ATGXYZ123"),
        ])
        errors = validate(tmp_path)
        assert any("invalid nucleotide" in e for e in errors)

    def test_duplicate_names(self, tmp_path):
        _write_manifest(tmp_path, [
            ["Mamu-A1*026:new01", "Mamu", "A1", "I", "coding",
             "Mamu/A1.fasta", "J. Karl", "2026-02-15", ""],
            ["Mamu-A1*026:new01", "Mamu", "A1", "I", "coding",
             "Mamu/A1.fasta", "J. Karl", "2026-02-15", ""],
        ])
        _write_fasta(tmp_path, "Mamu/A1.fasta", [
            ("Mamu-A1*026:new01", "ATGGCACCCTGAAC"),
        ])
        errors = validate(tmp_path)
        assert any("Duplicate" in e for e in errors)

    def test_duplicate_sequences(self, tmp_path):
        _write_manifest(tmp_path, [
            ["Mamu-A1*026:new01", "Mamu", "A1", "I", "coding",
             "Mamu/A1.fasta", "J. Karl", "2026-02-15", ""],
            ["Mamu-A1*026:new02", "Mamu", "A1", "I", "coding",
             "Mamu/A1.fasta", "J. Karl", "2026-02-15", ""],
        ])
        _write_fasta(tmp_path, "Mamu/A1.fasta", [
            ("Mamu-A1*026:new01", "ATGGCACCCTGAAC"),
            ("Mamu-A1*026:new02", "ATGGCACCCTGAAC"),  # Same sequence
        ])
        errors = validate(tmp_path)
        assert any("same sequence" in e for e in errors)


class TestValidateNoIpdCollisions:
    def test_no_collisions(self):
        manifest = [{"name": "Mamu-A1*026:new01"}]
        ipd_names = {"Mamu-A1*026:01:01:01", "Mamu-A1*026:02"}
        assert validate_no_ipd_collisions(manifest, ipd_names) == []

    def test_collision(self):
        manifest = [{"name": "Mamu-A1*026:01:01:01"}]
        ipd_names = {"Mamu-A1*026:01:01:01"}
        errors = validate_no_ipd_collisions(manifest, ipd_names)
        assert len(errors) == 1
        assert "collides" in errors[0]


class TestBuildMetadata:
    def test_builds_records(self):
        manifest = [{
            "name": "Mamu-A1*026:new01",
            "species": "Mamu",
            "locus": "A1",
            "class": "I",
            "seq_type": "coding",
            "submitter": "J. Karl",
            "date_added": "2026-02-15",
        }]
        sequences = {"Mamu-A1*026:new01": "ATGGCACCCTGAAC"}
        result = build_metadata(manifest, sequences)
        assert len(result) == 1
        r = result[0]
        assert r["a"].startswith("PROV")
        assert r["n"] == "Mamu-A1*026:new01"
        assert r["l"] == "A1"
        assert r["c"] == "I"
        assert r["s"] == "Mamu"
        assert r["p"] == "provisional"
        assert r["prov"] is True
        assert r["sub"] == "J. Karl"
        assert r["st"] == "coding"

    def test_skips_missing_sequence(self):
        manifest = [{"name": "Mamu-A1*026:new01", "species": "Mamu", "locus": "A1",
                      "class": "I", "seq_type": "coding", "submitter": "X", "date_added": "2026-01-01"}]
        result = build_metadata(manifest, {})
        assert result == []


class TestSequenceIndex:
    def test_build_from_genbank(self, tmp_path):
        """Build a sequence index from GenBank files on disk."""
        # Create a minimal GenBank file
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        gb_dir = tmp_path / "data" / "MHC" / "Mamu"
        gb_dir.mkdir(parents=True)
        rec = SeqRecord(Seq("ATGGCACCCTGAAC"), id="NHP01224", name="NHP01224",
                        description="Mamu-A1*026:01:01:01, A1 locus allele")
        rec.annotations["molecule_type"] = "DNA"
        from Bio import SeqIO
        with open(gb_dir / "A1.gb", "w") as fh:
            SeqIO.write([rec], fh, "genbank")

        index = build_sequence_index(tmp_path)
        h = sequence_hash("ATGGCACCCTGAAC")
        assert h in index
        assert index[h]["accession"] == "NHP01224"

    def test_saves_to_disk(self, tmp_path):
        gb_dir = tmp_path / "data" / "MHC" / "Mamu"
        gb_dir.mkdir(parents=True)
        # Empty dir, no GenBank files
        index = build_sequence_index(tmp_path)
        index_path = tmp_path / "provisional" / "ipd_sequences.json"
        assert index_path.exists()
        assert json.loads(index_path.read_text()) == {}


class TestCheckRetirements:
    def test_finds_match(self):
        manifest = [{"name": "Mamu-A1*026:new01"}]
        sequences = {"Mamu-A1*026:new01": "ATGGCACCCTGAAC"}
        h = sequence_hash("ATGGCACCCTGAAC")
        ipd_index = {h: {"accession": "NHP01224", "name": "Mamu-A1*026:05"}}

        retirements = check_retirements(manifest, sequences, ipd_index)
        assert len(retirements) == 1
        assert retirements[0]["provisional_name"] == "Mamu-A1*026:new01"
        assert retirements[0]["ipd_accession"] == "NHP01224"

    def test_no_match(self):
        manifest = [{"name": "Mamu-A1*026:new01"}]
        sequences = {"Mamu-A1*026:new01": "ATGGCACCCTGAAC"}
        retirements = check_retirements(manifest, sequences, {})
        assert retirements == []


class TestRetireAlleles:
    def test_retires_and_logs(self, tmp_path):
        _write_manifest(tmp_path, [
            ["Mamu-A1*026:new01", "Mamu", "A1", "I", "coding",
             "Mamu/A1.fasta", "J. Karl", "2026-02-15", ""],
            ["Mamu-B*017:new01", "Mamu", "B", "I", "coding",
             "Mamu/B.fasta", "J. Karl", "2026-02-15", ""],
        ])
        _write_fasta(tmp_path, "Mamu/A1.fasta", [
            ("Mamu-A1*026:new01", "ATGGCACCCTGAAC"),
        ])
        _write_fasta(tmp_path, "Mamu/B.fasta", [
            ("Mamu-B*017:new01", "ATGGCACCCTGAAT"),
        ])

        retire_alleles(tmp_path, [{
            "provisional_name": "Mamu-A1*026:new01",
            "ipd_accession": "NHP01224",
            "ipd_name": "Mamu-A1*026:05",
        }])

        # Check manifest
        remaining = load_manifest(tmp_path)
        assert len(remaining) == 1
        assert remaining[0]["name"] == "Mamu-B*017:new01"

        # Check retired.tsv
        retired_path = tmp_path / "provisional" / "retired" / "retired.tsv"
        with open(retired_path) as fh:
            lines = fh.readlines()
        # Header + original header + new entry
        assert any("Mamu-A1*026:new01" in line for line in lines)
        assert any("NHP01224" in line for line in lines)

        # Check FASTA file removed (was only allele)
        fasta_path = tmp_path / "provisional" / "sequences" / "Mamu" / "A1.fasta"
        assert not fasta_path.exists()


class TestAppendToManifest:
    def test_appends(self, tmp_path):
        _write_manifest(tmp_path, [])
        append_to_manifest(tmp_path, {
            "name": "Mamu-A1*026:new01",
            "species": "Mamu",
            "locus": "A1",
            "class": "I",
            "seq_type": "coding",
            "sequence_file": "Mamu/A1.fasta",
            "submitter": "J. Karl",
            "date_added": "2026-02-15",
            "notes": "",
        })
        result = load_manifest(tmp_path)
        assert len(result) == 1
        assert result[0]["name"] == "Mamu-A1*026:new01"


class TestAppendToFasta:
    def test_creates_and_appends(self, tmp_path):
        append_to_fasta(tmp_path, "Mamu/A1.fasta", "Mamu-A1*026:new01", "ATGGCACCCTGAAC")
        fasta_path = tmp_path / "provisional" / "sequences" / "Mamu" / "A1.fasta"
        assert fasta_path.exists()
        content = fasta_path.read_text()
        assert ">Mamu-A1*026:new01" in content
        assert "ATGGCACCCTGAAC" in content

    def test_wraps_long_sequence(self, tmp_path):
        long_seq = "A" * 150
        append_to_fasta(tmp_path, "Mamu/A1.fasta", "test", long_seq)
        fasta_path = tmp_path / "provisional" / "sequences" / "Mamu" / "A1.fasta"
        lines = fasta_path.read_text().strip().split("\n")
        assert lines[0] == ">test"
        assert len(lines[1]) == 70
        assert len(lines[2]) == 70
        assert len(lines[3]) == 10


class TestProvisionalToSeqrecord:
    def test_coding_sequence(self):
        sr = provisional_to_seqrecord(
            name="Mamu-A1*026:new01",
            sequence="ATGGCACCCTGAAC",
            species="Mamu",
            locus="A1",
            allele_class="I",
            seq_type="coding",
            species_info={"scientificName": "Macaca mulatta", "commonName": "Rhesus macaque", "taxon": 9544},
            date_added="2026-02-15",
            submitter="J. Karl",
            accession="PROVa3f2c",
        )
        assert sr.id == "PROVa3f2c"
        assert str(sr.seq) == "ATGGCACCCTGAAC"
        assert "provisional" in sr.description.lower()
        assert sr.annotations["molecule_type"] == "DNA"
        assert sr.annotations["organism"] == "Macaca mulatta"
        assert "Provisional" in sr.annotations["keywords"]
        assert "J. Karl" in sr.annotations["comment"]

        # Should have source, gene, and CDS features
        feature_types = [f.type for f in sr.features]
        assert "source" in feature_types
        assert "gene" in feature_types
        assert "CDS" in feature_types

        # Source should have mRNA mol_type for coding
        source = [f for f in sr.features if f.type == "source"][0]
        assert source.qualifiers["mol_type"] == ["mRNA"]

    def test_genomic_sequence(self):
        sr = provisional_to_seqrecord(
            name="Mamu-B*017:new01",
            sequence="ATGGCACCCTGAACGTACTG",
            species="Mamu",
            locus="B",
            allele_class="I",
            seq_type="genomic",
            species_info={"scientificName": "Macaca mulatta", "taxon": 9544},
            date_added="2026-02-15",
            submitter="J. Karl",
            accession="PROVb1c2d",
        )
        # Genomic should NOT have CDS feature
        feature_types = [f.type for f in sr.features]
        assert "source" in feature_types
        assert "gene" in feature_types
        assert "CDS" not in feature_types

        # Source should have genomic DNA mol_type
        source = [f for f in sr.features if f.type == "source"][0]
        assert source.qualifiers["mol_type"] == ["genomic DNA"]

    def test_date_annotation(self):
        sr = provisional_to_seqrecord(
            name="Test", sequence="ATGGC", species="Mamu", locus="A1",
            allele_class="I", seq_type="coding",
            species_info={}, date_added="2026-02-15", submitter="X",
        )
        assert sr.annotations["date"] == "15-FEB-2026"

    def test_taxon_dbxref(self):
        sr = provisional_to_seqrecord(
            name="Test", sequence="ATGGC", species="Mamu", locus="A1",
            allele_class="I", seq_type="coding",
            species_info={"scientificName": "Macaca mulatta", "taxon": 9544},
            date_added="2026-02-15", submitter="X",
        )
        source = [f for f in sr.features if f.type == "source"][0]
        assert "taxon:9544" in source.qualifiers["db_xref"][0]

    def test_mhc_class_keyword(self):
        sr = provisional_to_seqrecord(
            name="Test", sequence="ATGGC", species="Mamu", locus="A1",
            allele_class="II", seq_type="coding",
            species_info={}, date_added="2026-02-15", submitter="X",
        )
        assert "MHC class II" in sr.annotations["keywords"]


class TestBuildGenbank:
    def test_generates_genbank_file(self, tmp_path):
        manifest = [{
            "name": "Mamu-A1*026:new01",
            "species": "Mamu",
            "locus": "A1",
            "class": "I",
            "seq_type": "coding",
            "submitter": "J. Karl",
            "date_added": "2026-02-15",
        }]
        sequences = {"Mamu-A1*026:new01": "ATGGCACCCTGAAC"}
        species_map = {
            "Mamu": {"scientificName": "Macaca mulatta", "commonName": "Rhesus macaque", "taxon": 9544},
        }

        n = build_genbank(tmp_path, manifest, sequences, species_map)
        assert n == 1

        gb_path = tmp_path / "data" / "provisional" / "Mamu" / "A1.gb"
        assert gb_path.exists()

        # Parse the generated GenBank file
        records = list(SeqIO.parse(str(gb_path), "genbank"))
        assert len(records) == 1
        assert str(records[0].seq) == "ATGGCACCCTGAAC"

    def test_multiple_loci(self, tmp_path):
        manifest = [
            {"name": "Mamu-A1*026:new01", "species": "Mamu", "locus": "A1",
             "class": "I", "seq_type": "coding", "submitter": "X", "date_added": "2026-01-01"},
            {"name": "Mamu-B*017:new01", "species": "Mamu", "locus": "B",
             "class": "I", "seq_type": "coding", "submitter": "X", "date_added": "2026-01-01"},
        ]
        sequences = {
            "Mamu-A1*026:new01": "ATGGCACCCTGAAC",
            "Mamu-B*017:new01": "ATGGCACCCTGAAT",
        }
        species_map = {"Mamu": {"scientificName": "Macaca mulatta", "taxon": 9544}}

        n = build_genbank(tmp_path, manifest, sequences, species_map)
        assert n == 2

        assert (tmp_path / "data" / "provisional" / "Mamu" / "A1.gb").exists()
        assert (tmp_path / "data" / "provisional" / "Mamu" / "B.gb").exists()

    def test_skips_missing_sequence(self, tmp_path):
        manifest = [{"name": "Mamu-A1*026:new01", "species": "Mamu", "locus": "A1",
                      "class": "I", "seq_type": "coding", "submitter": "X", "date_added": "2026-01-01"}]
        n = build_genbank(tmp_path, manifest, {}, {})
        assert n == 0


class TestMetadataIndexProvisional:
    def test_includes_provisional_alleles(self, tmp_path):
        from pipeline.metadata_index import build_metadata_index

        prov = [{
            "a": "PROVa3f2c",
            "n": "Mamu-A1*026:new01",
            "l": "A1",
            "c": "I",
            "s": "Mamu",
            "p": "provisional",
            "da": "2026-02-15",
            "dm": "2026-02-15",
            "prov": True,
            "sub": "J. Karl",
            "st": "coding",
        }]

        output = tmp_path / "alleles.json"
        index = build_metadata_index([], [], output, provisional_alleles=prov)
        assert index["provisional_count"] == 1
        assert any(a["a"] == "PROVa3f2c" for a in index["alleles"])

    def test_no_provisional_count_when_empty(self, tmp_path):
        from pipeline.metadata_index import build_metadata_index

        output = tmp_path / "alleles.json"
        index = build_metadata_index([], [], output)
        assert "provisional_count" not in index


def _make_seqrecord(seq_id, sequence):
    """Helper to create a Bio.SeqRecord for testing."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    return SeqRecord(Seq(sequence), id=seq_id, description=seq_id)


def _setup_genbank(tmp_path, species, locus, alleles):
    """Helper to create a GenBank file with alleles for nearest-neighbor tests."""
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    gb_dir = tmp_path / "data" / "MHC" / species
    gb_dir.mkdir(parents=True, exist_ok=True)
    recs = []
    for allele_name, seq in alleles:
        rec = SeqRecord(Seq(seq), id="NHP00000", name="NHP00000",
                        description=f"{allele_name}, {locus} locus allele")
        rec.annotations["molecule_type"] = "DNA"
        recs.append(rec)
    with open(gb_dir / f"{locus}.gb", "w") as fh:
        SeqIO.write(recs, fh, "genbank")


class TestAddProvisionalAlleles:
    def test_single_record(self, tmp_path):
        _write_manifest(tmp_path, [])
        _setup_genbank(tmp_path, "Mamu", "A1",
                       [("Mamu-A1*026:01:01:01", "ATGGCACCCTGAAC")])

        records = [_make_seqrecord("seq1", "ATGGCACCCTGAAT")]
        names = add_provisional_alleles(
            tmp_path, records, "Mamu", "A1", "I", "coding", "J. Karl",
        )
        assert len(names) == 1
        # New namer detects synonymous protein → uses protein field in name
        assert names[0] == "Mamu-A1*026:01:new01"

        manifest = load_manifest(tmp_path)
        assert len(manifest) == 1
        assert manifest[0]["submitter"] == "J. Karl"

    def test_multiple_records_increment(self, tmp_path):
        _write_manifest(tmp_path, [])
        _setup_genbank(tmp_path, "Mamu", "A1",
                       [("Mamu-A1*026:01:01:01", "ATGGCACCCTGAAC")])

        records = [
            _make_seqrecord("seq1", "ATGGCACCCTGAAT"),
            _make_seqrecord("seq2", "ATGGCACCCTGAAG"),
        ]
        names = add_provisional_alleles(
            tmp_path, records, "Mamu", "A1", "I", "coding", "J. Karl",
        )
        assert len(names) == 2
        # New namer detects synonymous protein → uses protein field in name
        assert names[0] == "Mamu-A1*026:01:new01"
        assert names[1] == "Mamu-A1*026:01:new02"

        manifest = load_manifest(tmp_path)
        assert len(manifest) == 2

    def test_name_override_single(self, tmp_path):
        _write_manifest(tmp_path, [])

        records = [_make_seqrecord("seq1", "ATGGCACCCTGAAT")]
        names = add_provisional_alleles(
            tmp_path, records, "Mamu", "A1", "I", "coding", "J. Karl",
            name_override="Mamu-A1*999:new01",
        )
        assert names == ["Mamu-A1*999:new01"]

    def test_name_override_multi_raises(self, tmp_path):
        _write_manifest(tmp_path, [])

        records = [
            _make_seqrecord("seq1", "ATGGCACCCTGAAT"),
            _make_seqrecord("seq2", "ATGGCACCCTGAAG"),
        ]
        with pytest.raises(ValueError, match="cannot be used with multiple"):
            add_provisional_alleles(
                tmp_path, records, "Mamu", "A1", "I", "coding", "J. Karl",
                name_override="Mamu-A1*999:new01",
            )

    def test_duplicate_in_batch_raises(self, tmp_path):
        _write_manifest(tmp_path, [])

        records = [
            _make_seqrecord("seq1", "ATGGCACCCTGAAT"),
            _make_seqrecord("seq2", "ATGGCACCCTGAAT"),
        ]
        with pytest.raises(ValueError, match="same sequence"):
            add_provisional_alleles(
                tmp_path, records, "Mamu", "A1", "I", "coding", "J. Karl",
            )

    def test_duplicate_vs_existing_raises(self, tmp_path):
        _write_manifest(tmp_path, [
            ["Mamu-A1*026:new01", "Mamu", "A1", "I", "coding",
             "Mamu/A1.fasta", "J. Karl", "2026-02-15", ""],
        ])
        _write_fasta(tmp_path, "Mamu/A1.fasta", [
            ("Mamu-A1*026:new01", "ATGGCACCCTGAAT"),
        ])

        records = [_make_seqrecord("seq1", "ATGGCACCCTGAAT")]
        with pytest.raises(ValueError, match="existing provisional"):
            add_provisional_alleles(
                tmp_path, records, "Mamu", "A1", "I", "coding", "J. Karl",
            )

    def test_empty_sequence_raises(self, tmp_path):
        _write_manifest(tmp_path, [])

        records = [_make_seqrecord("seq1", "")]
        with pytest.raises(ValueError, match="empty sequence"):
            add_provisional_alleles(
                tmp_path, records, "Mamu", "A1", "I", "coding", "J. Karl",
            )

    def test_invalid_nucleotides_raises(self, tmp_path):
        _write_manifest(tmp_path, [])

        records = [_make_seqrecord("seq1", "ATGXYZ123")]
        with pytest.raises(ValueError, match="invalid nucleotide"):
            add_provisional_alleles(
                tmp_path, records, "Mamu", "A1", "I", "coding", "J. Karl",
            )

    def test_no_genbank_fallback(self, tmp_path):
        """When no GenBank file exists, names use *000:newNN fallback."""
        _write_manifest(tmp_path, [])

        records = [_make_seqrecord("seq1", "ATGGCACCCTGAAT")]
        names = add_provisional_alleles(
            tmp_path, records, "Mamu", "Z9", "I", "coding", "J. Karl",
        )
        assert names[0] == "Mamu-Z9*000:new01"
