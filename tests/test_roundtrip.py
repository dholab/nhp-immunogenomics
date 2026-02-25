"""End-to-end roundtrip tests: API JSON -> GenBank -> parse back."""

import io

from Bio import SeqIO

from pipeline.genbank_builder import allele_to_seqrecord


class TestRoundtrip:
    def _roundtrip(self, record, project):
        """Convert to GenBank and parse back."""
        sr = allele_to_seqrecord(record, project)
        buf = io.StringIO()
        SeqIO.write(sr, buf, "genbank")
        buf.seek(0)
        return SeqIO.read(buf, "genbank")

    def test_genomic_mhc_roundtrip(self, mhc_genomic_record):
        """Full genomic MHC record survives GenBank roundtrip."""
        parsed = self._roundtrip(mhc_genomic_record, "MHC")
        assert parsed.id == "NHP01224"
        assert len(parsed.seq) == 2874
        assert parsed.annotations["organism"] == "Macaca mulatta"

        # Verify CDS feature exists with correct structure
        cds_feats = [f for f in parsed.features if f.type == "CDS"]
        assert len(cds_feats) == 1
        cds = cds_feats[0]
        assert cds.qualifiers["gene"] == ["A1"]
        assert cds.qualifiers["allele"] == ["Mamu-A1*026:01:01:01"]
        assert "translation" in cds.qualifiers

    def test_coding_only_roundtrip(self, mhc_coding_only_record):
        """Coding-only MHC record survives GenBank roundtrip."""
        parsed = self._roundtrip(mhc_coding_only_record, "MHC")
        assert parsed.id == "NHP00001"
        assert len(parsed.seq) == 224

        source = [f for f in parsed.features if f.type == "source"][0]
        assert source.qualifiers["mol_type"] == ["mRNA"]

    def test_nhkir_roundtrip(self, nhkir_record):
        """NHKIR record survives GenBank roundtrip."""
        parsed = self._roundtrip(nhkir_record, "NHKIR")
        assert parsed.id == "NHP00001"

        exons = [f for f in parsed.features if f.type == "exon"]
        assert len(exons) == 9

    def test_cds_extracts_coding_sequence(self, mhc_genomic_record):
        """CDS join coordinates extract the correct coding sequence from genomic."""
        sr = allele_to_seqrecord(mhc_genomic_record, "MHC")
        cds = [f for f in sr.features if f.type == "CDS"][0]

        # Extract the CDS sequence from the genomic sequence
        extracted = cds.location.extract(sr.seq)
        coding_from_api = mhc_genomic_record["sequence"]["coding"]
        assert str(extracted).lower() == coding_from_api.lower()

    def test_cds_translation_matches_protein(self, mhc_genomic_record):
        """CDS /translation matches the protein from the API."""
        sr = allele_to_seqrecord(mhc_genomic_record, "MHC")
        cds = [f for f in sr.features if f.type == "CDS"][0]
        translation = cds.qualifiers["translation"][0]
        api_protein = mhc_genomic_record["sequence"]["protein"]
        assert translation == api_protein

    def test_multirecord_genbank(self, mhc_genomic_record, mhc_coding_only_record, nhkir_record):
        """Multiple records in one GenBank file parse correctly."""
        records = [
            allele_to_seqrecord(mhc_genomic_record, "MHC"),
            allele_to_seqrecord(mhc_coding_only_record, "MHC"),
            allele_to_seqrecord(nhkir_record, "NHKIR"),
        ]
        buf = io.StringIO()
        SeqIO.write(records, buf, "genbank")
        buf.seek(0)
        parsed = list(SeqIO.parse(buf, "genbank"))
        assert len(parsed) == 3
        assert parsed[0].id == "NHP01224"
        assert parsed[1].id == "NHP00001"
        assert parsed[2].id == "NHP00001"  # Same accession, different project
