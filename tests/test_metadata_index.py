"""Tests for metadata index building."""

import json

from pipeline.metadata_index import build_metadata_index


def _make_listing_entry(accession="NHP00001", name="Mamu-A1*001:01",
                        locus="A1", cls="I", species="Mamu",
                        sci_name="Macaca mulatta", common="rhesus monkey",
                        taxon=9544, previous=None):
    entry = {
        "accession": accession,
        "name": name,
        "locus": locus,
        "class": cls,
        "status": "public",
        "date_assigned": "2020-01-01",
        "date_modified": "2025-01-01",
        "organism.name": species,
        "organism.scientificName": sci_name,
        "organism.commonName": common,
        "organism.taxon": taxon,
    }
    if previous:
        entry["previous"] = previous
    return entry


class TestBuildMetadataIndex:
    def test_basic_structure(self, tmp_path):
        output = tmp_path / "alleles.json"
        mhc = [_make_listing_entry()]
        nhkir = [_make_listing_entry(
            "NHP00100", "Mamu-KIR3DL01*001", "KIR3DL01", "unknown",
            "Mamu", "Macaca mulatta", "rhesus monkey", 9544,
        )]
        index = build_metadata_index(mhc, nhkir, output, "3.16.0.0", "1.7.0.0")

        assert "generated" in index
        assert index["mhc_version"] == "3.16.0.0"
        assert index["nhkir_version"] == "1.7.0.0"
        assert len(index["alleles"]) == 2
        assert "Mamu" in index["species"]
        assert "A1" in index["loci"]["MHC"]
        assert "KIR3DL01" in index["loci"]["NHKIR"]

    def test_compact_keys(self, tmp_path):
        output = tmp_path / "alleles.json"
        mhc = [_make_listing_entry()]
        index = build_metadata_index(mhc, [], output)

        allele = index["alleles"][0]
        assert allele["a"] == "NHP00001"
        assert allele["n"] == "Mamu-A1*001:01"
        assert allele["l"] == "A1"
        assert allele["c"] == "I"
        assert allele["s"] == "Mamu"
        assert allele["p"] == "MHC"
        assert allele["da"] == "2020-01-01"
        assert allele["dm"] == "2025-01-01"

    def test_previous_designations(self, tmp_path):
        output = tmp_path / "alleles.json"
        mhc = [_make_listing_entry(previous=["Mamu-A*25", "Mamu-A1*001:01:old"])]
        index = build_metadata_index(mhc, [], output)

        assert index["alleles"][0]["prev"] == ["Mamu-A*25", "Mamu-A1*001:01:old"]

    def test_no_previous_omits_key(self, tmp_path):
        output = tmp_path / "alleles.json"
        mhc = [_make_listing_entry()]
        index = build_metadata_index(mhc, [], output)
        assert "prev" not in index["alleles"][0]

    def test_species_deduplication(self, tmp_path):
        output = tmp_path / "alleles.json"
        mhc = [
            _make_listing_entry("NHP00001"),
            _make_listing_entry("NHP00002", "Mamu-A1*002:01"),
            _make_listing_entry("NHP00003", "Mafa-B*001", "B", "I", "Mafa",
                               "Macaca fascicularis", "cynomolgus macaque", 9541),
        ]
        index = build_metadata_index(mhc, [], output)
        assert len(index["species"]) == 2
        assert "Mamu" in index["species"]
        assert "Mafa" in index["species"]
        assert index["species"]["Mafa"]["scientificName"] == "Macaca fascicularis"

    def test_loci_grouped_by_project(self, tmp_path):
        output = tmp_path / "alleles.json"
        mhc = [
            _make_listing_entry("NHP00001", locus="A1"),
            _make_listing_entry("NHP00002", locus="B"),
        ]
        nhkir = [
            _make_listing_entry("NHP00003", locus="KIR3DL01", cls="unknown"),
        ]
        index = build_metadata_index(mhc, nhkir, output)
        assert set(index["loci"]["MHC"]) == {"A1", "B"}
        assert set(index["loci"]["NHKIR"]) == {"KIR3DL01"}

    def test_output_is_minified_json(self, tmp_path):
        output = tmp_path / "alleles.json"
        build_metadata_index([_make_listing_entry()], [], output)
        content = output.read_text()
        # Minified JSON should not have spaces after separators
        assert ": " not in content
        assert ", " not in content
        parsed = json.loads(content)
        assert len(parsed["alleles"]) == 1

    def test_creates_parent_dirs(self, tmp_path):
        output = tmp_path / "a" / "b" / "alleles.json"
        build_metadata_index([], [], output)
        assert output.exists()
