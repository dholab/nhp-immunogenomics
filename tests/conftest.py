"""Shared test fixtures."""

import json
from pathlib import Path

import pytest

FIXTURES_DIR = Path(__file__).parent / "fixtures"


def load_fixture(name: str) -> dict:
    """Load a JSON fixture file."""
    return json.loads((FIXTURES_DIR / name).read_text())


@pytest.fixture
def mhc_genomic_record():
    """Full genomic MHC Class I allele (Mamu-A1*026:01:01:01, 8 exons, 7 introns)."""
    return load_fixture("allele_mhc_genomic.json")


@pytest.fixture
def mhc_coding_only_record():
    """CDS-only MHC Class II allele (Aona-DQA1*27:01, single exon, codon_start=3)."""
    return load_fixture("allele_mhc_coding_only.json")


@pytest.fixture
def mhc_full_class1_record():
    """Full MHC Class I allele (Mafa-B*028:04:01:01)."""
    return load_fixture("allele_mhc_full_class1.json")


@pytest.fixture
def nhkir_record():
    """Full NHKIR KIR record (Mamu-KIR3DL01*013)."""
    return load_fixture("allele_nhkir.json")


@pytest.fixture
def listing_response():
    """API listing response with 3 records and pagination cursor."""
    return load_fixture("listing_with_fields.json")


@pytest.fixture
def tmp_repo(tmp_path):
    """Create a temporary repo-like directory structure."""
    (tmp_path / "cache").mkdir()
    (tmp_path / "data").mkdir()
    (tmp_path / "docs").mkdir()
    return tmp_path
