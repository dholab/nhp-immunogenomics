"""Pipeline configuration constants."""

API_BASE = "https://www.ebi.ac.uk/cgi-bin/ipd/api"
PROJECTS = ("MHC", "NHKIR")

# API pagination limit (max 1000 per the API docs)
PAGE_LIMIT = 1000

# Rate limiting: seconds between API requests
REQUEST_DELAY = 0.25  # 4 requests/sec

# Retry configuration
MAX_RETRIES = 3
RETRY_BACKOFF = 2.0  # Exponential backoff base (seconds)

# Paths (relative to repo root)
DATA_DIR = "data"
CACHE_DIR = "cache"
DOCS_DIR = "docs"
VERSION_FILE = "version.json"

# Provisional allele paths
PROVISIONAL_DIR = "provisional"
PROVISIONAL_MANIFEST = "provisional/manifest.tsv"
PROVISIONAL_SEQUENCES_DIR = "provisional/sequences"
PROVISIONAL_DATA_DIR = "data/provisional"
PROVISIONAL_RETIRED = "provisional/retired/retired.tsv"
PROVISIONAL_SEQ_INDEX = "provisional/ipd_sequences.json"

# Fields to request in the listing endpoint (flat dot-notation keys)
LISTING_FIELDS = (
    "accession,name,locus,class,status,"
    "date_assigned,date_modified,previous,"
    "organism.name,organism.scientificName,"
    "organism.commonName,organism.taxon"
)
