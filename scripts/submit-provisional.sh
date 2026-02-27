#!/usr/bin/env bash
# Submit provisional alleles from a local FASTA file.
# This bypasses the GitHub Actions workflow_dispatch input size limit.
#
# Requirements: git, gh (GitHub CLI), uv
#
# Usage:
#   ./scripts/submit-provisional.sh \
#     --species Mamu --locus E --submitter "J. Karl" \
#     --fasta /path/to/sequences.fasta
set -euo pipefail

usage() {
    cat >&2 <<EOF
Usage: $0 --species SPECIES --locus LOCUS --submitter SUBMITTER --fasta FILE [OPTIONS]

Required:
  --species     Species prefix (e.g., Mamu, Mafa, Mane)
  --locus       Locus name (e.g., A1, B, DRB, KIR3DL01)
  --submitter   Submitter name (e.g., "J. Karl")
  --fasta       Path to FASTA file

Optional:
  --class       MHC class: I, II, or empty [default: I]
  --seq-type    Sequence type: coding or genomic [default: coding]
  --notes       Optional notes [default: ""]
EOF
    exit 1
}

SPECIES="" LOCUS="" SUBMITTER="" FASTA="" CLASS="I" SEQ_TYPE="coding" NOTES=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --species)   SPECIES="$2";   shift 2;;
        --locus)     LOCUS="$2";     shift 2;;
        --submitter) SUBMITTER="$2"; shift 2;;
        --fasta)     FASTA="$2";     shift 2;;
        --class)     CLASS="$2";     shift 2;;
        --seq-type)  SEQ_TYPE="$2";  shift 2;;
        --notes)     NOTES="$2";     shift 2;;
        -h|--help)   usage;;
        *)           echo "Unknown option: $1" >&2; usage;;
    esac
done

[[ -z "$SPECIES" || -z "$LOCUS" || -z "$SUBMITTER" || -z "$FASTA" ]] && usage
[[ -f "$FASTA" ]] || { echo "Error: FASTA file not found: $FASTA" >&2; exit 1; }

for cmd in git gh uv; do
    command -v "$cmd" >/dev/null || { echo "Error: $cmd is required but not found" >&2; exit 1; }
done

if [[ -n "$(git status --porcelain)" ]]; then
    echo "Error: working tree is not clean. Commit or stash changes first." >&2
    exit 1
fi

REPO_ROOT="$(git rev-parse --show-toplevel)"
ORIGINAL_BRANCH="$(git rev-parse --abbrev-ref HEAD)"
BRANCH="provisional/${SPECIES}-${LOCUS}-$(date +%s)"

echo "Creating branch: $BRANCH"
git checkout -b "$BRANCH"

echo "Adding provisional alleles..."
uv run nhp-pipeline provisional-add \
    --species "$SPECIES" \
    --locus "$LOCUS" \
    --class "$CLASS" \
    --seq-type "$SEQ_TYPE" \
    --submitter "$SUBMITTER" \
    --notes "$NOTES" \
    --fasta "$FASTA"

echo "Building provisional GenBank files..."
uv run nhp-pipeline provisional-build

added=$(git diff --unified=0 provisional/manifest.tsv | grep '^+[^+]' | grep -v '^+name' | wc -l | tr -d ' ')
names=$(git diff --unified=0 provisional/manifest.tsv | grep '^+[^+]' | grep -v '^+name' | cut -f1)

echo "Committing ${added} allele(s)..."
git add provisional/
git commit -m "Add ${added} provisional allele(s) at ${SPECIES}-${LOCUS}"

echo "Pushing..."
git push -u origin "$BRANCH"

echo "Creating pull request..."
gh pr create \
    --title "Add provisional: ${SPECIES}-${LOCUS} (${added} allele(s))" \
    --body "$(cat <<EOF
## Provisional Allele Submission

| Field | Value |
|-------|-------|
| **Species** | ${SPECIES} |
| **Locus** | ${LOCUS} |
| **Class** | ${CLASS} |
| **Seq type** | ${SEQ_TYPE} |
| **Submitter** | ${SUBMITTER} |
| **Count** | ${added} |
| **Notes** | ${NOTES} |

### Assigned names
\`\`\`
${names}
\`\`\`

Submitted via \`submit-provisional.sh\`.
This PR will be validated by the \`validate-provisional\` CI check.
EOF
)"

git checkout "$ORIGINAL_BRANCH"
echo "Done!"
