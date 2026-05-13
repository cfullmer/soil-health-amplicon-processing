#!/usr/bin/env bash
set -euo pipefail

# Manuscript-facing 16S post-merge workflow:
# taxonomy assignment, prokaryote-only filtering, and export.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ ! -f "${SCRIPT_DIR}/paths.env" ]]; then
  echo "Missing ${SCRIPT_DIR}/paths.env. Copy 00_paths_template.env to paths.env and edit it first." >&2
  exit 1
fi
source "${SCRIPT_DIR}/paths.env"

MERGED_SEQS="${MERGED_16S_SEQS:-${PROCESSED_DATA_DIR}/16S/rep_seqs/SH_16S_rep-seqs.qza}"
MERGED_TBL="${MERGED_16S_TABLE:-${PROCESSED_DATA_DIR}/16S/filtered_tables/SH_16S_filtered_table.qza}"
OUT_DIR="${POSTMERGE_16S_DIR:-${PROCESSED_DATA_DIR}/16S}"
TAX_DIR="${OUT_DIR}/taxonomy"
PROK_DIR="${OUT_DIR}/prok_only"
EXPORT_DIR="${OUT_DIR}/exported/non_rarefied"

mkdir -p \
  "${TAX_DIR}" \
  "${PROK_DIR}" \
  "${EXPORT_DIR}"

echo "== 16S taxonomy assignment =="
qiime feature-classifier classify-sklearn \
  --i-classifier "${CLASSIFIER_16S}" \
  --i-reads "${MERGED_SEQS}" \
  --o-classification "${TAX_DIR}/taxonomy-16S.qza"

qiime metadata tabulate \
  --m-input-file "${TAX_DIR}/taxonomy-16S.qza" \
  --o-visualization "${TAX_DIR}/taxonomy-16S.qzv"

qiime taxa barplot \
  --i-table "${MERGED_TBL}" \
  --i-taxonomy "${TAX_DIR}/taxonomy-16S.qza" \
  --m-metadata-file "${METADATA_16S}" \
  --o-visualization "${TAX_DIR}/taxa-bar-plots-16S.qzv"

echo "== 16S prokaryote-only filtering =="
qiime taxa filter-table \
  --i-table "${MERGED_TBL}" \
  --i-taxonomy "${TAX_DIR}/taxonomy-16S.qza" \
  --p-exclude "Unassigned,mitochondria,chloroplast,Eukaryota" \
  --o-filtered-table "${PROK_DIR}/table-prok-only.qza"

qiime feature-table filter-seqs \
  --i-data "${MERGED_SEQS}" \
  --i-table "${PROK_DIR}/table-prok-only.qza" \
  --o-filtered-data "${PROK_DIR}/filtered_rep_seqs_prok-only.qza"

echo "== 16S export =="
qiime tools export \
  --input-path "${TAX_DIR}/taxonomy-16S.qza" \
  --output-path "${TAX_DIR}/exported-taxonomy"

cp "${TAX_DIR}/exported-taxonomy/taxonomy.tsv" "${TAX_DIR}/exported-taxonomy/biom-taxonomy.tsv"
tr '\t' ',' < "${TAX_DIR}/exported-taxonomy/biom-taxonomy.tsv" > "${TAX_DIR}/exported-taxonomy/biom-taxonomy.csv"

python - <<'PY' "${TAX_DIR}/exported-taxonomy/biom-taxonomy.csv"
from pathlib import Path
import sys

csv_path = Path(sys.argv[1])
lines = csv_path.read_text().splitlines()
if lines:
    lines[0] = "#OTUID,taxonomy,confidence"
csv_path.write_text("\n".join(lines) + ("\n" if lines else ""))
PY

tr ',' '\t' < "${TAX_DIR}/exported-taxonomy/biom-taxonomy.csv" > "${TAX_DIR}/exported-taxonomy/biom-taxonomy.tsv"
rm -f "${TAX_DIR}/exported-taxonomy/biom-taxonomy.csv"

qiime tools export \
  --input-path "${PROK_DIR}/table-prok-only.qza" \
  --output-path "${EXPORT_DIR}"

qiime tools export \
  --input-path "${PROK_DIR}/filtered_rep_seqs_prok-only.qza" \
  --output-path "${EXPORT_DIR}/rep_seqs"

biom add-metadata \
  -i "${EXPORT_DIR}/feature-table.biom" \
  -o "${EXPORT_DIR}/feature-table-with-taxonomy.biom" \
  --observation-metadata-fp "${TAX_DIR}/exported-taxonomy/biom-taxonomy.tsv" \
  --sc-separated taxonomy

biom convert \
  -i "${EXPORT_DIR}/feature-table-with-taxonomy.biom" \
  -o "${EXPORT_DIR}/feature-table-with-taxonomy.tsv" \
  --to-tsv \
  --process-obs-metadata taxonomy \
  --header-key taxonomy
