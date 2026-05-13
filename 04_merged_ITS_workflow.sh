#!/usr/bin/env bash
set -euo pipefail

# Manuscript-facing merged ITS workflow:
# merge per-batch DADA2 artifacts, cluster at 97%, BLAST/QC against UNITE,
# assign taxonomy, filter to fungi, and export final representative sequences.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ ! -f "${SCRIPT_DIR}/paths.env" ]]; then
  echo "Missing ${SCRIPT_DIR}/paths.env. Copy 00_paths_template.env to paths.env and edit it first." >&2
  exit 1
fi
source "${SCRIPT_DIR}/paths.env"

OUT_DIR="${MERGED_ITS_WORKFLOW_DIR:-${PROCESSED_DATA_DIR}/ITS/merged_workflow}"
MERGED_DIR="${OUT_DIR}/merged"
CLUSTER_DIR="${OUT_DIR}/cluster97"
QC_DIR="${OUT_DIR}/qc_blast"
TAX_DIR="${OUT_DIR}/taxonomy"
FILTER_DIR="${OUT_DIR}/filtered_tables"
SEQ_DIR="${OUT_DIR}/filtered_seqs"
EXPORT_DIR="${OUT_DIR}/exported/non_rarefied"

mkdir -p \
  "${MERGED_DIR}" \
  "${CLUSTER_DIR}" \
  "${QC_DIR}" \
  "${TAX_DIR}" \
  "${FILTER_DIR}" \
  "${SEQ_DIR}" \
  "${EXPORT_DIR}/rep_seqs"

# Edit these batch-specific inputs as needed for the manuscript dataset.
BATCH_REP=(
  "${PROCESSED_DATA_DIR}/imported_seqs/imported-paired-end-seqs-NN21-ITS-trimmed-2024-08-23-rep-seqs.qza"
  "${PROCESSED_DATA_DIR}/imported_seqs/imported-paired-end-seqs-NN22-ITS-trimmed-2024-08-23-rep-seqs.qza"
  "${PROCESSED_DATA_DIR}/imported_seqs/imported-paired-end-seqs-NN24-ITS-trimmed-2024-08-23-rep-seqs.qza"
  "${PROCESSED_DATA_DIR}/archive/nn26/Fungal_ITS/dada2/rep-seqs-Fungal_ITS.qza"
  "${PROCESSED_DATA_DIR}/ITS/NN28/dada2/NN28-ITS-rep-seqs.qza"
  "${PROCESSED_DATA_DIR}/ITS/NN33/qiime2_out/dada2/rep-seqs-NN33-ITS.qza"
)

BATCH_TBL=(
  "${PROCESSED_DATA_DIR}/imported_seqs/imported-paired-end-seqs-NN21-ITS-trimmed-2024-08-23-table.qza"
  "${PROCESSED_DATA_DIR}/imported_seqs/imported-paired-end-seqs-NN22-ITS-trimmed-2024-08-23-table.qza"
  "${PROCESSED_DATA_DIR}/imported_seqs/imported-paired-end-seqs-NN24-ITS-trimmed-2024-08-23-table.qza"
  "${PROCESSED_DATA_DIR}/archive/nn26/Fungal_ITS/dada2/table-Fungal_ITS.qza"
  "${PROCESSED_DATA_DIR}/ITS/NN28/dada2/NN28-ITS-table.qza"
  "${PROCESSED_DATA_DIR}/ITS/NN33/qiime2_out/dada2/table-NN33-ITS.qza"
)

echo "== Merge per-batch ITS representative sequences =="
args=()
for f in "${BATCH_REP[@]}"; do args+=(--i-data "$f"); done
qiime feature-table merge-seqs "${args[@]}" \
  --o-merged-data "${MERGED_DIR}/merged-rep-seqs-ITS.qza"

echo "== Merge per-batch ITS feature tables =="
args=()
for f in "${BATCH_TBL[@]}"; do args+=(--i-tables "$f"); done
qiime feature-table merge "${args[@]}" \
  --p-overlap-method sum \
  --o-merged-table "${MERGED_DIR}/merged-feature-table-ITS.qza"

echo "== Cluster merged ITS features at 97% identity =="
qiime vsearch cluster-features-de-novo \
  --i-table "${MERGED_DIR}/merged-feature-table-ITS.qza" \
  --i-sequences "${MERGED_DIR}/merged-rep-seqs-ITS.qza" \
  --p-perc-identity 0.97 \
  --o-clustered-table "${CLUSTER_DIR}/merged-ITS-cluster97-table.qza" \
  --o-clustered-sequences "${CLUSTER_DIR}/merged-ITS-cluster97-rep-seqs.qza"

echo "== BLAST/QC merged ITS representative sequences against UNITE =="
qiime quality-control exclude-seqs \
  --i-query-sequences "${CLUSTER_DIR}/merged-ITS-cluster97-rep-seqs.qza" \
  --i-reference-sequences "${UNITE_ALL_SEQS}" \
  --p-method blast \
  --p-perc-identity 0.0 \
  --p-perc-query-aligned 0.83 \
  --o-sequence-hits "${QC_DIR}/quality-control-hits_UNITE_ALL_0.83.qza" \
  --o-sequence-misses "${QC_DIR}/quality-control-misses_UNITE_ALL_0.83.qza"

qiime feature-table filter-features \
  --i-table "${CLUSTER_DIR}/merged-ITS-cluster97-table.qza" \
  --m-metadata-file "${QC_DIR}/quality-control-hits_UNITE_ALL_0.83.qza" \
  --o-filtered-table "${QC_DIR}/merged-ITS-cluster97-table-qc-UNITE_ALL_0.83.qza"

qiime feature-table filter-seqs \
  --i-data "${CLUSTER_DIR}/merged-ITS-cluster97-rep-seqs.qza" \
  --m-metadata-file "${QC_DIR}/quality-control-hits_UNITE_ALL_0.83.qza" \
  --o-filtered-data "${QC_DIR}/merged-ITS-cluster97-rep-seqs-qc-UNITE_ALL_0.83.qza"

echo "== Assign taxonomy to QC-filtered merged ITS sequences =="
qiime feature-classifier classify-sklearn \
  --i-classifier "${CLASSIFIER_ITS}" \
  --i-reads "${QC_DIR}/merged-ITS-cluster97-rep-seqs-qc-UNITE_ALL_0.83.qza" \
  --o-classification "${TAX_DIR}/taxonomy-merged-ITS-cluster97-qc-UNITEALL083.qza"

qiime metadata tabulate \
  --m-input-file "${TAX_DIR}/taxonomy-merged-ITS-cluster97-qc-UNITEALL083.qza" \
  --o-visualization "${TAX_DIR}/taxonomy-merged-ITS-cluster97-qc-UNITEALL083.qzv"

echo "== Filter merged ITS table to fungi only =="
qiime taxa filter-table \
  --i-table "${QC_DIR}/merged-ITS-cluster97-table-qc-UNITE_ALL_0.83.qza" \
  --i-taxonomy "${TAX_DIR}/taxonomy-merged-ITS-cluster97-qc-UNITEALL083.qza" \
  --p-include "k__Fungi" \
  --p-exclude "Unassigned,Rhizaria,Protista,Metazoa" \
  --o-filtered-table "${FILTER_DIR}/merged-ITS-cluster97-fungi-only-table.qza"

qiime feature-table filter-seqs \
  --i-data "${QC_DIR}/merged-ITS-cluster97-rep-seqs-qc-UNITE_ALL_0.83.qza" \
  --i-table "${FILTER_DIR}/merged-ITS-cluster97-fungi-only-table.qza" \
  --o-filtered-data "${SEQ_DIR}/merged-ITS-cluster97-fungi-only-seqs.qza"

echo "== Export final merged ITS representative sequences =="
qiime tools export \
  --input-path "${SEQ_DIR}/merged-ITS-cluster97-fungi-only-seqs.qza" \
  --output-path "${EXPORT_DIR}/rep_seqs"
