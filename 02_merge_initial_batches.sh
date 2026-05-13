#!/usr/bin/env bash
set -euo pipefail

# Cleaned representation of the initial merged 16S / ITS workflow.
# Edit explicit batch artifact paths for the relevant project batches.

source "$(dirname "$0")/paths.env"

# Example 16S inputs from initial merged workflow.
REP_16S=(
  "${PROCESSED_DATA_DIR}/imported_seqs/imported-paired-end-seqs-NN18-16S-trimmed-2024-08-05-rep-seqs.qza"
  "${PROCESSED_DATA_DIR}/imported_seqs/imported-paired-end-seqs-NN19-16S-trimmed-2024-08-23-rep-seqs.qza"
  "${PROCESSED_DATA_DIR}/imported_seqs/imported-paired-end-seqs-NN21-16S-trimmed-2024-08-23-rep-seqs.qza"
  "${PROCESSED_DATA_DIR}/imported_seqs/imported-paired-end-seqs-NN25-16S-trimmed-2024-08-23-rep-seqs.qza"
)

TBL_16S=(
  "${PROCESSED_DATA_DIR}/imported_seqs/imported-paired-end-seqs-NN18-16S-trimmed-2024-08-05-table.qza"
  "${PROCESSED_DATA_DIR}/imported_seqs/imported-paired-end-seqs-NN19-16S-trimmed-2024-08-23-table.qza"
  "${PROCESSED_DATA_DIR}/imported_seqs/imported-paired-end-seqs-NN21-16S-trimmed-2024-08-23-table.qza"
  "${PROCESSED_DATA_DIR}/imported_seqs/imported-paired-end-seqs-NN25-16S-trimmed-2024-08-23-table.qza"
)

REP_ITS=(
  "${PROCESSED_DATA_DIR}/imported_seqs/imported-paired-end-seqs-NN21-ITS-trimmed-2024-08-23-rep-seqs.qza"
  "${PROCESSED_DATA_DIR}/imported_seqs/imported-paired-end-seqs-NN22-ITS-trimmed-2024-08-23-rep-seqs.qza"
  "${PROCESSED_DATA_DIR}/imported_seqs/imported-paired-end-seqs-NN24-ITS-trimmed-2024-08-23-rep-seqs.qza"
)

TBL_ITS=(
  "${PROCESSED_DATA_DIR}/imported_seqs/imported-paired-end-seqs-NN21-ITS-trimmed-2024-08-23-table.qza"
  "${PROCESSED_DATA_DIR}/imported_seqs/imported-paired-end-seqs-NN22-ITS-trimmed-2024-08-23-table.qza"
  "${PROCESSED_DATA_DIR}/imported_seqs/imported-paired-end-seqs-NN24-ITS-trimmed-2024-08-23-table.qza"
)

mkdir -p "${PROCESSED_DATA_DIR}/imported_seqs"

args=()
for f in "${REP_16S[@]}"; do args+=(--i-data "$f"); done
qiime feature-table merge-seqs "${args[@]}" \
  --o-merged-data "${PROCESSED_DATA_DIR}/imported_seqs/merged-rep-seqs-16S.qza"

args=()
for f in "${TBL_16S[@]}"; do args+=(--i-tables "$f"); done
qiime feature-table merge "${args[@]}" \
  --p-overlap-method sum \
  --o-merged-table "${PROCESSED_DATA_DIR}/imported_seqs/merged-feature-table-16S.qza"

args=()
for f in "${REP_ITS[@]}"; do args+=(--i-data "$f"); done
qiime feature-table merge-seqs "${args[@]}" \
  --o-merged-data "${PROCESSED_DATA_DIR}/imported_seqs/merged-rep-seqs-ITS.qza"

args=()
for f in "${TBL_ITS[@]}"; do args+=(--i-tables "$f"); done
qiime feature-table merge "${args[@]}" \
  --p-overlap-method sum \
  --o-merged-table "${PROCESSED_DATA_DIR}/imported_seqs/merged-feature-table-ITS.qza"

