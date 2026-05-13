#!/usr/bin/env bash
set -euo pipefail

# Manuscript-facing representation of the initial multi-batch processing workflow.
# This consolidates the relevant logic from the original manifest, import, summary,
# primer-trimming, and denoising-stat notebooks/scripts into one lean supplement file.
#
# Before running:
# 1. Copy `00_paths_template.env` to `paths.env`
# 2. Edit the paths and expected-sample files for your environment
# 3. Review the batch-specific primer settings below

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ ! -f "${SCRIPT_DIR}/paths.env" ]]; then
  echo "Missing ${SCRIPT_DIR}/paths.env. Copy 00_paths_template.env to paths.env and edit it first." >&2
  exit 1
fi
source "${SCRIPT_DIR}/paths.env"

MANIFEST_SEARCH_ROOT="${MANIFEST_SEARCH_ROOT:-${RAW_DATA_DIR}}"
MANIFEST_OUT_DIR="${MANIFEST_OUT_DIR:-${PROCESSED_DATA_DIR}/manifests}"
IMPORTED_DIR="${IMPORTED_DIR:-${PROCESSED_DATA_DIR}/imported_seqs}"
DEMUX_SUMMARY_DIR="${DEMUX_SUMMARY_DIR:-${PROCESSED_DATA_DIR}/demux_summaries}"
DENOISE_SUMMARY_DIR="${DENOISE_SUMMARY_DIR:-${PROCESSED_DATA_DIR}/denoising_stats}"
MASTER_DB="${MASTER_DB:-${PROCESSED_DATA_DIR}/HSH_MASTER_DATABASEallsamples - Master Database.csv}"
EXPECTED_16S="${EXPECTED_16S:-${PROCESSED_DATA_DIR}/ALL_SH_MICRO_SAMPLES - Entire16S.csv}"
EXPECTED_ITS="${EXPECTED_ITS:-${PROCESSED_DATA_DIR}/ALL_SH_MICRO_SAMPLES - EntireITS.csv}"

mkdir -p "${MANIFEST_OUT_DIR}" "${IMPORTED_DIR}" "${DEMUX_SUMMARY_DIR}" "${DENOISE_SUMMARY_DIR}"

# Primer patterns preserved from the original trimming notebooks.
# These are linked-primer patterns, not plain 5' primer sequences.
PRIMER_16S_F='^GTGYCAGCMGCCGCGGTAA...ATTAGAWACCCBNGTAGTCC'
PRIMER_16S_R='^GGACTACNVGGGTWTCTAAT...TTACCGCGGCKGCTGRCAC'
PRIMER_ITS_F='^CTTGGTCATTTAGAGGAAGTAA...GCATCGATGAAGAACGCAGC'
PRIMER_ITS_R='^GCTGCGTTCTTCATCGATGC...TTACTTCCTCTAAATGACCAAG'

echo "== Step 1: build updated manifests from demultiplexed FASTQs =="
MANIFEST_SEARCH_ROOT="${MANIFEST_SEARCH_ROOT}" \
MASTER_DB="${MASTER_DB}" \
MANIFEST_OUT_DIR="${MANIFEST_OUT_DIR}" \
Rscript - <<'RSCRIPT'
library(tidyverse)

manifests_dir_path <- Sys.getenv("MANIFEST_SEARCH_ROOT")
master_db_path <- Sys.getenv("MASTER_DB")
manifest_out_dir <- Sys.getenv("MANIFEST_OUT_DIR")

create_var_name <- function(path, base_dir) {
  relative_path <- str_replace(path, paste0("^", base_dir, "/"), "")
  subdirectory <- dirname(relative_path)
  file_name <- tools::file_path_sans_ext(basename(relative_path))

  paste0(subdirectory, "_", file_name) %>%
    str_replace_all("[/\\\\]", "_") %>%
    str_replace_all("[^[:alnum:]_]", "") %>%
    str_replace_all("^[^[:alpha:]]+", "")
}

csv_files <- list.files(
  manifests_dir_path,
  pattern = "(?i)(manifest).*\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

read_csv_headers_first_row <- function(file) {
  read_csv(file, show_col_types = FALSE)
}

read_csv_skip_first_row <- function(file) {
  read_csv(file, skip = 1, show_col_types = FALSE)
}

manifests <- list()
for (file in csv_files) {
  var_name <- create_var_name(file, manifests_dir_path)
  tibble_data <- read_csv_headers_first_row(file)
  colnames(tibble_data) <- tolower(colnames(tibble_data))

  if (all(c("sample-id", "absolute-filepath", "direction") %in% colnames(tibble_data))) {
    manifests[[var_name]] <- tibble_data
  }
}

master_db <- read_csv_skip_first_row(master_db_path)
colnames(master_db) <- tolower(colnames(master_db))

update_sample_ids <- function(manifest, master_db) {
  manifest %>%
    left_join(master_db, by = c("sample-id" = "old_sh_1")) %>%
    mutate(`sample-id` = ifelse(is.na(sh_1), `sample-id`, sh_1)) %>%
    mutate(`sample-id` = str_replace(`sample-id`, "^SH-", "")) %>%
    mutate(
      `sample-id` = ifelse(
        startsWith(`sample-id`, "20") & nchar(`sample-id`) > 9,
        substr(`sample-id`, 1, nchar(`sample-id`) - 4),
        `sample-id`
      )
    ) %>%
    select(`sample-id`, `absolute-filepath`, `direction`)
}

updated_manifests <- lapply(manifests, function(manifest) {
  update_sample_ids(manifest, master_db)
})

invisible(lapply(names(updated_manifests), function(name) {
  write_csv(updated_manifests[[name]], file.path(manifest_out_dir, paste0(name, "_updated.csv")))
}))

its_manifests <- bind_rows(lapply(names(updated_manifests), function(name) {
  if (str_detect(name, "ITS")) updated_manifests[[name]]
}))

s16_manifests <- bind_rows(lapply(names(updated_manifests), function(name) {
  if (str_detect(name, "16S")) updated_manifests[[name]]
}))

write_csv(its_manifests, file.path(manifest_out_dir, "ITS_combined_manifest.csv"))
write_csv(s16_manifests, file.path(manifest_out_dir, "16S_combined_manifest.csv"))
RSCRIPT

echo "== Step 2: QC manifests against expected sample lists =="
MANIFEST_OUT_DIR="${MANIFEST_OUT_DIR}" \
EXPECTED_16S="${EXPECTED_16S}" \
EXPECTED_ITS="${EXPECTED_ITS}" \
Rscript - <<'RSCRIPT'
library(tidyverse)

manifest_out_dir <- Sys.getenv("MANIFEST_OUT_DIR")
expected_16s <- Sys.getenv("EXPECTED_16S")
expected_its <- Sys.getenv("EXPECTED_ITS")

bind_csv_files_with_filename <- function(file_list) {
  file_list %>%
    map_dfr(~ read_csv(.x, show_col_types = FALSE) %>% mutate(filename = .x))
}

clean_expected <- function(path) {
  read.csv(path, header = FALSE)[, 1] %>%
    unique() %>%
    as_tibble() %>%
    rename(value = value) %>%
    mutate(
      value = str_remove(value, "&.*$"),
      value = str_replace(value, "/.*$", ""),
      value = str_remove(value, "SH-"),
      value = str_replace_all(value, " ", "-"),
      value = str_replace(value, "(\\d{4}-\\d{3})-\\d{3}", "\\1"),
      value = str_replace(value, "(2019-)([A-Z]+)", function(x) tolower(x))
    )
}

run_manifest_qc <- function(pattern, expected_path, prefix) {
  csv_files <- list.files(manifest_out_dir, pattern = pattern, full.names = TRUE)
  csv_files <- csv_files[!grepl("combined_manifest", basename(csv_files))]
  if (!length(csv_files)) return(invisible(NULL))

  all_samples <- bind_csv_files_with_filename(csv_files) %>%
    select("sample-id") %>%
    unique() %>%
    rename(value = "sample-id")

  expected_cleaned <- clean_expected(expected_path)

  only_in_manifest <- setdiff(all_samples, expected_cleaned)
  only_in_expected <- setdiff(expected_cleaned, all_samples)
  duplicates <- bind_csv_files_with_filename(csv_files) %>%
    group_by(`sample-id`, filename) %>%
    filter(n() > 2) %>%
    ungroup()

  write_tsv(only_in_manifest, file.path(manifest_out_dir, paste0(prefix, "_unexpected_samples.tsv")))
  write_tsv(only_in_expected, file.path(manifest_out_dir, paste0(prefix, "_missing_samples.tsv")))
  write_tsv(duplicates, file.path(manifest_out_dir, paste0(prefix, "_duplicate_rows.tsv")))
}

run_manifest_qc("16S.*updated.*\\.csv$", expected_16s, "16S")
run_manifest_qc("ITS.*updated.*\\.csv$", expected_its, "ITS")
RSCRIPT

echo "== Step 3: import updated manifests into QIIME2 =="
shopt -s nullglob globstar
for manifest in "${MANIFEST_OUT_DIR}"/*updated*.csv; do
  name="$(basename "${manifest}")"
  batch="$(echo "${name}" | grep -o 'NN[0-9][0-9]' | head -1 || true)"
  assay="$(echo "${name}" | grep -oE '16S|ITS' | head -1 || true)"

  [[ -n "${batch}" && -n "${assay}" ]] || continue

  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "${manifest}" \
    --input-format PairedEndFastqManifestPhred33 \
    --output-path "${IMPORTED_DIR}/imported-paired-end-seqs-${batch}-${assay}.qza"
done

echo "== Step 4: summarize imported demux artifacts =="
for input_file in "${IMPORTED_DIR}"/imported-paired-end-seqs-*.qza; do
  [[ -e "${input_file}" ]] || continue
  base="$(basename "${input_file}" .qza)"
  qiime demux summarize \
    --i-data "${input_file}" \
    --o-visualization "${DEMUX_SUMMARY_DIR}/${base}.qzv"
done

echo "== Step 5: trim primers =="
for input_file in "${IMPORTED_DIR}"/imported-paired-end-seqs-*.qza; do
  [[ -e "${input_file}" ]] || continue
  base="$(basename "${input_file}" .qza)"
  output_file="${IMPORTED_DIR}/${base}-trimmed.qza"

  if [[ "${base}" == *"16S"* ]]; then
    qiime cutadapt trim-paired \
      --i-demultiplexed-sequences "${input_file}" \
      --p-adapter-f "${PRIMER_16S_F}" \
      --p-adapter-r "${PRIMER_16S_R}" \
      --p-cores 4 \
      --p-discard-untrimmed \
      --p-overlap 10 \
      --p-minimum-length 50 \
      --p-error-rate 0.13 \
      --o-trimmed-sequences "${output_file}" \
      --verbose
  elif [[ "${base}" == *"ITS"* ]]; then
    qiime cutadapt trim-paired \
      --i-demultiplexed-sequences "${input_file}" \
      --p-adapter-f "${PRIMER_ITS_F}" \
      --p-adapter-r "${PRIMER_ITS_R}" \
      --p-cores 4 \
      --p-discard-untrimmed \
      --p-overlap 10 \
      --p-minimum-length 50 \
      --p-error-rate 0.13 \
      --o-trimmed-sequences "${output_file}" \
      --verbose
  fi
done

cat <<'EOF'
Special-case note:
Some 16S batches used a multiple-primer workflow rather than the single-primer loop above.
In the original notebooks, those batches were handled by importing the manifest, splitting or
grouping by primer set, and running `qiime cutadapt trim-paired` with primer-specific adapter
patterns for each subset before continuing into DADA2.
EOF

echo "== Step 6: run DADA2 on trimmed per-batch artifacts =="
# Edit these truncation lengths to match the manuscript runs for each batch/artifact.
# In the original workflow, these values were chosen after inspecting the trimmed
# demux summaries for each batch and then entered through a prompt-driven notebook.
# They are batch-specific because read quality profiles differed across runs and assays.
declare -A TRUNC_LEN_F=(
  ["imported-paired-end-seqs-NN18-16S-trimmed"]=222
  ["imported-paired-end-seqs-NN19-16S-trimmed"]=198
  ["imported-paired-end-seqs-NN21-16S-trimmed"]=231
  ["imported-paired-end-seqs-NN25-16S-trimmed"]=216
  ["imported-paired-end-seqs-NN21-ITS-trimmed"]=202
  ["imported-paired-end-seqs-NN22-ITS-trimmed"]=188
  ["imported-paired-end-seqs-NN24-ITS-trimmed"]=183
)

declare -A TRUNC_LEN_R=(
  ["imported-paired-end-seqs-NN18-16S-trimmed"]=176
  ["imported-paired-end-seqs-NN19-16S-trimmed"]=118
  ["imported-paired-end-seqs-NN21-16S-trimmed"]=176
  ["imported-paired-end-seqs-NN25-16S-trimmed"]=143
  ["imported-paired-end-seqs-NN21-ITS-trimmed"]=123
  ["imported-paired-end-seqs-NN22-ITS-trimmed"]=145
  ["imported-paired-end-seqs-NN24-ITS-trimmed"]=140
)

for trimmed_qza in "${IMPORTED_DIR}"/imported-paired-end-seqs-*-trimmed*.qza; do
  [[ -e "${trimmed_qza}" ]] || continue
  base="$(basename "${trimmed_qza}" .qza)"
  key="${base%-2024-*}"
  trunc_f="${TRUNC_LEN_F[$key]:-}"
  trunc_r="${TRUNC_LEN_R[$key]:-}"

  if [[ -z "${trunc_f}" || -z "${trunc_r}" ]]; then
    echo "Skipping ${trimmed_qza}: no truncation lengths configured for ${key}" >&2
    continue
  fi

  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "${trimmed_qza}" \
    --p-trim-left-f 0 \
    --p-trunc-len-f "${trunc_f}" \
    --p-trim-left-r 0 \
    --p-trunc-len-r "${trunc_r}" \
    --p-n-threads 0 \
    --p-n-reads-learn 1000000 \
    --o-representative-sequences "${IMPORTED_DIR}/${base}-rep-seqs.qza" \
    --o-table "${IMPORTED_DIR}/${base}-table.qza" \
    --o-denoising-stats "${IMPORTED_DIR}/${base}-denoising-stats.qza"
done

echo "== Step 7: summarize denoising statistics after DADA2 =="
for stats_qza in "${IMPORTED_DIR}"/**/*denoising-stats*.qza "${IMPORTED_DIR}"/*denoising-stats*.qza; do
  [[ -e "${stats_qza}" ]] || continue
  output_file="${DENOISE_SUMMARY_DIR}/$(basename "${stats_qza}" .qza).qzv"
  qiime metadata tabulate \
    --m-input-file "${stats_qza}" \
    --o-visualization "${output_file}"
done

cat <<'EOF'
Step 8: Generate per-batch DADA2 artifacts.
The original workflow then produced, for each batch and assay:
- feature table (`table.qza`)
- representative sequences (`rep-seqs.qza`)
- denoising stats (`denoising-stats.qza`)

Step 9: Hand off those per-batch artifacts to `02_merge_initial_batches.sh`, which performs:
- `qiime feature-table merge-seqs`
- `qiime feature-table merge`
for the initial 16S and ITS merged datasets.
EOF
