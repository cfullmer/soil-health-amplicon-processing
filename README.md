# Cleaned Manuscript Workflow Scripts

This directory contains the manuscript-facing bioinformatics workflow scripts for the DSFAS soil health microbiome analyses. The scripts are distilled from the original working notebooks and workflow notes into a smaller, reviewer-readable supplement.

These files are intended to document and reproduce the core QIIME2 processing logic without including exploratory notebooks, local scratch outputs, or large intermediate artifacts.

## Contents

| File | Purpose |
| --- | --- |
| `00_paths_template.env` | Template for local paths, metadata files, and classifier artifacts. Copy this to `paths.env` before running any script. |
| `01_initial_batch_processing.sh` | Builds and checks manifests, imports paired-end FASTQs into QIIME2, summarizes demultiplexed data, trims primers, runs per-batch DADA2, and tabulates denoising statistics. |
| `02_merge_initial_batches.sh` | Merges the initial 16S and ITS per-batch representative sequences and feature tables. |
| `03_postmerge_16S_workflow.sh` | Assigns 16S taxonomy, filters to prokaryote-only features, and exports the final non-rarefied 16S table and representative sequences. |
| `04_merged_ITS_workflow.sh` | Merges ITS batches, clusters features at 97% identity, performs UNITE-based sequence QC, assigns taxonomy, filters to fungi-only features, and exports final representative sequences. |

The parent directory includes `FILE_MANIFEST.tsv`, which maps each cleaned script back to the original source notebook or workflow document.

## Expected Environment

Run these scripts from an environment with the command-line tools used in the workflow:

- QIIME2 with the required plugins: `cutadapt`, `dada2`, `demux`, `feature-classifier`, `feature-table`, `quality-control`, `taxa`, `tools`, and `vsearch`
- R with `tidyverse`
- Python 3
- BIOM command-line tools
- Standard Unix shell utilities such as `cp`, `tr`, `grep`, and `mkdir`

The scripts use `bash` with `set -euo pipefail`, so missing variables, missing files, and command failures should stop execution.

## Setup

Copy the path template and edit it for the machine where the workflow will run:

```bash
cp 00_paths_template.env paths.env
```

At minimum, update these variables in `paths.env`:

- `PROJECT_ROOT`
- `RAW_DATA_DIR`
- `PROCESSED_DATA_DIR`
- `CLASSIFIER_16S`
- `CLASSIFIER_ITS`
- `UNITE_ALL_SEQS`
- `METADATA_16S`
- `METADATA_ITS`

Several scripts also define batch-specific input artifact arrays or truncation lengths inside the script body. Review those values before running, especially when applying the scripts outside the exact manuscript dataset.

## Recommended Run Order

Run the scripts from this directory after activating the appropriate QIIME2 environment:

```bash
cd /path/to/bioinformatics/scripts/manuscript_submission_package/cleaned_scripts
```

1. Prepare paths:

```bash
cp 00_paths_template.env paths.env
# edit paths.env
```

2. Process initial batches:

```bash
bash 01_initial_batch_processing.sh
```

This generates updated manifests, manifest QC reports, imported QIIME2 artifacts, demultiplexed summaries, trimmed sequence artifacts, per-batch DADA2 outputs, and denoising-stat visualizations.

3. Merge initial 16S and ITS batches:

```bash
bash 02_merge_initial_batches.sh
```

This creates initial merged representative-sequence and feature-table artifacts for 16S and ITS.

4. Run post-merge 16S processing:

```bash
bash 03_postmerge_16S_workflow.sh
```

This assigns SILVA-based taxonomy, removes non-target 16S features such as unassigned, mitochondrial, chloroplast, and eukaryotic features, then exports the non-rarefied prokaryote-only table and representative sequences.

5. Run merged ITS processing:

```bash
bash 04_merged_ITS_workflow.sh
```

This merges ITS batch artifacts, clusters at 97% identity, screens representative sequences against UNITE, assigns taxonomy, filters to fungi-only features, and exports final representative sequences.

## Key Inputs

The scripts expect the following project inputs to exist after `paths.env` is configured:

- Raw demultiplexed FASTQ files and manifest CSV files under `RAW_DATA_DIR`
- Master sample database used to update sample IDs
- Expected 16S and ITS sample lists for manifest QC
- Per-batch QIIME2 DADA2 artifacts for merge steps where they are not generated in the same run
- 16S classifier: SILVA 138 99% 515F/806R classifier
- ITS classifier: UNITE v10 dynamic classifier
- UNITE all-sequence reference artifact for ITS QC
- 16S and ITS metadata TSV files for taxonomy barplots and downstream summaries

## Main Outputs

Output locations are controlled by variables in `paths.env` and by script-level defaults. Important output groups include:

- Updated and combined manifests in `MANIFEST_OUT_DIR`
- Manifest QC reports: unexpected samples, missing samples, and duplicate rows
- Imported, trimmed, and denoised QIIME2 artifacts in `IMPORTED_DIR`
- Demultiplexed summaries in `DEMUX_SUMMARY_DIR`
- Denoising-stat visualizations in `DENOISE_SUMMARY_DIR`
- Merged 16S and ITS feature tables and representative sequences
- 16S taxonomy, prokaryote-only filtered tables, exported BIOM/TSV tables, and representative sequences
- ITS merged, cluster97, UNITE-QC, taxonomy, fungi-only filtered, and exported representative-sequence artifacts

## Notes for Reuse

These scripts are cleaned supplement scripts, not the full exploratory working history. They preserve the core analysis commands and parameters used for the manuscript workflow while replacing local machine paths with configurable variables.

Before reusing the scripts on a new dataset:

- Confirm primer patterns in `01_initial_batch_processing.sh`
- Confirm DADA2 truncation lengths in `01_initial_batch_processing.sh`
- Confirm batch artifact arrays in `02_merge_initial_batches.sh` and `04_merged_ITS_workflow.sh`
- Confirm taxonomy classifiers and metadata paths in `paths.env`
- Confirm output directories do not overwrite existing analysis results unintentionally

For provenance, use `../FILE_MANIFEST.tsv` to trace each cleaned script to its original notebook or workflow source.
