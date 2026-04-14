# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

`expressionIndex` is an R package for calculating transcriptional index scores from gene expression data (single-cell, spatial, and bulk RNA-seq). The core idea: an index score = mean(UP genes) / mean(DOWN genes) for a curated gene signature.

## Common Commands

```r
# Install dependencies and load
devtools::install_deps()
library(expressionIndex)

# Regenerate documentation from roxygen2 comments
devtools::document()

# Run package checks
devtools::check()

# Run tests (if added)
devtools::test()

# Install locally
devtools::install()
```

## Architecture

All functions live in two R files:

- **`R/index_calculations.R`** — all exported functions
- **`R/data.R`** — roxygen2 documentation stubs for the three built-in datasets

**Function hierarchy:**
- `calculate_index()` — core workhorse; operates on a genes × samples matrix, returns a data frame of index scores
- `calculate_index_seurat()` and `calculate_index_deseq()` — thin wrappers that extract a matrix from their respective objects, call `calculate_index()`, and write results back to object metadata
- `seurat_gene_var()` — computes per-cell gene variance using `sparseMatrixStats::colVars`, adds to Seurat metadata
- `rolling_variance_index()` — orders cells by index, computes rolling window averages of index and variance; requires both `calculate_index_seurat()` and `seurat_gene_var()` to have been run first
- `find_transition_genes()` — chunked Spearman rank correlation of each gene against an index column; chunk processing avoids memory overload
- `calculate_tpm.matrix()` / `calculate_tpm.DESeqDataSet()` — TPM calculation utilities (RPK → per-sample scaling)

**Built-in index datasets** (stored as `.rda` in `data/`):
- `mdi` — Microglia Development Index
- `lps_index` — LPS Inflammatory Microglia Index
- `fsi_index` — Fast-Spiking Interneuron Index

Each dataset is a data frame with columns `gene`, `direction` (`"UP"`/`"DOWN"`), and `valence`. Note: `fsi_index` uses `Gene` (capital G) as the gene column name, unlike `mdi` and `lps_index`.

## Key Dependencies

- `Seurat` — single-cell object access via `GetAssayData(..., layer = "data")`
- `DESeq2` / `SummarizedExperiment` — bulk RNA-seq object support; VST normalization via `DESeq2::vst(dds, blind = FALSE)`
- `sparseMatrixStats` — efficient sparse matrix column variance
- `matrixStats` — `rowRanks()` used in chunked correlation
- `Matrix` — sparse matrix row means for expression filtering
- `stringr` — `str_to_title()` to normalize gene name casing

## Documentation

Documentation is generated with roxygen2. After editing function docs in `R/index_calculations.R`, run `devtools::document()` to regenerate `man/*.Rd` and `NAMESPACE`. Do not manually edit `NAMESPACE` or files in `man/`.
