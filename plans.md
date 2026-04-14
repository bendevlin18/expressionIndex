# Plans

## Index Creation Function

Currently the package only supports calculating scores from pre-built indices (`mdi`, `lps_index`, `fsi_index`). The goal is to add a function that lets users **create their own index** from expression data.

The output will be a data frame with `gene`, `direction`, and `valence` columns — identical in format to the built-in index datasets — so it drops directly into `calculate_index()` and the rest of the existing pipeline.

### Reference Implementation

The Python version (`python/developmental_index.py`) defines the full workflow as a set of modular functions:

1. **`scale_expression()`** — min-max scales each gene's expression across all samples so every gene contributes equal weight
2. **`drop_unexpressed_genes()`** — removes genes with NaN values
3. **`identify_significant_genes()`** — runs a t-test between two sample groups (e.g. early vs. late), computes log2 fold change, and assigns `"UP"` / `"DOWN"` / `"N/A"` direction per gene (p < 0.05 threshold)
4. **`remove_insignificant_rows()`** — filters the data frame down to only significantly regulated genes
5. **`extract_regulated_genes()`** — returns the final index data frame with `gene`, `direction`, and `valence` (log2FC) columns, sorted by valence

### R Implementation Plan

Translate the Python pipeline into R, following the same modular structure. All functions go in `R/index_calculations.R`.

#### Helper functions (unexported)

- **`scale_expression(mat)`** — min-max scale each row (gene) across columns (samples). R equivalent: `(x - min(x)) / (max(x) - min(x))` applied per row via `apply(mat, 1, ...)`.
- **`drop_unexpressed_genes(mat)`** — remove rows where any value is `NA`, or where the row sum is 0.

#### Main exported function

```r
create_index(mat, group1, group2, p_threshold = 0.05, gene_column = "gene",
             direction_column = "direction", valence_column = "valence")
```

- `mat` — genes × samples expression matrix (pre-normalized: TPM, FPKM, VST, etc.)
- `group1` — character vector of column names representing the "low/early" state
- `group2` — character vector of column names representing the "high/late" state
- `p_threshold` — significance cutoff, default 0.05

**Steps inside `create_index()`:**

1. Call `drop_unexpressed_genes()` then `scale_expression()`
2. For each gene, run `stats::t.test()` between `group1` and `group2` columns; record p-value and t-statistic
3. Compute log2 fold change: `log2(mean(group2) / mean(group1))`
4. Assign direction: `"UP"` if log2FC > 0 and p < threshold, `"DOWN"` if log2FC < 0 and p < threshold, otherwise drop the gene
5. Return a data frame with columns `gene`, `direction`, `valence` (log2FC), sorted descending by `valence` — ready to pass to `calculate_index()`

#### Notes

- Use `stats::t.test()` (base R) — no new dependencies needed
- The `valence` column stores log2FC, matching the format of `mdi` and `lps_index`
- `fsi_index` uses a capital-G `Gene` column; the new function should output lowercase `gene` to match `mdi`/`lps_index` and the defaults in `calculate_index()`
