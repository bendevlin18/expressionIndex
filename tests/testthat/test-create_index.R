### Tests for create_index()
### Verifies parity with the Python pipeline in python/developmental_index.py:
###   scale_expression -> identify_significant_genes -> extract_regulated_genes
###
### Python uses scipy.stats.ttest_ind (equal_var=True by default).
### R uses stats::t.test(var.equal=TRUE) to match.
###
### Pre-computed expected values for dummy_mat (see helper-dummy_matrix.R):
###   After scaling, group means are:
###     GeneA: early=0.0625,  late=0.9375  -> log2FC = log2(15)   ~ 3.907
###     GeneB: early=0.9375,  late=0.0625  -> log2FC = log2(1/15) ~ -3.907
###     GeneC: early=0.5,     late=0.5333  -> small FC, not significant
###     GeneD: dropped (zero)
###     GeneE: early=0.125,   late=0.875   -> log2FC = log2(7)    ~ 2.807
###     GeneF: early=0.875,   late=0.125   -> log2FC = log2(1/7)  ~ -2.807

test_that("create_index returns a data frame with correct columns", {
  result <- create_index(dummy_mat, early_cols, late_cols)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("gene", "direction", "valence") %in% colnames(result)))
})

test_that("create_index identifies UP genes correctly", {
  result <- create_index(dummy_mat, early_cols, late_cols)
  up_genes <- result$gene[result$direction == "UP"]
  expect_true("GeneA" %in% up_genes)
  expect_true("GeneE" %in% up_genes)
})

test_that("create_index identifies DOWN genes correctly", {
  result <- create_index(dummy_mat, early_cols, late_cols)
  down_genes <- result$gene[result$direction == "DOWN"]
  expect_true("GeneB" %in% down_genes)
  expect_true("GeneF" %in% down_genes)
})

test_that("create_index excludes zero-expressed and non-significant genes", {
  result <- create_index(dummy_mat, early_cols, late_cols)
  expect_false("GeneD" %in% result$gene)   # zero-expressed, dropped
  expect_false("GeneC" %in% result$gene)   # not significant
})

test_that("create_index valence matches expected log2FC values", {
  result <- create_index(dummy_mat, early_cols, late_cols)

  geneA_val <- result$valence[result$gene == "GeneA"]
  geneB_val <- result$valence[result$gene == "GeneB"]
  geneE_val <- result$valence[result$gene == "GeneE"]
  geneF_val <- result$valence[result$gene == "GeneF"]

  expect_equal(geneA_val,  log2(15), tolerance = 1e-10)
  expect_equal(geneB_val, -log2(15), tolerance = 1e-10)
  expect_equal(geneE_val,  log2(7),  tolerance = 1e-10)
  expect_equal(geneF_val, -log2(7),  tolerance = 1e-10)
})

test_that("create_index output is sorted descending by valence", {
  result <- create_index(dummy_mat, early_cols, late_cols)
  expect_true(all(diff(result$valence) <= 0))
})

test_that("create_index direction column contains only UP or DOWN", {
  result <- create_index(dummy_mat, early_cols, late_cols)
  expect_true(all(result$direction %in% c("UP", "DOWN")))
})

test_that("create_index output plugs into calculate_index without error", {
  index_df <- create_index(dummy_mat, early_cols, late_cols)

  # Build a small target matrix using only the genes returned by create_index
  target_mat <- dummy_mat[index_df$gene, , drop = FALSE]

  expect_no_error(
    scores <- calculate_index(target_mat, index_df)
  )
  expect_s3_class(scores, "data.frame")
  expect_equal(colnames(scores), "index")
  expect_equal(nrow(scores), ncol(target_mat))
})

test_that("create_index respects custom column name arguments", {
  result <- create_index(dummy_mat, early_cols, late_cols,
                         gene_column      = "symbol",
                         direction_column = "reg",
                         valence_column   = "log2fc")
  expect_true(all(c("symbol", "reg", "log2fc") %in% colnames(result)))
  expect_false(any(c("gene", "direction", "valence") %in% colnames(result)))
})

test_that("create_index respects p_threshold", {
  # Very stringent threshold — only the most extreme genes should survive
  strict  <- create_index(dummy_mat, early_cols, late_cols, p_threshold = 0.001)
  lenient <- create_index(dummy_mat, early_cols, late_cols, p_threshold = 0.05)
  expect_lte(nrow(strict), nrow(lenient))
})
