### Tests for unexported helper functions: scale_expression, drop_unexpressed_genes
### These match the behaviour of scale_expression() and drop_unexpressed_genes()
### in python/developmental_index.py.

# ---- scale_expression -------------------------------------------------------

test_that("scale_expression min-max scales each gene row to [0, 1]", {
  scaled <- expressionIndex:::scale_expression(dummy_mat)

  # every row should have a minimum of 0 and maximum of 1
  row_mins <- apply(scaled, 1, min)
  row_maxs <- apply(scaled, 1, max)

  expect_true(all(abs(row_mins) < 1e-10),
              info = "all row minimums should be 0 after scaling")
  expect_true(all(abs(row_maxs - 1) < 1e-10),
              info = "all row maximums should be 1 after scaling")
})

test_that("scale_expression produces correct values for GeneA", {
  # GeneA: (1, 2, 1.5, 8, 9, 8.5) -> min=1, max=9, range=8
  expected <- c(0, 1/8, 0.5/8, 7/8, 1, 7.5/8)
  scaled   <- expressionIndex:::scale_expression(dummy_mat)

  expect_equal(as.numeric(scaled["GeneA", ]), expected, tolerance = 1e-10)
})

test_that("scale_expression produces correct values for GeneC", {
  # GeneC: (5, 5.5, 4.5, 5.2, 4.8, 5.1) -> min=4.5, max=5.5, range=1
  expected <- c(0.5, 1.0, 0.0, 0.7, 0.3, 0.6)
  scaled   <- expressionIndex:::scale_expression(dummy_mat)

  expect_equal(as.numeric(scaled["GeneC", ]), expected, tolerance = 1e-10)
})

test_that("scale_expression returns all-zero row for a constant gene", {
  const_mat <- matrix(c(5, 5, 5, 5), nrow = 1,
                      dimnames = list("ConstGene", NULL))
  result <- expressionIndex:::scale_expression(const_mat)
  expect_equal(as.numeric(result), c(0, 0, 0, 0))
})

test_that("scale_expression preserves row and column names", {
  scaled <- expressionIndex:::scale_expression(dummy_mat)
  expect_equal(rownames(scaled), rownames(dummy_mat))
  expect_equal(colnames(scaled), colnames(dummy_mat))
})

# ---- drop_unexpressed_genes -------------------------------------------------

test_that("drop_unexpressed_genes removes all-zero rows", {
  result <- expressionIndex:::drop_unexpressed_genes(dummy_mat)
  expect_false("GeneD" %in% rownames(result))
  expect_equal(nrow(result), 5)
})

test_that("drop_unexpressed_genes removes rows containing NA", {
  mat_with_na        <- dummy_mat
  mat_with_na["GeneC", "early1"] <- NA
  result <- expressionIndex:::drop_unexpressed_genes(mat_with_na)
  expect_false("GeneC" %in% rownames(result))
})

test_that("drop_unexpressed_genes keeps expressed genes intact", {
  result <- expressionIndex:::drop_unexpressed_genes(dummy_mat)
  expect_true(all(c("GeneA", "GeneB", "GeneC", "GeneE", "GeneF") %in% rownames(result)))
})
