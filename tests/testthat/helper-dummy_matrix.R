# Shared dummy RNAseq matrix used across all tests.
# 6 genes x 6 samples (3 early / 3 late).
# Values are chosen so expected outputs can be verified analytically
# and compared against the Python reference implementation.
#
# Expected behaviour after create_index(group1=early, group2=late):
#   GeneA -> UP,   valence = log2(15)   (high in late)
#   GeneB -> DOWN, valence = log2(1/15) (high in early)
#   GeneC -> (not significant, excluded)
#   GeneD -> (zero-expressed, dropped before testing)
#   GeneE -> UP,   valence = log2(7)
#   GeneF -> DOWN, valence = log2(1/7)

dummy_mat <- matrix(
  c(
    1.0, 2.0, 1.5,  8.0, 9.0, 8.5,   # GeneA
    8.0, 9.0, 8.5,  1.0, 2.0, 1.5,   # GeneB
    5.0, 5.5, 4.5,  5.2, 4.8, 5.1,   # GeneC
    0.0, 0.0, 0.0,  0.0, 0.0, 0.0,   # GeneD (zero)
    3.0, 4.0, 3.5,  6.0, 7.0, 6.5,   # GeneE
    7.0, 6.0, 6.5,  4.0, 3.0, 3.5    # GeneF
  ),
  nrow = 6, byrow = TRUE,
  dimnames = list(
    c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF"),
    c("early1", "early2", "early3", "late1", "late2", "late3")
  )
)

early_cols <- c("early1", "early2", "early3")
late_cols  <- c("late1",  "late2",  "late3")
