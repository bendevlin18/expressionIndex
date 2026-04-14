"""
Python parity script for expressionIndex.

Runs the same dummy matrix through the original Python pipeline
(python/developmental_index.py) and prints the expected outputs so they can
be compared against the R testthat results.

Run from the repo root:
    python tests/python_parity.py

Expected output (matches R test assertions in tests/testthat/test-create_index.R):
    GeneA  UP    log2(15)  ~  3.9069
    GeneE  UP    log2(7)   ~  2.8074
    GeneF  DOWN  log2(1/7) ~ -2.8074
    GeneB  DOWN  log2(1/15)~ -3.9069
"""

import sys
import os
import math
import numpy as np
import pandas as pd

# ── import the reference Python functions ────────────────────────────────────
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
from developmental_index import (
    scale_expression,
    drop_unexpressed_genes,
    identify_significant_genes,
    remove_insignificant_rows,
    extract_regulated_genes,
)

# ── build the same dummy matrix used in the R tests ──────────────────────────
data = {
    "early1": [1.0, 8.0, 5.0, 0.0, 3.0, 7.0],
    "early2": [2.0, 9.0, 5.5, 0.0, 4.0, 6.0],
    "early3": [1.5, 8.5, 4.5, 0.0, 3.5, 6.5],
    "late1":  [8.0, 1.0, 5.2, 0.0, 6.0, 4.0],
    "late2":  [9.0, 2.0, 4.8, 0.0, 7.0, 3.0],
    "late3":  [8.5, 1.5, 5.1, 0.0, 6.5, 3.5],
}
genes = ["GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF"]
df    = pd.DataFrame(data, index=genes)

early_cols = ["early1", "early2", "early3"]
late_cols  = ["late1",  "late2",  "late3"]

print("=" * 60)
print("STEP 1 — drop_unexpressed_genes")
print("=" * 60)
df = drop_unexpressed_genes(df)
print(f"Remaining genes: {list(df.index)}")
assert "GeneD" not in df.index, "GeneD should have been dropped"

print()
print("=" * 60)
print("STEP 2 — scale_expression")
print("=" * 60)
df = scale_expression(df)

# Verify known scaled values for GeneA
geneA_scaled = df.loc["GeneA"].values
expected_geneA = np.array([0, 1/8, 0.5/8, 7/8, 1, 7.5/8])
assert np.allclose(geneA_scaled, expected_geneA), \
    f"GeneA scaling mismatch:\n  got:      {geneA_scaled}\n  expected: {expected_geneA}"
print("GeneA scaled:", np.round(geneA_scaled, 6))

geneC_scaled = df.loc["GeneC"].values
expected_geneC = np.array([0.5, 1.0, 0.0, 0.7, 0.3, 0.6])
assert np.allclose(geneC_scaled, expected_geneC), \
    f"GeneC scaling mismatch:\n  got:      {geneC_scaled}\n  expected: {expected_geneC}"
print("GeneC scaled:", np.round(geneC_scaled, 6))

print()
print("=" * 60)
print("STEP 3 — identify_significant_genes")
print("=" * 60)
young = df[early_cols].copy()
old   = df[late_cols].copy()
df    = identify_significant_genes(df, young, old)
print(df[["pvals", "logdiff", "direction"]].to_string())

print()
print("=" * 60)
print("STEP 4 — remove_insignificant_rows + extract_regulated_genes")
print("=" * 60)
df_sig    = remove_insignificant_rows(df.copy())
regulated = extract_regulated_genes(df_sig)
print(regulated.to_string())

print()
print("=" * 60)
print("PARITY CHECK vs R expected values")
print("=" * 60)
expected = {
    "GeneA": ("UP",    math.log2(15)),
    "GeneB": ("DOWN", -math.log2(15)),
    "GeneE": ("UP",    math.log2(7)),
    "GeneF": ("DOWN", -math.log2(7)),
}

all_passed = True
for gene, (exp_dir, exp_val) in expected.items():
    row = regulated[regulated["gene"] == gene] if "gene" in regulated.columns \
          else regulated.loc[[gene]] if gene in regulated.index else None

    if row is None or len(row) == 0:
        print(f"  FAIL  {gene}: not found in output")
        all_passed = False
        continue

    got_dir = row["direction"].values[0]
    got_val = row["valence"].values[0]

    dir_ok = got_dir == exp_dir
    val_ok = abs(got_val - exp_val) < 1e-10

    status = "PASS" if dir_ok and val_ok else "FAIL"
    if not (dir_ok and val_ok):
        all_passed = False

    print(f"  {status}  {gene}: direction={got_dir} (expected {exp_dir}), "
          f"valence={got_val:.6f} (expected {exp_val:.6f})")

print()
if all_passed:
    print("All parity checks passed. R and Python outputs are identical.")
else:
    print("One or more parity checks FAILED. See output above.")
    sys.exit(1)
