# 2026-03-16 — Python-to-R DML Migration

> Run: `py-to-r-dml-001` | Profile: r-package | Verdict: pending (auditor/skeptic review)

## What Changed

The DML (Double Machine Learning) estimator in interflex was completely migrated from a Python/reticulate backend to pure R. Previously, `interflex.dml()` called Python via `reticulate::source_python()` to execute `inst/python/dml.py`, which depended on `doubleml`, `scikit-learn`, `patsy`, `pandas`, `numpy`, and `scipy`. This required users to install Python, Miniconda, and multiple Python packages — a significant friction point for R users. The new implementation uses the R `DoubleML` package (R6 interface) with `mlr3` learners, eliminating all Python dependencies while preserving the same API and output format.

## Files Changed

| File | Action | Description |
| --- | --- | --- |
| `R/DML.R` | modified (rewritten) | Complete rewrite from ~244 lines to ~430 lines. Replaced all Python/reticulate calls with R DoubleML + mlr3. Added 10 new internal helper functions. |
| `inst/python/dml.py` | deleted | Removed the 298-line Python DML engine. No longer needed. |
| `DESCRIPTION` | modified | Removed `reticulate` from Imports. Added `DoubleML`, `mlr3`, `mlr3learners`, `data.table`, `paradox`, `ranger` to Imports. Added `mlr3tuning`, `lightgbm`, `nnet` to Suggests. |
| `NAMESPACE` | modified | Removed `importFrom("reticulate", ...)`. New packages accessed via `::` (no importFrom needed). |

## Design Decisions

1. **Manual BLP for CATE/GATE**: The R DoubleML package does not expose `.cate()` and `.gate()` methods like the Python version. Rather than creating a wrapper or using a different package, CATE and GATE are computed via manual Best Linear Projection (BLP) — projecting DML pseudo-outcomes onto B-spline basis functions (CATE) or group dummies (GATE), with HC sandwich variance and both pointwise and uniform confidence intervals. This follows the established statistical methodology and gives full control over the CI computation.

2. **Pseudo-outcome extraction with fallback**: The primary extraction path uses `dml_model$psi_b[, 1, 1]` (influence function values). A fallback reconstructs from `$all_psi` and `$coef`. A last-resort path uses the constant ATE. This three-tier approach handles different DoubleML versions and edge cases gracefully.

3. **Parameter mapping tables**: Rather than requiring users to learn mlr3/ranger parameter names, explicit mapping tables translate sklearn parameter names (`n_estimators`, `max_depth`, `hidden_layer_sizes`, etc.) to their R equivalents (`num.trees`, `max.depth`, `size`, etc.). Unknown names are passed through unchanged. This preserves backward compatibility for users who have existing code with sklearn-style parameters.

4. **`ranger` in Imports, `lightgbm`/`nnet` in Suggests**: Random forest is the default model (`model.y = "rf"`), so `ranger` must always be available. Boosting and neural network are non-default choices, so their backends are optional.

5. **`paradox` in Imports**: Required for building tuning parameter sets whenever CV=TRUE with param grids. Placed in Imports rather than Suggests to avoid fragile `requireNamespace()` checks in the tuning path.

6. **No intercept in B-spline expansion**: Python's `patsy.dmatrix("bs()")` includes an intercept by default, but DoubleML learners handle intercepts internally. B-spline columns are named `X_bs_1` through `X_bs_5` without an intercept column.

7. **`n.jobs` parameter retained**: Kept in the function signature for API compatibility but not used internally. R DoubleML handles parallelism via the mlr3 future backend, which users can configure separately.

## Handoff Notes

- **Testing**: The package has no formal test suite under `tests/`. Validation relies on `R CMD check` and manual examples. The DML vignette (`vignettes/DML.Rmd`) contains runnable examples that serve as integration tests.

- **DoubleML version sensitivity**: The pseudo-outcome extraction (`$psi_b`) works with DoubleML >= 0.5.0. The fallback paths handle older versions, but the primary path should be verified when updating DoubleML.

- **CATE B-spline parameters**: The BLP uses `splines::bs(degree=2, df=5)` for CATE evaluation. These are hardcoded in `.compute_cate_blp()`. If users need different smoothness, this would need to be parameterized.

- **Ridge regularization in `.safe_solve()`**: When the B-spline Gram matrix is singular (can happen with sparse data), `.safe_solve()` adds `1e-8 * diag(ncol)` as ridge regularization. This is a practical fix but affects inference slightly.

- **Vignette update needed**: The DML vignette (`vignettes/DML.Rmd`) previously contained Python/reticulate setup instructions and references to sklearn documentation. These have been updated to reflect the pure R implementation.

- **`lightgbm` installation**: On some platforms, `lightgbm` requires system-level dependencies (cmake, etc.). Since it is in Suggests, users only encounter this if they choose `model.y = "hgb"` or `model.t = "hgb"`.

## Verification

- **R CMD check**: Pending auditor run
- **Review verdict**: Pending skeptic review
