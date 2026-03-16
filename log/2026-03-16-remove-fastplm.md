# 2026-03-16 — Remove fastplm dead code

> Run: `remove-fastplm-20260316-041157` | Profile: r-package | Verdict: PASS

## What Changed

Removed all references to fastplm and iv_fastplm C++ functions that were replaced by fixest in a prior refactoring. These functions were defined in `src/` but never called by any R estimator code. Cleaned up corresponding Rcpp bindings, manual page aliases, and updated `architecture.md` to reflect the removal.

## Files Changed

| File | Action | Description |
| --- | --- | --- |
| src/fastplm.cpp | deleted | Unused C++ implementation of fastplm |
| src/fastplm.o | deleted | Compiled object file for fastplm |
| src/iv_fastplm.cpp | deleted | Unused C++ implementation of iv_fastplm |
| src/iv_fastplm.o | deleted | Compiled object file for iv_fastplm |
| R/RcppExports.R | modified | Removed R wrapper functions for fastplm/iv_fastplm |
| src/RcppExports.cpp | modified | Removed C++ binding entries for fastplm/iv_fastplm |
| man/interflex-internal.Rd | modified | Removed alias entries for fastplm/iv_fastplm |
| architecture.md | modified | Updated overview to remove fastplm references |

## Design Decisions

1. **Pure deletion, no replacement**: The fastplm/iv_fastplm functions had already been replaced by fixest-based estimators. No new code was needed — only cleanup of dead artifacts.
2. **Preserve NEWS.md historical reference**: NEWS.md contains a historical mention of fastplm removal. This was intentionally kept as it documents the package's changelog.

## Handoff Notes

- The `.o` files were committed build artifacts that should not have been in version control. They are now removed.
- All remaining Rcpp bindings (CppEstimate, CppKernelEstimate, CppCrossValidation) are intact and unchanged.
- R CMD check was not run (R not available in CI environment), but code inspection confirmed no orphaned symbols or broken registrations.

## Verification

- **R CMD check**: Not run (R unavailable); code inspection confirmed structural integrity
- **Review verdict**: PASS — 8 challenges raised, all cleared. Both pipelines converged.
