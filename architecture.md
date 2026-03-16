# Architecture — interflex

> Updated by scribe for run `merge-bs-dml-20260316-032034` on 2026-03-16.
> Previous runs: `interflex-dml-refactor-20260315-212459` (2026-03-15 — utils.R refactoring), `merge-bs-dml-20260316-032034` (2026-03-16 — merge bs features into dml).

## Overview

**interflex** is an R package (v1.3.5) for diagnosing and visualizing multiplicative interaction models. It estimates non-linear marginal effects of a treatment (D) on an outcome (Y) across values of a moderator (X), supporting both discrete and continuous treatments. The package provides eight estimation strategies (linear, binning, kernel, GAM, raw, GRF, DML, lasso), unified behind a single `interflex()` entry point. Key external dependencies include ggplot2 (plotting), mgcv (GAM), grf (causal forests), glmnet (lasso/ridge), reticulate (Python DML integration), and Rcpp/RcppArmadillo (C++ linear algebra stubs).

---

## Module Structure

```mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
    subgraph API["API Layer"]
        IFX["interflex.R — router"]
        PLT["plot.R — S3 plot method"]
        PRD["predict.R — S3 predict"]
        TST["inter_test.R — t-tests"]
        TST2["ttest.R — t-tests (bs)"]
    end

    subgraph Estimators["Estimator Layer"]
        LIN["linear.R"]
        BIN["binning.R"]
        KER["kernel.R"]
        GAM["gam.R"]
        RAW["raw.R"]
        GRF["grf.R"]
        DML["DML.R"]
        LAS["lasso.R"]
        LSD["lasso_discrete.R"]
    end

    subgraph DMLSub["DML Sub-Estimators"]
        CME_I["estimate_cme_irm.R"]
        CME_P["estimate_cme_plr.R"]
        GTE_I["estimate_gte_irm.R"]
        GTE_P["estimate_gte_plr.R"]
    end

    subgraph Output["Output Layer"]
        PPL["plot_pool.R — pooled plots"]
    end

    subgraph Utils["Utilities"]
        UTL["utils.R — shared helpers"]
        UNI["uniform.R — uniform CI"]
        VCL["vcluster.R — cluster vcov"]
    end

    subgraph Python["Python Backend"]
        DMLPY["dml.py — DML engine"]
    end

    IFX --> LIN
    IFX --> BIN
    IFX --> KER
    IFX --> GAM
    IFX --> RAW
    IFX --> GRF
    IFX --> DML
    IFX --> LAS
    IFX --> LSD
    PLT --> PPL
    LAS --> CME_I
    LAS --> CME_P
    LAS --> GTE_P
    LSD --> CME_I
    LSD --> GTE_I
    LIN --> UTL
    BIN --> UTL
    KER --> UTL
    GRF --> UTL
    DML --> UTL
    LAS --> UTL
    LSD --> UTL
    RAW --> UTL
    BIN --> VCL
    LIN --> VCL
    KER --> VCL
    DML --> DMLPY

    style IFX fill:#1e90ff,stroke:#1565c0,color:#fff
    style TST2 fill:#1e90ff,stroke:#1565c0,color:#fff
    style LAS fill:#1e90ff,stroke:#1565c0,color:#fff
    style DMLPY fill:#1e90ff,stroke:#1565c0,color:#fff
```

### Module Reference

| Module / File | Layer | Purpose | Key Exports | Changed |
| --- | --- | --- | --- | --- |
| `R/interflex.R` | API | Main entry point; validates inputs, builds `treat.info`/`diff.info`, routes to estimator | `interflex()` | yes (defaults) |
| `R/plot.R` | API | S3 `plot.interflex()` method; renders marginal effect plots with density/histogram overlays | `plot.interflex()` | no |
| `R/predict.R` | API | S3 `predict.interflex()` method; computes predicted marginal effects at new X values | `predict.interflex()` | no |
| `R/inter_test.R` | API | Post-estimation t-test for difference in marginal effects (dml-style) | `inter.test()` | no |
| `R/ttest.R` | API | Post-estimation t-test for difference in marginal effects (bs-style copy) | `inter.test()` | **new** |
| `R/linear.R` | Estimator | Linear interaction model with delta/bootstrap/simulation variance | `interflex.linear()` | no |
| `R/binning.R` | Estimator | Binning estimator: splits X into bins, estimates within-bin effects | `interflex.binning()` | no |
| `R/kernel.R` | Estimator | Kernel estimator: local polynomial regression with bandwidth selection | `interflex.kernel()` | no |
| `R/gam.R` | Estimator | GAM estimator via `mgcv::gam()` with 3D visualization | `interflex.gam()` | no |
| `R/raw.R` | Estimator | Raw data scatter plots with LOESS smoothing | `interflex.raw()` | no |
| `R/grf.R` | Estimator | Generalized random forests via `grf::causal_forest()` | `interflex.grf()` | no |
| `R/DML.R` | Estimator | Double/debiased ML via Python (reticulate); supports cross-fitting | `interflex.dml()` | no |
| `R/lasso.R` | Estimator | Lasso/ridge DML for continuous moderators; calls CME/GTE sub-estimators | `interflex.lasso()` | yes (defaults) |
| `R/lasso_discrete.R` | Estimator | Lasso/ridge DML for discrete moderators (<5 unique X values) | `interflex.lasso_discrete()` | no |
| `R/estimate_cme_irm.R` | DML Sub | CME estimation via AIPW-Lasso (binary treatment, IRM) | `estimateCME_IRM()` | no |
| `R/estimate_cme_plr.R` | DML Sub | CME estimation via PO-Lasso (continuous treatment, PLRM) | `estimateCME_PLR()` | no |
| `R/estimate_gte_irm.R` | DML Sub | Group treatment effects via AIPW-Lasso (binary treatment, discrete X) | `estimateGTE_IRM()` | no |
| `R/estimate_gte_plr.R` | DML Sub | Group treatment effects via PO-Lasso (continuous treatment, discrete X) | `estimateGTE_PLR()` | no |
| `R/plot_pool.R` | Output | Pooled multi-treatment plot with overlaid CIs | `interflex.plot.pool()` | no |
| `R/utils.R` | Utils | Shared internal helpers: treat.info extraction, density, histograms | (internal: dot-prefixed) | no |
| `R/uniform.R` | Utils | Uniform confidence interval quantiles via bootstrap/delta method | `calculate_uniform_quantiles()`, `calculate_delta_uniformCI()` | no |
| `R/vcluster.R` | Utils | Cluster-robust variance-covariance matrix computation | `vcovCluster()` | no |
| `R/RcppExports.R` | Utils | Auto-generated Rcpp bindings (do not edit) | `rcpparma_hello_world()`, etc. | no |
| `inst/python/dml.py` | Python | DML estimation engine; B-spline expansion for linear/logistic models | `marginal_effect_for_treatment()` | yes (B-spline) |
| `DESCRIPTION` | Config | Package metadata; Imports, Depends, LinkingTo | N/A | no |
| `NAMESPACE` | Config | Export pattern, S3 methods, importFrom declarations | N/A | no |
| `.Rbuildignore` | Config | Patterns excluded from R CMD build | N/A | yes |
| `.gitignore` | Config | Git ignore patterns | N/A | **new** |
| `.lintr` | Config | Lintr configuration | N/A | **new** |
| `_pkgdown.yml` | Config | Pkgdown site configuration | N/A | **new** |
| `interflex.Rproj` | Config | RStudio project file | N/A | **new** |

---

## Function Call Graph

### Main Pipeline

```mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
    USR["User call"] --> IFX["interflex()"]
    IFX --> VAL["Input validation"]
    VAL --> TI["Build treat.info"]
    TI --> DI["Build diff.info"]
    DI --> ROUTE{{"estimator?"}}

    ROUTE -- linear --> LIN["interflex.linear()"]
    ROUTE -- binning --> BIN["interflex.binning()"]
    ROUTE -- kernel --> KER["interflex.kernel()"]
    ROUTE -- gam --> GAM["interflex.gam()"]
    ROUTE -- raw --> RAW["interflex.raw()"]
    ROUTE -- grf --> GRF["interflex.grf()"]
    ROUTE -- dml --> DML["interflex.dml()"]
    ROUTE -- lasso --> LSPLIT{{"X < 5 levels?"}}
    LSPLIT -- yes --> LSD["interflex.lasso_discrete()"]
    LSPLIT -- no --> LAS["interflex.lasso()"]

    LIN --> ETI[".extract_treat_info()"]
    BIN --> ETI
    KER --> ETI
    GRF --> ETI
    DML --> ETI
    LAS --> ETI
    LSD --> ETI
    RAW --> ETI

    LIN --> CD[".compute_density()"]
    BIN --> CD
    KER --> CD
    GRF --> CD
    DML --> CD
    LAS --> CD
    LSD --> CD

    LIN --> CH[".compute_histograms()"]
    BIN --> CH
    KER --> CH
    GRF --> CH
    DML --> CH
    LAS --> CH
    LSD --> CH

    LAS --> CME_P["estimateCME_PLR()"]
    LAS --> CME_I["estimateCME_IRM()"]
    LAS --> GTE_P["estimateGTE_PLR()"]
    LSD --> CME_I
    LSD --> GTE_I["estimateGTE_IRM()"]

    DML --> DMLPY["dml.py B-spline"]

    style IFX fill:#1e90ff,stroke:#1565c0,color:#fff
    style LAS fill:#1e90ff,stroke:#1565c0,color:#fff
    style DMLPY fill:#1e90ff,stroke:#1565c0,color:#fff
```

### Output Pipeline

```mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
    OUT["interflex output object"] --> PLOT["plot.interflex()"]
    OUT --> PRED["predict.interflex()"]
    OUT --> TEST["inter.test() (inter_test.R)"]
    OUT --> TEST2["inter.test() (ttest.R)"]
    PLOT --> POOL{{"pool = TRUE?"}}
    POOL -- yes --> PPL["interflex.plot.pool()"]
    POOL -- no --> GPLT["ggplot2 rendering"]
    PPL --> GPLT
    PRED --> GPLT

    style TEST2 fill:#1e90ff,stroke:#1565c0,color:#fff
```

### Function Reference

| Function | Defined In | Called By | Calls | Changed | Purpose |
| --- | --- | --- | --- | --- | --- |
| `interflex()` | `R/interflex.R` | user (exported) | all estimators | yes (defaults) | Validate inputs, build metadata, route to estimator |
| `plot.interflex()` | `R/plot.R` | user (S3 method) | `interflex.plot.pool()` | no | Render marginal effect plots |
| `predict.interflex()` | `R/predict.R` | user (S3 method) | ggplot2 | no | Compute and plot predicted marginal effects |
| `inter.test()` | `R/inter_test.R` | user (exported) | mgcv::gam | no | Test differences in marginal effects (dml style) |
| `inter.test()` | `R/ttest.R` | user (exported) | mgcv::gam | **new** | Test differences in marginal effects (bs copy) |
| `interflex.linear()` | `R/linear.R` | `interflex()` | `.extract_treat_info`, `.compute_density`, `.compute_histograms`, `vcovCluster` | no | Linear interaction model estimation |
| `interflex.binning()` | `R/binning.R` | `interflex()` | `.extract_treat_info`, `.compute_density`, `.compute_histograms`, `vcovCluster` | no | Binning estimator |
| `interflex.kernel()` | `R/kernel.R` | `interflex()` | `.extract_treat_info`, `.compute_density`, `.compute_histograms`, `vcovCluster` | no | Kernel estimator |
| `interflex.gam()` | `R/gam.R` | `interflex()` | `mgcv::gam` | no | GAM-based 3D surface estimation |
| `interflex.raw()` | `R/raw.R` | `interflex()` | `.extract_treat_info` | no | Raw scatter plots with LOESS |
| `interflex.grf()` | `R/grf.R` | `interflex()` | `.extract_treat_info`, `.compute_density`, `.compute_histograms`, `grf::causal_forest` | no | Causal forest estimation |
| `interflex.dml()` | `R/DML.R` | `interflex()` | `.extract_treat_info`, `.compute_density`, `.compute_histograms`, `reticulate::source_python` | no | Python-based DML estimation |
| `marginal_effect_for_treatment()` | `inst/python/dml.py` | `interflex.dml()` (via reticulate) | patsy, doubleml, sklearn | yes (B-spline) | Python DML with optional B-spline expansion |
| `interflex.lasso()` | `R/lasso.R` | `interflex()` | `.extract_treat_info`, `.compute_density`, `.compute_histograms`, `estimateCME_PLR`, `estimateCME_IRM`, `estimateGTE_PLR` | yes (defaults) | Lasso DML for continuous X |
| `interflex.lasso_discrete()` | `R/lasso_discrete.R` | `interflex()` | `.extract_treat_info`, `.compute_density`, `.compute_histograms`, `estimateCME_IRM`, `estimateGTE_IRM` | no | Lasso DML for discrete X |
| `estimateCME_IRM()` | `R/estimate_cme_irm.R` | `interflex.lasso`, `interflex.lasso_discrete` | `glmnet::cv.glmnet` | no | CME via AIPW-Lasso (binary D) |
| `estimateCME_PLR()` | `R/estimate_cme_plr.R` | `interflex.lasso` | `glmnet::cv.glmnet` | no | CME via PO-Lasso (continuous D) |
| `estimateGTE_IRM()` | `R/estimate_gte_irm.R` | `interflex.lasso_discrete` | `glmnet::cv.glmnet` | no | GTE via AIPW-Lasso (binary D, discrete X) |
| `estimateGTE_PLR()` | `R/estimate_gte_plr.R` | `interflex.lasso` | `glmnet::cv.glmnet` | no | GTE via PO-Lasso (continuous D, discrete X) |
| `.extract_treat_info()` | `R/utils.R` | 8 estimators + raw | — | no | Unpack treat.info list into local variables |
| `.compute_density()` | `R/utils.R` | 7 estimators | `stats::density` | no | Compute kernel density estimates |
| `.compute_histograms()` | `R/utils.R` | 7 estimators | `graphics::hist` | no | Compute histogram bin counts |
| `vcovCluster()` | `R/vcluster.R` | binning, kernel, linear | `sandwich::estfun`, `sandwich::bread` | no | Cluster-robust variance-covariance |
| `calculate_uniform_quantiles()` | `R/uniform.R` | binning, kernel, linear, lasso | — | no | Bootstrap uniform CI bands |
| `calculate_delta_uniformCI()` | `R/uniform.R` | linear, binning | `MASS::mvrnorm` | no | Delta-method uniform CI bands |
---

## Data Flow

```mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
    INPUT["User: interflex(estimator, data, Y, D, X, ...)"]
    INPUT --> VALIDATE["Validate inputs & coerce types"]
    VALIDATE --> TREAT["Build treat.info metadata"]
    TREAT --> DIFF["Build diff.info for contrasts"]
    DIFF --> ROUTE{{"Select estimator"}}

    ROUTE --> EST["Run estimator function"]
    EST --> ETI[".extract_treat_info()"]
    ETI --> MODEL["Fit statistical model"]
    MODEL --> ME["Compute marginal effects"]
    ME --> XDISTR{{"Xdistr setting?"}}
    XDISTR -- density --> DENS[".compute_density()"]
    XDISTR -- histogram --> HIST[".compute_histograms()"]
    XDISTR -- none --> SKIP["Skip distribution"]
    DENS --> BUILDOUT["Build output list"]
    HIST --> BUILDOUT
    SKIP --> BUILDOUT
    BUILDOUT --> FIGURE{{"figure = TRUE?"}}
    FIGURE -- yes --> PLOTINT["Generate ggplot"]
    FIGURE -- no --> RETURN["Return interflex object"]
    PLOTINT --> RETURN

    RETURN --> USER["User receives interflex object"]
    USER --> SPLOT["plot(out)"]
    USER --> SPRED["predict(out)"]
    USER --> STEST["inter.test(out)"]

    style ROUTE fill:#1e90ff,stroke:#1565c0,color:#fff
```

---

## Key Data Structures

### `treat.info` (built by `interflex()`, consumed by all estimators)

A named list containing treatment metadata. The new `.extract_treat_info()` utility unpacks this uniformly.

| Field | When Present | Content |
| --- | --- | --- |
| `treat.type` | always | `"discrete"` or `"continuous"` |
| `other.treat` | discrete | Named character vector of non-base treatment levels |
| `all.treat` | discrete | Named character vector of all treatment levels |
| `base` | discrete | Base treatment level (reference group) |
| `D.sample` | continuous | Named numeric vector of sampled treatment values |
| `ncols` | when set | Number of plot columns |

### `interflex` output object

A list of class `"interflex"` returned by each estimator, containing:

| Field | Content |
| --- | --- |
| `est.lin` / `est.bin` / `est.kernel` / etc. | Marginal effect estimates data frame |
| `diff.estimate` | Treatment contrast estimates |
| `figure` | ggplot object(s) |
| `hist.out`, `treat.hist`, `de`, `treat_den` | Distribution data for X-axis overlays |
| `treat.info`, `diff.info` | Metadata passed through |
| `model.coef`, `model.vcov` | Model coefficients and variance-covariance (linear, binning) |

---

## Estimator Architecture

| Estimator | Function | Treatment Type | Moderator Type | Method | Variance |
| --- | --- | --- | --- | --- | --- |
| `"linear"` | `interflex.linear()` | discrete or continuous | continuous | Parametric OLS/GLM with D*X interaction | delta, bootstrap, simulation |
| `"binning"` | `interflex.binning()` | discrete or continuous | continuous (binned) | Split X into bins, within-bin linear models | delta, bootstrap, simulation |
| `"kernel"` | `interflex.kernel()` | discrete or continuous | continuous | Local polynomial regression, CV bandwidth | bootstrap |
| `"gam"` | `interflex.gam()` | continuous only | continuous | `mgcv::gam()` smooth surface | GAM built-in |
| `"raw"` | `interflex.raw()` | discrete or continuous | continuous | Scatter + LOESS (no formal estimation) | none |
| `"grf"` | `interflex.grf()` | binary | continuous | `grf::causal_forest()` | forest-based |
| `"dml"` | `interflex.dml()` | discrete or continuous | continuous | Python DML via reticulate cross-fitting | cross-fit |
| `"lasso"` | `interflex.lasso()` | binary or continuous | continuous (>=5 levels) | PO-Lasso (PLRM) or AIPW-Lasso (IRM) | bootstrap |
| `"lasso"` | `interflex.lasso_discrete()` | binary or continuous | discrete (<5 levels) | GTE estimation via IRM or PLR | bootstrap |

### DML Two-Pipeline Architecture

The `"lasso"` estimator routes to sub-estimators based on treatment type and moderator cardinality:

- **Continuous D, continuous X**: `interflex.lasso()` calls `estimateCME_PLR()` (Partially Linear Regression Model)
- **Binary D, continuous X**: `interflex.lasso()` calls `estimateCME_IRM()` (Interactive Regression Model via AIPW)
- **Continuous D, discrete X** (<5 levels): `interflex.lasso_discrete()` calls `estimateGTE_PLR()`
- **Binary D, discrete X** (<5 levels): `interflex.lasso_discrete()` calls `estimateGTE_IRM()`

All four sub-estimators use `glmnet::cv.glmnet()` for regularized nuisance function estimation with basis expansion (polynomial, B-spline, or none).

### Python Integration (DML Estimator)

`interflex.dml()` uses `reticulate` to:
1. Validate the Python script path exists (`nzchar()` check)
2. Source the Python script via `reticulate::source_python()` wrapped in `tryCatch()`
3. Pass data and parameters to Python via `reticulate::dict()`
4. **B-spline expansion** (new in this run): When both `model_y` and `model_t` are linear, logistic, or regularization types, the Python script creates a B-spline basis expansion of X via `patsy.dmatrix("bs(x, degree=3, df=5)")` and appends the spline columns to covariates before fitting `DoubleMLData`. For non-linear models (random forest, etc.), the original data path is used unchanged.
5. Python performs cross-fitted DML estimation (scikit-learn, doubleml)
6. Results are returned to R for plotting

---

## Utility Functions (Added in Previous Refactoring Run)

Three internal utility functions in `R/utils.R` consolidate previously duplicated code. They are dot-prefixed (`.extract_treat_info`, `.compute_density`, `.compute_histograms`) to prevent auto-export via the `exportPattern("^[[:alpha:]]+")` rule in NAMESPACE.

| Function | Purpose | Used By | Lines Saved |
| --- | --- | --- | --- |
| `.extract_treat_info(treat.info)` | Unpacks the `treat.info` list into local variables (`treat.type`, `other.treat`, `all.treat`, `base`, `D.sample`, etc.) | 9 files: DML, binning, kernel, linear, grf, lasso, lasso_discrete, raw, inter_test | ~135 lines |
| `.compute_density(data, X, D, weights, treat.type, all.treat, all.treat.origin)` | Computes overall and per-treatment kernel density estimates for the X-axis distribution overlay | 7 files: DML, binning, kernel, linear, grf, lasso, lasso_discrete | ~175 lines |
| `.compute_histograms(data, X, D, weights, treat.type, all.treat, all.treat.origin)` | Computes overall and per-treatment histogram bin counts for the X-axis distribution overlay | 7 files: DML, binning, kernel, linear, grf, lasso, lasso_discrete | ~175 lines |

---

## Architectural Patterns

- **Router pattern**: `interflex()` is a monolithic router (~1400 lines) that validates all inputs, builds shared metadata (`treat.info`, `diff.info`), and dispatches to one of 9 estimator functions. Each estimator is a standalone function in its own file.

- **Shared metadata**: `treat.info` and `diff.info` are computed once by the router and passed to every estimator. The `.extract_treat_info()` utility provides uniform unpacking.

- **Inline plotting**: Each estimator builds its own ggplot figure internally rather than delegating to a separate plot function. The S3 `plot.interflex()` method re-renders from stored data.

- **Dot-prefix convention**: Internal helpers use `.` prefix to avoid export via the blanket `exportPattern("^[[:alpha:]]+")` rule, without requiring explicit `@keywords internal` tags.

- **Lasso moderator cardinality split**: The `"lasso"` estimator auto-selects `interflex.lasso_discrete()` when X has fewer than 5 unique values, switching from CME to GTE estimation.

---

## Notes

- **Previous run (refactor)**: -506 lines across 21 files. ~500 lines of duplicated code consolidated into `R/utils.R`. ~700 redundant boolean comparisons cleaned. 0 ERRORs, 5 WARNINGs (pre-existing), 2 NOTEs (pre-existing).
- **This run (bs merge)**: 5 files created, 4 files modified. Selectively integrated bs branch features (ttest.R, B-spline expansion, parameter defaults) into dml's refactored architecture without reverting the utils.R helper consolidation.
- **Duplicate ttest files**: Both `R/ttest.R` (bs copy, `== FALSE` style) and `R/inter_test.R` (dml copy, `!` style) exist and both export `inter.test()`. R will load one or the other depending on file ordering. A future cleanup should remove the duplicate.
- **Parameter default changes**: `spline.df` 5 to 4, `spline.degree` 3 to 2, `reduce.dimension` "bspline" to "kernel". These are intentional feature changes from the bs branch. Existing users relying on previous defaults will see different behavior.
- **show.uniform.CI bug prevented**: The bs branch passed `show.uniform.CI` to `interflex.lasso_discrete()` which does not accept it. This was caught during audit and removed before shipping.
- **No formal test suite**: The package does not have tests under `tests/`. Validation relies on R CMD check and manual examples.
