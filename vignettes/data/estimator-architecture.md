# interflex Estimator Architecture

## Conceptual Framework: Nuisance × Aggregation

Every CME/GATE estimator in **interflex** can be understood as combining two independent choices:

1. **Nuisance estimation**: How do we model the relationship between outcome, treatment, and covariates?
2. **Aggregation**: How do we summarize heterogeneous effects across values of the moderator X?

### Nuisance models

| Nuisance model | Description | interflex estimators |
|----------------|-------------|---------------------|
| Linear (OLS) | Parametric: Y ~ D*X + Z | `linear`, `binning` |
| Kernel | Local linear regression | `kernel` |
| Lasso / Ridge | Penalized regression with basis expansion | `lasso` |
| ML (DML) | mlr3 learners (random forest, neural net, xgboost, etc.) | `dml` |
| Causal Forest | Honest causal forest via `grf` | `grf` |

### Aggregation methods

| Aggregation | Description | When to use |
|-------------|-------------|-------------|
| **Smooth curve** (kernel/spline) | Project scores onto a smooth function of X | X is continuous — produces CME(x) curve |
| **Bin** (group average) | Average scores within each level of X | X is discrete — produces GATE(x) point estimates |

### The mapping

| | **Smooth** (CME curve) | **Bin** (GATE) |
|---|---|---|
| **Linear** | `estimator = 'linear'` | `estimator = 'linear', gate = TRUE` |
| **Lasso** | `estimator = 'lasso'` | `estimator = 'lasso', gate = TRUE` |
| **DML** | `estimator = 'dml'` | `estimator = 'dml', gate = TRUE` |
| **GRF** | `estimator = 'grf'` | `estimator = 'grf', gate = TRUE` |
| **Kernel** | `estimator = 'kernel'` | not supported (kernel is itself an aggregation) |
| **Binning** | `estimator = 'binning'` | not applicable (binning IS bin aggregation) |

Note: `kernel` and `binning` are "fused" estimators where nuisance and aggregation are inseparable. The semiparametric estimators (`lasso`, `dml`, `grf`) and `linear` (with `gate = TRUE`) cleanly separate the two steps.

## Signal types (binary treatment only)

When D is binary and `gate = TRUE`, the `lasso` estimator supports different signal constructions:

| Signal | Formula | Nuisance requirements |
|--------|---------|----------------------|
| `outcome` | μ̂₁(X,Z) − μ̂₀(X,Z) | outcome model only |
| `ipw` | DY/π̂ − (1−D)Y/(1−π̂) | propensity model only |
| `aipw` | (μ̂₁ − μ̂₀) + D(Y−μ̂₁)/π̂ − (1−D)(Y−μ̂₀)/(1−π̂) | both outcome + propensity (doubly robust) |

- `linear` with `gate = TRUE`: outcome signal only (parametric, no IPW/AIPW needed)
- `lasso` with `gate = TRUE`: all three signals available
- `dml` with `gate = TRUE`: uses DML orthogonal score (analogous to AIPW)
- `grf` with `gate = TRUE`: uses forest-based CATE estimates

For continuous D, the PLR (partially linear regression) framework is used instead — no signal choice needed.

## Current implementation status

| Estimator | `gate = TRUE` | `signal` | Status |
|-----------|--------------|----------|--------|
| `linear` | planned | outcome only | **not yet implemented** |
| `lasso` | works | outcome/ipw/aipw | implemented (auto-routed when X < 5 levels) |
| `dml` | works | DML score | implemented |
| `grf` | planned | forest CATE | **not yet implemented** |
| `kernel` | skip | — | n/a |
| `binning` | skip | — | n/a |

## Existing GATE infrastructure

These functions already support `outcome_model_type = "linear"` and are reusable:

- `estimateGTE()` / `bootstrapGTE()` in `estimate_gte_irm.R` — binary D, supports outcome/ipw/aipw signals
- `estimateGATE_PLR()` / `bootstrapGATE_PLR()` in `estimate_gte_plr.R` — continuous D, PLR framework
- `.compute_gate_blp()` in `DML.R` — DML-specific GATE via BLP on group dummies

## Unified output field

All estimators should populate `g.est` (not `g.est.dml`) when `gate = TRUE`. The `plot()` function checks `g.est` to enable `by.group = TRUE` plotting.

## Changes needed for `gate = TRUE` on `linear` and `grf`

1. Accept `gate` parameter in `interflex.linear()` and `interflex.grf()`
2. When `gate = TRUE`: call `bootstrapGTE` (binary D) or `bootstrapGATE_PLR` (continuous D) with appropriate `outcome_model_type`
3. Rename `g.est.dml` → `g.est` across codebase
4. Update `plot.R` `by.group` check from `g.est.dml` to `g.est`
5. Make `lasso` gate routing explicit via `gate = TRUE` instead of auto by X cardinality
