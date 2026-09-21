# Architecture -- interflex

> Updated by scriber for run `KBW-20260917` on 2026-09-20.
> Previous runs: `interflex-dml-refactor-20260315-212459` (2026-03-15 -- utils.R refactoring), `merge-bs-dml-20260316-032034` (2026-03-16 -- merge bs features into dml), `py-to-r-dml-001` (2026-03-16 -- Python-to-R DML migration), `REQ-20260403-080941` (2026-04-03 -- parallel RNG migration from doParallel to doFuture), `GATE-001` (2026-04-04 -- generalize GATE support across estimators), `BOOK-003` (2026-04-07 -- plot xlim/ylim coord_cartesian migration, defensive narrow-table helpers, ch6 fit/plot chunk split), `PAD-001` (2026-04-07 -- visible xlim padding for continuous-treatment plots via two-pass grid restriction + plot-time row filter), `PAD-002-discrete` (2026-04-07 -- extended PAD-001 PASS 2 row filter to binary and multi-arm discrete treatments; tightened panel window so scale limits match `coord_cartesian` exactly), `KBW-20260917` (2026-09-20 -- kernel adaptive-bandwidth normalizer fix: geometric mean at the data points instead of the density grid; drop/stop policy for degenerate local fits; 1.4.1 release prep).

## Padded-xlim invariant (PAD-001 / PAD-002)

For continuous-treatment plots, whenever the user supplies an explicit
`xlim = c(lo, hi)`, the visible curve (mean line, pointwise CI ribbon,
uniform CI bands) and the moderator distribution overlay (density ribbon
or histogram bars) end exactly at `lo`/`hi`, and `coord_cartesian(xlim =
.pad_xlim(xlim, mult = 0.04))` produces ~4% whitespace between those
endpoints and the panel edges. Two independent guarantees enforce this:

- **PASS 1 -- fit-time grid restriction.** When `xlim` is passed to
  `interflex(..., xlim = c(lo, hi))` and `treat.type == "continuous"`, the
  prediction grid is built as `seq(lo, hi, length.out = neval)` at every
  continuous-estimator construction site: `R/linear.R`, `R/kernel.R`,
  `R/binning.R`, `R/estimate_cme_plr.R` (lasso PLR path),
  `R/DML.R::.compute_cate_blp` (DML CATE path). A function-local
  `user_xlim_explicit` flag (default `FALSE`) is plumbed as a named
  argument from `interflex.R` through each dispatch -- the gate NEVER
  reads from `getOption(...)` inside helper files. Degenerate `xlim`
  (non-finite, `lo >= hi`) falls back to `NULL` with a warning at the
  `interflex.R` validation block. AIPW/IRM is binary-discrete only and
  out of scope; gate plots (`.compute_gate_blp`) and `R/gam.R` (which
  delegates to `mgcv::vis.gam` with its own internal grid) are untouched.

- **PASS 2 -- plot-time row filter.** When the user passes `xlim` to
  `plot.interflex(out, xlim = c(lo, hi))` on an `out` object that was
  fit WITHOUT `xlim`, the stored prediction tables still span the full
  data range. `plot.interflex` therefore builds a local filtered shadow
  `.out_filt` immediately after the `.user_xlim_in` snapshot at
  `R/plot.R:198`. The gate `.pad_xlim_gate` fires only when
  `treat.type == "continuous"` AND `.user_xlim_in` is finite, ordered
  length-2. When it fires, `.filter_xlim_rows` drops rows where
  `X < lo - eps | X > hi + eps` (with `eps = 1e-9 * max(1, |hi - lo|)`
  and an empty-keep guard) from every continuous `est.*` field
  (`est.lin`, `est.bin`, `est.kernel`, `est.dml`, `est.lasso`,
  `est.grf`). Downstream continuous branches rebind `est.* <- .out_filt$est.*`
  -- the filter therefore covers the mean line, pointwise CI ribbon,
  AND uniform CI bands in a single pass (all three read from the same
  `tempest` data frame). `out` itself is NEVER mutated; `.out_filt`
  is a fresh local.

- **PASS 2b -- density/histogram filter.** The Xdistr overlay bars and
  density ribbons extend into the padded whitespace if not filtered.
  Two additional helpers `.filter_xlim_density(dens, lo, hi)` and
  `.filter_xlim_histlike(h, lo, hi)` (also in `R/plot.R`) extend the
  `.out_filt` shadow to cover `out$de` (density: `$x`, `$y`) and
  `out$hist.out` (histogram: `$mids`, `$counts`, `$density`).
  The continuous branches in the `Xdistr == "density"` and
  `Xdistr %in% c("histogram","hist")` blocks of `plot.R` rebind
  `de <- .out_filt$de` and `hist.out <- .out_filt$hist.out` as the
  first statement inside their `if (treat.type == "continuous")`
  sub-block. Histogram bin-width `dist` has a fallback to the
  pre-filter spacing when the filtered hist has <2 mids. Discrete
  Xdistr fields (`out$de.tr`, `out$count.tr`) are intentionally NOT
  filtered.

- **Idempotence.** PASS 1 and PASS 2/2b compose cleanly. When `xlim`
  is supplied to BOTH `interflex()` and `plot()`, PASS 1 has already
  built the grid on `[lo, hi]`, so PASS 2's row filter is a no-op
  on a table whose rows are already inside the window. When `xlim`
  is supplied at plot time only, PASS 1 is silent and PASS 2 does
  the work. When `xlim` is absent at both layers, neither gate
  fires and the behavior is bit-identical to pre-PAD-001.

- **Forbidden paths (preserved).** `.pad_xlim(mult = 0.04)` is
  unchanged. `R/plot_pool.R`'s `coord_cartesian(.pad_xlim(...))`
  call sites are unchanged. No `scale_x_continuous(limits = ...)` or
  `oob = censor` is used anywhere -- the 4% whitespace is still
  delivered by `coord_cartesian`, not by scale clipping. `R/gam.R`
  (delegated to `mgcv::vis.gam`) is untouched.

- **PAD-002 -- discrete treatment coverage.** PAD-001 originally
  hard-gated PASS 2 on `treat.type == "continuous"`, so binary and
  multi-arm discrete fits bypassed the row filter entirely. PAD-002
  removes the `treat.type == "continuous"` clause from
  `.pad_xlim_gate` in `R/plot.R` so the gate fires for any fit type
  when `.user_xlim_in` is a finite ordered length-2 numeric. Two new
  helpers extend the `.out_filt` shadow to cover discrete-only
  fields: `.filter_xlim_count_tr(count_tr, mask, orig_len)` masks
  each per-treatment count vector in `out$count.tr` using the
  **unfiltered** `out$hist.out$mids` as the indexing basis (the
  mids snapshot MUST be captured before `.filter_xlim_histlike`
  mutates `.out_filt$hist.out$mids`), and `.filter_xlim_de_tr` is a
  thin wrapper applying `.filter_xlim_density` across the named list
  of per-arm density objects in `out$de.tr`. `g.est` and `g.est.dml`
  (GATE discrete outputs) are filtered via the existing
  `.filter_xlim_field` helper. Every new branch is guarded by a
  `%in% names(.out_filt)` presence check and degrades to a no-op
  when the field is absent or has unexpected shape.

- **PAD-002 -- rendering locals rebound under the gate.** The local
  bindings `de`, `de.tr`, `hist.out`, `count.tr` (at the top of the
  per-treatment rendering block in `plot.interflex`) are now
  assigned conditionally: when `.pad_xlim_gate` is TRUE, they read
  from `.out_filt$*`; when FALSE, they fall back to `out$*`
  byte-identically. The same gate-conditional rebinding is applied
  to every discrete `est.*` sub-branch -- binning (`est.lin`,
  `est.bin`), dml (`est.dml`), linear (`est.lin`), grf (`est.grf`),
  lasso (`est.lasso`) -- mirroring the pattern already present in
  the continuous sub-branches. Kernel's discrete sub-branch already
  reads `.out_filt$est.kernel` unconditionally and is untouched.

- **PAD-002 -- rect boundary clamp.** The discrete histogram overlay
  builds rect data frames from `hist.out$mids +/- dist/2`, so rect
  rows for boundary mids can extend up to half a bin width past the
  user `xlim`. Under the gate, `xmin <- pmax(xmin, .user_xlim_in[1])`
  and `xmax <- pmin(xmax, .user_xlim_in[2])` are applied at three
  sites in `plot.R`: the discrete control-arm `hist.col` overlay,
  the discrete per-treated-arm `hist.treat` overlay (inside the
  `for (char in other.treat)` loop), and the continuous `histX`
  overlay (for symmetry -- the continuous branch was already
  filtered but boundary mids can still overhang). The clamp
  shrinks the outermost bars rather than dropping them.

- **PAD-002 -- panel window matches coord limits exactly.** Even
  after the row filter and rect clamp cover all geom data, the
  rendered x-axis panel range still overshot `coord_cartesian`'s
  `.pad_xlim(xlim, 0.04)` window because ggplot2's default scale
  expansion adds another ~5% per side on top of `coord_cartesian`.
  Under the gate, `plot.interflex` now adds
  `ggplot2::scale_x_continuous(expand = ggplot2::expansion(0, 0))`
  immediately before the existing `coord_cartesian(xlim =
  final_xlim, ylim = final_ylim)` call at BOTH per-panel sites --
  the discrete per-arm panel loop (`for (char in other.treat)`)
  and the continuous per-label panel loop
  (`for (label in label.name)`). The rendered panel `x.range` is
  therefore exactly `.pad_xlim(xlim)` (verified to within `1e-6`
  in the PAD-002 test suite: `[0.22, 1.03]` for binary `xlim =
  c(0.25, 1)`, `[-2.16, 2.16]` for multi-arm `xlim = c(-2, 2)`).
  When the gate is FALSE the new `scale_x_continuous` line is
  never added, so the no-xlim path is byte-identical to the
  pre-PAD-002 behavior.

- **PAD-002 -- composition with PAD-001.** PAD-002 is additive: the
  PAD-001 continuous flow (row filter for `est.*`, density, and
  histogram) is unchanged in both code path and numerical result.
  The discrete extension reuses the same helpers and the same
  gate, so continuous and discrete fits now share one unified
  row-filter invariant. The gate remains a complete no-op when
  `xlim` is unset at plot time.

- **Encoding invariant.** All R source files under `R/` MUST remain
  ASCII-safe unless `DESCRIPTION` declares `Encoding: UTF-8`.
  `pkgload::source_one()` (used by `devtools::load_all` and therefore
  by `devtools::test`) calls `readLines()` without an explicit
  encoding; non-ASCII bytes in the native locale cause
  "unexpected end of input" parse errors even though plain
  `parse()` and `R CMD INSTALL` succeed. See PAD-001 process record.


## Kernel adaptive bandwidth (KBW-20260917)

The kernel estimator (`R/kernel.R`, `interflex.kernel()`) fits a local weighted
regression at every requested evaluation point x0 on the moderator. The local
Gaussian window is `h(x0) = bw * sqrt(g / f_eff(x0))`, where `f_eff` is a pilot
density of the moderator (`stats::density(X, weights = w)`, looked up at the
nearest grid point) and `g` is a normalizing constant (Abramson's square-root
rule; Silverman 1986, Sec. 5.3.1). Before this run, `g` was the geometric mean
of the density taken over the whole 512-point `density()` grid, including the
empty gaps and tails beyond the data range. On a moderator with gaps or a long
tail (e.g. the Malesky data on its raw scale), that made `g` far too small, so
windows near the data collapsed onto a handful of points, the local fit lost
coefficients (aliased), and the old code either zero-filled them or crashed
with "subscript out of bounds" / "$ operator is invalid for atomic vectors".

**The fix, in three parts:**

1. **Normalizer computed at the observations, not the grid.** `g` is now the
   sampling-weighted geometric mean of the pilot density AT THE OBSERVATIONS
   with positive weight and positive, finite density:
   `g = exp(sum(w_i * log f(X_i)) / sum(w_i))`. It is computed once per
   sample, immediately after the `density()` call that built the pilot
   density (the main fit, each CV training fold, each bootstrap replicate),
   and cached on the density object itself as `dens$adapt.g` /
   `dens$adapt.floor` (the smallest positive grid value, used as a floor when
   `f(x0)` is 0 or non-finite -- this can still happen at an EVALUATION point
   that falls in an empty gap, even though it essentially never happens at an
   observation).
2. **No degenerate fit is silently zero-filled or crashed on.** Every local
   fit now returns a `status` field (`"ok"` or one of five reason codes
   below) instead of writing 0 into unidentified coefficients and letting
   `vcovHC` silently drop them from the covariance matrix -- which is what
   produced the out-of-bounds index. An unusable point is dropped from every
   output table with one summary warning, or the call stops with an
   "Inappropriate bandwidth" error once 3 or fewer points remain usable.
3. **Consistency across all five sites.** `adaptive.bw.at()` (used by the
   support diagnostics) and all four `wls.*()` local-fit functions go through
   the same three helpers, so the window used to report ESS/support and the
   window used to fit are always the same number.

A fixed-`bw` fit under the new normalizer is mathematically identical to the
old grid-normalized fit at a rescaled `bw` (`bw * sqrt(g_grid / g)`), so
nothing about the model itself changed -- only what `bw` means. With a FIXED
`bw`, this makes windows wider than in 1.4.0 (typically 1.5-3x on
well-behaved moderators). With the DEFAULT cross-validated `bw`, CV
compensates by picking a smaller bandwidth, so results move only slightly.
Full before/after numbers: `runs/KBW-20260917/audit.md`.

### New internal helpers

All four are non-exported, defined at the top level of `R/kernel.R`, right
after `interflex.kernel()` closes (the same spot the pre-existing
`createFolds` helper already lived).

```mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
    KENT["interflex.kernel()"]
    VBW["Validate bw (Part G)"]
    DENS["density(X, w)"]
    PREP[".prepare_density()"]
    LOOP["Fit each X.eval"]
    LBW[".kernel_local_bw()"]
    NIDX[".nearest_index()"]
    WLSF["wls.*() local fit"]
    ENG["Engine glm/ivreg/feols"]
    VCU[".vcov_usable()"]
    STAT["status / aliased"]
    KEEP["Keep mask + warn/stop"]
    DIFF["gen.kernel.difference()"]
    CVERR["getError.CV()"]
    CVPREP["Prepare fold density"]
    BOOT["Bootstrap loop"]
    BPREP["Prepare boot density"]
    SDIAG["support.diagnostics()"]
    ADPT["adaptive.bw.at()"]

    KENT --> VBW --> DENS --> PREP --> LOOP
    LOOP --> WLSF
    WLSF --> LBW --> NIDX
    WLSF --> ENG
    WLSF --> VCU
    WLSF --> STAT --> KEEP
    KEEP --> DIFF
    KENT --> CVERR --> CVPREP --> WLSF
    KENT --> BOOT --> BPREP --> WLSF
    KENT --> SDIAG --> ADPT --> LBW

    style VBW fill:#1e90ff,stroke:#1565c0,color:#fff
    style PREP fill:#1e90ff,stroke:#1565c0,color:#fff
    style LBW fill:#1e90ff,stroke:#1565c0,color:#fff
    style NIDX fill:#1e90ff,stroke:#1565c0,color:#fff
    style WLSF fill:#1e90ff,stroke:#1565c0,color:#fff
    style VCU fill:#1e90ff,stroke:#1565c0,color:#fff
    style STAT fill:#1e90ff,stroke:#1565c0,color:#fff
    style KEEP fill:#1e90ff,stroke:#1565c0,color:#fff
    style DIFF fill:#1e90ff,stroke:#1565c0,color:#fff
    style CVPREP fill:#1e90ff,stroke:#1565c0,color:#fff
    style BPREP fill:#1e90ff,stroke:#1565c0,color:#fff
    style ADPT fill:#1e90ff,stroke:#1565c0,color:#fff
```

> Blue = new or changed in this run. `ENG` (the `glm`/`ivreg`/`feols` call
> itself) and `DENS` (the base-R `density()` call) are unchanged; everything
> that consumes their output is new or rewritten.

| Function | Purpose | Changed |
| --- | --- | --- |
| `.kernel_nearest_index(grid, v)` | Vectorised nearest-grid-point lookup, ties to the lower index (`findInterval`-based) | new |
| `.kernel_prepare_density(dens, x, w)` | Attaches `adapt.g` (weighted geometric mean at the observations) and `adapt.floor` to a `density()` object | new |
| `.kernel_local_bw(x0, bw, dens)` | `bw * sqrt(adapt.g / f_eff(x0))`; `NA_real_` when unusable (covers `bw` <= 0 too) | new |
| `.kernel_vcov_usable(V, needed)` | Shared usability check for a variance matrix (not NULL, has every needed name, every entry finite) | new |
| `wls.nofe()` / `wls.iv()` / `wls.fe()` / `wls.iv.fe()` | Local weighted fit at one evaluation point; now returns `status`/`aliased` alongside the result | **yes** |
| `adaptive.bw.at()` | Support-diagnostics helper; now a thin wrapper around `.kernel_local_bw()` | **yes** |
| `support.diagnostics.for.bw()` | ESS / `varD` diagnostics per evaluation point; logic unchanged, now consistent with the local fits | no |
| `getError.CV()` | Cross-validation loss per candidate bandwidth; failed candidates keep their own `bw`, fits skip the variance | **yes** |
| `gen.kernel.difference()` | Differences at `diff.values`; now takes a pre-computed `diff.fits` list instead of refitting | **yes** |

### Local-fit status codes (Part B)

Each `wls.*()` function returns
`list(result, model.vcov, model.df, data.touse, status, aliased)`. `result`
is always a full-length, correctly-named vector -- `NA_real_` when
`status != "ok"` -- never a vector with 0s spliced in for coefficients that
could not be estimated.

| `status` | Meaning | Effect on the point |
| --- | --- | --- |
| `"ok"` | usable local fit | kept |
| `"no usable kernel weights"` | `h(x0)` unusable, or no kernel weight is finite and > 0 | dropped |
| `"estimation failed"` | `glm`/`glm.nb`/`ivreg`/`feols` errored or returned a non-model | dropped |
| `"did not converge"` | `glm`/`glm.nb` reports `converged = FALSE` | dropped |
| `"coefficients not identified"` | a coefficient is NA/NaN/+-Inf, or (for `feols`) missing because fixest dropped it for collinearity | dropped; offending names recorded |
| `"variance not estimable"` | only when a variance is requested: the matrix call errors, is missing a name, has a non-finite entry, or (HC2 path) a hat value exceeds `LEVERAGE_MAX` | dropped |

**Drop/stop policy (main fit, Part C).** After fitting every requested
evaluation point, points with `status == "ok"` are kept; `coef.grid`,
`results`, `model.vcovs`, `model.dfs` and `X.eval` are all re-subset together
so they stay aligned (this alignment was broken before this run -- `gen.sd`
iterated over unfiltered `results`, including NA ones). If any point was
dropped, one warning names the first five dropped X values, a reason
breakdown, and (when relevant) the unidentified coefficient names. If 3 or
fewer points remain usable the call stops with an "Inappropriate bandwidth"
error; if half or fewer remain, a warning is raised instead of a stop.

### CV and bootstrap handling

- **Cross-validation (`getError.CV()`, Part E).** Each fold builds its own
  pilot density on the training rows only (`Xdensity.train`, prepared the
  same way as the main fit) and fits candidates with `vcov = FALSE` (CV never
  needs the variance; skipping it makes CV roughly 2x faster). A candidate
  whose fit errors now keeps its own `bw` value in the output row instead of
  losing it to `NA` -- before this run, that loss let the "no finite CV loss"
  fallback silently return `bw = -Inf`. After a bandwidth is selected it is
  checked to be a single finite positive number, or the call stops with a
  clear error.
- **Bootstrap (Part F).** Each replicate builds its own pilot density on the
  resampled rows (`Xdensity.boot`). Replicate fits at each `X.eval` point
  keep `vcov = FALSE` as before; unusable replicate points become NA and are
  excluded by the existing `na.rm` / uniform-band machinery. One warning
  reports how many replicates had at least one unusable point.
- **Differences (`gen.kernel.difference()`, Part D).** The local fits at
  `diff.values` are computed ONCE (`diff.fits`), right after the main fit,
  and reused -- including inside the bootstrap loop, which previously refit
  at `diff.values` on every replicate via a closure lookup that always read
  the MAIN sample's data (a pre-existing, silently-wrong-sample bug for
  `diff.estimate`'s bootstrap SE, fixed as a side effect of passing
  `diff.fits` explicitly per replicate). If any `diff.values` fit is
  unusable, all reported differences and their SEs are NA, with one warning.

### Data flow: one evaluation point

```mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
    SAMP["Sample: full/fold/boot"]
    DENS2["density(X, w)"]
    PREP2["Prepare adapt.g + floor"]
    GCHK{{"adapt.g finite > 0?"}}
    STOPALL["stop: density all-zero"]
    FE["For each X.eval point"]
    LBW2["h(x0)=bw*sqrt(g/f_eff)"]
    KW["Kernel weights k_i"]
    FR["Fit rows: k_i > 0"]
    FIT["Fit local model"]
    USE{{"status == ok?"}}
    KP["Keep: result + vcov"]
    DP["Drop: NA + reason"]
    MRG["Collect all points"]
    ENUF{{"usable > 3?"}}
    STOPFEW["stop: Inappropriate bw"]
    WARNC["warn if any dropped"]
    OUT["est/pred/link.kernel"]

    SAMP --> DENS2 --> PREP2 --> GCHK
    GCHK -- no --> STOPALL
    GCHK -- yes --> FE
    FE --> LBW2 --> KW --> FR --> FIT --> USE
    USE -- yes --> KP
    USE -- no --> DP
    KP --> MRG
    DP --> MRG
    MRG --> ENUF
    ENUF -- no --> STOPFEW
    ENUF -- yes --> WARNC
    WARNC --> OUT

    style PREP2 fill:#1e90ff,stroke:#1565c0,color:#fff
    style GCHK fill:#1e90ff,stroke:#1565c0,color:#fff
    style STOPALL fill:#1e90ff,stroke:#1565c0,color:#fff
    style LBW2 fill:#1e90ff,stroke:#1565c0,color:#fff
    style FR fill:#1e90ff,stroke:#1565c0,color:#fff
    style USE fill:#1e90ff,stroke:#1565c0,color:#fff
    style DP fill:#1e90ff,stroke:#1565c0,color:#fff
    style ENUF fill:#1e90ff,stroke:#1565c0,color:#fff
    style STOPFEW fill:#1e90ff,stroke:#1565c0,color:#fff
    style WARNC fill:#1e90ff,stroke:#1565c0,color:#fff
```

> Narrow vertical chain; the two decisions (`GCHK`, `USE`, `ENUF`) are the
> only branch points and each rejoins within one or two steps. `SAMP` is
> whichever sample built the density at that site: the full data for the
> main fit, the training fold for CV, or the resampled rows for a bootstrap
> replicate -- `PREP2` runs once per sample, not once per evaluation point.

Full validation evidence (per-scenario pass/fail, before/after tables):
`runs/KBW-20260917/audit.md`. Design rationale and rejected alternatives
(e.g. why the old zero-weight nudge had to be removed rather than kept):
`runs/KBW-20260917/spec.md` Part B.2 and `comprehension.md` Q4.

---

## Overview

**interflex** is an R package (v1.4.1) for diagnosing and visualizing multiplicative interaction models. It estimates non-linear marginal effects of a treatment (D) on an outcome (Y) across values of a moderator (X), supporting both discrete and continuous treatments. The package provides eight estimation strategies (linear, binning, kernel, GAM, raw, GRF, DML, lasso), unified behind a single `interflex()` entry point. Key external dependencies include ggplot2 (plotting), mgcv (GAM), grf (causal forests), glmnet (lasso/ridge), and DoubleML/mlr3 (DML estimation). The package is pure R (`NeedsCompilation: no`) -- no compiled code.

**GATE (Group Average Treatment Effects)**: When `gate = TRUE` is specified with a discrete moderator X, estimators compute group-level average treatment effects instead of (or in addition to) smooth conditional marginal effect curves. GATE is supported by `linear`, `grf`, `dml`, and `lasso` estimators. The unified output field `g.est` holds GATE results across all estimators, with standardized column names (`X`, `ME`, `sd`, `lower CI(95%)`, `upper CI(95%)`).

---

## Module Structure

```mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
    subgraph API["API Layer"]
        IFX["interflex.R -- router"]
        PLT["plot.R -- S3 plot method"]
        PRD["predict.R -- S3 predict"]
        TST["inter_test.R -- t-tests"]
        PRT["print.R -- S3 print"]
    end

    subgraph Estimators["Estimator Layer"]
        LIN["linear.R"]
        BIN["binning.R"]
        KER["kernel.R"]
        GAM["gam.R"]
        RAW["raw.R"]
        GRF["grf.R"]
        DML["DML.R -- pure R DML"]
        LAS["lasso.R"]
        LSD["lasso_discrete.R"]
    end

    subgraph DMLHelpers["DML Internal Helpers"]
        SDL[".set_dml_learner()"]
        MRF[".map_rf_params()"]
        MBP[".map_boosting_params()"]
        MNP[".map_nn_params()"]
        BTG[".build_dml_tuning_grid()"]
        CCB[".compute_cate_blp()"]
        CGB[".compute_gate_blp()"]
        EPS[".extract_psi()"]
        SSV[".safe_solve()"]
        RDE[".run_dml_estimation()"]
    end

    subgraph DMLSub["DML Sub-Estimators"]
        CME_I["estimate_cme_irm.R"]
        CME_P["estimate_cme_plr.R"]
        GTE_I["estimate_gte_irm.R"]
        GTE_P["estimate_gte_plr.R"]
    end

    subgraph GATEUtils["GATE Utilities"]
        GTU["gate_utils.R"]
    end

    subgraph Output["Output Layer"]
        PPL["plot_pool.R -- pooled plots"]
    end

    subgraph Utils["Utilities"]
        UTL["utils.R -- shared helpers"]
        UNI["uniform.R -- uniform CI"]
        VCL["vcluster.R -- cluster vcov"]
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
    DML --> SDL
    DML --> RDE
    RDE --> SDL
    RDE --> CCB
    RDE --> CGB
    RDE --> BTG
    CCB --> EPS
    CGB --> EPS
    CCB --> SSV
    CGB --> SSV
    SDL --> MRF
    SDL --> MBP
    SDL --> MNP
    DML --> UNI
    LIN --> GTE_I
    LIN --> GTE_P
    LIN --> GTU
    GRF --> GTU

    style IFX fill:#1e90ff,stroke:#1565c0,color:#fff
    style LIN fill:#1e90ff,stroke:#1565c0,color:#fff
    style GRF fill:#1e90ff,stroke:#1565c0,color:#fff
    style DML fill:#1e90ff,stroke:#1565c0,color:#fff
    style LSD fill:#1e90ff,stroke:#1565c0,color:#fff
    style PLT fill:#1e90ff,stroke:#1565c0,color:#fff
    style GTU fill:#1e90ff,stroke:#1565c0,color:#fff
```

> One unified diagram. Subgraph layers group related modules. Blue fill = modified in this run (GATE generalization).

### Module Reference

| Module / File | Layer | Purpose | Key Exports | Changed |
| --- | --- | --- | --- | --- |
| `R/interflex.R` | API | Main entry point; validates inputs, builds `treat.info`/`diff.info`, routes to estimator; validates `gate` param | `interflex()` | **yes** |
| `R/plot.R` | API | S3 `plot.interflex()` method; renders marginal effect plots with density/histogram overlays; supports `by.group` via unified `g.est` | `plot.interflex()` | **yes** |
| `R/predict.R` | API | S3 `predict.interflex()` method; computes predicted marginal effects at new X values | `predict.interflex()` | no |
| `R/inter_test.R` | API | Post-estimation t-test for difference in marginal effects (dml style) | `inter.test()` | no |
| `R/print.R` | API | S3 `print.interflex()` -- auto-prints the figure attached to an interflex object | `print.interflex()` | yes |
| `R/linear.R` | Estimator | Linear interaction model with delta/bootstrap/simulation variance; GATE via `bootstrapGTE`/`bootstrapGATE_PLR` | `interflex.linear()` | **yes** |
| `R/binning.R` | Estimator | Binning estimator: splits X into bins, estimates within-bin effects | `interflex.binning()` | no |
| `R/kernel.R` | Estimator | Kernel estimator: local polynomial regression; adaptive (Abramson) bandwidth normalized at the data points; degenerate local fits dropped/stopped, never zero-filled | `interflex.kernel()` | **yes** |
| `R/gam.R` | Estimator | GAM estimator via `mgcv::gam()` with 3D visualization | `interflex.gam()` | no |
| `R/raw.R` | Estimator | Raw data scatter plots with LOESS smoothing | `interflex.raw()` | no |
| `R/grf.R` | Estimator | Generalized random forests via `grf::causal_forest()`; GATE via aggregated forest CATEs | `interflex.grf()` | **yes** |
| `R/DML.R` | Estimator | Pure R DML via DoubleML + mlr3; outputs both `g.est` and `g.est.dml` (deprecated alias) | `interflex.dml()` | **yes** |
| `R/lasso.R` | Estimator | Lasso/ridge DML for continuous moderators; calls CME/GTE sub-estimators | `interflex.lasso()` | no |
| `R/lasso_discrete.R` | Estimator | Lasso/ridge DML for discrete moderators (<5 unique X values); adds `g.est` output field | `interflex.lasso_discrete()` | **yes** |
| `R/gate_utils.R` | GATE Utils | Column name standardization for bootstrap GATE returns to unified format | `.standardize_gate_columns()` | **new** |
| `R/estimate_cme_irm.R` | DML Sub | CME estimation via AIPW-Lasso (binary treatment, IRM) | `estimateCME()` | no |
| `R/estimate_cme_plr.R` | DML Sub | CME estimation via PO-Lasso (continuous treatment, PLRM) | `estimateCME_PLR()` | no |
| `R/estimate_gte_irm.R` | DML Sub | Group treatment effects via AIPW-Lasso (binary treatment, discrete X); `bootstrapGTE()` | `estimateGTE()`, `bootstrapGTE()` | no |
| `R/estimate_gte_plr.R` | DML Sub | Group treatment effects via PO-Lasso (continuous treatment, discrete X); `bootstrapGATE_PLR()` | `estimateGATE_PLR()`, `bootstrapGATE_PLR()` | no |
| `R/plot_pool.R` | Output | Pooled multi-treatment plot with overlaid CIs | `interflex.plot.pool()` | no |
| `R/utils.R` | Utils | Shared internal helpers: treat.info extraction, density, histograms | (internal: dot-prefixed) | no |
| `R/uniform.R` | Utils | Uniform confidence interval quantiles via bootstrap/delta method | `calculate_uniform_quantiles()`, `calculate_delta_uniformCI()` | no |
| `R/vcluster.R` | Utils | Cluster-robust variance-covariance matrix computation | `vcovCluster()` | no |
| `DESCRIPTION` | Config | Package metadata; Imports, Depends | N/A | no |
| `NAMESPACE` | Config | Export pattern, S3 methods, importFrom declarations | N/A | no |

---

## Function Call Graph

### Main Pipeline

```mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
    USR["User call"] --> IFX["interflex()"]
    IFX --> VAL["Input validation"]
    VAL --> GVAL["gate validation"]
    GVAL --> TI["Build treat.info"]
    TI --> DI["Build diff.info"]
    DI --> ROUTE{{"estimator?"}}

    ROUTE -- linear --> LIN["interflex.linear()"]
    ROUTE -- binning --> BIN["interflex.binning()"]
    ROUTE -- kernel --> KER["interflex.kernel()"]
    ROUTE -- gam --> GAM["interflex.gam()"]
    ROUTE -- raw --> RAW["interflex.raw()"]
    ROUTE -- grf --> GRF["interflex.grf()"]
    ROUTE -- dml --> DML["interflex.dml()"]
    ROUTE -- lasso --> LSPLIT{{"gate || X < 5?"}}
    LSPLIT -- yes --> LSD["interflex.lasso_discrete()"]
    LSPLIT -- no --> LAS["interflex.lasso()"]

    LIN --> LINGATE{{"gate?"}}
    LINGATE -- yes,binary D --> BGTE["bootstrapGTE()"]
    LINGATE -- yes,continuous D --> BGPLR["bootstrapGATE_PLR()"]
    LINGATE -- no --> LINCME["Linear CME output"]
    BGTE --> STDCOL[".standardize_gate_columns()"]
    BGPLR --> STDCOL

    GRF --> GRFGATE{{"gate?"}}
    GRFGATE -- yes --> GRFAGG["Aggregate forest CATEs"]
    GRFGATE -- no --> GRFCME["GRF CATE output"]

    DML --> RDE[".run_dml_estimation()"]
    RDE --> SDL[".set_dml_learner()"]
    RDE --> BTG[".build_dml_tuning_grid()"]
    RDE --> DMLFIT["DoubleML fit + tune"]
    DMLFIT --> CATE[".compute_cate_blp()"]
    DMLFIT --> GATE[".compute_gate_blp()"]
    CATE --> EPS[".extract_psi()"]
    GATE --> EPS
    CATE --> SSV[".safe_solve()"]
    GATE --> SSV
    DML --> UNI["calculate_delta_uniformCI()"]

    style IFX fill:#1e90ff,stroke:#1565c0,color:#fff
    style GVAL fill:#1e90ff,stroke:#1565c0,color:#fff
    style LIN fill:#1e90ff,stroke:#1565c0,color:#fff
    style LINGATE fill:#1e90ff,stroke:#1565c0,color:#fff
    style BGTE fill:#1e90ff,stroke:#1565c0,color:#fff
    style BGPLR fill:#1e90ff,stroke:#1565c0,color:#fff
    style STDCOL fill:#1e90ff,stroke:#1565c0,color:#fff
    style GRF fill:#1e90ff,stroke:#1565c0,color:#fff
    style GRFGATE fill:#1e90ff,stroke:#1565c0,color:#fff
    style GRFAGG fill:#1e90ff,stroke:#1565c0,color:#fff
    style LSPLIT fill:#1e90ff,stroke:#1565c0,color:#fff
```

### Output Pipeline

```mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
    OUT["interflex output object"] --> PLOT["plot.interflex()"]
    OUT --> PRED["predict.interflex()"]
    OUT --> TEST["inter.test() (inter_test.R)"]
    OUT --> PRT["print.interflex() (print.R)"]
    PLOT --> BYGRP{{"by.group?"}}
    BYGRP -- yes --> GEST["Read out$g.est"]
    BYGRP -- no --> POOL{{"pool = TRUE?"}}
    POOL -- yes --> PPL["interflex.plot.pool()"]
    POOL -- no --> GPLT["ggplot2 rendering"]
    GEST --> GPLT
    PPL --> GPLT
    PRED --> GPLT

    style PLOT fill:#1e90ff,stroke:#1565c0,color:#fff
    style BYGRP fill:#1e90ff,stroke:#1565c0,color:#fff
    style GEST fill:#1e90ff,stroke:#1565c0,color:#fff
```

### Function Reference

| Function | Defined In | Called By | Calls | Changed | Purpose |
| --- | --- | --- | --- | --- | --- |
| `interflex()` | `R/interflex.R` | user (exported) | all estimators | **yes** | Validate inputs (incl. gate), build metadata, route to estimator |
| `plot.interflex()` | `R/plot.R` | user (S3 method) | `interflex.plot.pool()` | **yes** | Render marginal effect or GATE plots; unified `g.est` check |
| `predict.interflex()` | `R/predict.R` | user (S3 method) | ggplot2 | no | Compute and plot predicted marginal effects |
| `inter.test()` | `R/inter_test.R` | user (exported) | mgcv::gam | no | Test differences in marginal effects |
| `interflex.linear()` | `R/linear.R` | `interflex()` | `bootstrapGTE`, `bootstrapGATE_PLR`, `.standardize_gate_columns` | **yes** | Linear interaction model; GATE via bootstrap GTE/GATE_PLR |
| `interflex.grf()` | `R/grf.R` | `interflex()` | `grf::causal_forest`, `predict` | **yes** | GRF CATE; GATE via aggregated forest predictions |
| `interflex.dml()` | `R/DML.R` | `interflex()` | `.run_dml_estimation`, `.compute_cate_blp`, `.compute_gate_blp` | **yes** | Pure R DML; output adds unified `g.est` field |
| `interflex.lasso_discrete()` | `R/lasso_discrete.R` | `interflex()` | `bootstrapGTE`, `bootstrapGATE_PLR` | **yes** | Lasso DML for discrete X; output adds `g.est` field |
| `.standardize_gate_columns()` | `R/gate_utils.R` | `interflex.linear()`, `interflex.grf()` | -- | **new** | Map bootstrap column names to unified GATE format |
| `.run_dml_estimation()` | `R/DML.R` | `interflex.dml()` | `.set_dml_learner`, `.build_dml_tuning_grid`, `.compute_cate_blp`, `.compute_gate_blp` | no | Core DML worker |
| `.compute_cate_blp()` | `R/DML.R` | `.run_dml_estimation` | `.extract_psi`, `.safe_solve`, `splines::bs` | no | CATE via BLP of pseudo-outcomes onto B-spline basis |
| `.compute_gate_blp()` | `R/DML.R` | `.run_dml_estimation` | `.extract_psi`, `.safe_solve` | no | GATE via BLP of pseudo-outcomes onto group dummies |
| `bootstrapGTE()` | `R/estimate_gte_irm.R` | `interflex.linear()`, `interflex.lasso_discrete()` | `estimateGTE`, `glmnet` | no | Bootstrap GTE for binary treatment (IRM framework) |
| `bootstrapGATE_PLR()` | `R/estimate_gte_plr.R` | `interflex.linear()`, `interflex.lasso_discrete()` | `estimateGATE_PLR`, `glmnet` | no | Bootstrap GATE for continuous treatment (PLR framework) |

---

## Data Flow

```mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
    INPUT["User: interflex(gate=TRUE, ...)"]
    INPUT --> VALIDATE["Validate inputs & coerce types"]
    VALIDATE --> GATEVAL{{"gate=TRUE?"}}
    GATEVAL -- yes --> GCHECK["Check: X discrete, estimator supported"]
    GCHECK --> TREAT["Build treat.info metadata"]
    GATEVAL -- no --> TREAT
    TREAT --> DIFF["Build diff.info for contrasts"]
    DIFF --> ROUTE{{"Select estimator"}}

    ROUTE --> LINENTRY["interflex.linear()"]
    LINENTRY --> LINCME["Compute smooth CME"]
    LINCME --> LINGATQ{{"gate?"}}
    LINGATQ -- yes,binary --> BGTE["bootstrapGTE(linear)"]
    LINGATQ -- yes,continuous --> BGPLR["bootstrapGATE_PLR(linear)"]
    LINGATQ -- no --> LINOUT["Return est.lin only"]
    BGTE --> STDCOL["Standardize columns"]
    BGPLR --> STDCOL
    STDCOL --> LINGEST["Attach g.est"]
    LINGEST --> RETURN["Return interflex object"]
    LINOUT --> RETURN

    ROUTE --> GRFENTRY["interflex.grf()"]
    GRFENTRY --> GRFCATE["Fit causal_forest"]
    GRFCATE --> GRFGATQ{{"gate?"}}
    GRFGATQ -- yes --> GRFAGG["Aggregate CATEs by X group"]
    GRFAGG --> GRFGEST["Attach g.est"]
    GRFGATQ -- no --> GRFOUT["Return est.grf only"]
    GRFGEST --> RETURN
    GRFOUT --> RETURN

    style GATEVAL fill:#1e90ff,stroke:#1565c0,color:#fff
    style GCHECK fill:#1e90ff,stroke:#1565c0,color:#fff
    style LINGATQ fill:#1e90ff,stroke:#1565c0,color:#fff
    style BGTE fill:#1e90ff,stroke:#1565c0,color:#fff
    style BGPLR fill:#1e90ff,stroke:#1565c0,color:#fff
    style STDCOL fill:#1e90ff,stroke:#1565c0,color:#fff
    style LINGEST fill:#1e90ff,stroke:#1565c0,color:#fff
    style GRFGATQ fill:#1e90ff,stroke:#1565c0,color:#fff
    style GRFAGG fill:#1e90ff,stroke:#1565c0,color:#fff
    style GRFGEST fill:#1e90ff,stroke:#1565c0,color:#fff
```

---

## Key Data Structures

### `treat.info` (built by `interflex()`, consumed by all estimators)

A named list containing treatment metadata. The `.extract_treat_info()` utility unpacks this uniformly.

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
| `est.lin` / `est.bin` / `est.kernel` / `est.dml` / etc. | Marginal effect estimates data frame |
| `g.est` | **Unified GATE estimates** (when `gate = TRUE`) -- named list of data.frames keyed by treatment arm |
| `g.est.dml` | Deprecated alias for `g.est` (DML only, backward compatibility) |
| `dml.models` | Tuned model info: `model.y`, `model.t` (DML only) |
| `dml.losses` | Nuisance losses from DML fit (DML only) |
| `diff.estimate` | Treatment contrast estimates |
| `figure` | ggplot object(s) |
| `hist.out`, `treat.hist`, `de`, `treat_den` | Distribution data for X-axis overlays |
| `treat.info`, `diff.info` | Metadata passed through |

### Unified GATE output columns (`g.est` data.frame)

Each element of `output$g.est` is a data.frame with standardized columns:

| Column | Type | Description |
| --- | --- | --- |
| `X` | numeric | Moderator level |
| `ME` | numeric | Group average treatment effect estimate |
| `sd` | numeric | Standard error |
| `lower CI(95%)` | numeric | Lower pointwise 95% CI |
| `upper CI(95%)` | numeric | Upper pointwise 95% CI |
| `lower uniform CI(95%)` | numeric | Lower uniform 95% CI (when available) |
| `upper uniform CI(95%)` | numeric | Upper uniform 95% CI (when available) |

When `CI = FALSE`, only `X` and `ME` columns are present.

---

## Estimator Architecture

| Estimator | Function | Treatment Type | Moderator Type | Method | Variance | GATE Support |
| --- | --- | --- | --- | --- | --- | --- |
| `"linear"` | `interflex.linear()` | discrete or continuous | continuous | Parametric OLS/GLM with D*X interaction | delta, bootstrap, simulation | **yes** (via bootstrapGTE/bootstrapGATE_PLR) |
| `"binning"` | `interflex.binning()` | discrete or continuous | continuous (binned) | Split X into bins, within-bin linear models | delta, bootstrap, simulation | no |
| `"kernel"` | `interflex.kernel()` | discrete or continuous | continuous | Local polynomial regression, adaptive (Abramson) bandwidth via CV or fixed `bw` | bootstrap | no |
| `"gam"` | `interflex.gam()` | continuous only | continuous | `mgcv::gam()` smooth surface | GAM built-in | no |
| `"raw"` | `interflex.raw()` | discrete or continuous | continuous | Scatter + LOESS (no formal estimation) | none | no |
| `"grf"` | `interflex.grf()` | binary | continuous | `grf::causal_forest()` | forest-based | **yes** (via aggregated CATEs) |
| `"dml"` | `interflex.dml()` | discrete or continuous | continuous | R DoubleML + mlr3 cross-fitting with BLP CATE/GATE | HC sandwich + uniform CI | **yes** (via `.compute_gate_blp()`) |
| `"lasso"` | `interflex.lasso()` / `interflex.lasso_discrete()` | binary or continuous | continuous/discrete | PO-Lasso (PLRM) or AIPW-Lasso (IRM) | bootstrap | **yes** (via bootstrapGTE/bootstrapGATE_PLR) |

### GATE Estimation Paths

| Estimator | Binary D | Continuous D | Method |
| --- | --- | --- | --- |
| `linear` | `bootstrapGTE()` with linear nuisance | `bootstrapGATE_PLR()` with linear nuisance | Bootstrap inference, column standardization via `.standardize_gate_columns()` |
| `grf` | Aggregate `causal_forest` CATEs by X group | N/A (GRF is binary-D only) | Forest variance: `SE = sqrt(mean(var_i) / n_group)` |
| `dml` | `.compute_gate_blp()` with IRM pseudo-outcomes | `.compute_gate_blp()` with PLR pseudo-outcomes | HC sandwich variance, uniform CIs |
| `lasso` | `bootstrapGTE()` via `lasso_discrete` | `bootstrapGATE_PLR()` via `lasso_discrete` | Bootstrap inference |

### DML Estimation Pipeline (Pure R)

- `R/DML.R` uses `DoubleML` R package (R6 classes) for core DML estimation
- `mlr3` + `mlr3learners` for ML model backends (ranger, glmnet, lightgbm, nnet)
- `DoubleML::DoubleMLIRM` for binary treatment, `DoubleML::DoubleMLPLR` for continuous
- CATE computed via manual BLP: project pseudo-outcomes onto B-spline basis
- GATE computed via manual BLP: project pseudo-outcomes onto group dummies
- Uniform CIs via existing `calculate_delta_uniformCI()` from `uniform.R`
- Parameter mapping tables translate sklearn names to mlr3/ranger/lightgbm equivalents

### DML Dependencies

| Package | Role | Import Type |
| --- | --- | --- |
| `DoubleML` | Core DML framework (R6 classes: `DoubleMLIRM`, `DoubleMLPLR`) | Imports |
| `mlr3` | ML task framework, `lrn()` for learner creation | Imports |
| `mlr3learners` | Standard learner implementations | Imports |
| `ranger` | Random forest backend (default model) | Imports |
| `data.table` | Required by DoubleML internals | Imports |
| `paradox` | Tuning parameter sets | Imports |
| `mlr3tuning` | Grid search tuning (when CV=TRUE) | Suggests |
| `lightgbm` | Gradient boosting backend (non-default) | Suggests |
| `nnet` | Neural network backend (non-default) | Suggests |

### Parallel RNG Dependencies

| Package | Role | Import Type |
| --- | --- | --- |
| `doFuture` | Future-based foreach backend; registers via `registerDoFuture()` | Imports |
| `doRNG` | Provides `%dorng%` for reproducible parallel RNG in bootstrap loops (used in `R/utils.R`) | Imports |
| `future` | Plan-based parallel execution (`plan(multisession)`) | Imports (pre-existing) |
| `parallelly` | `availableCores()` for robust core detection in legacy files | Imports (pre-existing) |

### Dependencies Removed (Previous Runs)

| Package | Reason |
| --- | --- |
| `reticulate` | No longer needed -- all Python code eliminated |
| `inst/python/dml.py` | Deleted -- Python DML engine replaced by R code |

---

## Utility Functions

| Function | Purpose | Used By |
| --- | --- | --- |
| `.extract_treat_info(treat.info)` | Unpacks the `treat.info` list into local variables | 9 files: DML, binning, kernel, linear, grf, lasso, lasso_discrete, raw, inter_test |
| `.compute_density(...)` | Computes kernel density estimates for X-axis distribution overlay | 7 files: DML, binning, kernel, linear, grf, lasso, lasso_discrete |
| `.compute_histograms(...)` | Computes histogram bin counts for X-axis distribution overlay | 7 files: DML, binning, kernel, linear, grf, lasso, lasso_discrete |
| `.standardize_gate_columns(nms, effect_col)` | Maps bootstrap GATE column names to unified format | linear.R, grf.R |

---

## Architectural Patterns

- **Router pattern**: `interflex()` is a monolithic router (~1400 lines) that validates all inputs, builds shared metadata (`treat.info`, `diff.info`), and dispatches to one of 9 estimator functions. Each estimator is a standalone function in its own file.

- **GATE as an orthogonal extension**: GATE support is layered on top of existing estimators via the `gate = TRUE` parameter. Each GATE-capable estimator first computes its normal CME output, then appends `g.est` if gate is requested. This avoids disrupting existing code paths.

- **Unified output contract**: All GATE-capable estimators produce `output$g.est` with identical column schema (`X`, `ME`, `sd`, CIs). The plotting code consumes `g.est` without needing to know which estimator produced it.

- **Column standardization bridge**: Bootstrap functions (`bootstrapGTE`, `bootstrapGATE_PLR`) return columns with different names (`GTE`/`GATE`, `SE`, `CI.lower`). The `.standardize_gate_columns()` helper maps these to the DML-originated format that the plotting code expects.

- **Shared metadata**: `treat.info` and `diff.info` are computed once by the router and passed to every estimator. The `.extract_treat_info()` utility provides uniform unpacking.

- **Inline plotting**: Each estimator builds its own ggplot figure internally rather than delegating to a separate plot function. The S3 `plot.interflex()` method re-renders from stored data.

- **Dot-prefix convention**: Internal helpers use `.` prefix to avoid export via the blanket `exportPattern("^[[:alpha:]]+")` rule.

- **Lasso moderator cardinality split**: The `"lasso"` estimator auto-selects `interflex.lasso_discrete()` when X has fewer than 5 unique values OR when `gate = TRUE` is explicitly set, switching from CME to GTE estimation.

- **Parameter compatibility layer**: The DML estimator accepts sklearn-style parameter names and maps them to R equivalents via `.map_*_params()` helpers.

- **Manual BLP for CATE/GATE**: Since the R DoubleML package lacks `.cate()` and `.gate()` methods, CATE and GATE are computed via Best Linear Projection of pseudo-outcomes onto basis functions, with HC sandwich variance estimation.

- **doFuture parallel backend**: All parallel `foreach` loops use `doFuture::registerDoFuture()` + `future::plan(future::multisession)` with `.options.future = list(seed = TRUE)` for reproducible L'Ecuyer-CMRG parallel RNG streams.

---

## Notes

- **Previous run (refactor, 2026-03-15)**: -506 lines across 21 files. ~500 lines of duplicated code consolidated into `R/utils.R`. ~700 redundant boolean comparisons cleaned.
- **Previous run (bs merge, 2026-03-16)**: 5 files created, 4 files modified. Selectively integrated bs branch features (ttest.R, B-spline expansion, parameter defaults).
- **Previous run (Python-to-R DML migration, 2026-03-16)**: `R/DML.R` completely rewritten (~244 lines old to ~430 lines new). `inst/python/dml.py` deleted (298 lines). Users no longer need Python to use the DML estimator.
- **Previous run (parallel RNG migration, 2026-04-03)**: 7 R files modified (8 parallel blocks). Replaced `doParallel` with `doFuture` for reproducible parallel RNG.
- **Previous run (GATE generalization, 2026-04-04)**: 7 R files modified, 1 new file (`gate_utils.R`), 1 new test file (`test-gate.R` with 49 tests). Generalized GATE support from DML-only to linear, grf, dml, and lasso estimators. Unified output field `g.est` with backward-compatible `g.est.dml` alias. Added input validation for `gate` parameter. Fixed 4 bugs during builder respawn (DML.R encoding, linear.R treatment label mismatch, grf.R column access, lasso_discrete.R type coercion).
- **This run (BOOK-003, 2026-04-07)**: Plot-layer cleanup and ch6 vignette restructure. Files modified: `R/plot.R`, `R/plot_pool.R`, `R/raw.R`, `R/predict.R`, `R/interflex.R`, `vignettes/06_discrete.qmd`, plus new test file `tests/testthat/test-plot-limits.R` (16 tests, all PASS). Three new internal helpers in `R/plot.R`: `.pad_xlim()` (2% symmetric padding for user xlim), `.append_yrange_ci()` (defensive yrange CI column accumulation), `.rename_est_ci()` (defensive colnames assignment for narrow estimator tables). Migrated all `xlim()`/`ylim()`/`scale_*_continuous(limits=)` calls in plot builders to `coord_cartesian()` so visual clipping no longer drops underlying ribbon data. Group-equalization loops in `plot.interflex` and `predict.interflex` collapsed to a single authoritative `coord_cartesian` per panel (one-coord rule). Added user-vs-default sentinel: `interflex()` sets `interflex.user_xlim_explicit` / `interflex.user_ylim_explicit` options at entry (8-line additive block in `R/interflex.R`); `plot.interflex()` reads them, snapshots `.user_xlim_in`/`.user_ylim_in`, and stamps them as attributes on the returned graph for cross-call recovery. This was required because the auto-trim feature added in commit d9b3075 makes raw `xlim` indistinguishable from user input inside `plot.interflex`. ch6 of the Quarto book split nine `dis_out_*` chunks into `ch6-*-fit` (`cache=TRUE`, fit only) + `ch6-*-plot` (`cache=FALSE`, plotting only) pairs, mirroring the ch2 pattern; second-render time drops from ~18 minutes to ~64 seconds (~17x speedup). Three builder respawns were needed: (1) xlim/ylim threading + write-surface expansion to `R/interflex.R`, (2) `.rename_est_ci` helper for narrow-DML colnames defect surfaced by Check 10b, (3) one-character em-dash -> `--` cleanup on `R/plot.R:52` to unblock `devtools::load_all` and clear a new `R CMD check` non-ASCII finding. Test-spec was revised mid-run from `ggplot_build($figure)` grob introspection (which inspects the wrapped canvas, not the inner ME plot) to a proper testthat unit-test file using `plot.interflex(out, show.all = TRUE)` to access the raw inner ggplot list `p.group` directly -- this is reusable regression coverage for any future plot-layer change.
- **No formal test suite prior to this run**: The package did not have tests under `tests/` before the DML migration. Test files have been added incrementally.
- **This run (KBW-20260917, 2026-09-20)**: Kernel estimator only (`R/kernel.R`, +389/-162 lines), plus small cookbook cleanups (`R/interflex.R`, `R/uniform.R`, `R/plot_pool.R`, `R/predict.R`, `R/raw.R`) and release-prep files (`DESCRIPTION` to 1.4.1, `.Rbuildignore`, `man/inter_test.Rd` `\value`, new `tests/testthat/test-kernel-adaptive-bw.R`). Fixed the adaptive-bandwidth normalizer (Abramson's rule): `g` is now the sampling-weighted geometric mean of the pilot density AT THE OBSERVATIONS, computed once per sample right after each `density()` call, instead of over the whole 512-point density grid (which included the empty tails and made windows collapse on moderators with gaps or long tails -- the "subscript out of bounds" / "$ operator is invalid" crashes on Malesky-style raw-X data). The four `wls.*()` local-fit functions and `adaptive.bw.at()` (support diagnostics) all now go through three new non-exported helpers (`.kernel_nearest_index()`, `.kernel_prepare_density()`, `.kernel_local_bw()`) so the normalizer and the density lookup have one definition. Every local fit now returns a `status` (`"ok"` or one of five reason codes) instead of silently zero-filling unidentified coefficients or crashing on an index error; unusable evaluation points are dropped with one summary warning, or the call stops when 3 or fewer points remain. See the "Kernel adaptive bandwidth (KBW-20260917)" section above for the full diagram and the status-code table. User-visible consequence: with a FIXED `bw`, local windows are now wider than in 1.4.0 (1.5-3x on well-behaved moderators; the vignette's `bw` sentence and the `man/interflex.Rd` `bw` item were rewritten accordingly); with the DEFAULT cross-validated `bw`, results move only slightly, since CV compensates by picking a smaller bandwidth. Full validation: `runs/KBW-20260917/audit.md` (AUDIT PASS -- R1 max rel. error 1.9e-14; R-MAL 7/7 previously-crashing Malesky raw-X calls now complete; R-SIM 0/200 failures across 5 DGPs; R-ID rescaling identity to 1.3e-13; R-CV worst case 0.12 CI half-widths / SE ratio within [0.99, 1.02]; R-EDGE E1-E9 as specified; R-OPT clean). Known pre-existing bug left unfixed by explicit spec decision: `wls.iv.fe` (kernel + IV + FE) still errors with `object 'excluded.iv' not found`.
