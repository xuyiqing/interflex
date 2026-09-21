# CRAN submission comments

## Submission

interflex 1.4.1, an update to interflex 1.4.0, which is currently on CRAN.

## Reason for this submission

This release fixes a bug in the kernel estimator (`estimator = "kernel"`). On
moderators with gaps or long tails, the adaptive bandwidth's normalizing
constant was computed over the whole 512-point density grid, including empty
regions beyond the data range. That could collapse local windows onto a
single observation and crash with "subscript out of bounds" or "$ operator
is invalid for atomic vectors". The normalizer is now computed at the
observations only (Abramson's square-root rule; Silverman 1986, Sec. 5.3.1).
Local fits that cannot be identified, fail to converge, or have no usable
kernel weights are now dropped with a warning, or the call stops with an
informative error, instead of being silently filled with zeros or crashing.

With a fixed `bw`, this makes kernel windows wider than in 1.4.0 (typically
1.5-3x on well-behaved moderators). Results under the default
cross-validated `bw` change only slightly, since cross-validation picks a
smaller bandwidth to compensate. See the "interflex 1.4.1" section of
NEWS.md for the full list of changes.

## Test environments

* local: macOS 26.6.2 (25G83), aarch64-apple-darwin20, R 4.5.3 (checked 2026-09-20)
* No GitHub Actions / win-builder / macbuilder run yet -- this submission
  has not been pushed; run those before the actual CRAN upload if broader
  platform coverage is wanted.

## R CMD check results

0 errors | 0 warnings | 0 notes

The package's full test suite (107 tests; all but three smoke tests are
skipped on CRAN with `skip_on_cran()`) also passes locally against the
built tarball: 0 failures, 0 errors.

The mlr3extralearners Suggests item appeared only as informational text
under "checking CRAN incoming feasibility", not as a NOTE, in the local
--as-cran run; the "Notes for the reviewer" paragraph about it below is
kept as-is, since it may still surface as a NOTE on CRAN's own submission
infrastructure.

## Reverse dependencies

None expected to be affected. This is a behavioral fix inside the kernel
estimator's internal bandwidth calculation: no exported function's
signature, arguments, or return structure changed.

## Notes for the reviewer

* The package Suggests 'mlr3extralearners', which is not on CRAN; it is
  declared via `Additional_repositories: https://mlr-org.r-universe.dev`
  and used only for optional double-machine-learning backends, guarded by
  `requireNamespace()`. This is unchanged from 1.4.0 and is expected to
  produce the standard "Suggests or Enhances not in mainstream
  repositories" NOTE.
