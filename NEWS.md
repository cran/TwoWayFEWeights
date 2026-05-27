# News

## 2.1.0

Performance / internals

- Re-implement core weight computations on top of `data.table` and Rcpp C++
  kernels, removing the runtime dependencies on `dplyr`, `magrittr`, and
  `rlang`.
- Add `cpp_weighted_mean`, `cpp_rev_cumsum_by_group`, `cpp_fdtr_wtilde2`, and
  `cpp_feS_delta` for the per-row group computations that previously used
  `dplyr::group_by()` + `dplyr::lead()` / `dplyr::lag()`.

## 2.0.3.99 (dev version)

Internals

- Add unit tests (#13 @grantmcdermott)
