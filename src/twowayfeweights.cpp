#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Inputs are assumed to be sorted by (g, t) where required.
// All kernels treat NA in the value vectors as 0 to mirror the
// na.rm = TRUE behaviour of the R code being replaced.

// [[Rcpp::export]]
double cpp_weighted_mean(NumericVector x, NumericVector w) {
  const R_xlen_t n = x.size();
  if (w.size() != n) stop("cpp_weighted_mean: x and w must have equal length");
  long double num = 0.0L, den = 0.0L;
  for (R_xlen_t i = 0; i < n; ++i) {
    const double xi = x[i], wi = w[i];
    if (NumericVector::is_na(xi) || NumericVector::is_na(wi)) continue;
    num += static_cast<long double>(xi) * static_cast<long double>(wi);
    den += static_cast<long double>(wi);
  }
  if (den == 0.0L) return NA_REAL;
  return static_cast<double>(num / den);
}

// Reverse cumulative sum within each contiguous group.
// Equivalent to: ave(x, g, FUN = function(v) rev(cumsum(rev(v))))
// when rows are already sorted so that group ids are contiguous.
// [[Rcpp::export]]
NumericVector cpp_rev_cumsum_by_group(IntegerVector g, NumericVector x) {
  const R_xlen_t n = g.size();
  if (x.size() != n) stop("cpp_rev_cumsum_by_group: g and x must have equal length");
  NumericVector out(n);
  if (n == 0) return out;

  R_xlen_t end = n - 1;
  while (end >= 0) {
    const int gid = g[end];
    R_xlen_t start = end;
    while (start - 1 >= 0 && g[start - 1] == gid) --start;
    long double acc = 0.0L;
    for (R_xlen_t i = start; i <= end; ++i) {
      const double xi = x[i];
      if (!NumericVector::is_na(xi)) acc += static_cast<long double>(xi);
    }
    long double running = acc;
    for (R_xlen_t i = start; i <= end; ++i) {
      out[i] = static_cast<double>(running);
      const double xi = x[i];
      if (!NumericVector::is_na(xi)) running -= static_cast<long double>(xi);
    }
    end = start - 1;
  }
  return out;
}

// Forward-looking residual reweighting used by the fdTR branch.
// For rows sorted by (g, t):
//   if next row in same group has TFactorNum + 1 == lead(TFactorNum):
//     w_tilde_2 = eps_2 - lead(eps_2) * (lead(P_gt) / P_gt)
//   else (or if result is NA / non-finite): w_tilde_2 = eps_2
// [[Rcpp::export]]
NumericVector cpp_fdtr_wtilde2(IntegerVector g,
                               NumericVector t,
                               NumericVector eps_2,
                               NumericVector P_gt) {
  const R_xlen_t n = g.size();
  if (t.size() != n || eps_2.size() != n || P_gt.size() != n)
    stop("cpp_fdtr_wtilde2: all inputs must have equal length");
  NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    const double eps_i = eps_2[i];
    double val = NA_REAL;
    const bool has_next = (i + 1 < n) && (g[i + 1] == g[i]);
    if (has_next) {
      const double t_i  = t[i];
      const double t_n  = t[i + 1];
      const double P_i  = P_gt[i];
      const double P_n  = P_gt[i + 1];
      const double e_n  = eps_2[i + 1];
      if (!NumericVector::is_na(t_i) && !NumericVector::is_na(t_n) &&
          (t_i + 1.0 == t_n) &&
          !NumericVector::is_na(P_i) && P_i != 0.0 &&
          !NumericVector::is_na(P_n) && !NumericVector::is_na(e_n) &&
          !NumericVector::is_na(eps_i)) {
        val = eps_i - e_n * (P_n / P_i);
      }
    }
    if (NumericVector::is_na(val) || !std::isfinite(val)) val = eps_i;
    out[i] = val;
  }
  return out;
}

// Lag-based delta for the feS branch.
// For rows sorted by (g, t):
//   if previous row in same group has TFactorNum - 1 == lag(TFactorNum):
//     delta_D    = D - lag(D)
//     keep       = TRUE
//     abs_delta  = |delta_D|
//     s_gt       = sign(delta_D)  (1 / -1 / 0)
//     nat_weight = P_gt * abs_delta
//   else: keep = FALSE; the other slots are filled with safe defaults so
//   that callers can subset by `keep` on the R side.
// [[Rcpp::export]]
List cpp_feS_delta(IntegerVector g,
                   NumericVector t,
                   NumericVector D,
                   NumericVector P_gt) {
  const R_xlen_t n = g.size();
  if (t.size() != n || D.size() != n || P_gt.size() != n)
    stop("cpp_feS_delta: all inputs must have equal length");

  NumericVector delta_D(n);
  IntegerVector s_gt(n);
  NumericVector abs_delta_D(n);
  NumericVector nat_weight(n);
  LogicalVector keep(n);

  for (R_xlen_t i = 0; i < n; ++i) {
    const bool has_prev = (i > 0) && (g[i - 1] == g[i]);
    bool ok = false;
    double dD = NA_REAL;
    if (has_prev) {
      const double t_i = t[i];
      const double t_p = t[i - 1];
      if (!NumericVector::is_na(t_i) && !NumericVector::is_na(t_p) &&
          (t_i - 1.0 == t_p) &&
          !NumericVector::is_na(D[i]) && !NumericVector::is_na(D[i - 1])) {
        dD = D[i] - D[i - 1];
        ok = !NumericVector::is_na(dD);
      }
    }

    if (ok) {
      delta_D[i]     = dD;
      const double a = std::fabs(dD);
      abs_delta_D[i] = a;
      s_gt[i]        = (dD > 0.0) ? 1 : ((dD < 0.0) ? -1 : 0);
      const double pgt = P_gt[i];
      nat_weight[i]  = NumericVector::is_na(pgt) ? NA_REAL : pgt * a;
      keep[i]        = true;
    } else {
      delta_D[i]     = NA_REAL;
      abs_delta_D[i] = NA_REAL;
      s_gt[i]        = 0;
      nat_weight[i]  = NA_REAL;
      keep[i]        = false;
    }
  }

  return List::create(
    _["delta_D"]     = delta_D,
    _["s_gt"]        = s_gt,
    _["abs_delta_D"] = abs_delta_D,
    _["nat_weight"]  = nat_weight,
    _["keep"]        = keep
  );
}
