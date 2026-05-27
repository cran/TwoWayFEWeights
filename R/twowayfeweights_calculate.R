#' Internal workhorse function for calculating the twoway FE weights
#' 
#' @param dat A data frame.
#' @param type Desired type of calculation.
#' @param controls A vector of controls.
#' @param treatments Additional treatments.
#' @returns A list.
#' @importFrom fixest feols
#' @importFrom data.table data.table setorderv set fifelse
#' @importFrom stats as.formula resid
#' @noRd
# =============================================================================
# twowayfeweights_calculate.R
# Optimized core calculation using data.table + Rcpp C++ kernels
# Replaces: twowayfeweights_calculate.R
# =============================================================================

twowayfeweights_calculate <- function(dt, type = c("feTR", "fdTR", "feS", "fdS"),
                                           controls = NULL, treatments = NULL) {
  
  type <- match.arg(type)
  if (!is.null(treatments) && type != "feTR") {
    stop("When `other_treatments` is specified, you need `type = 'feTR'`.")
  }
  
  type_TR <- type %in% c("feTR", "fdTR")
  type_fe <- type %in% c("feTR", "feS")
  
  # --- Weighted mean of D (or D0 for fdTR) ---
  if (type_TR) {
    DVAR <- if (type == "feTR") "D" else "D0"
    mean_D <- cpp_weighted_mean(dt[[DVAR]], dt$weights)
  }
  
  # --- P_gt: sum of weights by (G, T), then normalize ---
  obs <- sum(dt$weights)
  dt[, P_gt := sum(weights), by = .(G, T)]
  dt[, P_gt := P_gt / obs]
  
  if (type_TR) {
    dt[, nat_weight := P_gt * get(DVAR) / mean_D]
  }
  
  # --- Denominator regression (feols) ---
  if (is.null(controls)) controls_formula <- "1" else controls_formula <- controls
  fes <- "Tfactor"
  if (type_fe) fes <- c("G", fes)
  
  xvars <- c(controls_formula, treatments)
  fml_str <- paste0("D ~ ", paste(xvars, collapse = " + "),
                    " | ", paste(fes, collapse = " + "))
  fml <- as.formula(fml_str)
  
  if (type == "fdS") {
    sub_idx <- dt$weights != 0
    denom.lm <- fixest::feols(fml, data = dt[sub_idx], weights = dt$weights[sub_idx])
  } else {
    denom.lm <- fixest::feols(fml, data = dt, weights = dt$weights)
  }
  
  EPS_VAR <- if (type_fe) "eps_1" else "eps_2"
  
  if (type_fe || type == "fdS") {
    data.table::set(dt, j = EPS_VAR, value = stats::resid(denom.lm))
  } else if (type == "fdTR") {
    eps_vals <- stats::resid(denom.lm, na.rm = FALSE)
    eps_vals[is.na(eps_vals)] <- 0.0
    data.table::set(dt, j = EPS_VAR, value = eps_vals)
  }
  
  # --- Beta regression ---
  xvars_beta <- c("D", xvars)
  fml_beta_str <- paste0("Y ~ ", paste(xvars_beta, collapse = " + "),
                         " | ", paste(fes, collapse = " + "))
  fml_beta <- as.formula(fml_beta_str)
  
  if (type == "fdS") {
    sub_idx <- dt$weights != 0
    beta.lm <- fixest::feols(fml_beta, data = dt[sub_idx],
                             weights = dt$weights[sub_idx], only.coef = TRUE)
  } else {
    beta.lm <- fixest::feols(fml_beta, data = dt,
                             weights = dt$weights, only.coef = TRUE)
  }
  beta <- beta.lm[["D"]]
  
  # ===========================================================================
  # Type-specific weight calculations
  # ===========================================================================
  
  if (type == "feTR") {
    
    # W and weight_result (all in-place)
    eps_vec <- dt[[EPS_VAR]]
    D_vec <- dt[[DVAR]]
    
    # Stata (lines 119-124 and 611-612 of twowayfeweights.ado) always uses
    # weighted mean of eps_1 * D [aweight=weight_XX], regardless of whether
    # other_treatments are specified.  The previous unweighted-mean branch
    # caused Sum(weights) to miss the normalization factor when weights and
    # other_treatments were combined, so the main treatment's weights no longer
    # summed to 1 and the test_random_weights coefficient was off by that same
    # factor (correlation was unaffected because it is scale-invariant).
    denom_W <- cpp_weighted_mean(eps_vec * D_vec, dt$weights)
    
    dt[, W := get(EPS_VAR) * mean_D / denom_W]
    dt[, weight_result := W * nat_weight]
    
    # Other treatments
    if (!is.null(treatments)) {
      for (treatment in treatments) {
        varname <- fn_treatment_weight_rename(treatment)
        data.table::set(dt, j = varname,
                        value = dt$W * dt$P_gt * dt[[treatment]] / mean_D)
      }
    }
    
    # Remove temp columns
    dt[, c(EPS_VAR, "P_gt") := NULL]
    
    # Keep one obs per (G, Tfactor) cell
    dt <- dt[dt[, .I[1], by = .(G, Tfactor)]$V1]
    
  } else if (type == "feS") {
    
    # --- Reverse cumsum within groups (C++ kernel) ---
    data.table::setorderv(dt, c("G", "Tfactor"))
    
    eps_w <- dt[[EPS_VAR]] * dt$weights
    g_int <- as.integer(factor(dt$G))
    
    E_eps_1_g_ge_aux <- cpp_rev_cumsum_by_group(g_int, eps_w)
    weights_aux <- cpp_rev_cumsum_by_group(g_int, dt$weights)
    E_eps_1_g_ge <- E_eps_1_g_ge_aux / weights_aux
    
    data.table::set(dt, j = "E_eps_1_g_ge", value = E_eps_1_g_ge)
    
  } else if (type == "fdTR") {
    
    # eps_2 fallback (duplicate handling)
    dt[, eps_2 := ifelse(is.na(get(EPS_VAR)), 0, get(EPS_VAR))]
    
  }
  # fdS: nothing extra needed before beta reg
  
  # ===========================================================================
  # Post-beta calculations per type
  # ===========================================================================
  
  if (type == "fdTR") {
    
    # --- C++ kernel for w_tilde_2 ---
    data.table::setorderv(dt, c("G", "TFactorNum"))
    g_int <- as.integer(factor(dt$G))
    
    w_tilde_2 <- cpp_fdtr_wtilde2(g_int, dt$TFactorNum, dt$eps_2, dt$P_gt)
    data.table::set(dt, j = "w_tilde_2", value = w_tilde_2)
    
    dt[, w_tilde_2_E_D_gt := w_tilde_2 * D0]
    
    denom_W <- cpp_weighted_mean(dt$w_tilde_2_E_D_gt, dt$P_gt)
    dt[, W := w_tilde_2 * mean_D / denom_W]
    dt[, weight_result := W * nat_weight]
    
    # Cleanup temp columns
    dt[, c("eps_2", "P_gt", "w_tilde_2", "w_tilde_2_E_D_gt") := NULL]
    
  } else if (type == "feS") {
    
    # --- C++ kernel for delta calculations ---
    data.table::setorderv(dt, c("G", "TFactorNum"))
    g_int <- as.integer(factor(dt$G))
    
    delta_res <- cpp_feS_delta(g_int, dt$TFactorNum, dt$D, dt$P_gt)
    
    # Filter: keep only rows where delta_D is not NA
    keep <- delta_res$keep
    dt <- dt[keep]
    delta_D <- delta_res$delta_D[keep]
    s_gt <- delta_res$s_gt[keep]
    abs_delta_D <- delta_res$abs_delta_D[keep]
    nat_w <- delta_res$nat_weight[keep]
    
    data.table::set(dt, j = "delta_D", value = delta_D)
    data.table::set(dt, j = "s_gt", value = s_gt)
    data.table::set(dt, j = "abs_delta_D", value = abs_delta_D)
    data.table::set(dt, j = "nat_weight", value = nat_w)
    
    P_S <- sum(nat_w, na.rm = TRUE)
    dt[, nat_weight := nat_weight / P_S]
    dt[, om_tilde_1 := s_gt * E_eps_1_g_ge / P_gt]
    
    denom_W <- cpp_weighted_mean(dt$om_tilde_1, dt$nat_weight)
    dt[, W := om_tilde_1 / denom_W]
    dt[, weight_result := W * nat_weight]
    
    # Cleanup
    dt[, c("eps_1", "P_gt", "om_tilde_1", "E_eps_1_g_ge",
           "abs_delta_D", "delta_D") := NULL]
    
  } else if (type == "fdS") {
    
    dt[, s_gt := data.table::fifelse(D > 0, 1L,
                                     data.table::fifelse(D < 0, -1L, 0L))]
    dt[, abs_delta_D := abs(D)]
    dt[, nat_weight := P_gt * abs_delta_D]
    
    P_S <- sum(dt$nat_weight)
    dt[, nat_weight := nat_weight / P_S]
    dt[, W := s_gt * eps_2]
    
    denom_W <- cpp_weighted_mean(dt$W, dt$nat_weight)
    dt[, W := W / denom_W]
    dt[, weight_result := W * nat_weight]
    
    # Cleanup
    dt[, c("eps_2", "P_gt", "abs_delta_D") := NULL]
  }
  
  return(list(dat = dt, beta = beta))
}
