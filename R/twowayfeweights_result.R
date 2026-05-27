#' Internal workhorse function for creating the return object of a
#' `twowayfeweights()` call.
#'
#' @param dat A data.table, as per the return object from
#'   `twowayfeweights_calculate()`.
#' @param beta Coefficient value of the treatment variable ("D"), again as per
#'   the return object of `twowayfeweights_calculate()`.
#' @param random_weights A vector indicating the column names of random weights.
#' @param treatments A vector indicating the column names of other treatments.
#' @returns A list.
#' @details This function is normally run directly after
#'   `twowayfeweights_calculate()`.
#' @importFrom data.table as.data.table is.data.table setnames setorder shift
#' @importFrom stats weighted.mean
#' @noRd
twowayfeweights_result <- function(dat, beta, random_weights, treatments = NULL) {

  if (!data.table::is.data.table(dat)) {
    dat <- data.table::as.data.table(dat)
  }

  limit_sensitivity <- 1e-10

  zero_below_eps <- function(x) {
    data.table::fifelse(!is.na(x) & abs(x) < limit_sensitivity, 0, x)
  }

  dat[, weight_result := zero_below_eps(weight_result)]

  if (is.null(treatments)) {

    ret <- twowayfeweights_summarize_weights(dat, "weight_result")

    W_mean <- stats::weighted.mean(dat$W, dat$nat_weight, na.rm = TRUE)
    M <- sum(dat$nat_weight != 0, na.rm = TRUE)
    W_sd <- sqrt(sum(dat$nat_weight * (dat$W - W_mean)^2, na.rm = TRUE)) *
            sqrt(M / (M - 1))
    sensibility <- abs(beta) / W_sd

    dat_result <- dat[, list(T, G, weight = weight_result)]

    ret$dat_result  <- dat_result
    ret$beta        <- beta
    ret$sensibility <- sensibility

    if (length(random_weights) > 0) {
      ret$mat <- twowayfeweights_test_random_weights(dat, random_weights)
    }

    if (ret$sum_minus < 0) {
      dat_sens <- dat[weight_result != 0]
      data.table::setorder(dat_sens, W, -G, -T)
      dat_sens[, Wsq := nat_weight * W^2]
      dat_sens[, P_k := cumsum(nat_weight)]
      dat_sens[, S_k := cumsum(weight_result)]
      dat_sens[, T_k := cumsum(Wsq)]
      data.table::setorder(dat_sens, -W, G, T)

      N <- nrow(dat_sens)
      dat_sens[, sens_measure2 := abs(beta) /
                  sqrt(T_k + S_k^2 / (1 - P_k))]
      dat_sens[, indicator := as.numeric(W < -S_k / (1 - P_k))]
      dat_sens[1L, indicator := 0]
      dat_sens[, indicator_l := data.table::shift(indicator, n = 1L,
                                                  type = "lag", fill = -1)]
      dat_sens[, indicator := pmax(indicator, indicator_l)]
      total_indicator <- sum(dat_sens$indicator)
      ret$sensibility2 <- dat_sens$sens_measure2[N - total_indicator + 1L]
    }

    ret$tot_cells <- sum(dat$nat_weight != 0, na.rm = TRUE)

  } else {

    for (v in c("result", treatments)) {
      colname <- paste0("weight_", v)
      data.table::set(dat, j = colname,
                      value = zero_below_eps(dat[[colname]]))
    }

    columns <- c("T", "G", "weight_result")
    ret <- twowayfeweights_summarize_weights(dat, "weight_result")
    ret$tot_cells <- sum(dat$nat_weight != 0, na.rm = TRUE)

    if (length(random_weights) > 0) {
      ret$mat <- twowayfeweights_test_random_weights(dat, random_weights)
    }

    for (treatment in treatments) {
      varname <- fn_treatment_weight_rename(treatment)
      columns <- c(columns, varname)
      ret2 <- twowayfeweights_summarize_weights(dat, varname)
      ret[[treatment]] <- ret2
      ret[[treatment]]$tot_cells <- sum(dat[[treatment]] != 0, na.rm = TRUE)
    }

    dat_result <- dat[, ..columns]
    data.table::setnames(dat_result, "weight_result", "weight")

    ret$beta       <- beta
    ret$dat_result <- dat_result
  }

  ret
}
