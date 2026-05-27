#' Internal function for normalizing a (g, t)-varying variable.
#'
#' Replaces the values of `varname` with their `(G, T)` cell mean whenever the
#' variable shows any within-cell variation.
#'
#' @param df A data frame or data.table.
#' @param varname Name of the column to normalize.
#' @returns A list with `retcode` (TRUE if normalization was performed) and
#'   `df` (the possibly updated data.table).
#' @importFrom data.table as.data.table is.data.table
#' @noRd
twowayfeweights_normalize_var <- function(df, varname) {

  if (!data.table::is.data.table(df)) {
    df <- data.table::as.data.table(df)
  }

  cell_stats <- df[, list(tmp_mean_gt = mean(get(varname), na.rm = TRUE),
                          tmp_sd_gt   = stats::sd(get(varname), na.rm = TRUE)),
                   by = c("G", "T")]

  any_within_var <- sum(cell_stats$tmp_sd_gt, na.rm = TRUE) > 0

  if (any_within_var) {
    df[cell_stats, on = c("G", "T"),
       (varname) := i.tmp_mean_gt]
  }

  list(retcode = any_within_var, df = df)
}
