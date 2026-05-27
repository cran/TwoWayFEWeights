# Silence R CMD check NOTEs for data.table non-standard evaluation.
utils::globalVariables(c(
  ":=", ".", ".I", ".SD",
  "..columns", "..original_names", "..random_weights",
  "D", "D0", "G", "P_gt", "P_k", "S_k", "T", "TFactorNum",
  "T_k", "Tfactor", "W", "Wsq", "eps_2", "i.tmp_mean_gt",
  "indicator", "indicator_l", "nat_weight", "om_tilde_1",
  "sens_measure2", "tmp_mean_gt", "tmp_sd_gt",
  "w_tilde_2_E_D_gt", "weight_result", "weights"
))
