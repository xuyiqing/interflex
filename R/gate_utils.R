#' @title Standardize GATE Column Names
#' @description Internal helper to rename bootstrap GATE columns to the
#'   unified format used by the plotting functions.
#' @keywords internal
.standardize_gate_columns <- function(nms, effect_col) {
    nms[nms == effect_col]         <- "ME"
    nms[nms == "SE"]               <- "sd"
    nms[nms == "CI.lower"]         <- "lower CI(95%)"
    nms[nms == "CI.upper"]         <- "upper CI(95%)"
    nms[nms == "CI.lower.uniform"] <- "lower uniform CI(95%)"
    nms[nms == "CI.upper.uniform"] <- "upper uniform CI(95%)"
    nms
}
