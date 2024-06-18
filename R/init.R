#' @useDynLib GFusionMed
#' @importFrom Rcpp sourceCpp
#' @importFrom rms matinv
#' @importFrom Matrix bdiag
#' @import RcppArmadillo
#' @import RcppEigen
NULL

.onLoad <- function(libname, pkgname) {
  suppressPackageStartupMessages({
    if (!requireNamespace("RcppArmadillo", quietly = TRUE)) {
      stop("RcppArmadillo package is required but not installed.")
    }
    if (!requireNamespace("RcppEigen", quietly = TRUE)) {
      stop("RcppEigen package is required but not installed.")
    }
  })
}