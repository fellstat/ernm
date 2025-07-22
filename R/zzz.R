#' @useDynLib ernm, .registration = TRUE
#' @exportPattern ^[^\.]
#' @import methods
#' @import Rcpp
#' @import BH
#' @import trust
#' @import moments
NULL

.onLoad <- function(libname, pkgname){
    # Load the Rcpp module "ernm"
    Rcpp::loadModule("ernm", TRUE)
    # Skip initialization when running devtools::document()
    if (nzchar(Sys.getenv("ROXYGEN_PKG"))) {
        return()
    }
    .Call("_ernm_initStats")
    .Call("_ernm_initToggles")
    register_rcpp_net_methods()
}
