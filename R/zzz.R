#' @useDynLib ernm, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import BH
#' @import trust
#' @import moments
loadModule("ernm",TRUE)
.onLoad <- function(libname, pkgname){
    # Skip initialization when running devtools::document()
    if (nzchar(Sys.getenv("ROXYGEN_PKG"))) {
        return()
    }
    .Call("_ernm_initStats")
    .Call("_ernm_initToggles")
}
.onUnload <- function(libpath) {}
