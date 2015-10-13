


#' An ernm plug-in for easy C++ prototyping and access
#' @param ... plug-in arguments
#' @examples
#' library(inline)
#' registerPlugin("ernm",inlineErnmPlugin)
#' src <- "
#'		Rcpp::IntegerMatrix tmp(0,2);
#'		ernm::BinaryNet<ernm::Directed> net(tmp,Rcpp::as<int>(n));
#'		return net;
#'		"
#' emptyNetwork <- cxxfunction(signature(n="integer"), src, plugin="ernm")
#' net <- emptyNetwork(10)
#' net[1:10,1:10]
inlineErnmPlugin <- Rcpp:::Rcpp.plugin.maker(
		include.before = "#include <ernm.h>", 
		libs           = "", 
		package        = "ernm"
)
