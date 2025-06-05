# Do imports from other packages 
#' @importFrom stats na.omit acf cov.wt density pnorm quantile scatter.smooth sd update
#' @importFrom graphics par
#' @importFrom network list.vertex.attributes `%v%` `%v%<-` network network.vertex.names<- as.network `%n%` is.directed
#' @importFrom rlang .data
NULL

#' Register Statistics
#' @name registerDirectedStatistic
#' @return no return value
#' @aliases registerDirectedStatistic registerUndirectedStatistic  
#' registerDirectedOffset 
#' registerUndirectedOffset
#' @usage registerDirectedStatistic
NULL

#' Metropolis Samplers
#' @name ErnmSamplers
#' @docType class
#' @aliases 
#' Rcpp_DirectedMetropolisHastings-class
#' Rcpp_UndirectedMetropolisHastings-class
#' Rcpp_UndirectedCdSampler-class
#' Rcpp_DirectedCdSampler-class
#' Rcpp_UndirectedGibbsCdSampler-class
#' Rcpp_DirectedGibbsCdSampler-class
#' Rcpp_UndirectedGibbsCdSampler2-class
#' Rcpp_DirectedGibbsCdSampler2-class
#' DirectedMetropolisHastings
#' UndirectedMetropolisHastings
#' UndirectedCdSampler
#' DirectedCdSampler
#' UndirectedGibbsCdSampler
#' DirectedGibbsCdSampler
#' UndirectedGibbsCdSampler2
#' DirectedGibbsCdSampler2
NULL

#' Models
#' @name ErnmModels
#' @docType class
#' @aliases DirectedModel UndirectedModel 
#' Rcpp_DirectedModel-class Rcpp_UndirectedModel-class
#' Rcpp_DirectedTaperedModel-class Rcpp_UndirectedTaperedModel-class
#' DirectedTaperedModel UndirectedTaperedModel
NULL

#' BinaryNet
#' @name BinaryNet
#' @docType class
#' @aliases DirectedNet UndirectedNet Rcpp_DirectedNet-class Rcpp_UndirectedNet-class
NULL

#' Internal Symbols
#' @name call-symbols
#' @description Internal symbols used to access compiles code.
#' @docType methods
#' @aliases _rcpp_module_boot_ernm _ernm_initToggles _ernm_initStats initLatent
NULL

#' runErnmCppTests
#' @name runErnmCppTests
#' @description Runs the internal C++ tests for the ernm package.
#' @return A logical value indicating whether all tests passed.
#' @examples
#' runErnmCppTests()
#'
NULL

#' DirectedNet Class
#'
#' An S4 (old-style) class representing a directed network.
#'
#' @name DirectedNet-class
#' @docType class
#' @exportClass Rcpp_DirectedNet
#' @aliases DirectedNet
NULL

#' UndirectedNet Class
#'
#' An S4 (old-style) class representing an undirected network.
#'
#' @name UndirectedNet-class
#' @docType class
#' @exportClass Rcpp_UndirectedNet
#' @aliases UndirectedNet
NULL


