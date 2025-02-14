#' The inline plug-in for ernm
#' @name inlineErnmPlugin
#' @aliases inlineErnmPlugin
#' @usage inlineErnmPlugin
NULL

#' Register Statistics
#' @name registerDirectedStatistic
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
#' 
NULL

#' Models
#' @name ErnmModels
#' @docType class
#' @aliases DirectedModel UndirectedModel 
#' Rcpp_DirectedModel-class Rcpp_UndirectedModel-class
#' Rcpp_DirectedTaperedModel-class Rcpp_UndirectedTaperedModel-class
NULL

#' BinaryNet
#' @name BinaryNet
#' @docType class
#' @aliases DirectedNet UndirectedNet Rcpp_DirectedNet-class Rcpp_UndirectedNet-class
NULL

#' DirectedMetropolisHastings Class
#'
#' A class implementing the Metropolis-Hastings algorithm for directed networks.
#' @name DirectedMetropolisHastings
#' @docType class
#' @aliases Rcpp_DirectedMetropolisHastings
#' @section Methods:
#' \describe{
#'   \item{\code{setModel(model)}}{Sets the model to be used.}
#'   \item{\code{getModel()}}{Returns the current model.}
#'   \item{\code{setDyadToggleType(type)}}{Sets the dyad toggle type.}
#'   \item{\code{generateSample()}}{Generates a network sample.}
#'   \item{\code{generateSampleStatistics()}}{Computes sample statistics.}
#' }
#'
#' @examples
#' \dontrun{
#'   dmh <- new(DirectedMetropolisHastings)
#'   dmh$setModel(myDirectedModel)
#'   sample <- dmh$generateSample()
#' }
#'
#' @export
"DirectedMetropolisHastings"

#' DirectedTaperedModel Class
#'
#' A class representing a tapered model for directed networks.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{setTau(tau)}}{Sets the tapering parameter tau.}
#'   \item{\code{tau()}}{Returns the current tau parameters.}
#'   \item{\code{setCenters(centers)}}{Sets the center parameters for tapering.}
#'   \item{\code{centers()}}{Retrieves the center parameters.}
#' }
#'
#'
#' @export
"DirectedTaperedModel"

#' UndirectedCdSampler Class
#'
#' A class implementing a CdSampler for undirected networks.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{setModel(model)}}{Sets the model to be sampled.}
#'   \item{\code{getModel()}}{Returns the current model.}
#'   \item{\code{setDyadToggleType(type)}}{Sets the dyad toggle type.}
#'   \item{\code{setVertexToggleType(type)}}{Sets the vertex toggle type.}
#'   \item{\code{setDyadProbability(prob)}}{Sets the probability for dyad toggling.}
#'   \item{\code{generateSample()}}{Generates a network sample.}
#'   \item{\code{generateSampleStatistics()}}{Computes sample statistics.}
#' }
#'
#' @export
"UndirectedCdSampler"

#' UndirectedGibbsCdSampler Class
#'
#' A class implementing a Gibbs-based CdSampler for undirected networks.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{setModel(model)}}{Sets the model to be sampled.}
#'   \item{\code{getModel()}}{Returns the current model.}
#'   \item{\code{setDyadToggleType(type)}}{Sets the dyad toggle type.}
#'   \item{\code{setVertexToggleType(type)}}{Sets the vertex toggle type.}
#'   \item{\code{setDyadProbability(prob)}}{Sets the probability for dyad toggling.}
#'   \item{\code{generateSample()}}{Generates a network sample.}
#'   \item{\code{generateSampleStatistics()}}{Computes sample statistics.}
#' }
#'
#' @export
"UndirectedGibbsCdSampler"

#' UndirectedGibbsCdSampler2 Class
#'
#' A variant of the Gibbs-based CdSampler for undirected networks.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{setModel(model)}}{Sets the model to be sampled.}
#'   \item{\code{getModel()}}{Returns the current model.}
#'   \item{\code{setDyadToggleType(type)}}{Sets the dyad toggle type.}
#'   \item{\code{setVertexToggleType(type)}}{Sets the vertex toggle type.}
#'   \item{\code{setDyadProbability(prob)}}{Sets the probability for dyad toggling.}
#'   \item{\code{generateSample()}}{Generates a network sample.}
#'   \item{\code{generateSampleStatistics()}}{Computes sample statistics.}
#' }
#'
#'
#' @export
"UndirectedGibbsCdSampler2"

#' UndirectedMetropolisHastings Class
#'
#' A class implementing the Metropolis-Hastings algorithm for undirected networks.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{setModel(model)}}{Sets the model for the sampler.}
#'   \item{\code{getModel()}}{Returns the current model.}
#'   \item{\code{setDyadToggleType(type)}}{Sets the dyad toggle type.}
#'   \item{\code{setVertexToggleType(type)}}{Sets the vertex toggle type.}
#'   \item{\code{setDyadProbability(prob)}}{Sets the probability for dyad toggling.}
#'   \item{\code{generateSample()}}{Generates a network sample.}
#'   \item{\code{generateSampleStatistics()}}{Computes sample statistics.}
#' }
#'
#'
#' @export
"UndirectedMetropolisHastings"

#' UndirectedTaperedModel Class
#'
#' A class representing a tapered model for undirected networks.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{setTau(tau)}}{Sets the tapering parameter.}
#'   \item{\code{tau()}}{Retrieves the current tau parameter(s).}
#'   \item{\code{setCenters(centers)}}{Sets the center parameters.}
#'   \item{\code{centers()}}{Retrieves the center parameters.}
#' }
#'
#'
#' @export
"UndirectedTaperedModel"

#' Internal Symbols
#' @name call-symbols
#' @description Internal symbols used to access compiles code.
#' @docType methods
#' @aliases _rcpp_module_boot_ernm _ernm_initToggles _ernm_initStats initLatent
NULL

#' runErnmCppTests
#'
#' Runs the internal C++ tests for the ernm package.
#'
#' @return A logical value indicating whether all tests passed.
#'
#' @examples
#' \dontrun{
#'   runErnmCppTests()
#' }
#'
#' @export
"runErnmCppTests"



