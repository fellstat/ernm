# Do imports from other packages
#' @importFrom stats na.omit acf cov.wt density pnorm quantile scatter.smooth sd update
#' @importFrom graphics par
#' @importFrom network list.vertex.attributes `%v%` `%v%<-` network network.vertex.names<- as.network `%n%` is.directed
#' @importFrom rlang .data
NULL

#' Register statistics
#' @name registerDirectedStatistic
#' @return no return value
#' @aliases registerDirectedStatistic registerUndirectedStatistic
#' registerDirectedOffset
#' registerUndirectedOffset
#' @usage registerDirectedStatistic
NULL

#' Metropolis samplers
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

#' Internal symbols
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

#' DirectedNet class
#'
#' An S4 (old-style) class representing a directed network.
#'
#' @name DirectedNet-class
#' @docType class
#' @exportClass Rcpp_DirectedNet
NULL

#' UndirectedNet class
#'
#' An S4 (old-style) class representing an undirected network.
#'
#' @name UndirectedNet-class
#' @docType class
#' @exportClass Rcpp_UndirectedNet
NULL


#' ERNM model terms
#' @name ernm-terms
#' @docType methods
#' @section Statistic Descriptions:
#' \describe{
#' \item{\code{edges} (directed)  (undirected)}{
#' \emph{Edges:} This term adds one network statistic equal to the number of edges
#' (i.e. nonzero values) in the network. }
#'
#' \item{\code{ reciprocity() } (directed)}{ A count of the number of pairs of actors
#' \eqn{i} and \eqn{j} for which \eqn{(i{\rightarrow}j)}{(i,j)} and \eqn{(j{\rightarrow}i)}{(j,i)}
#' both exist.
#' }
#'
#' \item{\code{ star(k, direction="in")  }  (directed)  (undirected)}{
#' The \code{k} argument is a vector of distinct integers.
#' This term adds one network statistic to the model for each element in \code{k}.
#' The \eqn{i}th such statistic counts the number of distinct \code{k[i]}-stars in the network,
#' where a \eqn{k}-star is defined to be a node \eqn{N} and a set of \eqn{k} different nodes
#' \eqn{\{O_1, \dots, O_k\}} such that the ties \eqn{\{N, O_i\}} exist for \eqn{i=1, \dots, k}.
#' For directed networks, direction indicates whether the count is of in-stars (direction="in")
#' or out-stars (direction="out")}
#'
#' \item{\code{triangles()} (directed)  (undirected)}{
#' This term adds one statistic to the model equal to the number of triangles
#' in the network. For an undirected network, a triangle is defined to be any
#' set \eqn{\{(i,j), (j,k), (k,i)\}} of three edges. For a directed network, a
#' triangle is defined as any set of three edges \eqn{(i{\rightarrow}j)}{(i,j)}
#' and \eqn{(j{\rightarrow}k)}{(j,k)} and either \eqn{(k{\rightarrow}i)}{(k,i)}
#' or \eqn{(k{\leftarrow}i)}{(i,k)}. }
#'
#'
#'
#' \item{\code{transitivity()}  (undirected)}{
#' The Soffer-Vazquez   transitivity. This is clustering metric that adjusts for large degree
#' differences and is described by C in Equation 6 of #' \url{https://pubmed.ncbi.nlm.nih.gov/16089694/}. Note
#' The approximation of the number of possible shared neighbors between node i and j of min(d_i,d_j) - 1
#' in this implementation.
#'  }
#'
#' \item{\code{ nodeMatch(name) } (directed)  (undirected)}{
#'  For categorical network nodal variable 'name,' the number of edges between nodes with the same
#'  variable value.
#'  }
#' \item{\code{ nodeMix(name) } (directed)  (undirected)}{
#' For categorical network nodal variable 'name,' adds one statistic for each combination of levels of the
#' variable equal to the count of edges between those levels.
#'   }
#' \item{\code{ homophily(name) } (directed)  (undirected)}{
#' A degeneracy robust homophily term for use in untapered models. See Fellows (2012).
#'   }
#'
#' \item{\code{ degree(d, direction="undirected", lessThanOrEqual=FALSE) } (directed)
#' (undirected)}{
#' The \code{d} argument is a vector of distinct integers. This term adds one
#' network statistic to the model for each element in \code{d}; the \eqn{i}th
#' such statistic equals the number of nodes in the network of degree
#' \code{d[i]}, i.e. with exactly \code{d[i]} edges.
#' For directed networks if direction="undirected"
#' degree is counted as the sum of the in and out degrees of a node. If direction="in" then in-degrees are
#' used and direction="out" indicates out-degrees.
#'
#' If lessThanOrEqual=TRUE, then the count is the number of nodes with degree less than or equal to d.
#'   }
#'
#' \item{\code{ nodeCov(name, direction="undirected") } (directed)  (undirected)}{
#' The \code{name} argument is a character string giving the name of a
#' numeric attribute in the network's vertex attribute list.
#' This term adds a single network statistic to the model equaling the sum of
#' \code{name(i)} and \code{name(j)} for all edges \eqn{(i,j)} in the
#' network. For categorical variables, levels are coded as 1,..,nlevels`.
#' If direction="in", only in-edges are counted. If direction="out" only
#' out-edges are counted.
#'   }
#'
#' \item{\code{gwesp(alpha)}   (directed)  (undirected)}{
#' This term is just like \code{gwdsp} except it adds a statistic equal to the
#' geometrically weighted \emph{edgewise} (not dyadwise) shared partner
#' distribution with decay parameter
#' \code{alpha} parameter, which should be non-negative.
#'   }
#'
#' \item{\code{ gwdegree(alpha, direction="undirected") }   (directed)  (undirected)}{
#' This term adds one network statistic to the model equal to the weighted
#' degree distribution with decay controlled by the \code{decay} parameter.
#' The \code{alpha} parameter is the same as theta_s in equation (14) in Hunter (2007).
#'
#' For directed networks if direction="undirected" degree is counted as the sum of the in and
#' out degrees of a node. If direction="in" then in-degrees are used ans direction="out"
#' indicates out-degrees.
#'
#'  }
#'
#' \item{\code{ gwdsp(alpha) }  (directed)  (undirected)}{
#'
#'  This term adds one network statistic to the model equal to the geometrically
#'  weighted dyadwise shared partner distribution with decay parameter
#'  \code{decay} parameter, which should be non-negative.
#'
#'  }
#'
#' \item{\code{ esp(d, type=2) } (directed)  (undirected)}{
#'
#' This term adds one network
#' statistic to the model for each element in \code{d} where the \eqn{i}th such
#' statistic equals the number of \emph{edges} (rather than dyads) in the
#' network with exactly \code{d[i]} shared partners. This term can be used with
#' directed and undirected networks. For directed networks the count depends on type:
#'
#' type = 1     :   from -> to -> nbr -> from
#'
#' type = 2     :   from -> to <- nbr <- from (homogeneous)
#'
#' type = 3     :   either type 1 or 2
#'
#' type = 4     :   all combinations of from -> to <-> nbr <-> from
#'
#'  }
#'
#' \item{\code{ geoDist(long, lat, distCuts=Inf) } (undirected)}{
#'
#' given nodal variables for longitude and latitude, calculates the sum of the
#' great circle distance between connected nodes. distCuts splits this into
#' separate statistics that count the sum of the minimum of the cut point and the
#' distance.
#'
#' }
#'
#'
#'
#'
#' }
#' @references
#' Fellows, Ian Edward. Exponential family random network models. University of California, Los Angeles, 2012.
NULL

# TODO: replace "@examplesIf FALSE" with something that will run everywhere except CRAN
# TODO: Add a reference to this documentation in every @param statement that takes an ERNM formula
# TODO: Check title formatting. are things like "Creates an ERNM likelihood model with missing data" allowed, or does it need to be "Creates an ERNM Likelihood Model with Missing Data"

#' ERNM formula
#' @name ernm-formula
#' @docType methods
#' @section Formula:
#'
#' ERNM models are specified using a formula interface. The formula expresses what statistics are in the model,
#' along with any parameters that the statistics take. It also defines what nodal covariates are considered random
#' and if the edges are considered random. Finally, it specifies the network that the model is defined on.
#'
#' The format of the formula follows the pattern
#'
#' ..network.. ~ ..statistics.. | ..random..
#'
#' The ..network.. place should be a single object that can be coerced into a native ernm network using as.BinaryNet.
#'
#' The ..statistics.. place may contain any statistic supported by the package (see \code{\link{ernm-terms}}). Multiple statistics may be added
#' by separating them with "+".
#'
#' The ..random.. place defines what is considered random within the network. By default the dyads are considered random,
#' however, this may be changed by including noDyad in the random place. Additionally, the names of vertex attributes
#' can be added to make them stochastic.
#'
#' For example, a simple erdos-renyi random graph model can be specified on the samplike dataset (see data(samplike)) by
#'
#' samplike ~ edges()
#'
#' A simple multinomial ALAAM model for the group variable can be specified by setting the dyads to be fixed and group to be random.
#'
#' samplike ~ nodeCount("group") | noDyad + group
#'
#' A full ernm model with terms for edges, triangles, and homopily looks like
#'
#' samplike ~ edges() + triangles() + nodeCount("group") + homophily("group")
#'
#'
#'
#'
NULL




