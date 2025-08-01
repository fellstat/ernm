% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roxygen.R
\docType{methods}
\name{ernm-terms}
\alias{ernm-terms}
\title{ERNM model terms}
\description{
ERNM model terms
}
\section{Statistic Descriptions}{

\describe{
\item{\code{edges} (directed)  (undirected)}{
\emph{Edges:} This term adds one network statistic equal to the number of edges
(i.e. nonzero values) in the network. }

\item{\code{ reciprocity() } (directed)}{ A count of the number of pairs of actors
\eqn{i} and \eqn{j} for which \eqn{(i{\rightarrow}j)}{(i,j)} and \eqn{(j{\rightarrow}i)}{(j,i)}
both exist.
}

\item{\code{ star(k, direction="in")  }  (directed)  (undirected)}{
The \code{k} argument is a vector of distinct integers.
This term adds one network statistic to the model for each element in \code{k}.
The \eqn{i}th such statistic counts the number of distinct \code{k[i]}-stars in the network,
where a \eqn{k}-star is defined to be a node \eqn{N} and a set of \eqn{k} different nodes
\eqn{\{O_1, \dots, O_k\}} such that the ties \eqn{\{N, O_i\}} exist for \eqn{i=1, \dots, k}.
For directed networks, direction indicates whether the count is of in-stars (direction="in")
or out-stars (direction="out")}

\item{\code{triangles()} (directed)  (undirected)}{
This term adds one statistic to the model equal to the number of triangles
in the network. For an undirected network, a triangle is defined to be any
set \eqn{\{(i,j), (j,k), (k,i)\}} of three edges. For a directed network, a
triangle is defined as any set of three edges \eqn{(i{\rightarrow}j)}{(i,j)}
and \eqn{(j{\rightarrow}k)}{(j,k)} and either \eqn{(k{\rightarrow}i)}{(k,i)}
or \eqn{(k{\leftarrow}i)}{(i,k)}. }



\item{\code{transitivity()}  (undirected)}{
The Soffer-Vazquez   transitivity. This is clustering metric that adjusts for large degree
differences and is described by C in Equation 6 of #' \url{https://pubmed.ncbi.nlm.nih.gov/16089694/}. Note
The approximation of the number of possible shared neighbors between node i and j of min(d_i,d_j) - 1
in this implementation.
 }

\item{\code{ nodeMatch(name) } (directed)  (undirected)}{
 For categorical network nodal variable 'name,' the number of edges between nodes with the same
 variable value.
 }
\item{\code{ nodeMix(name) } (directed)  (undirected)}{
For categorical network nodal variable 'name,' adds one statistic for each combination of levels of the
variable equal to the count of edges between those levels.
  }
\item{\code{ homophily(name) } (directed)  (undirected)}{
A degeneracy robust homophily term for use in untapered models. See Fellows (2012).
  }

\item{\code{ degree(d, direction="undirected", lessThanOrEqual=FALSE) } (directed)
(undirected)}{
The \code{d} argument is a vector of distinct integers. This term adds one
network statistic to the model for each element in \code{d}; the \eqn{i}th
such statistic equals the number of nodes in the network of degree
\code{d[i]}, i.e. with exactly \code{d[i]} edges.
For directed networks if direction="undirected"
degree is counted as the sum of the in and out degrees of a node. If direction="in" then in-degrees are
used and direction="out" indicates out-degrees.

If lessThanOrEqual=TRUE, then the count is the number of nodes with degree less than or equal to d.
  }

\item{\code{ nodeCov(name, direction="undirected") } (directed)  (undirected)}{
The \code{name} argument is a character string giving the name of a
numeric attribute in the network's vertex attribute list.
This term adds a single network statistic to the model equaling the sum of
\code{name(i)} and \code{name(j)} for all edges \eqn{(i,j)} in the
network. For categorical variables, levels are coded as 1,..,nlevels`.
If direction="in", only in-edges are counted. If direction="out" only
out-edges are counted.
  }

\item{\code{gwesp(alpha)}   (directed)  (undirected)}{
This term is just like \code{gwdsp} except it adds a statistic equal to the
geometrically weighted \emph{edgewise} (not dyadwise) shared partner
distribution with decay parameter
\code{alpha} parameter, which should be non-negative.
  }

\item{\code{ gwdegree(alpha, direction="undirected") }   (directed)  (undirected)}{
This term adds one network statistic to the model equal to the weighted
degree distribution with decay controlled by the \code{decay} parameter.
The \code{alpha} parameter is the same as theta_s in equation (14) in Hunter (2007).

For directed networks if direction="undirected" degree is counted as the sum of the in and
out degrees of a node. If direction="in" then in-degrees are used ans direction="out"
indicates out-degrees.

 }

\item{\code{ gwdsp(alpha) }  (directed)  (undirected)}{

 This term adds one network statistic to the model equal to the geometrically
 weighted dyadwise shared partner distribution with decay parameter
 \code{decay} parameter, which should be non-negative.

 }

\item{\code{ esp(d, type=2) } (directed)  (undirected)}{

This term adds one network
statistic to the model for each element in \code{d} where the \eqn{i}th such
statistic equals the number of \emph{edges} (rather than dyads) in the
network with exactly \code{d[i]} shared partners. This term can be used with
directed and undirected networks. For directed networks the count depends on type:

type = 1     :   from -> to -> nbr -> from

type = 2     :   from -> to <- nbr <- from (homogeneous)

type = 3     :   either type 1 or 2

type = 4     :   all combinations of from -> to <-> nbr <-> from

 }

\item{\code{ geoDist(long, lat, distCuts=Inf) } (undirected)}{

given nodal variables for longitude and latitude, calculates the sum of the
great circle distance between connected nodes. distCuts splits this into
separate statistics that count the sum of the minimum of the cut point and the
distance.

}




}
}

\references{
Fellows, Ian Edward. Exponential family random network models. University of California, Los Angeles, 2012.
}
