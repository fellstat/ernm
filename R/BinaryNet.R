

#' Convert an Rcpp_UndirectedNet to a network object
#'
#' Converts a native ernm network into a network object
#' from the network package.
#'
#' @param x the object
#' @param ... unused
#' @return an undirected network object
#' @examples
#' edge_list <- matrix(c(1,3),ncol=2)
#' # create a network with an edge from 1 -> 3
#' ernm_net <- new(UndirectedNet,edge_list,5)
#'
#' # convert to a network object from the network package
#' network_net <- as.network(net)
#' network_net
#' @export
as.network.Rcpp_UndirectedNet <- function(x,...){
	el <- x$edges()
	attr(el,"n") <- n <- x$size()
	nw <- network(el,directed=FALSE)

	for(i in which(x$nMissing(1:n)>0)){
		nas <- which(is.na(x[i,1:n]))
		nw[i,nas] <- NA
	}

	vn <- x$variableNames(TRUE)
	if(length(vn)>0){
		for(i in 1:length(vn)){
			vals <- x[[vn[i]]]
			if(vn[i] == "vertex.names"){
				network.vertex.names(nw) <- as.character(vals)
			}else if(vn[i] == "na"){

			}else{
				nw %v% vn[i] <- if(is.factor(vals)) as.character(vals) else as.vector(vals)
			}
		}
	}
	nw
}

#' Convert an Rcpp_DirectedNet to a network object
#'
#' Converts a native ernm network into a network object
#' from the network package.
#'
#' @param x the object
#' @param ... unused
#' @return a directed network object
#' @examples
#' edge_list <- matrix(c(1,3),ncol=2)
#' # create a network with an edge from 1 -> 3
#' ernm_net <- new(DirectedNet,edge_list,5)
#'
#' # convert to a network object from the network package
#' network_net <- as.network(net)
#' network_net
#' @export
as.network.Rcpp_DirectedNet <- function(x,...){
	el <- x$edges()
	attr(el,"n") <- n <- x$size()
	nw <- network(el,directed=TRUE)

	for(i in which(x$nMissing(1:n)>0)){
		nas <- which(is.na(x[i,1:n]))
		nw[i,nas] <- NA
	}

	vn <- x$variableNames(TRUE)
	if(length(vn)>0){
		for(i in 1:length(vn)){
			vals <- x[[vn[i]]]
			if(vn[i] == "vertex.names"){
				network.vertex.names(nw) <- as.character(vals)
			}else if(vn[i] == "na"){

			}else{
				nw %v% vn[i] <- if(is.factor(vals)) as.character(vals) else as.vector(vals)
			}
		}
	}
	nw
}
#' Plot an Rcpp_DirectedNet object
#' @param x the object
#' @param ... additional parameters for plot.network
#' @return No return value, invisibly NULL
#' @method plot Rcpp_DirectedNet
#' @examples
#'
#' # create ring network with 5 vertices
#' edge_list <- matrix(c(1,2,2,3,3,4,4,5,5,1),ncol=2, byrow=TRUE)
#' ernm_net <- new(DirectedNet,edge_list,5)
#'
#' # basic plot
#' plot(ernm_net)
#'
#' # change vertex point size (see plot.network)
#' plot(ernm_net, vertex.cex=.5)
#'
#' @export
plot.Rcpp_DirectedNet <- function(x,...){
	x <- as.network(x)
	plot(x,...)
}

#' Plot an UndirectedNet object
#' @param x the object
#' @param ... additional parameters for plot.network
#' @return No return value, invisibly NULL
#' @method plot Rcpp_UndirectedNet
#' @examples
#'
#' # create ring network with 5 vertices
#' edge_list <- matrix(c(1,2,2,3,3,4,4,5,5,1),ncol=2, byrow=TRUE)
#' ernm_net <- new(UndirectedNet,edge_list,5)
#'
#' # basic plot
#' plot(ernm_net)
#'
#' # change vertex point size (see plot.network)
#' plot(ernm_net, vertex.cex=.5)
#'
#' @export
plot.Rcpp_UndirectedNet <- function(x,...){
	x <- as.network(x)
	plot(x,...)
}

#' Convert a network to either an Rcpp_UndirectedNet or Rcpp_DirectedNet object
#'
#'
#'
#' @param x the object
#' @param ... unused
#' @export
#' @return either an Rcpp_UndirectedNet or Rcpp_DirectedNet object
#' @examples
#'
#' data(samplike)
#'
#' # convert Sampson's monks into a native ernm network
#' net <- as.BinaryNet(samplike)
#' net[["group"]]
#' net[1:5,1:5]
#'
as.BinaryNet <- function(x,...){
	if(inherits(x,c("UndirectedNet","Rcpp_UndirectedNet")))
		return(x)
	if(inherits(x,c("DirectedNet","Rcpp_DirectedNet")))
		return(x)
	if(!inherits(x,"network"))
		stop("x must be a BinaryNet or network object")
	directed <- network::is.directed(x)
	el <- as.matrix(x,matrix.type="edgelist")
	n <- attr(el,"n")
	if(directed)
		net <- new(DirectedNet,el,n)
	else
		net <- new(UndirectedNet,el,n)
	vn <- list.vertex.attributes(x)
	for(v in vn)
		net[[v]] <- x %v% v

	net
}

#' Subsetting and assignment for ernm network objects
#'
#' These methods allow standard subsetting (`[`) and assignment (`[<-`) for `Rcpp_DirectedNet` and `Rcpp_UndirectedNet` objects.
#'
#' @param x an `Rcpp_DirectedNet` or `Rcpp_UndirectedNet` object.
#' @param i,j index vectors.
#' @param ... currently unused.
#' @param maskMissing Logical. Should missing values be masked by NA?
#' @param drop Ignored (present for compatibility).
#' @param value Values to assign (for `[<-` only).
#'
#' @return A modified object or extracted submatrix depending on the method.
#' @name extract-methods
#' @rdname extract-methods
#' @docType methods
#' @aliases [,DirectedNet-method [,UndirectedNet-method [<-,DirectedNet-method [<-,UndirectedNet-method
#' @examples
#'
#' # convert the Sampson's monks network into a native ernm network
#' data(samplike)
#' sampnet <- as.BinaryNet(samplike)
#' sampnet
#'
#' # get the number of nodes and edges in the network
#' sampnet$size()
#' sampnet$nEdges()
#'
#' # Extract and assign vertex attributes
#' sampnet[["group"]]
#' sampnet[["newvar"]] <- rnorm(18)
#' sampnet[["newvar"]]
#'
#' # get the edge matrix between the first 5 vertices
#' sampnet[1:5,1:5]
#'
#' # add an edge 2 --> 3
#' sampnet[2,3] <- TRUE
#'
#' # Make the dyad 4 --> 1 missing
#' sampnet[4,1] <- NA
#' sampnet[1:5,1:5]
#'
#' # get the in and out degrees for each vertex
#' sampnet$inDegree(1:18)
#' sampnet$outDegree(1:18)
#'
NULL

#' indexing
#' @rdname extract-methods
setMethod("[", c("Rcpp_DirectedNet"),
		function(x, i, j, ..., maskMissing=TRUE, drop=TRUE)
		{
			x$`[`(i,j,maskMissing)
		})

#' indexing
#' @rdname extract-methods
setMethod("[", c("Rcpp_UndirectedNet"),
		function(x, i, j, ..., maskMissing=TRUE, drop=TRUE)
		{
			x$`[`(i,j,maskMissing)
		})

#' indexing
#' @rdname extract-methods
setMethod("[<-", c("Rcpp_DirectedNet"),
		function(x, i, j, ..., value)
		{
			if(is.vector(value)){
				if(length(value)==length(i) && length(j)==1)
					value <- as.matrix(as.logical(value))
				else if(length(value)==length(j) && length(i)==1)
					value <- t(as.matrix(as.logical(value)))
				else
					stop("invalid assignment")
			}
			x$`[<-`(i,j,value)
			x
		})

#' indexing
#' @rdname extract-methods
setMethod("[<-", c("Rcpp_UndirectedNet"),
		function(x, i, j, ..., value)
		{
			if(is.vector(value)){
				if(length(value)==length(i) && length(j)==1)
					value <- as.matrix(as.logical(value))
				else if(length(value)==length(j) && length(i)==1)
					value <- t(as.matrix(as.logical(value)))
				else
					stop("invalid assignment")
			}
			x$`[<-`(i,j,value)
			x
		})

