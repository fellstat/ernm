

#' convert and UndirectedNet to a network object
#' @param x the object
#' @param ... unused
#' @return a undirected network object
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

#' convert and Rcpp_DirectedNet to a network object
#' @param x the object
#' @param ... unused
#' @return a directed network object
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

# incase binary net from cpp gets passed
#' @title Conversion to network object
#' @name as.network.UndirectedNet
#' @param x the object
#' @param ... unused
#' @return a undirected network object
#' @export
as.network.UndirectedNet <- function(x, ...) {
  as.network.Rcpp_UndirectedNet(x, ...)
}

# incase binary net from cpp gets passed
#' @title Conversion to network object
#' @name as.network.DirectedNet
#' @param x the object
#' @param ... unused
#' @return a directed network object
#' @export
as.network.DirectedNet <- function(x, ...) {
  as.network.Rcpp_DirectedNet(x, ...)
}

#' plot an Rcpp_UnirectedNet object
#' @param x the object
#' @param ... additional parameters for plot.network
#' @return No return value, invisibly NULL
#' @method plot Rcpp_DirectedNet
#' @export
plot.Rcpp_DirectedNet <- function(x,...){
	x <- as.network(x)
	plot(x,...)
}

#' plot an UndirectedNet object
#' @param x the object
#' @param ... additional parameters for plot.network
#' @return No return value, invisibly NULL
#' @method plot Rcpp_UndirectedNet
#' @export
plot.Rcpp_UndirectedNet <- function(x,...){
	x <- as.network(x)
	plot(x,...)
}

#' convert and network to either an UndirectedNet or DirectedNet object
#' @param x the object
#' @param ... unused
#' @export
#' @return a BinaryNet object
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

#' Subsetting and assignment for Net objects
#'
#' These methods allow standard subsetting (`[`) and assignment (`[<-`) for `DirectedNet` and `UndirectedNet` objects.
#'
#' @param x A `DirectedNet` or `UndirectedNet` object.
#' @param i,j Index vectors.
#' @param ... Currently unused.
#' @param maskMissing Logical. Should missing values be masked by NA?
#' @param drop Ignored (present for compatibility).
#' @param value Values to assign (for `[<-` only).
#'
#' @return A modified object or extracted submatrix depending on the method.
#' @name extract-methods
#' @rdname extract-methods
#' @docType methods
#' @aliases [,DirectedNet-method [,UndirectedNet-method [<-,DirectedNet-method [<-,UndirectedNet-method
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

