#' Register setter methods for Rcpp net objeccts
#' @name register_rcpp_net_methods
#' @return no return value
register_rcpp_net_methods <- function() {
  if (!methods::isClass("Rcpp_DirectedNet")) return()
  if (!methods::isClass("Rcpp_UndirectedNet")) return()
  
  setMethod("[", c("Rcpp_DirectedNet"),
            function(x, i, j, ..., maskMissing=TRUE, drop=TRUE)
            {
              x$`[`(i,j,maskMissing)
            })
  
  setMethod("[", c("Rcpp_UndirectedNet"),
            function(x, i, j, ..., maskMissing=TRUE, drop=TRUE)
            {
              x$`[`(i,j,maskMissing)
            })
  
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
}
