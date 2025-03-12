
#' MCMC Standard Error by Batch
#'
#' Computes the MCMC standard error from a statistic vector using a batching method.
#'
#' @param x A numeric vector of statistics.
#' @param expon A numeric value controlling the batch size; default is 0.5.
#' @return A numeric value representing the estimated standard error.
#' @export
mcmcse <- function(x, expon=.5){
	
	n <- length(x)
	b <- floor(n^expon)
	a <- floor(n/b)
	m <- mean(x)
	err <- rep(NA,a)
	for(i in 1:a){
		err[i] <- ( m - mean(x[((i-1)*b+1):(i*b)]) )^2
	}
	sqrt( (b/(a-1))*sum(err) ) / sqrt(n)
}


#' MCMC Effective Sample Size
#'
#' Computes the effective sample size from a statistic vector.
#'
#' @param x A numeric vector.
#' @return A numeric value representing the effective sample size.
#' @references
#' Kass, R. E., Carlin, B. P., Gelman, A., & Neal, R. M. (1998). 
#' "Markov Chain Monte Carlo in Practice: A Roundtable Discussion."
#' *The American Statistician*, 52(2), 93-100. 
#' DOI: \doi{10.2307/2685466}
#' @export
mcmcEss <- function(x){
	rho <- acf(x,plot=FALSE)$acf[2]
	length(x) * (1-rho) / (1+rho)
}



#' Create an ERNM Package Skeleton
#'
#' Creates a skeleton for a package extending the ernm package by copying an example package.
#'
#' @param path A character string specifying the directory where the package skeleton will be created.
#' @return A logical value indicating whether the copy was successful.
#' @export
ernmPackageSkeleton <- function( path = "."){
	pkgPath <- find.package("ernm")
	p <- file.path(pkgPath,"examplePackage","ErnmExtension")
	file.copy(p,path,recursive=TRUE)
}

#' Used to indicate a required parameter
#' @noRd
.required <- function(){
	r <- NA
	class(r) <- ".requiredParam"
	r
}

#' Parses evaulated parameters as if in a function call.
#' using names and positional matching
#' @param lis a list of parameter values
#' @param params a named list of default parameters
#' @importFrom stats na.omit
#' @noRd
.matchParams <- function(lis,params){
  n <- length(lis)
  nm <- names(lis)
  np <- names(params)
  if(is.null(nm)){
    params[1:n] <- lis
    names(params) <- np
    return(params)
  }
  result <- params
  mch <- pmatch(nm,np)
  result[na.omit(mch)] <- lis[!is.na(mch)]
  j <- 1
  for(i in 1:n){
    if(!is.na(mch)[i]){
      next
    }
    while(j %in% na.omit(mch))
      j <- j + 1
    if(nm[i]=="")
      result[[j]] <- lis[[i]]
    else
      stop(paste("unknown named parameter ",nm[i]))
    j <- j + 1
  }
  for(i in 1:length(result)){
  	if(inherits(result[[i]],".requiredParam"))
  		stop(paste("parameter",names(result)[i],"required but not present"))
  }
  result
}