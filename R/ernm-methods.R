#' Print ernm object
#' @param x x
#' @param ... unused
#' @return No return value, prints summary
#' @export
#' @method print ernm
print.ernm <- function(x,...){
  cat("                          ",x$m$name(),"\n")
  cat("Domain:\n Random graph =",x$m$randomGraph())
  rv <- x$m$randomVariables()
  if(length(rv)>0)
    cat(", Random variables =", rv)
  cat("\n\n")
  r <- rbind()
  #cat("Parameters:\n")
  #print(x$theta)
  samp <- x$sample
  if(is.list(samp))
    samp <- samp[[1]]
  #cat("Mean values:\n")
  #print(colMeans(samp))
  r <- rbind(x$theta,colMeans(samp))
  rownames(r) <- c("Parameters","Mean Values")
  print(r)
  offset <- attr(samp,"offset")
  if(!is.null(offset) && ncol(offset)>0){
    cat("\nOffset Mean Values:\n")
    print(colMeans(offset))
  }
}

# TODO: I've removed AIC/BIC from summary. summary should not have anything in it that takes significant time. We need to change the paper to reflect the different output.
#' Summary for ernm object
#' @param object object
#' @param ... unused
#' @return a data frame summary of the model
#' @export
#' @method summary ernm
summary.ernm <- function(object, ...){

  theta <- object$theta
  cv <- solve(object$info)
  se <- sqrt(diag(cv))
  z <- theta/se
  p.value <- 2*pnorm(abs(z),lower.tail=FALSE)
  d <- data.frame(theta,se,z,p.value)
  rownames(d) <- make.unique(names(theta))
  d
}



#' MCMC approximate log-likelihood
#'
#' @param object an ernm object
#' @param ... unused
#' @export
#' @method logLik ernm
#' @examplesIf FALSE
#' fit <- ernm(samplike ~ edges() + nodeCount("group") + nodeMatch("group") | group)
#' logLik(fit)
#' AIC(fit)
#' BIC(fit)
#'
logLik.ernm <- function(object, ...){
  theta <- object$theta
  # TODO: why are we generating these samples? Don't we get them for free in object$sample?
  # If they are needed. MCMC options should be passed as parameters.
  if(!is.null(object$m$missSamp)){
    n_sim <- dim(object$sample$unconditional)[1]
    samples <- object$m$generateSampleStatistics(10000,100,n_sim*10)$unconditional
  }else{
    n_sim <- dim(object$sample)[1]
    samples <- object$m$generateSampleStatistics(10000,100,n_sim*10)
  }

  sample_calc <- apply(samples,1,function(x){sum(theta*x)})
  max_term <- max(sample_calc)
  const_approx <- log(mean(exp(sample_calc - max_term))) + max_term
  logLik <- sum(theta*ernm::calculateStatistics(object$formula)) - const_approx # TODO: I don't think this is correct for a fit with missing data.
  net <- as.BinaryNet(eval(object$formula[[2]],envir=environment(object$formula)))
  n_verts <- net$size()
  n_dyads <- n_verts*(n_verts-1)*(1 - 0.5*(!net$isDirected())) #
  attr(logLik, "df") <- length(theta)
  attr(logLik, "nobs") <- n_verts*(length(object$m$randomVariables())!=0) + n_dyads*object$m$randomGraph()
  logLik
}


#' Parameter covariance matrix
#' @param object object
#' @param ... unused
#' @return covariance matrix
#' @export
#' @method vcov ernm
vcov.ernm <- function(object,...){
  solve(object$info)
}

#' Access ERNM parameters
#' @param object object
#' @param ... unused
#' @return parameter vector
#' @export
#' @method coef ernm
coef.ernm <- function(object,...){
  object$theta
}

#' Plot an ernm object
#' @param  x the object
#' @param ... unused
#' @return No return value, plots the likelihood history
#' @export
#' @method plot ernm
plot.ernm <- function(x,...){
  plot(x$likelihoodHistory,main="Likelihood convergence",
       ylab="Change in log-likelihood",xlab="iteration")
}

