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

#' Summary for ernm object
#' @param object object
#' @param ... unused
#' @return an ErnmSummary object
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
  if(!is.null(object$m$missSamp)){
    criteria <- NULL
  }else{
    criteria <- data.frame(AIC=AIC(object), BIC=BIC(object))
    rownames(criteria) <- ""
    attr(d, "criteria") <- criteria
  }
  class(d) <- c("ErnmSummary", "data.frame")
  d
}

#' Print a ERNM summary object
#' @param x the object
#' @param ... parameters passed to print.data.frame
#' @return x invisibly
#' @export
#' @method print ErnmSummary
print.ErnmSummary <- function(x, ...){
  # rename the columns:
  colnames(x) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  stats::printCoefmat(as.matrix(x),
                      digits = 3,
                      signif.stars = TRUE,
                      P.values = TRUE,
                      has.Pvalue = TRUE)
  crit <- attr(x,"criteria")
  if(!is.null(crit)){
    cat("\nInformation Criteria:\n")
    print(crit, ...)
  }
  invisible(x)
}


#' MCMC approximate log-likelihood
#'
#' @param object an ernm object
#' @param ... unused
#' @importFrom stats AIC BIC
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
  if(!is.null(object$m$missSamp)){
    stop("Missing data models are not supported yet")
  }
  samples <- object$sample

  sample_calc <- apply(samples,1,function(x){sum(theta*x)})
  max_term <- max(sample_calc)
  const_approx <- log(mean(exp(sample_calc - max_term))) + max_term
  logLik <- sum(theta*ernm::calculateStatistics(object$formula)) - const_approx
  net <- as.BinaryNet(eval(object$formula[[2]],envir=environment(object$formula)))
  n_verts <- net$size()
  n_dyads <- n_verts*(n_verts-1)*(1 - 0.5*(!net$isDirected())) #
  attr(logLik, "df") <- length(theta)
  attr(logLik, "nobs") <- n_verts*(length(object$m$randomVariables())!=0) + n_dyads*object$m$randomGraph()
  class(logLik) <- "logLik"
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

