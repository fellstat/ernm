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
#' @param include_AIC whether to include AIC in the summary, will slow down
#' @return a data frame summary of the model
#' @export
#' @method summary ernm
summary.ernm <- function(object,include_AIC = TRUE,...){

  theta <- object$theta
  cv <- solve(object$info)
  se <- sqrt(diag(cv))
  z <- theta/se
  p.value <- 2*pnorm(abs(z),lower.tail=FALSE)
  d <- data.frame(theta,se,z,p.value)
  rownames(d) <- make.unique(names(theta))

  # Compute AIC and BIC with latest sample - no bridge sampling for now
  # generate new bigger sample
  if(include_AIC){
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
    logLik <- sum(theta*ernm::calculateStatistics(object$formula)) - const_approx
    net <- eval(object$formula[[2]],envir=environment(object$formula))
    n_verts <- net %n% 'n'
    n_dyads <- n_verts*(n_verts-1)*(1 - 0.5*(!is.directed(net)))
    BIC <- -2*logLik + length(theta)*log(n_verts*(length(object$m$randomVariables())!=0) + n_dyads*object$m$randomGraph())
    AIC <- -2*logLik + 2*length(theta)
    BIC <- round(BIC,2)
    AIC <- round(AIC,2)

    print(round(d,4))
    cat(paste("\nBIC:",BIC,"AIC:",AIC, "(lower is better)\n"))
  }

  # Return the data frame invisibly
  invisible(d)
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

