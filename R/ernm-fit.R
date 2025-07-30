
#' Fit an ernm
#'
#' This is a lower level MCMC-MLE fitting function for ERNM. Users should generally use
#' the ernm() function instead.
#'
#' @param sampler the ErnmModel
#' @param theta0 initial starting values
#' @param mcmcBurnIn MCMC burn in
#' @param mcmcInterval MCMC interval
#' @param mcmcSampleSize MCMC sample size
#' @param minIter minimum number of MCMC-MLE iterations
#' @param maxIter maximum number of MCMC-MLE iterations
#' @param objectiveTolerance convergence criteria on change in log likelihood ratio
#' @param gradTolerance convergence criteria on scaled gradient
#' @param meanStats optional target statistics for the mean value parameters
#' @param verbose level of verbosity 0, 1, or 2
#' @param method the optimization method to use. "bounded" uses trust regions around the MCMC sample and is generally preferable. See Fellows (2012) for details.
#' @export
#' @return an ernm object
#' @references
#' Fellows, Ian Edward. Exponential family random network models. University of California, Los Angeles, 2012.
ernmFit <- function(sampler,
                    theta0,
                    mcmcBurnIn=10000,
                    mcmcInterval=100,
                    mcmcSampleSize=10000,
                    minIter=3,
                    maxIter=40,
                    objectiveTolerance=.5,
                    gradTolerance=.25,
                    meanStats,
                    verbose=1,
                    method=c("bounded","newton")){
	method <- match.arg(method)
	sampler$initialize()
	stats <- sampler$modelStatistics()
	if(!missing(meanStats))
		stats <- meanStats
	logLikelihoodFun <- sampler$logLikelihood

	if(!missing(theta0)){
		sampler$setThetas(theta0)
	}
	theta0 <- sampler$thetas()

	iter<-1
	converged <- FALSE
	likHistory <- c()
	gradHistory <- list()
	trace <- list()
	while(iter<maxIter){
		sample <- sampler$generateSampleStatistics(mcmcBurnIn,mcmcInterval,mcmcSampleSize)
		sampStats <- sapply(sampler$statistics(sample),function(x)colMeans(as.matrix(x)))
		sampSd <- sapply(sampler$statistics(sample),function(x)apply(x,2,sd))
		if(verbose>0){
			cat("iteration:",iter,"\n")
			cat("       parameters:\n")
			print(theta0)
			cat("sample statistics:\n")
			cat("       means:\n")
			sampStats <- t(sampStats)
			if(nrow(sampStats)==1){
				sampStats <- rbind(sampStats,stats)
				rownames(sampStats) <- c("simulated","observed")
			}
			print(sampStats)
			cat("       std:\n")
			print(t(sampSd))
		}
		if(verbose>1){
			nr <- mcmcSampleSize
			nc <- ncol(sampler$statistics(sample)[[1]])
			# make sure to fix pr on exit
			oldpar <- par(no.readonly = TRUE)
			on.exit(par(oldpar))
			par(mfrow=c(ceiling(nc/3),3))
			for(i in 1:nc)
				try(scatter.smooth(1:nr,sampler$statistics(sample)[[1]][,i],col="red"))
		}
		scaledGrad <- sampler$scaledGradient(theta0,sample,theta0,stats)
		maxGrad <- max(abs(scaledGrad))
		lastTheta <- theta0
		if(method == "bounded"){
			tty <- try(trustRes <- trust(objfun=logLikelihoodFun,
				parinit=theta0,
				rinit=1,
				rmax=100,
				parscale=rep(1,length(theta0)), minimize=FALSE,
				sample=sample,
				theta0=theta0,
				stats=stats))
			if(inherits(tty,"try-error"))
			if(!trustRes$converged)
				warning("Trust: convergance not met")
			theta0 <- trustRes$argument
		}else if(method == "newton"){
			llk <- sampler$logLikelihood(theta0,sample,theta0,stats)
			theta0 <- theta0 - drop(solve(llk$hessian, llk$grad))
			#theta0 <- theta0 - drop(qr.solve(llk$hessian) %*% llk$grad)
		}
		llik <- sampler$logLikelihood(theta0,sample,lastTheta,stats)$value
		if(verbose>0){
			cat("\nlog likelihood improved by: ",llik,"\n",
					"maximum scaled gradient: ",maxGrad,"\n")
		}
		likHistory <- c(likHistory,llik)
		gradHistory[[length(gradHistory) + 1]] <- scaledGrad
		trace[[length(trace) + 1]] <- lastTheta

		sampler$setThetas(theta0)
		if(iter>minIter && llik<objectiveTolerance && maxGrad<gradTolerance){
			sampler$setThetas(lastTheta)
			converged<-TRUE
			break
		}
		iter<-iter+1
	}
	if(verbose>0)
		cat("\n")
	if(verbose>0)
		cat("\n")
	result <- list(theta=lastTheta,
	               converged=converged,
	               iter=iter,
	               info=sampler$info(sample),
	               objectiveDiff=llik,
	               maxScaledGradiant=maxGrad,
	               sample=sample,
	               likelihoodHistory=likHistory,
	               gradientHistory=gradHistory,
	               trace=trace,
	               m=sampler)
	class(result) <- "ernm"
	result
}

