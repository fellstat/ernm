
#' fit an ernm model
#' @param sampler the ErnmModel
#' @param theta0 initial starting values
#' @param mcmcBurnIn burn in
#' @param mcmcInterval interval
#' @param mcmcSampleSize sample size
#' @param minIter minimum number of iterations
#' @param maxIter maximum number of iterations
#' @param objectiveTolerance convergance criteria on change in log likelihood ratio
#' @param gradTolerance convergance criteria on scaled gradient
#' @param meanStats if non-missing, these are the target statistics
#' @param verbose level of verbosity 0, 1, or 2
#' @param method the optimization method to use
#' @param stepScale scale the optimization step
ernmFit <- function(sampler,theta0,
		mcmcBurnIn=10000, mcmcInterval=100, mcmcSampleSize=10000,
		minIter=3, maxIter=40, objectiveTolerance=.5, gradTolerance=.25,
		meanStats,verbose=1,method=c("bounded","newton"),stepScale =1){
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
	if(verbose==0){
		cat("calculating: ")
	}
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
		}else{
			cat(".")
		}
		
		if(verbose>1){
			nr <- mcmcSampleSize
			nc <- ncol(sampler$statistics(sample)[[1]])
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
				browser()
			if(!trustRes$converged)
				warning("Trust: convergance not met")
			theta0 <- trustRes$argument
		}else if(method == "newton"){
			llk <- sampler$logLikelihood(theta0,sample,theta0,stats)
			theta0 <- theta0 - stepScale*drop(solve(llk$hessian, llk$grad))
			#theta0 <- theta0 - drop(qr.solve(llk$hessian) %*% llk$grad)
		}
		llik <- sampler$logLikelihood(theta0,sample,lastTheta,stats)$value
		#browser()
		#grad <- sampler$logLikelihood(lastTheta,sample,lastTheta,stats)$gradient
		#scaledGrad <- grad/sampSd
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
	if(verbose>=0)
		cat("\n")
	if(verbose>0)
		cat("\n")
	result <- list(theta=lastTheta,converged=converged,iter=iter,info=sampler$info(sample),
			objectiveDiff=llik, maxScaledGradiant=maxGrad,sample=sample,
			likelihoodHistory=likHistory,gradientHistory=gradHistory,trace=trace,m=sampler)
	class(result) <- "ernm"
	result
}

#' print
#' @param x x
#' @param ... unused
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

#' summary
#' @param object object
#' @param ... unused
#' @method summary ernm
summary.ernm <- function(object,...){
	
	theta <- object$theta
	cv <- solve(object$info)
	se <- sqrt(diag(cv))
	z <- theta/se
	p.value <- 2*pnorm(abs(z),lower.tail=FALSE)
	d <- data.frame(theta,se,z,p.value)
	rownames(d) <- make.unique(names(theta))
	d
}

#' parameter covariance matrix
#' @param object object
#' @param ... unused
#' @method vcov ernm
vcov.ernm <- function(object,...){
	solve(object$info)
}


#' plot an ernm object
#' @param  x the object
#' @param ... unused
#' @method plot ernm
plot.ernm <- function(x,...){
	plot(x$likelihoodHistory,main="Likelihood convergance",
			ylab="Change in log-likelihood",xlab="iteration")
}

