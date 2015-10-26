

#' Creates a probability model for a latent ordered network model
createLatentOrderLikelihood <- function(formula, theta=NULL){
	env <- environment(formula)
	net <- as.BinaryNet(eval(formula[[2]],envir=env))
	model <- createCppModel(formula)
	clss <- class(net)
	networkEngine <- substring(clss,6,nchar(clss)-3)
	LikType <- eval(parse(text=paste(networkEngine,"LatentOrderLikelihood",sep="")))
	lik <- new(LikType, model)
	if(!is.null(theta)){
		lik$setThetas(theta)
	}
	lik
}


.latentOrderObjective <- function(params, lolik, seed, nReplicates){
	lolik$setThetas(params)
	nparams = length(params)
	set.seed(seed)
	samples <- lolik$fullLogLik(as.integer(nReplicates), 0.0)
	lliks <- sapply(samples, function(x) x$logLik) - sapply(samples, function(x) x$logPartition)
	medLik <- median(lliks)
	combinedLik <- log(sum(exp(lliks - medLik)) / length(lliks)) + medLik
	
	derivDenom <- sapply(lliks, function(x) mean(exp(lliks - x)))
	condDerivs <- sapply(samples, function(x) x$derivative)
	deriv <- rowMeans(sweep(condDerivs, 2, derivDenom,FUN="/"))
	
	condHessians <- lapply(samples, function(x) x$hessian)
	hessian <- matrix(0, nrow=nparams, ncol=nparams)
	for(i in 1:nparams){
		for(j in 1:nparams){
			ddij <- sapply(samples, function(x) x$hessian[i,j])
			di <- sapply(samples, function(x) x$derivative[i])
			dj <- sapply(samples, function(x) x$derivative[j])
			hessian[i,j] <- mean(ddij / derivDenom) + mean(di * dj / derivDenom) - mean(di / derivDenom) * mean(dj / derivDenom)
		}
	}
	list(value = combinedLik, gradient=deriv, hessian=hessian)
}

.latentOrderUnweightedObjective <- function(params, lolik, seed, nReplicates){
	lolik$setThetas(params)
	nparams <- length(params)
	set.seed(seed)
	samples <- lolik$fullLogLik(as.integer(nReplicates), 0.0)
	lliks <- sapply(samples, function(x) x$logLik) - sapply(samples, function(x) x$logPartition)
	combinedLik <- mean(lliks)
	
	condDerivs <- sapply(samples, function(x) x$derivative)
	deriv <- rowMeans(condDerivs)
	print(combinedLik)
	print(deriv)
	condHessians <- lapply(samples, function(x) x$hessian)
	hessian <- matrix(0, nrow=nparams, ncol=nparams)
	for(i in 1:nparams){
		for(j in 1:nparams){
			ddij <- sapply(samples, function(x) x$hessian[i,j])
			hessian[i,j] <- mean(ddij)
		}
	}
	list(value = combinedLik, gradient=deriv, hessian=hessian)
}

newtonOptim <- function(fn, par, tol= 10e-5, nHalfSteps=0, maxIter=100, ...){
	lastPar <- NULL
	lastObs <- list(value=-Inf)
	hsCount <- 0
	iter <- 0
	while(iter < maxIter){
		iter <- iter + 1
		obs <- fn(par)
		if(hsCount < nHalfSteps && obs$value < lastObs$value && !is.null(lastPar)){
			cat("Half step back\n")
			par <- (lastPar + par) / 2
			hsCount <- hsCount + 1
			next
		}else{
			hsCount <- 0
		}
		lastPar <- par
		#inv <- solve(obs$hessian)
		#inv <- -chol2inv(chol(-obs$hessian))
		#par <- par -  inv %*% obs$gradient
		#browser()
		par <- par -  solve(obs$hessian, obs$gradient)
		if(abs(lastObs$value - obs$value) < tol){
			break
		}else if(iter < maxIter){
			lastObs <- obs
		}
	}
	obs$argument <- drop(par)
	obs$converged <- abs(lastObs$value - obs$value) < tol
	obs
}

#' Fits a latent ordered network model using maximum likelihood
latentOrderFit <- function(formula, nReplicates=500L, theta=NULL,
		seed = floor(runif(1, 0,.Machine$integer.max)), variational=FALSE ){
	
	lolik <- createLatentOrderLikelihood(formula, theta=theta)
	thetaInit <- lolik$getModel()$thetas()
	if(variational)
		method <- "variational"
	else
		method <- "ml"
	calculateObjective <- function(params){
		cat(".")
		if(variational)
			.latentOrderUnweightedObjective(params, lolik, seed, nReplicates)
		else	
			.latentOrderObjective(params, lolik, seed, nReplicates)
	}
	#browser()
	params <- trust(calculateObjective, thetaInit, 1, 1e+20, minimize=FALSE, blather=TRUE)
	#params <- newtonOptim(calculateObjective, thetaInit)
	#browser()
	if(!params$converged){
		print(params)
		stop("Convergance failure")
	}
	lolik$setThetas(params$argument)
	set.seed(seed)
	samples <- lolik$fullLogLik(nReplicates, 0)
	result <- list(logLikelihood=params$value, theta=params$argument, 
			gradient=params$gradient, hessian=params$hessian, method=method, 
			likelihoodModel=lolik, randomSeed=seed, nReplicates=nReplicates, samples=samples)
	class(result) <- "latentOrderFit"
	result
}

.latentOrderMMObjective <- function(params, lolik, seed, nReplicates, observedStatistics){
	lolik$setThetas(params)
	nparams <- length(params)
	set.seed(seed)
	samples <- lapply(1:nReplicates,function(x) lolik$generateNetwork())
	#browser()
	momentConditions = rowMeans(sapply(samples,function(samp) observedStatistics - samp$expectedStats - samp$emptyNetworkStats))
	gradient <- matrix(0, nrow=nparams, ncol=nparams)
	for(i in 1:nparams){
		for(j in 1:nparams){
			ddij <- sapply(samples, function(x) x$gradient1[i,j])
			gradient[i,j] <- mean(ddij)
		}
	}
	#browser()
	list(value=sum(momentConditions^2), gradient=momentConditions, hessian=gradient)
}

latentOrderFitMM <- function(formula, nReplicates=500L, theta=NULL,
		seed = floor(runif(1, 0,.Machine$integer.max)), variational=FALSE ){
	
	lolik <- createLatentOrderLikelihood(formula, theta=theta)
	thetaInit <- lolik$getModel()$thetas()
	observedStatistics <- calculateStatistics(formula=formula)
	calculateObjective <- function(params){
		cat(".")
		cat("params")
		print(params)
		result <- .latentOrderMMObjective(params, lolik, seed, nReplicates, observedStatistics)
		#browser()
		print(result)
		result
	}
	browser()
	optim(thetaInit, calculateObjective)
	newtonOptim(calculateObjective, thetaInit, nHalfSteps=0)
}
#' Fits a latent ordered network model using generalized method of moments (i.e. the variational approximation)
latentOrderGmmFit <- function(formula, nReplicates=500L, theta=NULL,
		seed = floor(runif(1, 0,.Machine$integer.max)), approxSampleSize=100000){
	set.seed(seed)
	lolik <- createLatentOrderLikelihood(formula, theta=theta)
	
	n <- lolik$getModel()$getNetwork()$size()
	if(lolik$getModel()$getNetwork()$isDirected()){
		pSample <- approxSampleSize / ( nReplicates * n * (n - 1)) 
	}else{
		pSample <- 2 * approxSampleSize / ( nReplicates * n * (n - 1)) 
	}
	samples <- lolik$fullLogLik(nReplicates, pSample)
	predictors <- do.call(rbind, lapply(samples,function (x) t(sapply(x$samples, function(y) y$changeStats))))
	outcome <- unlist(lapply(samples,function (x) sapply(x$samples, function(y) y$hasEdge)))
	theta <- glm.fit(predictors, outcome, family=binomial(),intercept=FALSE)$coef
	obj <- .latentOrderObjective(theta, lolik, seed, nReplicates)
	lolik$setThetas(theta)
	set.seed(seed)
	samples <- lolik$fullLogLik(nReplicates, pSample)
	
	result <- list(logLikelihood=obj$value, theta=theta, 
			gradient=obj$gradient, hessian=obj$hessian, method="gmm", 
			likelihoodModel=lolik, randomSeed=seed, nReplicates=nReplicates, 
			frame=list(y=outcome,X=predictors), samples=samples)	
	class(result) <- "latentOrderFit"
	result
}


#' parameter covariance matrix
#' @param object object
#' @param ... unused
#' @method vcov latentOrderFit
vcov.latentOrderFit <- function(object,...){
	solve(-object$hessian)
}
