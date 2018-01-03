

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

.createLatentOrderLikelihoodFromTerms <- function(terms, net, theta=NULL){
  net <- as.BinaryNet(net)
  model <- .makeCppModelFromTerms(terms, net, theta)
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
		if(hsCount < nHalfSteps && obs$value > lastObs$value && !is.null(lastPar)){
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

.latentOrderMMObjective2 <- function(params, lolik, seed, nReplicates, observedStatistics){
	lolik$setThetas(params)
	nparams <- length(params)
	set.seed(seed)
	samples <- lapply(1:nReplicates,function(x) lolik$generateNetwork())
	momentConditions = rowMeans(sapply(samples,function(samp) observedStatistics - samp$stats))
	cprod <- lapply(samples, function(x) outer(x$stats - x$expectedStats, x$stats))
	gradient <- -Reduce(`+`,cprod) / nReplicates
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
		result <- .latentOrderMMObjective2(params, lolik, seed, nReplicates, observedStatistics)
		#browser()
		print(result)
		result
	}
	browser()
	#optim(thetaInit, function(x) calculateObjective(x)$value)
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



elogFit <- function(formula, theta, nsamp=1000, hotellingTTol= .1, nHalfSteps=10, maxIter=100,
		startingStepSize=maxStepSize, maxStepSize=.5, order=NULL){
	
	lolik <- createLatentOrderLikelihood(formula, theta=theta)
	if(!is.null(order)){
		lolik$setOrder(as.integer(rank(order, ties.method = "min")))
	}
	obsStats <- lolik$getModel()$statistics()
	stepSize <- startingStepSize
	lastTheta <- NULL
	hsCount <- 0
	iter <- 0
	while(iter < maxIter){
		iter <- iter + 1
		
		#generate networks
		lolik$setThetas(theta)
		stats <- matrix(0,ncol=length(theta),nrow=nsamp)
		estats <- matrix(0,ncol=length(theta),nrow=nsamp)
		for(i in 1:nsamp){
			cat(".")
			samp <- lolik$generateNetwork()
			stats[i,] <- samp$stats + samp$emptyNetworkStats
			estats[i,] <- samp$expectedStats + samp$emptyNetworkStats
		}
		cat("\n")
		
		momentCondition <- obsStats - colMeans(stats)
		
		#calculate gradient of moment conditions
		grad <- matrix(0,ncol=length(theta),nrow=length(theta))
		for(i in 1:length(theta)){
			for(j in 1:length(theta)){
				#grad[i,j] <- -mean(stats[,i] * (stats[,j] - estats[,j]))
				grad[i,j] <- -(cov(stats[,i], stats[,j]) - cov(stats[,i], estats[,j]))
			}
		}
		
		
		cat("Moment Conditions:\n")
		print(momentCondition)
		
		
		#calculate inverse of gradient
		invFailed <- inherits(try(gradInv <- solve(grad),silent = TRUE),"try-error")
		#invFailed <- inherits(try(gradInv <- solve(-var(stats)),silent = TRUE),"try-error")
		pairs(stats)
		#browser()
		if(hsCount < nHalfSteps && invFailed && !is.null(lastTheta)){
			cat("Half step back\n")
			theta <- (lastTheta + theta) / 2
			hsCount <- hsCount + 1
			stepSize <- stepSize / 2
			next
		}else{
			stepSize <- min(maxStepSize, stepSize * 1.1)
			hsCount <- 0
		}
		lastTheta <- theta
		theta <- theta - stepSize * gradInv %*% momentCondition
		
		#Hotelling's T^2 test
		hotT <- momentCondition %*% solve(var(stats)/nrow(stats)) %*% momentCondition
		pvalue <- pchisq(hotT,df=length(theta), lower.tail = FALSE)
		cat("Hotelling's T2 p-value: ",pvalue,"\n")
		cat("Theta:\n")
		print(theta)
		if(pvalue > hotellingTTol){
			break
		}else if(iter < maxIter){
			
		}
	}
	vcov <- gradInv %*% var(stats) %*% t(gradInv)
	
	result <- list(theta=lastTheta,
			stats=stats,
			estats=estats, 
			net=samp$network,
			grad=grad, 
			vcov=vcov, 
			likelihoodModel=lolik)
	class(result) <- c("elog","list")
	result
}

summary.elog <- function(x, ...){
	theta <- fit$theta
	se <- sqrt(diag(fit$vcov))
	pvalue <- 2 * pnorm(abs(theta / se),lower.tail = FALSE)
	stats <- fit$likelihoodModel$getModel()$statistics()
	result <- data.frame(observed_statistics=stats, theta=theta, se=se, pvalue=round(pvalue,4))
	rownames(result) <- names(stats)
	result
}


elogGmmFit <- function(formula, auxFormula, theta, nsamp=1000, hotellingTTol= .1, nHalfSteps=10, maxIter=100,
		startingStepSize=.1, maxStepSize=.5, order=NULL, cluster=NULL){
	
	lolik <- createLatentOrderLikelihood(formula, theta=theta)
	if(!is.null(order)){
		lolik$setOrder(as.integer(rank(order, ties.method = "min")))
	}
	terms <- .prepModelTerms(formula)
	auxTerms <- .prepModelTerms(auxFormula)
	auxModel <- createCppModel(auxFormula)
	#browser()
	auxModel$setNetwork(lolik$getModel()$getNetwork())
	auxModel$calculate()
	obsStats <- auxModel$statistics()
	#obsStats <- lolik$getModel()$statistics()
	stepSize <- startingStepSize
	lastTheta <- NULL
	lastObjective <- Inf
	hsCount <- 0
	iter <- 0
	if(!is.null(cluster)){
	  clusterEvalQ(cluster, {
	    library(ernm)
	    library(network)
	  })
	  network <- as.network(lolik$getModel()$getNetwork())
	  clusterExport(cluster, "terms", envir = environment())
	  clusterExport(cluster, "auxTerms", envir = environment())
	  clusterExport(cluster, "network", envir = environment())
	}
	while(iter < maxIter){
		iter <- iter + 1
		
		#generate networks
		lolik$setThetas(theta)
		stats <- matrix(0,ncol=length(theta),nrow=nsamp)
		estats <- matrix(0,ncol=length(theta),nrow=nsamp)
		auxStats <- matrix(0,ncol=length(obsStats),nrow=nsamp)
		sumCov <- array(0, dim=c(nsamp,length(theta),length(theta)))
		if(is.null(cluster)){
			for(i in 1:nsamp){
				cat(".")
				samp <- lolik$generateNetwork()
				auxModel$setNetwork(samp$network)
				auxModel$calculate()
				auxStats[i,] <-  auxModel$statistics()
				stats[i,] <- samp$stats + samp$emptyNetworkStats
				estats[i,] <- samp$expectedStats + samp$emptyNetworkStats
				sumCov[i,,] <- samp$sumCov
			}
			cat("\n")
		}else{
		  workingNetwork <- as.network(lolik$getModel()$getNetwork())
			worker <- function(i, theta){
			  network <- as.BinaryNet(network)
			  lolik <- ernm:::.createLatentOrderLikelihoodFromTerms(terms, network, theta)
			  auxModel <- ernm:::.makeCppModelFromTerms(auxTerms, network)
			  samp <- lolik$generateNetwork()
			  auxModel$setNetwork(samp$network)
			  auxModel$calculate()
			  list(stats=samp$stats + samp$emptyNetworkStats,
			    estats = samp$expectedStats + samp$emptyNetworkStats,
			    sumCov = samp$sumCov)
			}
			worker(1, theta)

			results <- parallel::parLapply(cluster, 1:4, worker, theta=theta)
			browser()
		}
		
		#calculate gradient of moment conditions
		grad <- matrix(0,ncol=length(theta),nrow=length(obsStats))
		for(i in 1:length(obsStats)){
			for(j in 1:length(theta)){
				#cat(i," ",j," ", cov(auxStats[,i], stats[,j]), " ", cov(auxStats[,i], estats[,j]),"\n")
				#grad[i,j] <- -mean(stats[,i] * (stats[,j] - estats[,j]))
				grad[i,j] <- -(cov(auxStats[,i], stats[,j]) - cov(auxStats[,i], estats[,j]))
			}
		}
		
		W <- diag( 1 / (diag(var(auxStats)) + 1) )
		#W <- diag( 1 / (obsStats + 1) )
		hess <- matrix(0,ncol=length(theta),nrow=length(theta))
		mh <- colMeans(auxStats)
		for(i in 1:length(theta)){
			mgi <- mean(stats[,i])
			mGi <- mean(estats[,i])
			for(j in 1:length(theta)){
				mgj <- mean(stats[,j])
				mGj <- mean(estats[,j])
				d2m <- -(colMeans(auxStats * stats[,i] * stats[,j]) - mh * mgi * mgj -
							colMeans(auxStats * stats[,i] * estats[,j]) + mh * mgi * mGj -
							colMeans(auxStats * estats[,i] * stats[,j]) + mh * mGi * mgj +
							colMeans(auxStats * estats[,i] * estats[,j]) - mh * mGi * mGj -
							colMeans(auxStats * stats[,i]) + colMeans(auxStats * sumCov[,i,j]))
				hess[i,j] <- t(grad[,i,drop=FALSE]) %*% W %*% grad[,j,drop=FALSE] + d2m %*% W %*% (obsStats - mh)
				#cat(i," ",j," ", t(grad[,i,drop=FALSE]) %*% W %*% grad[,j,drop=FALSE] , " ", 
				#		d2m %*% W %*% (obsStats - mh),"\n")
			}
		}    
		#browser()
		
		#W <- diag(nrow=length(obsStats),ncol=length(obsStats))
		diffs <- -sweep(auxStats, 2, obsStats)
		transformedDiffs <- t(t(grad) %*% W %*% t(diffs))
		momentCondition <- colMeans(transformedDiffs)
		
		objective <- colMeans(diffs) %*% W %*% colMeans(diffs)
		cat("Objective:\n")
		print(objective)
		objCrit <- max(-1000000, objective - lastObjective) / (lastObjective + 1)
		
		cat("Moment Conditions:\n")
		print(momentCondition)
		
		
		#calculate inverse of gradient
		invFailed <- inherits(try(gradInv <- solve(t(grad) %*% W %*% grad),silent = TRUE),"try-error")
		#invFailed <- inherits(try(gradInv <- solve(-var(stats)),silent = TRUE),"try-error")
		pairs(stats)
		#browser()
		if(hsCount < nHalfSteps && !is.null(lastTheta) && (invFailed || objCrit > .3)){
			cat("Half step back\n")
			theta <- (lastTheta + theta) / 2
			hsCount <- hsCount + 1
			stepSize <- stepSize / 2
			cat("Theta:\n")
			print(theta)
			next
		}else{
			stepSize <- min(maxStepSize, stepSize * 1.25)
			hsCount <- 0
		}
		print(stepSize)
		lastTheta <- theta
		theta <- theta - stepSize * gradInv %*% momentCondition
		#theta <- theta - stepSize * solve(hess) %*% momentCondition
		#browser()
		lastObjective <- objective
		
		print("auxStat Diffs:")
		print(colMeans(diffs) / sqrt(diag(var(diffs))))
		
		#Hotelling's T^2 test
		hotT <- momentCondition %*% solve(var(transformedDiffs)/nrow(transformedDiffs)) %*% momentCondition
		pvalue <- pchisq(hotT,df=length(theta), lower.tail = FALSE)
		cat("Hotelling's T2 p-value: ",pvalue,"\n")
		cat("Theta:\n")
		print(theta)
		#browser()
		if(pvalue > hotellingTTol){
			break
		}else if(iter < maxIter){
			
		}
	}
	omega <- var(auxStats)
	vcov <- solve(t(grad) %*% W %*% grad) %*% 
			t(grad) %*% W %*% omega %*% t(W) %*% grad %*% 
			solve(t(grad) %*% t(W) %*% grad) 
	#vcov <- gradInv %*% var(stats) %*% t(gradInv)
	
	result <- list(theta=lastTheta,
			stats=stats,
			estats=estats, 
			auxStats=auxStats,
			obsStats=obsStats,
			net=samp$network,
			grad=grad, 
			vcov=vcov, 
			likelihoodModel=lolik)
	class(result) <- c("elog","list")
	result
}