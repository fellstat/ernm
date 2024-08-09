
#' creates an ERNM likelihood model
#' @param observedSampler a sampler
#' @param unobservedSampler a sampler conditional upon the observed values
#' @param ... additional parameters for the log likelihood
MissingErnmModel <- function(observedSampler,
                             unobservedSampler,
                             ...){
	res <- list(obsSamp = observedSampler, missSamp = unobservedSampler)
	
	res$generateSampleStatistics <- function(burnin,interval,size){
		
		obs <- observedSampler$generateSampleStatistics(burnin,interval,size)
		miss <-unobservedSampler$generateSampleStatistics(burnin,interval,size)
		list(uncondiational=obs,conditional=miss)
	}
	
	res$name <- function(){
		"ERNM Model (with missing values)"
	}

	res$randomGraph <- function(){
		observedSampler$getModel()$hasRandomGraph()
	}
	
	res$randomVariables <- function(){
		observedSampler$getModel()$getRandomVariables()
	}
	
	res$initialize <- function(){
		res$obsSamp$getModel()$calculate()
		res$missSamp$getModel()$calculate()
	}
	
	res$thetas <- function(){
		res$obsSamp$getModel()$thetas()
	}
	
	res$setThetas <- function(newThetas){
		res$obsSamp$getModel()$setThetas(newThetas)
		res$missSamp$getModel()$setThetas(newThetas)
	}
	
	res$info <- function(sample){
		cov(sample[[1]]) - cov(sample[[2]])
	}
	
	res$modelStatistics <- function(full=TRUE){
		if(full)
			as.matrix(res$obsSamp$getModel()$statistics())
		else
			as.matrix(res$missSamp$getModel()$statistics())
	}
	
	res$statistics <- function(sample){
		sample[1:2]
	}
	
	res$logLikelihood <- function(theta,sample,theta0,stats){
		marErnmLikelihood(theta=theta, sample=sample, theta0=theta0, stats=stats, ...)
	}
	
	res$scaledGradient <- function(theta,sample,theta0,stats){
		likes <- res$logLikelihood(theta=theta, sample=sample,
				theta0=theta0, stats=stats)
		grad <- likes$gradient
		#vars <- diag(likes$hess)
		vars <- diag(var(sample[[1]]))
		sd <- sqrt(abs(vars))
		scl <- ifelse(sd<.Machine$double.eps,1 + abs(theta),sd)
		grad / scl
	}
	
	class(res) <- c("MarModel","ErnmModel")
	res
}



#' creates an ERNM likelihood model
#' @param sampler a sampler
#' @param logLik a log likelihood function (optional)
#' @param ... additional parameters for the log likelihood
FullErnmModel <- function(sampler, logLik, ...){
	
	sampler <- sampler
	model <- sampler$getModel()
	if(class(model)[1] == 'Rcpp_UndirectedModel' | class(model)[1] == 'Rcpp_DirectedModel'){
	    tapered <- FALSE
	}else{
	    tapered <- TRUE
	}
	
	res <- list(sampler=sampler)
	
	res$generateSampleStatistics <- function(burnin,interval,size){
		sampler$generateSampleStatistics(burnin,interval,size)
	}
	
	res$name <- function(){
		"ERNM Model"
	}
	
	res$randomGraph <- function(){
		sampler$getModel()$hasRandomGraph()
	}
	
	res$randomVariables <- function(){
		sampler$getModel()$getRandomVariables()
	}
	
	res$initialize <- function(){
		sampler$getModel()$calculate()
	}
	
	res$thetas <- function(){
		sampler$getModel()$thetas()
	}
	
	res$setThetas <- function(newThetas){
		sampler$getModel()$setThetas(newThetas)
	}
	
	res$modelStatistics <- function(){
		sampler$getModel()$statistics()
	}
	
	res$statistics <- function(sample){
		list(sample)
	}
	# if not tapered use regular ERNM likelihood
	if(tapered){
	    res$info <- function(sample){
	        theta_dep <- FALSE
	        tau <- sampler$getModel()$tau()
	        theta <- sampler$getModel()$thetas()
	        mu_hat <- apply(sample,2,mean)
	        np <- length(theta)
	        
	        # Adjusted estimated info based on tapering
	        covar <- cov(as.matrix(sample),as.matrix(sample))
	        covar_2 <- cov(as.matrix(sample),as.matrix(sample**2))
	        
	        # USe sweep as in the tapered.ERGM package    
	        left_mat = sweep(covar,2,2*tau,'*')
	        if(theta_dep){
	            right_mat = covar + (-2)*sweep(covar,2,mu_hat/(tau**2),FUN = "*") + sweep(covar_2,2,1/(tau**2),FUN = "*")
	        }else{
	            right_mat = covar
	        }
	        dmu = solve(left_mat + diag(1,nrow(left_mat))) %*% right_mat
	        
	        # note direct COPY from ergm.tapered code
	        # second derivative of log likelihoods
	        ddll <- diag(rep(0,np))
	        dimnames(ddll) <- list(colnames(covar), colnames(covar))
	        for(i in 1:np){
	            for(j in 1:np){
	                ddll[i,j] <- -dmu[i,j] - sum(2*tau*dmu[,i]*dmu[,j]) 
	            }
	        }
	        -ddll
	    }
        if(!missing(logLik))
            res$logLikelihood <- logLik
        else
            res$logLikelihood <- function(theta,sample,theta0,stats){
                mod <- sampler$getModel()
                centers <- mod$centers()
                tau <- mod$tau()
                taperedErnmLikelihood(theta=theta,
                                      centers=centers,
                                      tau=tau,
                                      sample=sample,
                                      theta0=theta0,
                                      stats=stats, ...)
            }
	}else{
	    res$info <- function(sample){
	        cov(as.matrix(sample))
	    }
	    if(!missing(logLik))
	        res$logLikelihood <- logLik
	    else
	        res$logLikelihood <- function(theta,sample,theta0,stats){
	            fullErnmLikelihood(theta=theta,
	                               sample=sample,
	                               theta0=theta0,
	                               stats=stats,
	                               ...)
	        }
	}
	
	res$scaledGradient <- function(theta,sample,theta0,stats){
		grad <- res$logLikelihood(theta=theta, sample=sample,
				theta0=theta0, stats=stats)$gradient
		vars <- diag(var(sample))
		sd <- sqrt(abs(vars))
		scl <- ifelse(sd<.Machine$double.eps,1 + abs(theta),sd)
		grad / scl
	}
	
	class(res) <- c("FullyObservedModel","ErnmModel")
	res
}