
#' creates an ERNM likelihood model
#' @param observedSampler a sampler
#' @param unobservedSampler a sampler conditional upon the observed values
#' @param ... additional parameters for the log likelihood
MissingErnmModel <- function(observedSampler, unobservedSampler, ...){
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
	res$info <- function(sample){
		cov(as.matrix(sample))
	}
	
	if(!missing(logLik))
		res$logLikelihood <- logLik
	else
		res$logLikelihood <- function(theta,sample,theta0,stats){
			fullErnmLikelihood(theta=theta, sample=sample, theta0=theta0, stats=stats, ...)
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



#' creates an ERNM likelihood model
#' @param sampler a sampler
#' @param ... additional parameters for the log likelihood
ReGibbsModel <- function(sampler, ...){
	
	sampler <- sampler
	
	res <- list(sampler=sampler)
	
	res$generateSampleStatistics <- function(burnin,interval,size){
		sampler$generateSampleStatistics(burnin,interval,size)
	}
	
	res$name <- function(){
		"Reduced Entropy Gibbs Model"
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
	
	
	res$thetaVar <- function(sample){
		error()
	}
	
	res$info <- function(sample){
		cov(as.matrix(sample))
	}
	
	res$logLikelihood <- function(theta,sample,theta0,stats){
		mod <- sampler$getModel()
		centers <- mod$centers()
		betas <- mod$betas()
		isThetaDependent <- mod$isThetaDependent()
		GmmObjective(theta=theta, centers=centers, betas=betas, isThetaDependent=isThetaDependent, 
				sample=sample, theta0=theta0, stats=stats, ...)
	}
	
	res$scaledGradient <- function(theta,sample,theta0,stats){
		grad <- res$logLikelihood(theta=theta, sample=sample,
				theta0=theta0, stats=stats)$gradient
		vars <- diag(var(sample))
		sd <- sqrt(abs(vars))
		scl <- ifelse(sd<.Machine$double.eps,1 + abs(theta),sd)
		grad / scl
	}
	
	class(res) <- c("ReGibbsModel")
	res
}