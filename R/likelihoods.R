

#' likelihood for a fully observed ernm
#' @param theta parameters
#' @param sample mcmc sample
#' @param theta0 parameter values which generated sample
#' @param stats observed statistics
#' @param minEss minimum effective sample size
#' @param damping a damping parameter
#' @param method the method of partition function approximation to use
#' @param order the order of the cumulant approximation
#' @export
#' @return a list with value, gradient, and hessian
fullErnmLikelihood <- function(theta,
                               sample,
                               theta0,
                               stats,
                               minEss=5,
                               damping=.05,
                               method=c("cumulant","sample"),
                               order=3){
	method <- match.arg(method)

	llik <- function(){
		thetadiff <- theta-theta0
		expon <- mat %*% thetadiff

		#if thetaExtended is (almost) outside convex hull don't trust theta
		thetaExtended <- theta + (theta-theta0)*damping
		exponExtended <- mat %*% (thetaExtended - theta0)
		n <- nrow(mat)
		wts <- exp(exponExtended) / sum(exp(exponExtended))
		wtdSe <- apply(mat,2,function(x) mcmcse(wts*x))
		se <- apply(mat,2,function(x) sd(x/n)/sqrt(n))
		ess <- min(n * se^2 / wtdSe^2)
		if(is.finite(ess) && ess<minEss){
			if(all(theta==theta0))
				stop("Insufficient MCMC variation")
			return(-Inf)
		}

		if(method=="sample"){
			mLnorm <- mean(expon)
			vLnorm <- var(drop(expon))
			part <- log(mean(exp(expon)))
			llr <- sum(stats*thetadiff) - part
		}else if(method=="cumulant"){
			x <- drop(expon)
			mom <- all.moments(x,order.max=order)
			cum <- all.cumulants(mom)[-1]
			part <- sum(cum/factorial(1:order))
			llr <- sum(stats*thetadiff) - part
		}
		llr
	}

	dl <- function(){
		expon <- mat %*% (theta-theta0)
		wts <- exp(expon)/sum(exp(expon),na.rm=TRUE)
		if(any(is.nan(wts)))
			wts <- rep(1/length(wts),length(wts))
		diff <- stats - apply(mat,2,function(x) sum(wts*x,na.rm=TRUE))
		diff[!is.finite(diff)] <- 0
		diff
	}
	hessian <- function(){
		expon <- mat %*% (theta-theta0)
		wts <- drop(exp(expon)/sum(exp(expon)))
		if(any(is.nan(wts)))
			wts <- rep(1/length(wts),length(wts))
		wmat <- mat*sqrt(drop(wts))
		-crossprod(wmat)

	}

	#center statistics
	mat <- sample
	mns <- colMeans(mat)
	for(i in 1:ncol(mat))
		mat[,i] <- mat[,i] - mns[i]
	stats <- stats - mns


	value<-llik()
	grad<-dl()
	hess<-hessian()
	hess[upper.tri(hess)]<-t(hess)[upper.tri(hess)]
	list(value=value,gradient=as.vector(grad),hessian=hess)
}


#' likelihood for an ernm with missing data
#' @param theta parameters
#' @param sample mcmc sample
#' @param theta0 parameter values which generated sample
#' @param stats observed statistics
#' @param minEss minimum effective sample size
#' @param damping a damping parameter
#' @export
#' @return a list with value, gradient, and hessian
marErnmLikelihood <- function(theta,sample,theta0,stats,minEss=5, damping=.1){

	llik <- function(){
		#thetadiff <- theta-theta0
		fullExpon <- full %*% (theta-theta0)
		missExpon <- miss %*% (theta-theta0)
		n <- nrow(full)
		#########
		#if thetaExtended is (almost) outside convex hull don't trust theta
		thetaExtended <- theta + (theta-theta0)*damping

		#Full
		exponExtended <- full %*% (thetaExtended - theta0)
		wts <- exp(exponExtended) / sum(exp(exponExtended))
		fullWtdSe <- apply(full,2,function(x) mcmcse(wts*x))
		fullSe <- apply(full,2,function(x) sd(x/n)/sqrt(n))

		#Conditional
		exponExtended <- miss %*% (thetaExtended - theta0)
		wts <- exp(exponExtended) / sum(exp(exponExtended))
		missWtdSe <- apply(miss,2,function(x) mcmcse(wts*x))
		missSe <- apply(miss,2,function(x) sd(x/n)/sqrt(n))

		ess <- min(n * (missSe^2 + fullSe^2) / (missWtdSe^2 + fullWtdSe^2))
		if(!is.finite(ess) || ess<minEss){
			if(all(theta==theta0))
				stop("Insufficient MCMC variation")
			return(-Inf)
		}
		#
		#########


		#calculate unconditional term
		x <- drop(fullExpon)
		n <- length(x)
		mn <- mean(x)
		std <- sd(x)
		skew <- sqrt(n) * sum(x^3)/(sum(x^2)^(3/2)) * ((1 - 1/n))^(3/2)
		if(!is.finite(skew))
			skew <- 0
		fullPart <- mn+std^2/2 + (std^3/6)*skew

		#calculate conditional term
		x <- drop(missExpon)
		n <- length(x)
		mn <- mean(x)
		std <- sd(x)
		skew <- sqrt(n) * sum(x^3)/(sum(x^2)^(3/2)) * ((1 - 1/n))^(3/2)
		if(!is.finite(skew))
			skew <- 0
		missPart <- mn+std^2/2 + (std^3/6)*skew
		llr <- missPart - fullPart

		llr
	}
	dl <- function(){

		fullExpon <- full %*% (theta-theta0)
		fullWts <- exp(fullExpon)/sum(exp(fullExpon))

		missExpon <- miss %*% (theta-theta0)
		missWts <- exp(missExpon)/sum(exp(missExpon))

		if(all(abs(theta-theta0)<.Machine$double.eps))
			fullWts <- missWts <- rep(1 / length(missExpon),length(missExpon))

		diff <- apply(miss,2,function(x) sum(missWts*x)) - apply(full,2,function(x) sum(fullWts*x))
		diff
	}
	hessian <- function(){
		fullExpon <- drop(full %*% (theta-theta0))
		fullWts <- exp(fullExpon)/sum(exp(fullExpon))
		missExpon <- drop(miss %*% (theta-theta0))
		missWts <- exp(missExpon)/sum(exp(missExpon))

		if(all(abs(theta-theta0)<.Machine$double.eps))
			fullWts <- missWts <- rep(1,length(missExpon))

		if(any(is.na(fullWts)))
			fullWts <- rep(1,length(fullWts))
		if(any(is.na(missWts)))
			missWts <- rep(1,length(missWts))

		fullCov <- cov.wt(full,drop(fullWts))$cov
		missCov <- cov.wt(miss,drop(missWts))$cov
		missCov-fullCov
	}

	#center statistics
	#miss : only missing values are toggled
	#full : whole model toggleing
	full <- sample[[1]]
	fullMns <- colMeans(full)
	for(i in 1:ncol(full))
		full[,i] <- full[,i] - fullMns[i]

	miss <- sample[[2]]
	for(i in 1:ncol(miss))
		miss[,i] <- miss[,i] - fullMns[i]
	value<-llik()
	grad<-dl()
	hess<-hessian()
	hess[upper.tri(hess)]<-t(hess)[upper.tri(hess)]
	if(any(!is.finite(grad)))
		value <- -Inf
	list(value=value,gradient=as.vector(grad),hessian=hess)
}





#' (E(g(X)) - g(x_o)^2 for TaperedModel
#' @param theta parameters
#' @param centers center of statistics
#' @param tau tapering parameter
#' @param sample mcmc sample
#' @param theta0 parameter values which generated sample
#' @param stats observed statistics
#' @param minEss minimum effective sample size
#' @param damping a damping parameter
#' @export
#' @return a list with value, gradient, and hessian
taperedErnmLikelihood <- function(theta,
                                  centers,
                                  tau,
                                  sample,
                                  theta0,
                                  stats,
                                  minEss=5,
                                  damping=.05){

	b <- function(t){tau}
	db <- function(t){rep(0,length(t))}


	llik <- function(){
		thetadiff <- theta-theta0

		expon <- mat %*% thetadiff - penmat %*% (b(theta) - b(theta0))

		#if thetaExtended is (almost) outside convex hull don't trust theta
		thetaExtended <- theta + (theta-theta0)*damping
		exponExtended <- mat %*% (thetaExtended - theta0) - penmat %*% (b(thetaExtended) - b(theta0))
		n <- nrow(mat)
		wts <- exp(exponExtended) / sum(exp(exponExtended))
		wtdSe <- apply(mat,2,function(x) mcmcse(wts*x))
		se <- apply(mat,2,function(x) sd(x/n)/sqrt(n))
		ess <- min(n * se^2 / wtdSe^2)
		if(is.finite(ess) && ess<minEss){
			if(all(theta==theta0))
				stop("Insufficient MCMC variation")
			return(-Inf)
		}
		diff <- stats - apply(mat,2,function(x) sum(wts*x,na.rm=TRUE))
		- .5 * sum(diff^2)
	}

	dl <- function(){
		expon <- mat %*% (theta-theta0) - penmat %*% (b(theta) - b(theta0))
		wts <- exp(expon)/sum(exp(expon),na.rm=TRUE)
		if(any(is.nan(wts)))
			wts <- rep(1/length(wts),length(wts))
		diff <- stats - apply(mat,2,function(x) sum(wts*x,na.rm=TRUE))
		diff[!is.finite(diff)] <- 0
		diff
	}

	hessian <- function(){

		expon <- mat %*% (theta-theta0) - penmat %*% (b(theta) - b(theta0))
		wts <- drop(exp(expon)/sum(exp(expon)))
		if(any(is.nan(wts)))
			wts <- rep(1/length(wts),length(wts))
		wmat <- mat*sqrt(drop(wts))
		-crossprod(wmat)

	}

	#center statistics
	mat <- sample
	penmat <- sample
	mns <- colMeans(mat)
	for(i in 1:ncol(mat)){
		mat[,i] <- mat[,i] - mns[i]
		penmat[,i] <- penmat[,i] - centers[i]
	}
	stats <- stats - mns

	#penalty term
	#pen <- sum(mns - centers)^2

	value<-llik()
	grad<-dl()
	hess<-hessian()
	hess[upper.tri(hess)]<-t(hess)[upper.tri(hess)]
	list(value=value,gradient=as.vector(grad),hessian=hess)
}
