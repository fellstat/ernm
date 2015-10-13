

initLatent <- function(name, levels, lower=NULL,upper=NULL){
	cont <- TRUE
	if(!missing(levels)){
		cont <- FALSE
	}
	if(cont)
		return(list(name=name,type='continuous',lower=lower,upper=upper))
	else
		return(list(name=name,type='discrete',levels=levels))
}

#' creates a model
#' @param formula the model formula
#' @param ignoreMnar ignore missing not at random offsets
#' @param cloneNet should the network be cloned
#' @param theta the model parameters.
createCppModel <- function(formula,ignoreMnar=TRUE,cloneNet=TRUE,theta=NULL){
	form <- formula
	env <- environment(form)
	net <- as.BinaryNet(eval(form[[2]],envir=env))
	if(cloneNet)
		net <- net$clone()
	noDyad <- FALSE
	randomVars <- character()
	if(!is.symbol(formula[[3]]) && as.character(formula[[3]][[1]])=="|"){
		lastTerm <- FALSE
		tmp <- formula[[3]][[3]]
		while(!lastTerm){
			if(is.symbol(tmp)){
				lastTerm <- TRUE
				term <- tmp
			} else if(as.character(tmp[[1]])=="+"){
				term <- tmp[[3]]
			}else{
				lastTerm <- TRUE
				term <- tmp
			}
			if(length(term)==1){
				val <- as.character(term)
				if(val == "noDyad")
					noDyad <- TRUE
				else
					randomVars[length(randomVars)+1] <- val
			}else if(as.character(term[[1]]) == "latent"){
				term[[1]] <- as.name("initLatent")
				latent <- eval(term,envir=env)
				if(latent$type=="discrete")
					var <- factor(rep(NA,net$size()), levels=latent$levels)
				else{
					var <- as.numeric(rep(NA,net$size()))
					attr(var,"lowerBound") <- latent$lower
					attr(var,"upperBound") <- latent$upper
				}
				net[[latent$name]] <- var
				randomVars[length(randomVars)+1] <- latent$name
			}
			if(!lastTerm)
				tmp <- tmp[[2]]
		}
		form[[3]] <- formula[[3]][[2]]
	}
	clss <- class(net)
	networkEngine <- substring(clss,6,nchar(clss)-3)
	ModelType <- eval(parse(text=paste(networkEngine,"Model",sep="")))
	
	model <- new(ModelType)
	model$setNetwork(net)
	
	tmp <- form[[3]]
	lastTerm <- FALSE
	stats <- list()
	offsets <- list()
	while(!lastTerm){
		ls <- length(stats)
		lo <- length(offsets)
		term <- if(is.symbol(tmp)){
			lastTerm <- TRUE
			term <- tmp
		} else if(as.character(tmp[[1]])=="+"){
			tmp[[3]]
		}else{
			lastTerm <- TRUE
			tmp
		}
		
		name <- if(is.symbol(term)) as.character(term) else as.character(term[[1]])
		args <- NULL
		if(name=="offset" || name=="constraint"){
			term <- term[[2]]
			name <- if(is.symbol(term)) as.character(term) else as.character(term[[1]])
			if(length(term)>1){
				term[[1]] <- as.name("list")
				args <- eval(term,envir=env)
			}else{
				args <- list()
			}
			offsets[[lo+1]] <- args
			names(offsets)[lo+1] <- name
		}else if(name=="mnar"){
			if(!ignoreMnar){
				term <- term[[2]]
				name <- if(is.symbol(term)) as.character(term) else as.character(term[[1]])
				if(length(term)>1){
					term[[1]] <- as.name("list")
					args <- eval(term,envir=env)
				}else{
					args <- list()
				}
				offsets[[lo+1]] <- args
				names(offsets)[lo+1] <- name		
			}
		}else{
			if(length(term)>1){
				term[[1]] <- as.name("list")
				args <- eval(term,envir=env)
			}else{
				args <- list()
			}
			stats[[ls+1]] <- args
			names(stats)[ls+1] <- name
		}
		if(!lastTerm)
			tmp <- tmp[[2]]
	}
	offsets <- rev(offsets)
	stats <- rev(stats)
	if(length(stats)>0)
		for(i in 1:length(stats)){
			t <- try(model$addStatistic(names(stats)[i],stats[[i]]), silent=TRUE)
			if(inherits(t,"try-error")){
				to <- try(model$addOffset(names(offsets)[i],offsets[[i]]), silent=TRUE)
				if(inherits(to,"try-error"))
					stop(t)
			}
		}
	if(length(offsets)>0)
		for(i in 1:length(offsets))
			model$addOffset(names(offsets)[i],offsets[[i]])
	
	model$setRandomVariables(randomVars)
	model$setRandomGraph(!noDyad)
	if(!is.null(theta))
		model$setThetas(theta)
	model
}

#' create a sampler
#' @param formula the model formula
#' @param dyadToggle the method of sampling to use. Defaults to alternating between nodal-tie-dyad and neighborhood toggling.
#' @param vertexToggle the method of vertex attribuate sampling to use.
#' @param nodeSamplingPercentage how often the nodes should be toggled
#' @param ignoreMnar ignore missing not at random offsets
#' @param theta parameter values
#' library(network)
#' data(flo)
#' nflo<-network(flo,directed=FALSE) 
#' nw <- as.BinaryNet(nflo)
#' mod <- createCppModel(nw ~ edges() + triangles(),theta=c(-.5,0))
#' samp <- createCppSampler(nw ~ edges() + triangles(),theta=c(-.5,0))
#' samp$generateSampleStatistics(100,100,100)
createCppSampler <- function(formula, dyadToggle = NULL, dyadArgs=list(), vertexToggle = NULL,
		vertexArgs=list(), nodeSamplingPercentage=0.2, ignoreMnar=TRUE, theta=NULL){
	cppModel <- createCppModel(formula, ignoreMnar=ignoreMnar, theta=theta)
	net <- cppModel$getNetwork()
	clss <- class(net)
	networkEngine <- substring(clss,6,nchar(clss)-3)	
	SamplerClass <- eval(parse(text=paste0(networkEngine,"MetropolisHastings")))
	cppSampler <- new(SamplerClass,cppModel)
	if(!is.null(dyadToggle))
		cppSampler$setDyadToggleType(dyadToggle,dyadArgs)
	if(!is.null(vertexToggle))
		cppSampler$setVertexToggleType(vertexToggle, vertexArgs)
	cppSampler$setDyadProbability(1-nodeSamplingPercentage[1])
	cppSampler
}



#' simulate statistics
#' @param formula the model formula
#' @param theta model parameters
#' @param nodeSamplingPercentage how often the nodes should be toggled
#' @param mcmcBurnIn burn in
#' @param mcmcInterval interval
#' @param mcmcSampleSize sample size
#' @param ignoreMnar ignore missing not at random offsets
#' @param ... additional arguments to createCppSampler
simulateStatistics <- function(formula, theta, 
		nodeSamplingPercentage=0.2,
		mcmcBurnIn=10000, mcmcInterval=100, mcmcSampleSize=100,
		ignoreMnar=TRUE, ...){
	s <- createCppSampler(formula,ignoreMnar=ignoreMnar,theta=theta,...)
	s$generateSampleStatistics(mcmcBurnIn,mcmcInterval,mcmcSampleSize)
}

#'calculate model statistics from a formula
#' @param formula An ernm formula
calculateStatistics <- function(formula){
	createCppModel(formula,clone=FALSE,ignoreMnar=FALSE)$statistics()
}


#' fits an ERNM model
#' @param formula model formula
#' @param modelType the likelihood model type to use
#' @param modelArgs additional arguments for ModelType
#' @param fullToggles a character vector of length 2 indicating the dyad and vertex toggle types for the unconditional simulations
#' @param missingToggles a character vector of length 2 indicating the dyad and vertex toggle types for the conditional simulations
#' @param nodeSamplingPercentage how often are nodal variates toggled
#' @param ... additional parameters for ernmFit
#' library(network)
#' data(flo)
#' nflo<-network(flo,directed=FALSE) 
#' nw <- as.BinaryNet(nflo)
#' set.seed(1)
#' ernm(nw ~ edges() + triangles())
ernm <- function(formula, modelArgs=list(), 
		fullToggles=c("Compound_NodeTieDyad_Neighborhood","DefaultVertex"),
		missingToggles=c("Compound_NodeTieDyadMissing_NeighborhoodMissing",
				"VertexMissing"),
		nodeSamplingPercentage=0.2,
		..., modelType=NULL){

	fullCppSampler <- createCppSampler(formula,
			dyadToggle=fullToggles[1],
			vertexToggle=fullToggles[2],
			nodeSamplingPercentage=nodeSamplingPercentage[1])	
	net <- fullCppSampler$getModel()$getNetwork()
	
	isMissDyads <- sum(net$nMissing(1:net$size()))!=0
	vars <- fullCppSampler$getModel()$getRandomVariables( )
	isMissVars <- any(sapply(vars,function(x)any(is.na(net[[x]]))))
	if(is.null(modelType)){
		if(!isMissVars && !isMissDyads)
			modelType <- FullErnmModel
		else
			modelType <- MissingErnmModel
	}
	if(length(nodeSamplingPercentage)==1)
		nodeSamplingPercentage <- c(nodeSamplingPercentage,nodeSamplingPercentage)
	if(identical(modelType,FullErnmModel)){
		model <- do.call(modelType,c(fullCppSampler,modelArgs))
	}else{
		missCppSampler <- createCppSampler(formula,
				dyadToggle=missingToggles[1],
				vertexToggle=missingToggles[2],
				nodeSamplingPercentage=nodeSamplingPercentage[2])
		if(!isMissVars)
			missCppSampler$setDyadProbability(1)
		if(!isMissDyads)
			missCppSampler$setDyadProbability(0)
		if(!isMissDyads && !isMissVars)
			stop("The model type is for missing data, but none present.")
		model <- do.call(modelType,c(fullCppSampler,missCppSampler,modelArgs))
	}
	fit <- ernmFit(model,...)
	fit$formula <- formula
	fit
}



