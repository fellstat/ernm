

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

#' Creates a C++ representation of an ERNM model
#' @param formula the model formula
#' @param ignoreMnar ignore missing not at random offsets
#' @param cloneNet should the network be cloned
#' @param theta the model parameters.
#' @param modelArgs additional arguments for the model, e.g. tapering parameters
#' @export
#' @return a Model object
#' @examples
#'
#' edge_list <- matrix(numeric(),ncol=2)
#' net <- new(UndirectedNet,edge_list,5)
#'
#' rcpp_model <- createCppModel(net ~ edges(), theta = 0)
#'
#' rcpp_sampler <- new(UndirectedMetropolisHastings, rcpp_model)
#'
#' # Run MCMC to generate 30 networks with burnin=10 and an interval of 20 steps between each network
#' networks <- rcpp_sampler$generateSample(10,20,30)
#'
#' sapply(networks, function(net) net$nEdges()) # number of edges in each network
#'
createCppModel <- function(formula,
                           ignoreMnar=TRUE,
                           cloneNet=TRUE,
                           theta=NULL,
                           modelArgs = list(modelClass='Model')){
	form <- formula
	env <- environment(form)
	net_raw <- eval(form[[2]], envir = env)
	if(network::is.network(net_raw)){
	  # delete the "na" vertex attribute used for node activity - since can be nuisance in CPP
	  net_raw <- network::delete.vertex.attribute(eval(form[[2]], envir = env), "na")
	}
	net <- as.BinaryNet(net_raw)
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
	ModelType <- eval(parse(text=paste(networkEngine,modelArgs$modelClass,sep="")))

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
	if(modelArgs$modelClass == "TaperedModel"){
        if(is.null(modelArgs$centers)){
            stop("centers must be specified in modelArgs for a tapered model")
        }
	    if(is.null(modelArgs$tau)){
	        warning("no tau must be specified in modelArgs for a tapered model \n will set it to 1/(2*mu)")
	        modelArgs$tau <- 1/(2*modelArgs$centers)
	    }
	    model$setTau(modelArgs$tau)
	    model$setCenters(modelArgs$centers)
	    if(!is.null(modelArgs$thetaDependent)){
	        model$thetaDependent(modelArgs$thetaDependent)
	    }
	}
	model
}

#' Create a C++ MCMC sampler
#' @param formula the model formula
#' @param modelArgs additional arguments for the model, e.g. tapering parameters
#' @param dyadArgs list of args for dyad
#' @param dyadToggle the method of sampling to use. Defaults to alternating between nodal-tie-dyad and neighborhood toggling.
#' @param vertexToggle the method of vertex attribute sampling to use.
#' @param vertexArgs list of args for vertex
#' @param nodeSamplingPercentage how often the nodes should be toggled
#' @param ignoreMnar ignore missing not at random offsets
#' @param theta parameter values
#' @param ... additional parameters to be passed to createCppModel
#' @details
#' Available dyad toggles are:
#'
#' 'RandomDyad' : chooses a dyad randomly to toggle
#'
#' 'TieDyad' : chooses a dyad randomly with probability .5 and a random edge with probability .5
#'
#' 'NodeTieDyad' : chooses a random vertex and then chooses an out dyad from that vertex with probability .5 and an out edge with probability .5
#'
#' 'Neighborhood' : The neighborhood proposal starts by choosing a random vertex (a) and then
#' selecting two random neighbors (b and c) Then a random non-a neighbor of b is also selected
#' (d). The neighborhood toggle selects the dyad b- for toggling with probability 50\% and the dyad d-c otherwise.
#'
#' 'Compound_NodeTieDyad_Neighborhood' : chooses 'NodeTieDyad' with probability .5 and 'Neighborhood' otherwise. This is the default proposal.
#'
#' 'Tetrad' : A toggle that preserves network degrees. This is useful when degrees are considered fixed.
#'
#' 'RandomMissingDyad' : RandomDyad, but restricted to missing dyads.
#'
#' 'NodeTieDyadMissing' : NodeTieDyad, but restricted to missing dyads.
#'
#' 'NeighborhoodMissing' : Neighborhood, but restricted to missing dyads.
#'
#' 'Compound_NodeTieDyadMissing_NeighborhoodMissing' : Compound_NodeTieDyad_Neighborhood, but restricted to missing dyads.
#'
#' Available vertex toggles are:
#'
#' 'DefaultVertex' : The default toggle.
#'
#' 'DefaultVertexMissing' : DefaultVertex, but restricted to missing values.
#'
#'
#'
#' @export
#' @return a MetropolisHastings object
#' @examples
#'
#' edge_list <- matrix(numeric(),ncol=2)
#' net <- new(UndirectedNet,edge_list,5)
#' net[["group"]] <- c("a","b","a","b","a")
#'
#' # create a simple ernm model sampler
#' sampler <- createCppSampler(net ~ edges() + nodeCount("group") | group, theta = c(0,0))
#'
#' # generate network statistics for 10 networks with a burn in of 100 steps and 200 steps between them
#' sampler$generateSampleStatistics(100,200,10)
#'
#' # generate a network a further 100 steps later
#' sampler$generateSample(0,100,1)[[1]]
#'
#' # make a sampler using Tie-Dyad proposals for the graph
#' td_sampler <- createCppSampler(net ~ edges() + nodeCount("group") | group,
#'   theta = c(0,0),
#'   dyadToggle = "TieDyad")
#' td_sampler$generateSampleStatistics(100,200,10)
#'
createCppSampler <- function(formula,
                             modelArgs = list(modelClass='Model'),
                             dyadToggle = NULL,
                             dyadArgs=list(),
                             vertexToggle = NULL,
                             vertexArgs=list(),
                             nodeSamplingPercentage=0.2,
                             ignoreMnar=TRUE,
                             theta=NULL,
                             ...){
	cppModel <- createCppModel(formula,
	                           ignoreMnar=ignoreMnar,
	                           theta=theta,
	                           modelArgs=modelArgs,
	                           ...)
	net <- cppModel$getNetwork()
	clss <- class(net)
	networkEngine <- substring(clss,6,nchar(clss)-3)
	SamplerClass <- eval(parse(text=paste0(networkEngine,"MetropolisHastings")))
	cppSampler <- new(SamplerClass)
	cppSampler$setModel(cppModel)
	if(!is.null(dyadToggle))
		cppSampler$setDyadToggleType(dyadToggle,dyadArgs)
	if(!is.null(vertexToggle))
		cppSampler$setVertexToggleType(vertexToggle, vertexArgs)
	cppSampler$setDyadProbability(1-nodeSamplingPercentage[1])
	cppSampler
}



#' Simulate statistics
#'
#' Generates a MCMC chain for an ernm model and returns the sample statistics
#'
#' @param formula the model formula
#' @param theta model parameters
#' @param nodeSamplingPercentage how often the nodes should be toggled
#' @param mcmcBurnIn burn in
#' @param mcmcInterval interval
#' @param mcmcSampleSize sample size
#' @param ignoreMnar ignore missing not at random offsets
#' @param modelArgs additional arguments for the model, e.g. tapering parameters
#' @param ... additional arguments to createCppSampler
#' @export
#' @return a matrix of statistics
#' @examples
#' edge_list <- matrix(numeric(),ncol=2)
#' net <- new(UndirectedNet,edge_list,5)
#' net[["group"]] <- c("a","b","a","b","a")
#'
#' # generate sample statistics from a simple ernm model with positive homophily
#' stats <- simulateStatistics(net ~ edges() + nodeCount("group") + homophily("group") | group,
#'   theta = c(0,0,2),
#'   mcmcBurnIn=10)
#' colMeans(stats)
simulateStatistics <- function(formula,
                               theta,
                               nodeSamplingPercentage=0.2,
                               mcmcBurnIn=10000,
                               mcmcInterval=100,
                               mcmcSampleSize=100,
                               ignoreMnar=TRUE,
                               modelArgs=list(modelClass='Model'),
                               ...){
	s <- createCppSampler(formula,ignoreMnar=ignoreMnar,theta=theta,...)
	s$generateSampleStatistics(mcmcBurnIn,mcmcInterval,mcmcSampleSize)
}

#' Calculate network statistics
#'
#' Calculates all network statistics
#'
#' @param formula An ernm formula
#' @export
#' @return a named vector of statistics
#' @examples
#' data(samplike)
#' calculateStatistics(samplike ~ edges() + nodeCount("group") + nodeMatch("group") + homophily("group") + triangles())
calculateStatistics <- function(formula){
	createCppModel(formula,cloneNet=FALSE,ignoreMnar=FALSE)$statistics()
}


#' Fits an ERNM model
#'
#' ernm() fits exponential family random network models, which is an extension of exponential family random graph models, where
#' nodal covariates may be considered stochastic. Additionally, the function may be used to fit models where the nodal covariates
#' are random, but the graph is fixed (i.e. ALAAMs )
#'
#' @param formula an ernm model formula
#' @param tapered should the model be tapered
#' @param tapering_r the tapering parameter (tau = 1/(tapering_r^2 +5))
#' @param modelArgs additional arguments for the model, e.g. tapering parameters that override the defaults
#' @param nodeSamplingPercentage how often are nodal variates toggled
#' @param modelType either FullErnmModel or MissingErnmModel if NULL will check for missingness
#' @param likelihoodArgs additional arguments for the ernmLikelihood
#' @param fullToggles a character vector of length 2 indicating the dyad and vertex toggle types for the unconditional simulations
#' @param missingToggles a character vector of length 2 indicating the dyad and vertex toggle types for the conditional simulations
#' @param ... additional parameters for ernmFit
#' @export
#' @return a fitted model
#' @examplesIf FALSE
#'
#' data(samplike)
#'
#' # fit a tapered model to the samplike dataset where group is considered random
#' fit <- ernm(samplike ~ edges() + nodeCount("group") + nodeMatch("group") | group)
#' summary(fit)
#'
#' # fit an untapered model. The homophily term is a degeneracy robust version
#' # of nodeMatch, which should be used instead of nodeMatch when tapering is not
#' # present. See Fellows (2012)
#' fit2 <- ernm(samplike ~ edges() + nodeCount("group") + homophily("group") | group, tapered=FALSE)
#' summary(fit2)
#'
#' # standard ergms may be fit within ernm
#' library(network)
#' data(flo)
#' flomarriage <- network(flo,directed=FALSE)
#' fit_flo <- ernm(flomarriage ~ edges() + star(2) + triangles(), tapered=FALSE)
#' summary(fit_flo)
#'
#' # ALAAMs can be fit by specifying that edges are considered fixed using noDyad
#' fit3 <- ernm(samplike ~ nodeCount("group") + nodeMatch("group") | group + noDyad)
#' summary(fit3)
#'
#'
#' @references
#' Fellows, Ian Edward. Exponential family random network models. University of California, Los Angeles, 2012.
#'
ernm <- function(formula,
                 tapered = TRUE,
                 tapering_r = 3,
                 modelArgs=list(),
                 nodeSamplingPercentage=0.2,
                 modelType=NULL,
                 likelihoodArgs = list(),
                 fullToggles=c("Compound_NodeTieDyad_Neighborhood","DefaultVertex"),
                 missingToggles=c("Compound_NodeTieDyadMissing_NeighborhoodMissing","VertexMissing"),
                 ...){

	# if we don't get given things in modelArgs set them to reasonable defaults
	if(tapered){
	    stats <- calculateStatistics(formula)
	    if(is.null(modelArgs$tau)){
            tau <- 1 / (tapering_r^2 * (stats + 5))
            tau[stats<0] <- Inf
            modelArgs$tau <- tau
	    }
	    if(is.null(modelArgs$centers)){
	        modelArgs$centers <- stats
	    }
        modelArgs$modelClass <- "TaperedModel"
	}else{
	    modelArgs$modelClass <- "Model"
	}

    # make the full sampler
    fullCppSampler <- createCppSampler(formula,
                                       dyadToggle=fullToggles[1],
                                       vertexToggle=fullToggles[2],
                                       nodeSamplingPercentage=nodeSamplingPercentage[1],
                                       modelArgs = modelArgs)
    net <- fullCppSampler$getModel()$getNetwork()
    if(net$size() <4){
        stop("ERNM does not currently support networks with fewer than 4 nodes")
    }
	isMissDyads <- sum(net$nMissing(1:net$size()))!=0
	vars <- fullCppSampler$getModel()$getRandomVariables( )
	isMissVars <- any(sapply(vars,function(x)any(is.na(net[[x]]))))
	if(is.null(modelType)){
	    if(!isMissVars && !isMissDyads)
			modelType <- FullErnmModel
	    else{
	        if(tapered==TRUE){
	            stop("tapering is not supported for missing data yet")}
	        modelType <- MissingErnmModel
	    }
	}else{
	    if(!(modelType %in% c('FullErnmModel','MissingErnmModel'))){
	        stop("modelType must be either FullErnmModel or MissingErnmModel")
	    }
	}
	if(length(nodeSamplingPercentage)==1)
		nodeSamplingPercentage <- c(nodeSamplingPercentage,nodeSamplingPercentage)
	if(identical(modelType,FullErnmModel)){
		model <- do.call(modelType,c(fullCppSampler,likelihoodArgs))
	}else{
	    if(modelArgs$modelClass == "TaperedModel"){
	        stop("tapering is not supported for missing data yet")
	    }
		missCppSampler <- createCppSampler(formula,
				dyadToggle=missingToggles[1],
				vertexToggle=missingToggles[2],
				nodeSamplingPercentage=nodeSamplingPercentage[2],
				modelArgs = modelArgs)
		if(!isMissVars)
			missCppSampler$setDyadProbability(1)
		if(!isMissDyads)
			missCppSampler$setDyadProbability(0)
		if(!isMissDyads && !isMissVars)
			stop("The model type is for missing data, but none present.")
		model <- do.call(modelType,c(fullCppSampler,missCppSampler,likelihoodArgs))
	}
	fit <- ernmFit(model,...)
	fit$formula <- formula
	fit
}



