

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
#' @param theta the model parameters.
#' @param modelClass The trailing name of the model class. For now, either Model or ReModel
createCppModel <- function(formula,cloneNet=TRUE,theta=NULL, modelClass="Model"){
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
	ModelType <- eval(parse(text=paste(networkEngine,modelClass,sep="")))
	
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




#'calculate model statistics from a formula
#' @param formula An ernm formula
calculateStatistics <- function(formula){
	createCppModel(formula,clone=FALSE,ignoreMnar=FALSE)$statistics()
}

