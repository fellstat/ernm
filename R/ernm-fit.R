
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
#' @export
#' @return ernm object
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
	result <- list(theta=lastTheta,converged=converged,iter=iter,info=sampler$info(sample),
			objectiveDiff=llik, maxScaledGradiant=maxGrad,sample=sample,
			likelihoodHistory=likHistory,gradientHistory=gradHistory,trace=trace,m=sampler)
	class(result) <- "ernm"
	result
}

#' print
#' @param x x
#' @param ... unused
#' @export
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
#' @export
#' @method summary ernm
summary.ernm <- function(object,...){
	
	theta <- object$theta
	cv <- solve(object$info)
	se <- sqrt(diag(cv))
	z <- theta/se
	p.value <- 2*pnorm(abs(z),lower.tail=FALSE)
	d <- data.frame(theta,se,z,p.value)
	rownames(d) <- make.unique(names(theta))
	
	# Compute AIC and BIC with latest sample - no bridge sampling for now
	# generate new bigger sample
	if(!is.null(object$m$missSamp)){
	    n_sim <- dim(object$sample$unconditional)[1]
	    samples <- object$m$generateSampleStatistics(10000,100,n_sim*10)$unconditional
	}else{
	    n_sim <- dim(object$sample)[1]
	    samples <- object$m$generateSampleStatistics(10000,100,n_sim*10)
	}
	
	sample_calc <- apply(samples,1,function(x){sum(theta*x)})
	max_term <- max(sample_calc)
	const_approx <- log(mean(exp(sample_calc - max_term))) + max_term
    logLik <- sum(theta*ernm::calculateStatistics(object$formula)) - const_approx
    net <- eval(object$formula[[2]],envir=environment(object$formula))
    n_verts <- net %n% 'n'
    n_dyads <- n_verts*(n_verts-1)*(1 - 0.5*(!is.directed(net)))
    BIC <- -2*logLik + length(theta)*log(n_verts*(length(object$m$randomVariables())!=0) + n_dyads*object$m$randomGraph())
    AIC <- -2*logLik + 2*length(theta)
    BIC <- round(BIC,2)
    AIC <- round(AIC,2)
    
	print(round(d,4))
	cat(paste("\nBIC:",BIC,"AIC:",AIC, "(lower is better)\n"))
	# Return the data frame invisibly
	invisible(d)
}

#' parameter covariance matrix
#' @param object object
#' @param ... unused
#' @export
#' @method vcov ernm
vcov.ernm <- function(object,...){
	solve(object$info)
}


#' plot an ernm object
#' @param  x the object
#' @param ... unused
#' @export
#' @method plot ernm
plot.ernm <- function(x,...){
	plot(x$likelihoodHistory,main="Likelihood convergance",
			ylab="Change in log-likelihood",xlab="iteration")
}

#' print
#' @param models named list of ernm models to be to be compared (can be length 1
#' @param observed_network the observed network
#' @param stats_formula the formula for the statistics
#' @param style the style of the plot, either 'histogram' or 'boxplot'
#' @param scales the scales of the plot, either 'fixed' or 'free'
#' @param print whether to print the plot
#' @param n_sim the number of simulations to run
#' @param burnin the burnin for the MCMC simulation
#' @param interval the samplling interval for MCMC simualtion
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @return A list containing goodness-of-fit plots and simulated statistics
#' @export
#' @description Goodness of fit plot for ERNM models, particularly suited for comparing models
ernm_gof <- function(models,
                observed_network = NULL,
                stats_formula,
                style = "histogram",
                scales = "fixed",
                print = TRUE,
                n_sim = 10000,
                burnin = 10000,
                interval = 100){
        # Helper function to simulate networks and calculate statistics
        calculate_gof_stats <- function(model, name) {
            # Simulate networks
            sims <- model$m$sampler$generateSample(burnin,interval,n_sim)
            
            # Convert simulations to network objects and calculate statistics
            stats <- lapply(sims, function(sim) {
                if(sim$isDirected()){
                    sim <- ernm::as.network.DirectedNet(sim)
                } else {
                    sim <- ernm::as.network.UndirectedNet(sim)
                }
                new_formula <- update(stats_formula, sim ~ .)
                environment(new_formula) <- environment()
                ernm::calculateStatistics(new_formula)
            })
            
            # Combine statistics into a data frame
            stats_df <- as.data.frame(do.call(rbind, stats))
            stats_df$model <- name
            return(stats_df)
        }
        
        # Ensure models is a named list
        if (!is.list(models) || is.null(names(models))) {
            stop("The `models` argument must be a named list of models.")
        }
        
        # Calculate statistics for each model
        all_sim_stats <- do.call(rbind,lapply(names(models), function(name) {
            calculate_gof_stats(models[[name]], name)
        }))
        
        # If observed network is provided, calculate observed statistics
        if (!is.null(observed_network)) {
            new_formula <- update(stats_formula, observed_network ~ .)
            environment(new_formula) <- environment()
            observed_stats <- ernm::calculateStatistics(new_formula)
            observed_stats <- as.data.frame(t(observed_stats))
            observed_stats$model <- "observed"
        } else {
            observed_stats <- NULL
        }
        
        # Combine observed and simulated stats if observed is provided
        combined_stats <- if (!is.null(observed_stats)) {
            rbind(all_sim_stats, observed_stats)
        } else {
            all_sim_stats
        }
        
        # Pivot to long format for plotting
        long_stats <- combined_stats %>%
            tidyr::pivot_longer(
                cols = -.data$model, 
                names_to = "statistic", 
                values_to = "value"
            )
        
        # Calculate means of simulated statistics for plotting
        means <- long_stats %>%
            filter(.data$model != "observed") %>%
            group_by(.data$model, .data$statistic) %>%
            summarize(value = mean(.data$value), .groups = "drop")
        
        # Get unique statistics from the data
        unique_stats <- unique(long_stats$statistic)
        
        # Loop over each statistic and create a plot
        if(style == 'histogram'){
            plots <- list()
            for (stat_name in unique_stats) {
                stat_plot <- ggplot(long_stats %>% filter(.data$model != "observed", .data$statistic == stat_name), aes(x = .data$value, fill = .data$model)) +
                    geom_histogram(aes(y = after_stat(density)),alpha = 0.6, position = 'identity') +
                    geom_vline(
                        data = means %>% filter(.data$model != "observed",.data$statistic == stat_name),
                        aes(xintercept = .data$value, linetype = "Mean"),
                        color = "black", size = 0.8
                    ) +
                    geom_vline(data = long_stats %>% filter(.data$model == "observed", .data$statistic == stat_name) %>% select(.data$value),
                               aes(xintercept = .data$value, linetype = "observed"),
                               color = "red", linewidth = 0.8) +
                    facet_wrap(~.data$model,nrow =length(models), scales = scales) +
                    labs(
                        title = paste("Goodness-of-Fit: Histogram of", stat_name),
                        x = "Value",
                        y = "Frequency",
                        fill = "Model"
                    )
    
                # Save the plot to the list
                plots[[stat_name]] <- stat_plot
                
                # Print the plot if desired
                if(print){
                    print(stat_plot)
                }
            }
        }
        
        if(style == 'boxplot'){
            observed_data <- long_stats %>%
                filter(.data$model == "observed") %>%
                select(-.data$model)
            observed_data <- do.call(rbind,lapply(names(models),function(m){
                observed_data$model <- m
                observed_data
            }))
            
            stat_plot <-  ggplot(long_stats %>% filter(.data$model != "observed"), aes(x = .data$statistic, y = .data$value, fill = .data$model)) +
                geom_boxplot(alpha = 0.6, outlier.shape = NA, show.legend = TRUE) +
                facet_wrap(~.data$model, nrow = length(models), scales = scales) +
                geom_point(
                    data = observed_data,
                    aes(x = .data$statistic, y = .data$value, color = "Observed"),
                    size = 3,
                    show.legend = TRUE
                ) +
                # Define the point color legend
                scale_color_manual(
                    name = "Observation",
                    values = c("Observed" = "red")
                ) +
                # Separate the guides
                guides(
                    fill = guide_legend(order = 1),
                    color = guide_legend(order = 2)
                ) +
                coord_cartesian(ylim = c(0, quantile(long_stats$value,0.98))) +
                # Labels and theme
                labs(
                    title = "Goodness of fit boxplot",
                    x = "Statistic",
                    y = "Value"
                )
            
            # Print the plot if desired
            if(print){
                print(stat_plot)
            }
            plots <- list(stat_plot)
        }

        # Return the simulated statistics as a data frame
        return(list(stat = combined_stats,plots = plots))
    }
