#' Print ernm object
#' @param x x
#' @param ... unused
#' @return No return value, prints summary
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

#' Summary for ernm object
#' @param object object
#' @param ... unused
#' @param include_AIC whether to include AIC in the summary, will slow down
#' @return a data frame summary of the model
#' @export
#' @method summary ernm
summary.ernm <- function(object,include_AIC = TRUE,...){
  
  theta <- object$theta
  cv <- solve(object$info)
  se <- sqrt(diag(cv))
  z <- theta/se
  p.value <- 2*pnorm(abs(z),lower.tail=FALSE)
  d <- data.frame(theta,se,z,p.value)
  rownames(d) <- make.unique(names(theta))
  
  # Compute AIC and BIC with latest sample - no bridge sampling for now
  # generate new bigger sample
  if(include_AIC){
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
  }

  # Return the data frame invisibly
  invisible(d)
}

#' Parameter covariance matrix
#' @param object object
#' @param ... unused
#' @return covariance matrix
#' @export
#' @method vcov ernm
vcov.ernm <- function(object,...){
  solve(object$info)
}

#' Access ERNM parameters
#' @param object object
#' @param ... unused
#' @return parameter vector
#' @export
#' @method coef ernm
coef.ernm <- function(object,...){
  object$theta
}

#' Plot an ernm object
#' @param  x the object
#' @param ... unused
#' @return No return value, plots the likelihood history
#' @export
#' @method plot ernm
plot.ernm <- function(x,...){
  plot(x$likelihoodHistory,main="Likelihood convergence",
       ylab="Change in log-likelihood",xlab="iteration")
}

#' Goodness of fit for ERNM models
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
