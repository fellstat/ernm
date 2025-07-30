#' Goodness of fit for ERNM model
#'
#' Goodness of fit in ERNM is done by comparing simulated networks from the ernm model to
#' the observed network. If the observed network is typical of the simulated networks
#' it is considered to be well fit.
#'
#' @param models named list of ernm models to be to be compared (can be length 1
#' @param observed_network the observed network
#' @param stats_formula the formula for the statistics
#' @param style the style of the plot, either 'histogram' or 'boxplot'
#' @param scales the scales of the plot, either 'fixed' or 'free'
#' @param print whether to print the plot
#' @param n_sim the number of simulations to run
#' @param burnin the burnin for the MCMC simulation
#' @param interval the sampling interval for MCMC simulation
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @return A list containing goodness-of-fit plots and simulated statistics
#' @export
#' @description Goodness of fit plot for ERNM models, particularly suited for comparing models
#' @examplesIf FALSE
#' data(samplike)
#' fit_basic <- ernm(samplike ~ edges() + nodeCount("group") + nodeMatch("group") | group)
#' fit_tri <- ernm(samplike ~ edges() + nodeCount("group") + nodeMatch("group") + triangles() | group)
#'
#' # how well is the triangle term fit?
#' gof <- ernm_gof(
#'   list(
#'     basic = fit_basic,
#'     with_triangles = fit_tri
#'   ),
#'   observed_network = samplike,
#'   stats_formula = samplike ~ triangles(),
#'   n_sim = 100
#' )
#'
#' # look at the fit over all edgewise shared partners
#' gof <- ernm_gof(
#'   list(
#'     basic = fit_basic,
#'     with_triangles = fit_tri
#'   ),
#'   style="boxplot",
#'   observed_network = samplike,
#'   stats_formula = samplike ~ esp(1:10),
#'   n_sim = 100
#' )
#'
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
    calc <- function(sim) {
      new_formula <- update(stats_formula, sim ~ .)
      newenv <- new.env(parent = environment(stats_formula))
      assign("sim", sim, envir = newenv)
      environment(new_formula) <- newenv
      ernm::calculateStatistics(new_formula)
    }
    stats <- list()
    # burn in
    model$m$sampler$generateSample(burnin,0,1)
    # run sims
    for(i in seq_len(n_sim)){
      stats[[i]] <- calc(model$m$sampler$generateSample(interval,0,1)[[1]])
    }

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
    newenv <- new.env(parent = environment(stats_formula))
    assign("observed_network", observed_network, envir = newenv)
    environment(new_formula) <- newenv
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
  # TODO: the boxplot should preserve the statistic ordering. doing esp(1:10) yields alphabetic ordering.
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

  # TODO: This should return an "ErnmSummary" object which has a nice s3 print method.

  # Return the simulated statistics as a data frame
  return(list(stat = combined_stats,plots = plots))
}
