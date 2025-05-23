% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ernm-fit.R
\name{ernm_gof}
\alias{ernm_gof}
\title{print}
\usage{
ernm_gof(
  models,
  observed_network = NULL,
  stats_formula,
  style = "histogram",
  scales = "fixed",
  print = TRUE,
  n_sim = 10000,
  burnin = 10000,
  interval = 100
)
}
\arguments{
\item{models}{named list of ernm models to be to be compared (can be length 1}

\item{observed_network}{the observed network}

\item{stats_formula}{the formula for the statistics}

\item{style}{the style of the plot, either 'histogram' or 'boxplot'}

\item{scales}{the scales of the plot, either 'fixed' or 'free'}

\item{print}{whether to print the plot}

\item{n_sim}{the number of simulations to run}

\item{burnin}{the burnin for the MCMC simulation}

\item{interval}{the samplling interval for MCMC simualtion}
}
\value{
A list containing goodness-of-fit plots and simulated statistics
}
\description{
Goodness of fit plot for ERNM models, particularly suited for comparing models
}
