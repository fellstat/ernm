# Since we define S4 methods for these classes from cpp - need to declare them
# So that we don't get complaints during roxygen2, needs to be sourced first
# hence the 00- prefix
try(setOldClass("Rcpp_DirectedNet"),silent = TRUE)
try(setOldClass("Rcpp_UndirectedNet"),silent = TRUE)
