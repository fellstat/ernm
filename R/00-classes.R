# Since we define S4 methods for these classes from cpp - need to declare them
# So that we don't get complaints during roxygen2, needs to be sourced first 
# hence the 00- prefix
setOldClass("Rcpp_DirectedNet")
setOldClass("Rcpp_UndirectedNet")