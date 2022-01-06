


loadModule("ernm",TRUE)
.onLoad <- function(libname, pkgname){
	.Call("initStats")
	.Call("initToggles")
}

.onUnload <- function(libpath) {}



#setwd("/Users/ianfellows/Documents/Eclipse_workspace/")
#library(roxygen2)
#roxygenize('ernm',roclets='rd')
