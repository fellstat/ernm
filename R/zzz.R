loadModule("ernm",TRUE)
.onLoad <- function(libname, pkgname){
	.C("initStats")
	.C("initToggles")
}

.onUnload <- function(libpath) {}

