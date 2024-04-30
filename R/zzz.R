loadModule("ernm",TRUE)
.onLoad <- function(libname, pkgname){
	.Call("_ernm_initStats")
	.Call("_ernm_initToggles")
}
.onUnload <- function(libpath) {}
