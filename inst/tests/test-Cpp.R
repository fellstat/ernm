# TODO: Add comment
# 
# Author: ianfellows
###############################################################################


library(ernm)
library(network)

context("C++ tests")

test_that("C++",{
	.C("runErnmTests")	
	NULL
})