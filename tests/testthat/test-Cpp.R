
context("C++ tests")
library(testthat)

test_that("C++", {
  runErnmCppTests()
  expect_true(TRUE)
})
