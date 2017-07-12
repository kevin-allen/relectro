library(relectro)
library(testthat)
context("ElectroProject")
test_that("ElectroProject", {
 
  ## Generate new ElectroProject
  ep <- new("ElectroProject")
          
  ## Do tests
  expect_equal(ep@directory, "")
  rm(ep)

})
