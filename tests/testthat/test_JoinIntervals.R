library(relectro)
library(testthat)
context("JoinIntervals")
test_that("timeWithinIntervals",{
  expect_equal(sum(timeWithinIntervals(x=c(0,1,2,3),s=1,e=2)),0)
  expect_equal(timeWithinIntervals(x=c(0,1,2,3),s=1,e=3),c(F,F,T,F))
})