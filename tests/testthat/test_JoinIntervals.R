library(relectro)
library(testthat)
context("JoinIntervals")
test_that("timeWithinIntervals",{
  expect_equal(sum(timeWithinIntervals(x=c(0,1,2,3),s=1,e=2)),0)
  expect_equal(timeWithinIntervals(x=c(0,1,2,3),s=1,e=3),c(F,F,T,F))
  expect_equal(timeWithinIntervals(x=c(0,1,2,3,4,5),s=c(1,4),e=c(2,6)),c(F,F,F,F,F,T))
})

test_that("removeTimeOutsideIntervalsFromTimeStamps",{

  expect_error(removeTimeOutsideIntervalsFromTimeStamps(x=c(1),s=c(2),e=c(4)),"Timestamps outside intervals")
  expect_equal(removeTimeOutsideIntervalsFromTimeStamps(x=c(3,9),s=c(2,6),e=c(4,10)),c(1,5))
  expect_equal(removeTimeOutsideIntervalsFromTimeStamps(x=c(3,7,9),s=c(2,6),e=c(4,10)),c(1,3,5))
  
})
  


test_that("equallySpacedTimePointsWithinIntervals",{
  expect_equal(equallySpacedTimePointsWithinIntervals(min = 10,max=100,by=10,s = c(40),e=c(80)),c(40,50,60,70,80))
  expect_equal(equallySpacedTimePointsWithinIntervals(min = 50,max=100,by=10,s = c(40),e=c(80)),c(50,60,70,80))
  expect_equal(equallySpacedTimePointsWithinIntervals(min = 35,max=119,by=10,s = c(45,100),e=c(80,120)),c(45,55,65,75,105,115))
  expect_equal(equallySpacedTimePointsWithinIntervals(min = 30,max=200,by=10,s = c(40,50),e=c(45,65)),c(40,55,65))
})
