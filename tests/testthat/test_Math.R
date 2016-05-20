library(relectro)
library(testthat)
context("Math")
test_that("Math",{
  x<-matrix(data=c(rep(0,14),1,1),ncol=4)
  cm<-centerOfMass(x)
  expect_equal(cm,c(3.5,4.0))
  
  ## arrays out of bound, 11 samples, 10 is the last index
  expect_error(datFilesGetOneChannel(df,channelNo=0,firstSample=0,lastSample=11))
  
  rm(x)
})