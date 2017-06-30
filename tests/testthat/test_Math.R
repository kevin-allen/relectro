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

test_that("Smoothing degrees",{
  # this is a straight line looping 3 times
  x<-rep(c(seq(0,360,1)),2)
  ## there are some artefacts at the begining and end of the signal
  ## when all kernel does not fit 
  sd=3
  kernelSize=sd*3
  diff<-x-smoothGaussian(x,sd=sd,type = "degrees")
  plot(x)
  lines(smoothGaussian(x,sd=sd,type = "degrees"),col="red")
  lines(diff*100)
  
  ## apart from these artefacts, accept a difference of less than 1 degree
  ## other artifact occurs at the 0-360 transition
  expect_lt(max(diff[kernelSize:(length(diff)-kernelSize)]),1)
  rm(x,sd,kernelSize,diff)
})

test_that("smooth_double_gaussian_circular",{
  x<-seq(0,360,10)
  sd=2
  kernelSize=sd*3
  y<-smoothGaussian(x,sd=sd,type="circular")
  ## the kernel size for smoothing has a size of +-3 sd,
  expect_equal(sum(abs(x-y)>0.01),kernelSize*2)
  expect_equal((y[1]+y[length(y)])/2,180)
  expect_equal(sum(x),sum(y))
  rm(x,y,sd,kernelSize)
})
