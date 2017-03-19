library(relectro)
library(testthat)
context("RecSession")
test_that("RecSession",{
 
  #######################################
  ## set a RecObject                   ##
  #######################################
  rs<-new("RecSession")
  setRecSession(rs,session="test-31012017-0103",
                path="/data/processing/test/test-31012017-0103",
                samplingRate = 20000,
                nChannels = 49,
                nTrials = 3,
                nElectrodes = 12,
                trialNames=c("test-31012017_01","test-31012017_02","test-31012017_03"),
                channelsTetrode = matrix(0:47,ncol = 4,byrow = TRUE),
                env=c("sqr70","rest","sqr70"),
                stim=c("none","none","none"),
                electrodeLocation = rep("ca1",12))
          
  ## not much testing done here
  expect_equal(mean(pt@x),25.5)
  expect_equal(mean(pt@y),25.5)

  rm(pt,maxx,minx,maxy,miny,x,y,hd,sp1)
})