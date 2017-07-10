library(relectro)
library(testthat)
context("RecSession")
test_that("RecSession",{
 
  #######################################
  ## set a RecObject                   ##
  #######################################
  rs<-new("RecSession")
  rs<-setRecSession(rs,
                    session="test-31012017-0103",
                path="/data/processing/test/test-31012017-0103",
                samplingRate = 20000,
                nChannels = 49,
                nTrials = 3,
                nElectrodes = 12,
                trialNames=c("test-31012017_01","test-31012017_02","test-31012017_03"),
                channelsTetrode = matrix(0:47,ncol = 4,byrow = TRUE),
                environment=c("sqr70","rest","sqr70"),
                stimulation=c("train_fake","none","train_50"),
                setup=c("A","A","A"),
                environmentFamiliarity=c("fam","fam","fam"),
                electrodeLocation = rep("ca1",12),
                pxPerCm = 10,
                resSamplesPerWhdSample=400)
          
  ## not much testing done here
  expect_equal(rs@session,"test-31012017-0103")
  expect_equal(rs@pxPerCm,10)
  expect_equal(containsEnvironment(rs,environment="sqr70"),TRUE)
  expect_equal(containsEnvironment(rs,environment="rest"),TRUE)
  expect_equal(containsEnvironment(rs,environment="typo"),FALSE)
  expect_equal(rs@stimulation,c("train_fake","none","train_50"))
  expect_equal(containsStimulation(rs,stimulation="train_fake"),TRUE)
  expect_equal(containsStimulation(rs,stimulation="none"),TRUE)
  expect_equal(containsStimulation(rs,stimulation="train_50"),TRUE)
  expect_equal(containsStimulation(rs,stimulation="typo"),FALSE)
  expect_equal(containsElectrodeLocation(rs,location="typo"),FALSE)
  expect_equal(containsElectrodeLocation(rs,location="ca1"),TRUE)
  rm(rs)
})
