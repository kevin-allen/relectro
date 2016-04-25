library(relectro)
context("SpikeTimeAutocorrelation")

test_that("Sum of spikes a bin",{
  st<-new("SpikeTrain")
  st<-setSpikeTrain(st,res=c(20000,20001),clu=c(1,1),samplingRate=20000)
  st<-spikeTimeAutocorrelation(st,binSizeMs = 100,windowSizeMs = 1000,probability = F)
  expect_equal(sum(st@auto),2) # test that all spikes are considered
  rm(st)
})