library(relectro)
context("SpikeTimeAutocorrelation")
test_that("Sum of spikes in spike-time autocorrelation",{
  st<-new("SpikeTrain")
  st<-setSpikeTrain(st,res=c(20000,20001),clu=c(1,1),samplingRate=20000)
  st<-spikeTimeAutocorrelation(st,binSizeMs = 100,windowSizeMs = 1000,probability = F)
  expect_equal(sum(st@auto),2) # test that all spikes are considered
  st<-setSpikeTrain(st,res=c(20000,39999),clu=c(1,1),samplingRate=20000)
  st<-spikeTimeAutocorrelation(st,binSizeMs = 100,windowSizeMs = 1000,probability = F)
  expect_equal(sum(st@auto),2) # test that all spikes are considered
  st<-setSpikeTrain(st,res=c(20000,40000),clu=c(1,1),samplingRate=20000)
  st<-spikeTimeAutocorrelation(st,binSizeMs = 100,windowSizeMs = 1000,probability = F)
  expect_equal(sum(st@auto),0) # test that all spikes are considered
  st<-setSpikeTrain(st,res=c(20000,20001),clu=c(1,1),samplingRate=20000)
  st<-spikeTimeAutocorrelation(st,binSizeMs = 100,windowSizeMs = 1000,probability = T)
  expect_equal(sum(st@auto),1) # test that all spikes are considered
  rm(st)
})

context("set intervals")
test_that("set intervals",{
  st<-new("SpikeTrain")
  st<-setSpikeTrain(st,res=c(1,20000),clu=c(1,1),samplingRate=20000)
  expect_equal((sum(st@endInterval-st@startInterval)/st@samplingRate),1)
  st<-setIntervals(st,s=c(0),e=c(40000))
  expect_equal((sum(st@endInterval-st@startInterval)/st@samplingRate),2)
})

context("Firing rate")
test_that("Sum of spikes in spike-time autocorrelation",{
  st<-new("SpikeTrain")
  st<-setSpikeTrain(st,res=c(0,20000),clu=c(1,1),samplingRate=20000)
  st<-meanFiringRate(st)  
  expect_equal(st@meanFiringRate,2)
  st<-setIntervals(st,s=c(0),e=c(20000))
  st<-meanFiringRate(st)  
  expect_equal(st@meanFiringRate,0)
  st<-setSpikeTrain(st,res=c(1,19999),clu=c(1,1),samplingRate=20000)
  st<-setIntervals(st,s=c(0),e=c(20000))
  st<-meanFiringRate(st)  
  expect_equal(st@meanFiringRate,2)
  st<-setSpikeTrain(st,res=c(1,19999),clu=c(1,1),samplingRate=20000)
  st<-setIntervals(st,s=c(0),e=c(40000))
  st<-meanFiringRate(st)
  expect_equal(st@meanFiringRate,1)
})