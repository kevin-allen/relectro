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
  expect_equal(st@startResIndexc,0)
  expect_equal(st@endResIndexc,2)
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


context("Instantaneous firing rate")
test_that("Sum of ifr when know spikes",{
  st<-new("SpikeTrain")
  st<-setSpikeTrain(st,res=c(10000,20000),clu=c(1,1),samplingRate=20000)
  st<-setIntervals(st,s=c(0),e=c(30000))
  st<-ifr(st) 
  ## sum of ifr should give use 2 spikes.
  expect_equal(sum(st@ifr[1,]/(1000/st@ifrWindowSizeMs)),2,tolerance = .000001)
  st<-setSpikeTrain(st,res=c(10000,30000),clu=c(1,1),samplingRate=20000)
  st<-setIntervals(st,s=c(20000),e=c(40000))
  st<-ifr(st) 
  sum(st@ifr[1,]/(1000/st@ifrWindowSizeMs))
  expect_equal(sum(st@ifr[1,]/(1000/st@ifrWindowSizeMs)),1,tolerance = .000001)
  
  st@ifrWindowSizeMs=1
  st<-setIntervals(st,s=c(0),e=c(10001))
  st<-ifr(st)
  plot(st@ifrTime,st@ifr,type='l')
  ## is the kernel centered on the spikes
  expect_equal(sum(st@ifr[1,]/(1000/st@ifrWindowSizeMs)),0.5,tolerance = 0.1)
})