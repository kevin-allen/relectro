library(relectro)
library(testthat)
context("set intervals")
test_that("set intervals",{
  st<-new("SpikeTrain")
  st<-setSpikeTrain(st,res=c(1,19999),clu=c(1,1),samplingRate=20000)
  expect_equal((sum(st@endInterval-st@startInterval)/st@samplingRate),1)
  st<-setIntervals(st,s=c(0),e=c(40000))
  
  ## spikes all valid
  expect_equal(st@startResIndexc,0)
  expect_equal(st@endResIndexc,1)
  expect_equal((sum(st@endInterval-st@startInterval)/st@samplingRate),2)
  
  ## all invalid spikes, out of bound indices are returned
  st<-setSpikeTrain(st,res=c(10000,20000),clu=c(1,1),samplingRate=20000)
  st<-setIntervals(st,s=c(20000),e=c(40000))
  expect_equal(st@startResIndexc,2)
  expect_equal(st@endResIndexc,2)
  rm(st)
})

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
  
  ## all spikes outside intervals
  st<-setSpikeTrain(st,res=c(20000,20001),clu=c(1,1),samplingRate=20000)
  st<-setIntervals(st,s=c(40000),e=c(50000))
  st<-spikeTimeAutocorrelation(st,binSizeMs = 100,windowSizeMs = 1000,probability = T)
  expect_equal(sum(st@auto),0)
  st<-spikeTimeAutocorrelation(st,binSizeMs = 100,windowSizeMs = 1000,probability = F)
  expect_equal(sum(st@auto),0)
  rm(st)
})


context("Firing rate")
test_that("Sum of spikes in spike-time autocorrelation",{
  st<-new("SpikeTrain")
  st<-setSpikeTrain(st,res=c(0,19999),clu=c(1,1),samplingRate=20000)
  st<-meanFiringRate(st)  
  expect_equal(st@meanFiringRate,2)
  st<-setIntervals(st,s=c(0),e=c(19999))
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
  ## no spikes in intervals
  st<-setSpikeTrain(st,res=c(1,19999),clu=c(1,1),samplingRate=20000)
  st<-setIntervals(st,s=c(20000),e=c(40000))
  st<-meanFiringRate(st)
  expect_equal(st@meanFiringRate,0)
  rm(st)
})

context("Instantaneous firing rate")
test_that("Sum of ifr when know spikes",{
  st<-new("SpikeTrain")
  st<-setSpikeTrain(st,res=c(10000,20000),clu=c(1,1),samplingRate=20000)
  st<-setIntervals(st,s=c(0),e=c(30000))
  st<-ifr(st,windowSizeMs=50, spikeBinMs=1, kernelSdMs=50)
  
  ## sum of ifr should give use 2 spikes.
  expect_equal(sum(st@ifr[1,]/(1000/st@ifrWindowSizeMs)),2,tolerance = .000001)
  st<-setSpikeTrain(st,res=c(10000,30000),clu=c(1,1),samplingRate=20000)
  st<-setIntervals(st,s=c(20000),e=c(40000))
  st<-ifr(st,windowSizeMs=50, spikeBinMs=1, kernelSdMs=50)
  expect_equal(sum(st@ifr[1,]/(1000/st@ifrWindowSizeMs)),1,tolerance = .000001)
  st<-setIntervals(st,s=c(0),e=c(10001))
  st<-ifr(st,windowSizeMs=1, spikeBinMs=1, kernelSdMs=50)
  ## is the kernel centered on the spikes
  expect_equal(sum(st@ifr[1,]/(1000/st@ifrWindowSizeMs)),0.5,tolerance = 0.1)
  
  ## intervals after last spikes
  st<-new("SpikeTrain")
  st<-setSpikeTrain(st,res=c(10000,20000),clu=c(1,1),samplingRate=20000)
  st<-setIntervals(st,s=c(30000),e=c(40001))
  st<-ifr(st,windowSizeMs=50, spikeBinMs=1, kernelSdMs=50)
  expect_equal(sum(st@ifr[1,]),0)
  expect_equal(length(st@ifr[1,]),10)
  
  ## intervals before and after last spikes
  st<-new("SpikeTrain")
  st<-setSpikeTrain(st,res=c(10000,20000),clu=c(1,1),samplingRate=20000)
  st<-setIntervals(st,s=c(0,40000),e=c(30001,50001))
  st<-ifr(st,windowSizeMs=50, spikeBinMs=1, kernelSdMs=50)
  expect_equal(length(st@ifr[1,]),40)
  expect_equal(sum(st@ifr[1,]/(1000/st@ifrWindowSizeMs)),2,tolerance = .000001)
  rm(st)
})

context("Spike-time crosscorrelation to  events")
test_that("Crosscorrelation to events",{
  st<-new("SpikeTrain")
  st<-setSpikeTrain(st,res=c(10000,20000),clu=c(1,1),samplingRate=20000)
  st<-setEvents(st,events=c(10020,20020))
  st<-spikeTimeCrosscorrelationEvents(st,binSizeMs = 1,windowSizeMs = 100,probability = F) 
  expect_equal(sum(st@crossEvents),2)
  df<-spikeTimeCrosscorrelationEventsAsDataFrame(st)
  expect_equal(df[which(df$time==-0.5),"count"],2)
  
  ## all spikes outside the intervals
  st<-setSpikeTrain(st,res=c(10000,20000),clu=c(1,1),samplingRate=20000)
  st<-setIntervals(st,s=c(20000),e=c(40000))
  st<-setEvents(st,events=c(10020,20020))
  st<-spikeTimeCrosscorrelationEvents(st,binSizeMs = 1,windowSizeMs = 100,probability = F) 
  expect_equal(sum(st@crossEvents),0)
  rm(st,df)
})

context("Spike-time crosscorrelation between neurons")
test_that("Crosscorrelation between neurons",{
  st<-new("SpikeTrain")
  st<-setSpikeTrain(st,res=c(20000,20100),clu=c(1,2),samplingRate=20000)
  st<-spikeTimeCrosscorrelation(st,binSizeMs=1,windowSizeMs=50,probability=F)
  expect_equal(length(as.numeric(st@cross)),1*50*2)
  expect_equal(sum(st@cross),1)
  df<-spikeTimeCrosscorrelationAsDataFrame(st)
  expect_equal(df[which(df$time==5.5),"count"],1)

  st<-setSpikeTrain(st,res=c(20000,20100),clu=c(1,2),samplingRate=20000)
  st<-setIntervals(st,s=c(30000),e=c(40000)) # interval after spikes
  st<-spikeTimeCrosscorrelation(st,binSizeMs=1,windowSizeMs=50,probability=F)
  expect_equal(sum(st@cross),0)
  

  st<-setSpikeTrain(st,res=c(20000,20100),clu=c(1,2),samplingRate=20000)
  st<-setIntervals(st,s=c(20000),e=c(40000)) # interval after spikes
  st<-spikeTimeCrosscorrelation(st,binSizeMs=1,windowSizeMs=50,probability=F)
  expect_equal(sum(st@cross),0)

  st<-setSpikeTrain(st,res=c(20000,20100),clu=c(1,2),samplingRate=20000)
  st<-setIntervals(st,s=c(19999),e=c(20101)) # interval after spikes
  st<-spikeTimeCrosscorrelation(st,binSizeMs=1,windowSizeMs=50,probability=F)
  expect_equal(sum(st@cross),1)
  
})

