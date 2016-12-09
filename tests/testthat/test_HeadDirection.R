library(relectro)
library(testthat)
context("test HeadDirection object")
test_that("spike head direction",{
  ## get spike position
  resPerWhlSamples<-400
  nSpikes=4
  spikeTimes<-as.integer(seq(resPerWhlSamples,by=resPerWhlSamples,length=nSpikes))
  hd<-as.numeric(1:10)
  
  ## get spike head direction
  results<-.Call("spike_head_direction_cwrap",
                 hd,
                 length(hd),
                 as.integer(spikeTimes),
                 as.integer(nSpikes),
                 as.integer(resPerWhlSamples),
                 as.integer(0),
                 as.integer(20000),
                 as.integer(1))
  expect_equal(length(spikeTimes),length(results))
  expect_equal(sum(hd[1:length(spikeTimes)]),sum(results))
  
  ## interpolation
  resPerWhlSamples<-400
  nSpikes=4
  spikeTimes<-as.integer(seq(resPerWhlSamples,by=resPerWhlSamples,length=nSpikes)+resPerWhlSamples/2)
  hd<-as.numeric(1:10)
  results<-.Call("spike_head_direction_cwrap",
                 hd,
                 length(hd),
                 as.integer(spikeTimes),
                 as.integer(nSpikes),
                 as.integer(resPerWhlSamples),
                 as.integer(0),
                 as.integer(20000),
                 as.integer(1))
  expect_equal(sum(results),sum(hd[1:length(spikeTimes)]+0.5))
  
  ## wrapping around 360 and 0 degrees
  resPerWhlSamples<-400
  nSpikes=4
  spikeTimes<-as.integer(seq(resPerWhlSamples,by=resPerWhlSamples,length=nSpikes)+resPerWhlSamples/2)
  hd<-as.numeric(c(357:359,1:4))
  results<-.Call("spike_head_direction_cwrap",
                 hd,
                 length(hd),
                 as.integer(spikeTimes),
                 as.integer(nSpikes),
                 as.integer(resPerWhlSamples),
                 as.integer(0),
                 as.integer(20000),
                 as.integer(1))
  expect_equal(results[3],0)
  
  hd<-as.numeric(c(357:359,10:4))
  results<-.Call("spike_head_direction_cwrap",
                 hd,
                 length(hd),
                 as.integer(spikeTimes),
                 as.integer(nSpikes),
                 as.integer(resPerWhlSamples),
                 as.integer(0),
                 as.integer(20000),
                 as.integer(1))
  expect_equal(results[3],4.5)

  hd<-as.numeric(c(3:0,359:357))
  results<-.Call("spike_head_direction_cwrap",
                 hd,
                 length(hd),
                 as.integer(spikeTimes),
                 as.integer(nSpikes),
                 as.integer(resPerWhlSamples),
                 as.integer(0),
                 as.integer(20000),
                 as.integer(1))
  results
  expect_equal(results[4],359.5)
  rm(hd,spikeTimes,resPerWhlSamples,nSpikes,results)
})

test_that("occupancy histograms",{
  ### animal is at all location only once, from 1 to 50 in a 1d vector
  resSamplesPerWhlSample=400
  samplingRateDat=20000
  pt<-new("Positrack")
  hD<-as.numeric(0:359)
  x=rnorm(n=length(hD),mean = 40,sd=5)
  y=rnorm(n=length(hD),mean = 40,sd=5)
  pt<-setPositrack(pt, pxX=x, pxY=y, hd=hD, 
                   resSamplesPerWhlSample=resSamplesPerWhlSample,samplingRateDat = samplingRateDat,pxPerCm = 1)
  hd<-new("HeadDirection")
  results<-.Call("occupancy_histogram_cwrap",
                      hd@nBinHisto,
                      hd@degPerBin,
                      pt@hd,
                      length(pt@hd),
                      pt@resSamplesPerWhlSample/pt@samplingRateDat*1000, ## ms per whl samples
                      as.integer(0),
                      as.integer(400000000),
                      as.integer(1),
                      as.integer(pt@resSamplesPerWhlSample),
                      hd@histoRepetitions+1)
  ## animal was 10 times in each bin, resSamplesPerWhlSample/samplingRateDat*1000*10 = 200
  expect_false(any(results!=resSamplesPerWhlSample/samplingRateDat*1000*(length(hD)/hd@nBinHisto)))
  expect_equal(length(results),hd@nBinHisto)
  rm(x,y,hd,pt,results,samplingRateDat,resSamplesPerWhlSample)
})

test_that("hd rate histograms",{
  resSamplesPerWhlSample=400
  samplingRateDat=20000
  
  pt<-new("Positrack")
  hD<-as.numeric(0:359)
  x=rnorm(n=length(hD),mean = 40,sd=5)
  y=rnorm(n=length(hD),mean = 40,sd=5)
  pt<-setPositrack(pt, pxX=x, pxY=y, hd=hD, 
                   resSamplesPerWhlSample=resSamplesPerWhlSample,samplingRateDat = samplingRateDat,pxPerCm = 1)

  st<-new("SpikeTrain")
  st<-setSpikeTrain(st,res=c(400),clu=c(1),samplingRate=samplingRateDat)
  st<-setIntervals(st,s=c(0),e=c(4000000))
  
  hd<-new("HeadDirection")
  hd@smoothOccupancySd=0
  hd@smoothRateHistoSd=0
  hd<-headDirectionHisto(hd = hd,st=st,pt = pt)
  expect_false(any(hd@occupancy!=resSamplesPerWhlSample/samplingRateDat*1000*(length(hD)/hd@nBinHisto)))
  expect_equal(hd@histo[1],1*1000/(resSamplesPerWhlSample/samplingRateDat*1000*(length(hD)/hd@nBinHisto)))
  
  st<-setSpikeTrain(st,res=c(400,401),clu=c(1,1),samplingRate=samplingRateDat)
  st<-setIntervals(st,s=c(0),e=c(4000000))
  hd<-headDirectionHisto(hd = hd,st=st,pt = pt)
  expect_equal(hd@histo[1],2*1000/(resSamplesPerWhlSample/samplingRateDat*1000*(length(hD)/hd@nBinHisto)))
  expect_true(any(hd@histo[-1]==0))
  
  rm(pt,hD,x,y,st,samplingRateDat,resSamplesPerWhlSample)
})


test_that("spike-triggered head-direction occupancy histograms",{
  resSamplesPerWhlSample=400
  samplingRateDat=20000
  
  pt<-new("Positrack")
  hD<-as.numeric(rep(0:359,50)) # gives us a second per bins
  x=rnorm(n=length(hD),mean = 40,sd=5)
  y=rnorm(n=length(hD),mean = 40,sd=5)
  pt<-setPositrack(pt, pxX=x, pxY=y, hd=hD, 
                   resSamplesPerWhlSample=resSamplesPerWhlSample,samplingRateDat = samplingRateDat,pxPerCm = 1)
  
  st<-new("SpikeTrain")
  spikeTimes=seq(400,400*50,400)
  length(spikeTimes)
  st<-setSpikeTrain(st,res=spikeTimes,clu=rep(1,length(spikeTimes)),samplingRate=samplingRateDat)
  st<-setIntervals(st,s=c(0),e=c(4000000))
  
  hd<-new("HeadDirection")
  hd@smoothOccupancySd=0
  hd@smoothRateHistoSd=0
  minIsiMs=0
  maxIsiMs=1000
  
  results<- .Call("spike_triggered_head_direction_occupancy_histo_cwrap",
                  as.integer(hd@nBinHisto),
                  hd@degPerBin,
                  as.integer(st@cellList),
                  length(st@cellList),
                  pt@hd,
                  length(pt@hd),
                  as.integer(st@res),
                  as.integer(st@clu),
                  st@nSpikes,
                  as.integer(st@startInterval),
                  as.integer(st@endInterval),
                  length(st@startInterval),
                  pt@resSamplesPerWhlSample/pt@samplingRateDat*1000,
                  as.integer(pt@resSamplesPerWhlSample),
                  hd@smoothOccupancySd,
                  hd@smoothRateHistoSd,
                  minIsiMs,
                  maxIsiMs,
                  as.integer(st@samplingRate))
  
  # each data point in the pt object add 20 ms in a bin.
  # it moves at 1 deg per 20ms so their will be 10x20ms per bin for each spikes
  # we have 50 spikes
  expect_equal(max(results),hd@degPerBin*resSamplesPerWhlSample/samplingRateDat*1000*st@nSpikes)
  
  # animal rotates only in one direction so first half of histo not visited
  expect_false(any(results[1:(hd@nBinHisto/2)]!=-1))
  
  # there are 50 hd samples per seconds, so max of 50 degree range for each spike
  # so only 5 bins not at -1
  expect_equal(sum(results!=-1),5)
  
  rm(pt,hD,x,y,st,samplingRateDat,resSamplesPerWhlSample,results)
})


test_that("spike-triggered head-direction rate histograms",{
  resSamplesPerWhlSample=400
  samplingRateDat=20000
  
  pt<-new("Positrack")
  hD<-as.numeric(rep(0:359,50)) # gives us a second per bins
  x=rnorm(n=length(hD),mean = 40,sd=5)
  y=rnorm(n=length(hD),mean = 40,sd=5)
  pt<-setPositrack(pt, pxX=x, pxY=y, hd=hD, 
                   resSamplesPerWhlSample=resSamplesPerWhlSample,samplingRateDat = samplingRateDat,pxPerCm = 1)
  
  st<-new("SpikeTrain")
  spikeTimes=seq(400,400*10,400)
  length(spikeTimes)
  st<-setSpikeTrain(st,res=spikeTimes,clu=rep(1,length(spikeTimes)),samplingRate=samplingRateDat)
  st<-setIntervals(st,s=c(0),e=c(4000000))
  
  hd<-new("HeadDirection")
  hd@smoothOccupancySd=0
  hd@smoothRateHistoSd=0
  minIsiMs=0
  maxIsiMs=1000
  
  results<- .Call("spike_triggered_head_direction_occupancy_histo_cwrap",
                  as.integer(hd@nBinHisto),
                  hd@degPerBin,
                  as.integer(st@cellList),
                  length(st@cellList),
                  pt@hd,
                  length(pt@hd),
                  as.integer(st@res),
                  as.integer(st@clu),
                  st@nSpikes,
                  as.integer(st@startInterval),
                  as.integer(st@endInterval),
                  length(st@startInterval),
                  pt@resSamplesPerWhlSample/pt@samplingRateDat*1000,
                  as.integer(pt@resSamplesPerWhlSample),
                  hd@smoothOccupancySd,
                  hd@smoothRateHistoSd,
                  minIsiMs,
                  maxIsiMs,
                  as.integer(st@samplingRate))
  expect_equal(max(results),hd@degPerBin*resSamplesPerWhlSample/samplingRateDat*1000*st@nSpikes)
  
  hd<-spikeTriggeredHeadDirectionHisto(hd,st,pt,minIsiMs,maxIsiMs)

  ## all 10 spikes are within 10 degrees and within 1 seconds
  ## the peak number of spike in first bin is equal to number of spike pairs n*(n-1)/2
  expect_equal(max(hd@histo),((st@nSpikes*(st@nSpikes-1))/2)/(max(results)/1000))
  
  ## peak should be at bin just after 180 deg
  expect_equal(which.max(hd@histo),length(hd@histo)/2+1)
  
  ## need to test the intervals ##
  
  ## need to test with head direction cell##

  rm(pt,hD,x,y,st,samplingRateDat,resSamplesPerWhlSample,results)
})

