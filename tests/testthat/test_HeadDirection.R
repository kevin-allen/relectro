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
  nSpikes=50
  spikeTimes=seq(400,400*nSpikes,400)
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
  # it moves at 1 deg per 20ms so their will be 10x20ms per bin for each spikes (200 per bin)
  # we have 50 spikes 200x50 = 10000
  expect_equal(max(results),hd@degPerBin*resSamplesPerWhlSample/samplingRateDat*1000*st@nSpikes)
  
  # animal rotates only in one direction so first half of histo not visited
  expect_false(any(results[1:(hd@nBinHisto/2)]!=-1))
  
  # there are 50 hd samples per seconds, so max of 50 degree range for each spike
  # so only 5 bins not at -1
  expect_equal(sum(results!=-1),5)
  
  ## there should be a sum occupancy time of 1000 ms x number of spikes
  expect_equal(sum(results[which(results!=-1.0)]),(maxIsiMs-minIsiMs)*length(spikeTimes))
  
  
  ############################
  # same test with a different direction
  # try a different turning direction 
  pt<-new("Positrack")
  hD<-as.numeric(rep(359:0,50)) # gives us a second per bins
  x=rnorm(n=length(hD),mean = 40,sd=5)
  y=rnorm(n=length(hD),mean = 40,sd=5)
  pt<-setPositrack(pt, pxX=x, pxY=y, hd=hD, 
                   resSamplesPerWhlSample=resSamplesPerWhlSample,samplingRateDat = samplingRateDat,pxPerCm = 1)
  
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
  # it moves at 1 deg per 20ms so their will be 10x20ms per bin for each spikes (200 per bin)
  # we have 50 spikes 200x50 = 10000
  expect_equal(max(results),hd@degPerBin*resSamplesPerWhlSample/samplingRateDat*1000*st@nSpikes)
  
  # animal rotates only in one direction so first half of histo not visited
  # here the bin at 0-10 degree has some time, so we need to go +2
  expect_false(any(results[(hd@nBinHisto/2+2):hd@nBinHisto]!=-1))
  
  # there are 50 hd samples per seconds, so max of 50 degree range for each spike
  # But in the counter clockwise direction we also get value in 0-10 bin, so 6 bins
  expect_equal(sum(results!=-1),6)
  
  ## there should be a sum occupancy time of 1000 ms x number of spikes
  expect_equal(sum(results[which(results!=-1.0)]),(maxIsiMs-minIsiMs)*length(spikeTimes))
  
  ##### test what would happen if the intervals are long and it wraps up
  #### 
  #### one degree per hd datapoint, wrap up will occur every 360 datapoint
  #### test this in the clockwise direction
  pt<-new("Positrack")
  hD<-as.numeric(rep(0:359,50)) # gives us a second per bins
  x=rnorm(n=length(hD),mean = 40,sd=5)
  y=rnorm(n=length(hD),mean = 40,sd=5)
  pt<-setPositrack(pt, pxX=x, pxY=y, hd=hD, 
                   resSamplesPerWhlSample=resSamplesPerWhlSample,samplingRateDat = samplingRateDat,pxPerCm = 1)
  
  
  ### maximum time window without wrapping
  minIsiMs=0
  maxIsiMs=resSamplesPerWhlSample*360/samplingRateDat*1000
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
  
  ### all bins should have the same time
  expect_true(all(results==hd@degPerBin*resSamplesPerWhlSample/samplingRateDat*1000*st@nSpikes))
  expect_equal(sum(results[which(results!=-1.0)]),(maxIsiMs-minIsiMs)*length(spikeTimes))
  
  ### now with some wrapping around 
  minIsiMs=0
  maxIsiMs=resSamplesPerWhlSample*360/samplingRateDat*1000*1.5
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
  
  expect_equal(sum(results[which(results!=-1.0)]),(maxIsiMs-minIsiMs)*length(spikeTimes))
  expect_equal(results[length(results)]/2,results[1])
  
  
  ######################################################
  #### test whether it works fine with intervals   #####
  ######################################################
  
  minIsiMs=0
  maxIsiMs=2000
  st<-new("SpikeTrain")
  nSpikes=100
  spikeTimes=seq(400,400*nSpikes,400)
  st<-setSpikeTrain(st,res=spikeTimes,clu=rep(1,length(spikeTimes)),samplingRate=samplingRateDat)
  
  
  ## this interval should have no effect on occ map
  st<-setIntervals(st,s=c(0),e=c(st@res[nSpikes]+maxIsiMs*samplingRateDat/1000))
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
  
  expect_equal(sum(results[which(results!=-1.0)]),(maxIsiMs-minIsiMs)*length(spikeTimes))
  
  ## 
  minIsiMs=0
  maxIsiMs=1000
  st<-setIntervals(st,s=c(0),e=c(st@res[nSpikes]))
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
  #### this interval should affect only the spikes of the last 1000 ms of recording,
  #### on average, affected spikes only have half of the IsiMs time
  affected_spikes=length(which(st@res>st@res[length(st@res)]-maxIsiMs*samplingRateDat/1000))
  ## the -10 is probably half of a whd sample
  expect_equal((maxIsiMs-minIsiMs)*(length(spikeTimes)-affected_spikes)+(maxIsiMs-minIsiMs)/2*affected_spikes - 10 ,sum(results[which(results!=-1.0)]))
  
  #############################
  ## try with two intervals  ##
  #############################
  minIsiMs=0
  maxIsiMs=2000
  st<-new("SpikeTrain")
  nSpikes=50
  spikeTimes=seq(400,400*nSpikes,400)
  st<-setSpikeTrain(st,res=spikeTimes,clu=rep(1,length(spikeTimes)),samplingRate=samplingRateDat)
  st<-setIntervals(st,s=c(0,st@res[nSpikes]/2+1),
                   e=c(st@res[nSpikes]/2+1,st@res[nSpikes]+maxIsiMs*samplingRateDat/1000))
  
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
  
  expect_equal(sum(results[which(results!=-1.0)]),(maxIsiMs-minIsiMs)*length(spikeTimes))
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
  nSpikes=5
  spikeTimes=seq(400,400*nSpikes,400)
  st<-setSpikeTrain(st,res=spikeTimes,clu=rep(1,length(spikeTimes)),samplingRate=samplingRateDat)
  st<-setIntervals(st,s=c(0),e=c(4000000))
  
  hd<-new("HeadDirection")
  hd@smoothOccupancySd=0
  hd@smoothRateHistoSd=0
  minIsiMs=0
  maxIsiMs=1000
  hd<-spikeTriggeredHeadDirectionHisto(hd = hd,st = st,pt = pt,minIsiMs = minIsiMs,maxIsiMs = maxIsiMs)
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
  
  ## 5 spikes within the same bin, occupancy has 10 datapoint in the same bin * 20 ms per spike, so 1000 ms in the bin
  total_spikes=sum(seq(nSpikes-1,1))
  total_time=max(results)
  expect_equal(max(hd@histo),total_spikes/total_time*1000)
  
  #############################################################
  ##### test the phase difference between pairs of spikes  ####
  #############################################################
  st<-new("SpikeTrain")
  nSpikes=20
  spikeTimes=seq(400,400*20*nSpikes,400*20)
  minIsiMs=0
  maxIsiMs=2000
  st<-setSpikeTrain(st,res=spikeTimes,clu=rep(1,length(spikeTimes)),samplingRate=samplingRateDat)
  st<-setIntervals(st,s=c(0),e=c(4000000))
  hd<-spikeTriggeredHeadDirectionHisto(hd = hd,st = st,pt = pt,minIsiMs = minIsiMs,maxIsiMs = maxIsiMs)
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
  total_time=max(results)
  hd@histo[which(hd@histo==-1.0)]<-NA
  ## there should be 19 spikes at bin 25 degree
  expect_equal(hd@histo  [which(hd@histoDegree==25)]*total_time/1000,19)
  expect_equal(hd@histo  [which(hd@histoDegree==45)]*total_time/1000,18)
  expect_equal(hd@histo  [which(hd@histoDegree==65)]*total_time/1000,17)
  
  rm(pt,hD,x,y,st,samplingRateDat,resSamplesPerWhlSample,results)
})


test_that("spike-triggered head-direction rate crosscorrelation histograms",{
  resSamplesPerWhlSample=400
  samplingRateDat=20000
  
  pt<-new("Positrack")
  hD<-as.numeric(rep(0:359,50)) # gives us a second per bins
  x=rnorm(n=length(hD),mean = 40,sd=5)
  y=rnorm(n=length(hD),mean = 40,sd=5)
  pt<-setPositrack(pt, pxX=x, pxY=y, hd=hD, 
                   resSamplesPerWhlSample=resSamplesPerWhlSample,samplingRateDat = samplingRateDat,pxPerCm = 1)
  
  st<-new("SpikeTrain")
  nSpikes=10 ## only even numbers
  degreeDiff=41 ## each sample of the pt object change 1 degree
  spikeTimes=seq(400,400*nSpikes*degreeDiff,400*degreeDiff)
  st<-setSpikeTrain(st,res=spikeTimes,clu=rep(c(1,2),length(spikeTimes)/2),samplingRate=samplingRateDat)
  st<-setIntervals(st,s=c(0),e=c(4000000))
  
  hd<-new("HeadDirection")
  hd@smoothOccupancySd=0
  hd@smoothRateHistoSd=0
  minIsiMs=0
  maxIsiMs=1000
  hd<-spikeTriggeredHeadDirectionCrossHisto(hd,st,pt,minIsiMs = minIsiMs,maxIsiMs = maxIsiMs)
  #plot(hd@histoDegree,hd@crossHisto,xlim=c(-50,50))
  ## there are 5 spikes of cell 1, so total time in occupancy is 5 sec divided by 5 bins so 1 sec per bin
  ## There should be a bin at 5 Hz and the other 4 at 0.
  expect_equal(hd@crossHisto[which.min(abs(hd@histoDegree-degreeDiff))],5)
  expect_equal(sum(hd@crossHisto[which(hd@crossHisto!=-1.0)]),5)
  expect_equal(length(which(hd@crossHisto!=-1.0)),5)

  
  #####
  ### try the same thing with head direction changes being counter clockwise
  pt<-new("Positrack")
  hD<-as.numeric(rep(359:0,50)) # gives us a second per bins
  x=rnorm(n=length(hD),mean = 40,sd=5)
  y=rnorm(n=length(hD),mean = 40,sd=5)
  pt<-setPositrack(pt, pxX=x, pxY=y, hd=hD, 
                   resSamplesPerWhlSample=resSamplesPerWhlSample,samplingRateDat = samplingRateDat,pxPerCm = 1)
  st<-new("SpikeTrain")
  nSpikes=50 ## only even numbers
  degreeDiff=39 ## each sample of the pt object change 1 degree, with few spikes will fail
  spikeTimes=seq(400,400*nSpikes*degreeDiff,400*degreeDiff)
  st<-setSpikeTrain(st,res=spikeTimes,clu=rep(c(1,2),length(spikeTimes)/2),samplingRate=samplingRateDat)
  st<-setIntervals(st,s=c(0),e=c(4000000))
  
  minIsiMs=0
  maxIsiMs=1000
  hd<-spikeTriggeredHeadDirectionCrossHisto(hd,st,pt,minIsiMs = minIsiMs,maxIsiMs = maxIsiMs)
  
  ## there are 25 spikes of cell 1, so total time in occupancy is 25 sec divided by 5 bins so 5 sec per bin
  ## There should be a 25 spikes in one bin divided by 5 sec so 5 Hz. Other 4 should be 0 Hz.
  expect_equal(hd@crossHisto[which(hd@histoDegree==-35)],5)
  expect_equal(sum(hd@crossHisto[which(hd@crossHisto!=-1.0)]),5)
  expect_equal(length(which(hd@crossHisto!=-1.0)),5)
  
  #####
  ## try with longer time intervals
  st<-new("SpikeTrain")
  nSpikes=50 ## only even numbers
  degreeDiff=10 ## each sample of the pt object change 1 degree, with few spikes will fail
  spikeTimes=seq(400,400*nSpikes*degreeDiff,400*degreeDiff)
  st<-setSpikeTrain(st,res=spikeTimes,clu=rep(c(1,2),length(spikeTimes)/2),samplingRate=samplingRateDat)
  st<-setIntervals(st,s=c(0),e=c(4000000))
  minIsiMs=0
  maxIsiMs=7000
  hd<-spikeTriggeredHeadDirectionCrossHisto(hd,st,pt,minIsiMs = minIsiMs,maxIsiMs = maxIsiMs)
  #plot(hd@histoDegree,hd@crossHisto,type='l') output should look like a saw
  expect_equal(max(hd@crossHisto),5)
  val<-hd@crossHisto[1:(length(hd@crossHisto)/2)]
  expect_equal(sum(val[c(TRUE,FALSE)]),0)
  expect_true(all(val[c(FALSE,TRUE)]!=0))
  
  pt<-new("Positrack")
  hD<-as.numeric(rep(0:359,50)) # gives us a second per bins
  x=rnorm(n=length(hD),mean = 40,sd=5)
  y=rnorm(n=length(hD),mean = 40,sd=5)
  pt<-setPositrack(pt, pxX=x, pxY=y, hd=hD, 
                   resSamplesPerWhlSample=resSamplesPerWhlSample,samplingRateDat = samplingRateDat,pxPerCm = 1)
  
  minIsiMs=0
  maxIsiMs=7200
  hd<-spikeTriggeredHeadDirectionCrossHisto(hd,st,pt,minIsiMs = minIsiMs,maxIsiMs = maxIsiMs)
  #plot(hd@histoDegree,hd@crossHisto,type='l')
  val<-hd@crossHisto[which(hd@histoDegree>0)]
  expect_equal(sum(val[c(TRUE,FALSE)]),0)
  expect_true(all(val[c(FALSE,TRUE)]!=0))
  
})



test_that("bidirectionality score",{
  ## test the function that calculate the bidirectionality score
library(circular)
a<-as.numeric(rvonmises(100, circular(0),15,control.circular=list(units="degrees")))
b<-as.numeric(rvonmises(100, circular(pi/3),15,control.circular=list(units="degrees")))
c<-c(a,b)
d<-hist(c,breaks = seq(0,360,10))

## we want a one at the beginning and end of acf vector
HdHisto<-d$counts
HdBidirectionalityScore(d$counts)

HdBidirectionalityScore<-function(HdHisto){
  ac<-acf.circ(HdHisto) # get a circular autocorrelation
  ac<-c(ac,ac[1]) # add a 1 at the end
  plot(ac)
  peak<-which.peaks(ac,partial=FALSE,decreasing=FALSE) # find peaks that are not at beginning or end 
  if(length(peak)==0) return(0) # if no peak return 0 
  peak<-peak[which.max(ac[peak])] # get highest peak
  
  trough<-which.peaks(ac,partial=FALSE,decreasing=TRUE)
  
  # trough before and after the trough
  trb<-trough[which(trough<peak)]
  tra<-trough[which(trough>peak)]              
  # get the lowest trough before and after the trough
  trb<-trb[which.min(ac[trb])]
  tra<-tra[which.min(ac[tra])]
  # peak - (mean of peak before and after)
  return(ac[peak]-(mean(c(ac[trb],ac[tra]))))
}


# need to find a peak that is after the drop from 1  and before the rise to 1








})
