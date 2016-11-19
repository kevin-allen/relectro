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

## test rate histo
## test spike trigerred rate histo


