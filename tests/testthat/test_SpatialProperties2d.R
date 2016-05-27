library(relectro)
library(testthat)
context("test firing rate map 2d")
test_that("spike position",{
  
  
  ## get spike position
  resPerWhlSamples<-400
  nSpikes=4
  spikeTimes<-as.integer(seq(resPerWhlSamples,by=resPerWhlSamples,length=nSpikes))
  x<-as.numeric(1:10)
  y<-as.numeric(21:30)
  results<-.Call("spike_position_cwrap",
        x,y,length(x),
        spikeTimes,length(spikeTimes),
                 as.integer(resPerWhlSamples),
                 as.integer(0), # begin interval 1
                 as.integer(20000), # end interval 1
                 as.integer(1)) # number intervals
  expect_equal(length(spikeTimes),length(results[1,]))
  expect_equal(sum(x[1:length(spikeTimes)]),sum(results[1,]))
  
  
  resPerWhlSamples<-400
  nSpikes=4
  spikeTimes<-as.integer(seq(resPerWhlSamples,by=resPerWhlSamples,length=nSpikes)+resPerWhlSamples/2)
  x<-as.numeric(1:10)
  y<-as.numeric(21:30)
  results<-.Call("spike_position_cwrap",
                 x,y,length(x),
                 spikeTimes,length(spikeTimes),
                 as.integer(resPerWhlSamples),
                 as.integer(0), # begin interval 1
                 as.integer(20000), # end interval 1
                 as.integer(1)) # number intervals
  expect_equal(sum(results[1,]),sum(x[1:length(spikeTimes)]+0.5))
  rm(x,y,resPerWhlSamples,nSpikes,results)
})
test_that("occupancy maps",{
  
  sp<-new("SpatialProperties2d")
  ### animal is at all location only once, from 1 to 50 in a 2d matrix
  pt<-new("Positrack")
  maxx=50
  minx=1
  maxy=50
  miny=1
  x=rep(rep(c(seq(minx,maxx),seq(maxx,minx)),(maxy-miny+1)/2),3)-0.5
  y=rep(rep(seq(miny,maxy),each=maxx-minx+1),3)-0.5
  hd<-(sin(cumsum(rnorm(mean=0, sd=0.3, n=length(x))))+1)/2*360
  pt@defaultXYSmoothing=0
  pt<-setPositrack(pt, pxX=x, pxY=y, hd=hd, 
                   resSamplesPerWhlSample=400,samplingRateDat = 20000,pxPerCm = 1)
  
  
  sp@nRowMap=as.integer(((max(pt@x)+1)/sp@cmPerBin)+1) # x in R is a row
  sp@nColMap=as.integer(((max(pt@y)+1)/sp@cmPerBin)+1) # y in R is a col
  
  results<-.Call("occupancy_map_cwrap",
          sp@nRowMap,
          sp@nColMap,
          1, # cm per bin
          1,
          pt@x,
          pt@y,
          length(pt@x),
          pt@resSamplesPerWhlSample/pt@samplingRateDat*1000, ## ms per whl samples
          as.integer(0),
          as.integer(length(pt@x)*pt@resSamplesPerWhlSample+500),
          as.integer(1),
          as.integer(pt@resSamplesPerWhlSample))
  
  ## the animal was 3 times in each bin, resSamplesPerWhlSample/samplingRateDat*1000*3 = 60
  expect_false(any(results!=60))
  
  rm(maxx,minx,maxy,miny,x,y,hd,pt,st,sp)
})


test_that("firing rate maps",{
  
  sp<-new("SpatialProperties2d")
  ### animal is at all location only once, from 1 to 50 in a 2d matrix
  pt<-new("Positrack")
  maxx=50
  minx=1
  maxy=50
  miny=1
  x=rep(rep(c(seq(minx,maxx),seq(maxx,minx)),(maxy-miny+1)/2),3)-0.5
  y=rep(rep(seq(miny,maxy),each=maxx-minx+1),3)-0.5
  hd<-(sin(cumsum(rnorm(mean=0, sd=0.3, n=length(x))))+1)/2*360
  pt@defaultXYSmoothing=0
  pt<-setPositrack(pt, pxX=x, pxY=y, hd=hd, 
                   resSamplesPerWhlSample=400,samplingRateDat = 20000,pxPerCm = 1)
  
  
  st<-new("SpikeTrain")
  ## set the spike trains in the object
  st<-setSpikeTrain(st=st,res=spikeTimes,clu=rep(1,nSpikes),samplingRate=20000)
  st<-setIntervals(st,s=0,e=length(x)*pt@resSamplesPerWhlSample+pt@resSamplesPerWhlSample)

  sp@cmPerBin=1
  sp@smoothRateMapSd=0
  sp@smoothOccupancySd=0
  sp<-firingRateMap2d(sp,st,pt)   
  #firingRateMapsPlot(maps=sp@maps,names(sp@cellList))
  ## occ maps the animal was 3 times in each bin, resSamplesPerWhlSample/samplingRateDat*1000*3 = 60
  expect_equal(max(sp@occupancy),60)
  ## one spike in 60 ms time 
  expect_equal(max(sp@maps),1/60*1000)
  max(sp@maps)
  
  
  ## make sure the filter does not affect the sum of the firing rates 
  
  sp@cmPerBin=1
  sp@smoothRateMapSd=3
  sp@smoothOccupancySd=3
  nSpikes=maxx
  spikeTimes<-as.integer(seq(resPerWhlSamples,by=resPerWhlSamples,length=nSpikes))
  st<-setSpikeTrain(st=st,res=spikeTimes,clu=rep(1,nSpikes),samplingRate=20000)
  st<-setIntervals(st,s=0,e=length(x)*pt@resSamplesPerWhlSample+pt@resSamplesPerWhlSample)
  
  sp<-firingRateMap2d(sp,st,pt)   
  firingRateMapsPlot(maps=sp@maps,names(sp@cellList))
  getMapStats(sp,st,pt)

  
})