library(relectro)
library(testthat)
context("test firing rate map 2d")
test_that("firing rate map",{
  
  
  ## get spike position
  spikeTimes<-as.integer(c(400,800,1200,1600))
  x<-as.numeric(1:10)
  y<-as.numeric(21:30)
  .Call("spike_position_cwrap",
        x,y,length(x),
        spikeTimes,length(spikeTimes),
                 as.integer(400), # res samples per whl samples
                 as.integer(0), # begin interval 1
                 as.integer(20000), # end interval 1
                 as.integer(1)) # number intervals
  
  
  
  
  
  
  
  
  
  
  
  
  sp<-new("SpatialProperties2d")
  ### animal is at all location only once, from 1 to 50 in a 2d matrix
  pt<-new("Positrack")
  maxx=50
  minx=1
  maxy=50
  miny=1
  x=rep(c(seq(minx,maxx),seq(maxx,minx)),(maxy-miny+1)/2)
  y=rep(seq(miny,maxy),each=maxx-minx+1)
  hd<-(sin(cumsum(rnorm(mean=0, sd=0.3, n=length(x))))+1)/2*360
  pt@defaultXYSmoothing=0
  pt<-setPositrack(pt, pxX=x, pxY=y, hd=hd, 
                   resSamplesPerWhlSample=400,samplingRateDat = 20000,pxPerCm = 1)
  
  
  ### set some spikes
  st<-new("SpikeTrain")
  st<-setSpikeTrain(st=st,res=c(400,800),clu=c(1,1),samplingRate=20000)
  
  st
  sp<-getMapStats(sp,st,pt) 
  
  sp@xSpikes
  sp@ySpikes
  
  
  sp@maps
  
  
  rm(maxx,minx,maxy,miny,x,y,hd,pt,st,sp)
})
