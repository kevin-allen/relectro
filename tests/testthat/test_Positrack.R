library(relectro)
library(testthat)
context("Positrack")
test_that("Positrack",{
  pt<-new("Positrack")
  ##########################################
  ## make a fake pt object                ##
  ## constant speed, filling all position ##
  ##########################################
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
  
  
  ## not much testing done here
  expect_equal(mean(pt@x),25.5)
  expect_equal(mean(pt@y),25.5)
  
  ## change pxPerCm
  pt<-setPositrack(pt, pxX=x, pxY=y, hd=hd, 
                   resSamplesPerWhlSample=400,samplingRateDat = 20000,pxPerCm = 1)
  expect_equal(max(pt@x),50)
  expect_equal(max(pt@y),50)
  pt<-setPositrack(pt, pxX=x, pxY=y, hd=hd, 
                   resSamplesPerWhlSample=400,samplingRateDat = 20000,pxPerCm = 2)
  expect_equal(max(pt@x),25)
  expect_equal(max(pt@y),25)
  
  ##############################
  ## test the speed function ###
  ##############################
  ## speed should be of 1 cm per sample, at 50 Hz, so 50 cm/sec
  pt<-setPositrack(pt, pxX=x, pxY=y, hd=hd, 
                   resSamplesPerWhlSample=400,samplingRateDat = 20000,pxPerCm = 1)
  sp1<- .Call("speed_from_whl_cwrap",
              pt@x,
              pt@y,
              length(pt@x),
              1.0,
              pt@samplingRateDat, 
              pt@resSamplesPerWhlSample)
  sp1[which(sp1==-1.0)]<-NA
  expect_equal(mean(sp1,na.rm=T),50)
  expect_equal(mean(sp1,na.rm=T),
               mean(sqrt(diff(pt@x)^2+diff(pt@y)^2)/pt@pxPerCm*(pt@samplingRateDat/pt@resSamplesPerWhlSample))) # speed calculated in R
  expect_equal(length(pt@speed),length(pt@x))
  pt@speed
  expect_equal(pt@speed[which(pt@speed!=-1.0)],sp1[which(!is.na(sp1))])
  
  
  rm(pt,maxx,minx,maxy,miny,x,y,hd)
})