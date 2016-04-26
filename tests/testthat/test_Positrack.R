library(relectro)
context("Positrack")
test_that("Positrack",{
  pt<-new("Positrack")
  
  ## make a fake pt object
  maxx=50
  minx=1
  maxy=50
  miny=1
  x=rep(seq(minx,maxx),maxy-miny+1)
  y=rep(seq(miny,maxy),each=maxx-minx+1)
  hd<-(sin(cumsum(rnorm(mean=0, sd=0.3, n=length(x))))+1)/2*360
  pt@defaultXYSmoothing=0
  pt<-setPositrack(pt, x=x, y=y, hd=hd, 
               resSamplesPerWhlSample=400,samplingRateDat = 20000,pxPerCm = 1)
  
  
  ## not much testing done here
  expect_equal(mean(pt@x),25.5)
  expect_equal(mean(pt@y),25.5)
  pt<-setPositrack(pt, x=x, y=y, hd=hd, 
                   resSamplesPerWhlSample=400,samplingRateDat = 20000,pxPerCm = 2)
  expect_equal(max(pt@x),25)
  expect_equal(max(pt@y),25)
  
  rm(pt,maxx,minx,maxy,miny,x,y,hd)
})