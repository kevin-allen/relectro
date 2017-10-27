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
  sp<-new("SpatialProperties2d")
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
  rm(maxx,minx,maxy,miny,x,y,hd,pt,sp)
})

test_that("occupancy3D",{
  ### animal is at all location only once, from 1 to 50 in a 2d matrix
  pt<-new("Positrack")
  maxx=50
  minx=1
  maxy=50
  miny=1
  minhd=1
  maxhd=1
  degPerBin=10
  nHdBins=as.integer(ceiling(360/degPerBin))
  repetitions=2
  x=rep(rep(c(seq(minx,maxx),seq(maxx,minx)),(maxy-miny+1)/2),repetitions)-0.5
  y=rep(rep(seq(miny,maxy),each=maxx-minx+1),repetitions)-0.5
  hd<-rep(1,length(x))

  pt@defaultXYSmoothing=0
  pt<-setPositrack(pt, pxX=x, pxY=y, hd=hd,
                   resSamplesPerWhlSample=400,samplingRateDat = 20000,pxPerCm = 1)
  sp<-new("SpatialProperties2d")
  sp@nRowMap=as.integer(((max(pt@x)+1)/1)+1) # x in R is a row
  sp@nColMap=as.integer(((max(pt@y)+1)/1)+1) # y in R is a col
  results<-.Call("occupancy3D_cwrap",
                 sp@nRowMap,
                 sp@nColMap,
                 nHdBins,
                 1, # cm per bin
                 1, # cm per bin
                 degPerBin,
                 pt@x,
                 pt@y,
                 pt@hd,
                 length(pt@x),
                 pt@resSamplesPerWhlSample/pt@samplingRateDat*1000, ## ms per whl samples
                 as.integer(0),
                 as.integer(length(pt@x)*pt@resSamplesPerWhlSample+500),
                 as.integer(1),
                 as.integer(pt@resSamplesPerWhlSample))
  ## check that the sum of time in occupancy3D array is correct
  expect_equal(sum(results)/length(x),400/20000*1000)


  pt<-new("Positrack")
  maxx=5
  minx=1
  maxy=10
  miny=1
  minhd=21
  maxhd=21
  degPerBin=10
  nHdBins=as.integer(ceiling(360/degPerBin))
  repetitions=2
  x=rep(rep(c(seq(minx,maxx),seq(maxx,minx)),(maxy-miny+1)/2),repetitions)-0.5
  y=head(rep(rep(seq(miny,maxy),each=maxx-minx+1),repetitions)-0.5,length(x))
  hd<-as.numeric(rep(seq(minhd,maxhd),length(x)))
  pt@defaultXYSmoothing=0
  pt<-setPositrack(pt, pxX=x, pxY=y, hd=hd,
                   resSamplesPerWhlSample=400,samplingRateDat = 20000,pxPerCm = 1)
  sp<-new("SpatialProperties2d")
  sp@nRowMap=as.integer(((max(pt@x)+1)/1)+1) # x in R is a row
  sp@nColMap=as.integer(((max(pt@y)+1)/1)+1) # y in R is a col
  results<-.Call("occupancy3D_cwrap",
                 sp@nRowMap,
                 sp@nColMap,
                 nHdBins,
                 1, # cm per bin
                 1, # cm per bin
                 degPerBin,
                 pt@x,
                 pt@y,
                 pt@hd,
                 length(pt@x),
                 pt@resSamplesPerWhlSample/pt@samplingRateDat*1000, ## ms per whl samples
                 as.integer(0),
                 as.integer(length(pt@x)*pt@resSamplesPerWhlSample+500),
                 as.integer(1),
                 as.integer(pt@resSamplesPerWhlSample))
  ## check that data are ordered correctly, here the vector from c function is transformed to 3D array with x, y , hd
  a<-array(results,dim = c(nHdBins,sp@nColMap,sp@nRowMap))
  at<-aperm(a,c(3,2,1))
  expect_equal(sum(at[,,3])/length(x),400/20000*1000)
  expect_equal(sum(at[,,1])/length(x),0)



  #### use occupancyMap3d() directly instead of c function
  pt<-new("Positrack")
  maxx=5
  minx=1
  maxy=10
  miny=1
  minhd=21
  maxhd=21
  degPerBin=10
  nHdBins=as.integer(ceiling(360/degPerBin))
  repetitions=2
  x=rep(rep(c(seq(minx,maxx),seq(maxx,minx)),(maxy-miny+1)/2),repetitions)-0.5
  y=head(rep(rep(seq(miny,maxy),each=maxx-minx+1),repetitions)-0.5,length(x))
  hd<-as.numeric(rep(seq(minhd,maxhd),length(x)))
  pt@defaultXYSmoothing=0
  pt<-setPositrack(pt, pxX=x, pxY=y, hd=hd,
                   resSamplesPerWhlSample=400,samplingRateDat = 20000,pxPerCm = 1)
  sp<-new("SpatialProperties2d")
  st<-new("SpikeTrain")
  ## set the spike trains in the object
  st<-setSpikeTrain(st=st,res=c(1,length(x)*400+401),clu=c(1,1),samplingRate=20000)

  sp<-occupancyMap3d(sp,st,pt)

  expect_equal(sum(sp@occupancy3D[,,-3]),0)
  expect_equal(sum(sp@occupancy3D[,,3])/length(x),400/20000*1000)

  rm(maxx,minx,maxy,miny,minhd,maxhd,a,at,x,y,hd,pt,sp)
})

test_that("distributive_ratio",{
  ### try when the position explain all the apparent HD activity
  pt<-new("Positrack")
  sp<-new("SpatialProperties2d")
  st<-new("SpikeTrain")
  hd<-new("HeadDirection")

  hd@degPerBin=10
  nHdBins=as.integer(ceiling(360/hd@degPerBin))
  ##
  x<-(rep(sin(seq(0,2*pi,0.05)),10)+1)*20
  y<-(rep(sin(seq(0,2*pi,0.05)+pi/2),10)+1)*20
  HD<-(rep(sin(seq(0,2*pi,0.05)),10)+1)*180
  pt@defaultXYSmoothing=0
  pt<-setPositrack(pt, pxX=x, pxY=y, hd=HD,
                   resSamplesPerWhlSample=400,samplingRateDat = 20000,pxPerCm = 1)
    ## set the spike trains in the object
  res<-seq(401,400*length(x),by=400*length(x)/10)
  st<-setSpikeTrain(st=st,res=res,clu=rep(1,length(res)),samplingRate=20000)
  ## get the hd histo
  hd@smoothOccupancySd=0
  hd@smoothRateHistoSd=0
  hd<-headDirectionHisto(hd,st,pt) # observed histo
  ## get the firing rate maps
  sp@smoothOccupancySd=0
  sp@smoothRateMapSd=0
  sp<-firingRateMap2d(sp,st,pt)
  firingRateMapPlot(sp@maps[,,1])
  DR<-distributiveRatioFromHdHisto(sp,st,pt,hd,nRowMap=NA,nColMap=NA)
  ## should be close to 0 because all hd selectivity comes from the correlation between position and HD
  expect_equal(DR,0,tolerance=0.1)


  ### try when the position does not explain HD at all
  pt<-new("Positrack")
  sp<-new("SpatialProperties2d")
  st<-new("SpikeTrain")
  hd<-new("HeadDirection")
  hd@degPerBin=10
  nHdBins=as.integer(ceiling(360/hd@degPerBin))
  ##
  x<-(rep(sin(seq(0,2*pi,0.1)),100)+1)*20
  y<-(rep(sin(seq(0,2*pi,0.1)+pi/2),100)+1)*20
  HD<-(rep(sin(seq(0,2*pi,0.066)),200)+1)*180
  HD<-head(HD,length(x))
  pt@defaultXYSmoothing=0
  pt<-setPositrack(pt, pxX=x, pxY=y, hd=HD,
                   resSamplesPerWhlSample=400,samplingRateDat = 20000,pxPerCm = 1)

  ## set the spike trains in the object
  res<-which(HD>175&HD<185)*400+1 ## fire when HD is near 180
  st<-setSpikeTrain(st=st,res=res,clu=rep(1,length(res)),samplingRate=20000)

  ## get the hd histo
  hd@smoothOccupancySd=0
  hd@smoothRateHistoSd=0
  hd<-headDirectionHisto(hd,st,pt) # observed histo
  ## get the firing rate maps
  sp@smoothOccupancySd=0
  sp@smoothRateMapSd=0
  sp<-firingRateMap2d(sp,st,pt)
  plot(hd@histo,type='l')
  firingRateMapPlot(sp@maps[,,1])
  DR<-distributiveRatioFromHdHisto(sp,st,pt,hd,nRowMap=NA,nColMap=NA)
  ## should be close to 0 because all hd selectivity comes from the correlation between position and HD
  expect_equal(DR,1,tolerance=0.1)

  rm(HD,x,y,hd,pt,sp,st,res,DR)
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
  nSpikes=1
  spikeTimes<-(length(x)*pt@resSamplesPerWhlSample/2)-(maxx/2*pt@resSamplesPerWhlSample)
  st<-setSpikeTrain(st=st,res=spikeTimes,clu=rep(1,nSpikes),samplingRate=20000)
  st<-setIntervals(st,s=0,e=length(x)*pt@resSamplesPerWhlSample+pt@resSamplesPerWhlSample)

  sp@cmPerBin=1
  sp@smoothRateMapSd=0
  sp@smoothOccupancySd=0
  sp<-firingRateMap2d(sp,st,pt)
 # firingRateMapsPlot(maps=sp@maps,names(sp@cellList))
  ## occ maps the animal was 3 times in each bin, resSamplesPerWhlSample/samplingRateDat*1000*3 = 60
  expect_equal(max(sp@occupancy),60)
  ## one spike in 60 ms time
  expect_equal(max(sp@maps),1/60*1000)

  ## make sure the filter does not affect the sum of the firing rates
  sumSmoothZero<-sum(sp@maps[which(sp@maps!=-1.0)])
  sp@cmPerBin=1
  sp@smoothRateMapSd=3
  sp@smoothOccupancySd=3
  sp<-firingRateMap2d(sp,st,pt)
#  firingRateMapsPlot(maps=sp@maps,names(sp@cellList))
  sumSmoothThree<-sum(sp@maps[which(sp@maps!=-1.0)])
  expect_equal(sumSmoothZero,sumSmoothThree)


  ## create a map with a specific size
  sp<-firingRateMap2d(sp,st,pt,nRowMap = 100, nColMap = 101)
  expect_equal(dim(sp@maps),c(100,101,1))
})

test_that("border score, CM and DM in rectangular environments",{

  ### animal is at all location only once, from 1 to 50 in a 2d matrix
  pt<-new("Positrack")
  maxx=50
  minx=1
  maxy=40
  miny=1
  x=rep(rep(c(seq(minx,maxx),seq(maxx,minx)),(maxy-miny+1)/2),3)-0.5
  y=rep(rep(seq(miny,maxy),each=maxx-minx+1),3)-0.5
  hd<-(sin(cumsum(rnorm(mean=0, sd=0.3, n=length(x))))+1)/2*360
  pt@defaultXYSmoothing=0
  pt<-setPositrack(pt, pxX=x, pxY=y, hd=hd,
                   resSamplesPerWhlSample=400,samplingRateDat = 20000,pxPerCm = 1)

  st<-new("SpikeTrain")
  ## set the spike trains in the object
  spikeTimes<- (1:(maxx))*pt@resSamplesPerWhlSample
  st<-setSpikeTrain(st=st,res=spikeTimes,clu=rep(1,length(spikeTimes)),samplingRate=20000)
  st<-setIntervals(st,s=0,e=length(x)*pt@resSamplesPerWhlSample+pt@resSamplesPerWhlSample)

  sp<-new("SpatialProperties2d")
  sp@cmPerBin=1
  sp@smoothRateMapSd=0
  sp@smoothOccupancySd=0
  sp<-firingRateMap2d(sp,st,pt)
  #firingRateMapPlot(m=sp@maps[,,1])


  b<-borderDetection(sp,border="rectangular")

 # firingRateMapPlot(m=b)

  # returned matrix should be same size as occ map
  expect_equal(dim(b),dim(sp@occupancy))

  # number of pixels part of the border should be as follows
  expect_equal(length(b[b!=0]),(maxx-minx)*2+(maxy-miny)*2)

  ## get the bottom column of the map with spikes
  sp<-getMapStats(sp,st,pt,border="rectangular")
  m<-sp@maps[,2,1]
  m<-m[which(!is.na(m))]
  l1<-length(m[which(m!=0.0)])
  CM<-l1/length(m)
  expect_equal(CM,sp@borderCM)


  ## have a cell with firing covering only half of the wall
  expectedCM<-0.5
  spikeTimes<-(1:(maxx*expectedCM))*pt@resSamplesPerWhlSample
  st<-setSpikeTrain(st=st,res=spikeTimes,clu=rep(1,length(spikeTimes)),samplingRate=20000)
  st<-setIntervals(st,s=0,e=length(x)*pt@resSamplesPerWhlSample+pt@resSamplesPerWhlSample)
  sp<-getMapStats(sp,st,pt,border="rectangular")
  expect_equal(sp@borderCM,expectedCM)

  ## expect DM to be 0 because firing is limited to border pixels
  expect_equal(sp@borderDM,0)


  nSpikes=1
  spikeTimes<-(length(x)*pt@resSamplesPerWhlSample/2)-(maxx/2*pt@resSamplesPerWhlSample)
  st<-setSpikeTrain(st=st,res=spikeTimes,clu=rep(1,nSpikes),samplingRate=20000)
  st<-setIntervals(st,s=0,e=length(x)*pt@resSamplesPerWhlSample+pt@resSamplesPerWhlSample)
  sp<-firingRateMap2d(sp,st,pt)
 # firingRateMapPlot(m=sp@maps[,,1])
  sp<-getMapStats(sp,st,pt,border="rectangular")
  expect_equal(sp@borderCM,NaN) # because there is no field detected


  sp@smoothRateMapSd=1.5
  sp<-firingRateMap2d(sp,st,pt)
  #firingRateMapPlot(m=sp@maps[,,1])
  sp<-getMapStats(sp,st,pt,border="rectangular")
  expect_equal(sp@borderCM,0)
  expect_gt(sp@borderDM,0.95) ## all firing rate are in the middle of the map


  rm(maxx,minx,maxy,miny,x,y,st,sp,pt,spikeTimes,expectedCM,l1,CM,m,b,hd)
})

test_that("border score, CM and DM in circular environments",
          {

          ### animal is at all location only once, from 1 to 50 in a 2d matrix
          pt<-new("Positrack")
          maxx=50
          minx=1
          maxy=40
          miny=1
          x=rep(rep(c(seq(minx,maxx),seq(maxx,minx)),(maxy-miny+1)/2),3)-0.5
          y=rep(rep(seq(miny,maxy),each=maxx-minx+1),3)-0.5
          hd<-(sin(cumsum(rnorm(mean=0, sd=0.3, n=length(x))))+1)/2*360
          pt@defaultXYSmoothing=0
          pt<-setPositrack(pt, pxX=x, pxY=y, hd=hd,
          resSamplesPerWhlSample=400,samplingRateDat = 20000,pxPerCm = 1)

          st<-new("SpikeTrain")
          ## set the spike trains in the object
          spikeTimes<- (1:(maxx))*pt@resSamplesPerWhlSample
          st<-setSpikeTrain(st=st,res=spikeTimes,clu=rep(1,length(spikeTimes)),samplingRate=20000)
          st<-setIntervals(st,s=0,e=length(x)*pt@resSamplesPerWhlSample+pt@resSamplesPerWhlSample)

          sp<-new("SpatialProperties2d")
          sp@cmPerBin=1
          sp@smoothRateMapSd=0
          sp@smoothOccupancySd=0
          sp<-firingRateMap2d(sp,st,pt)
          #firingRateMapPlot(m=sp@maps[,,1])

          ## number of bins at the border
          b<-borderDetection(sp,border="circular")
          #firingRateMapPlot(m=b)
          expect_equal(length(b[b!=0]),(maxx-minx)*2+(maxy-miny)*2)

          ## all bins at border so DM = 0
          sp<-getMapStats(sp,st,pt,border="circular")
          expect_equal(sp@borderDM,0)

          ## firing rate now cover one x lenght.
          #firingRateMapPlot(m=sp@maps[,,1])
          m<-sp@maps[,,1]
          expect_equal(length(m[which(m!=-1.0&m!=0)])/length(b[which(b==2)]),sp@borderCM)



          sp@smoothRateMapSd=2
          spikeTimes<- (1:2)*pt@resSamplesPerWhlSample
          st<-setSpikeTrain(st=st,res=spikeTimes,clu=rep(1,length(spikeTimes)),samplingRate=20000)
          st<-setIntervals(st,s=0,e=length(x)*pt@resSamplesPerWhlSample+pt@resSamplesPerWhlSample)
          sp<-firingRateMap2d(sp,st,pt)
          sp<-getMapStats(sp,st,pt,border="circular")
         #firingRateMapPlot(m=sp@maps[,,1])
        #  firingRateMapPlot(m=b)
          expect_equal(sp@mapPolarity,1.0,tolerance=0.01)


          ## let's try with a 45 deg rotated box
          radius=10
          r<- -radius:radius
          x<-rep(unlist(sapply(r,function(x){12+(-(11-(abs(x)))):(11-(abs(x))) })),3)
          y<-rep(unlist(sapply(r,function(x){11+rep(x,length(12+(-(11-(abs(x)))):(11-(abs(x))))) })),3)
          hd<-(sin(cumsum(rnorm(mean=0, sd=0.3, n=length(x))))+1)/2*360
          pt@defaultXYSmoothing=0
          pt<-setPositrack(pt, pxX=x, pxY=y, hd=hd,
                           resSamplesPerWhlSample=400,samplingRateDat = 20000,pxPerCm = 1)

          st<-new("SpikeTrain")
          ## set the spike trains in the object
          spikeTimes<- (1:3)*pt@resSamplesPerWhlSample
          st<-setSpikeTrain(st=st,res=spikeTimes,clu=rep(1,length(spikeTimes)),samplingRate=20000)
          st<-setIntervals(st,s=0,e=length(x)*pt@resSamplesPerWhlSample+pt@resSamplesPerWhlSample)

          sp<-new("SpatialProperties2d")
          sp@cmPerBin=1
          sp@smoothRateMapSd=2
          sp<-firingRateMap2d(sp,st,pt)
          #firingRateMapPlot(m=sp@maps[,,1])

          b<-borderDetection(sp,border="circular")
          #firingRateMapPlot(m=b)
          sp<-getMapStats(sp,st,pt,border="circular")

          ## predicted length of border
          expect_equal((length(r)-2)*2+(2*3),length(b[which(b==2)]))

          rm(maxx,minx,maxy,miny,x,y,hd,st,sp,pt,spikeTimes,b,r,m)
})

test_that("spike triggered firing maps",
          {

            ### animal is at all location only once, from 1 to 50 in a 2d matrix
            pt<-new("Positrack")
            maxx=50
            minx=1
            maxy=40
            miny=1
            x=rep(rep(c(seq(minx,maxx),seq(maxx,minx)),(maxy-miny+1)/2),3)-0.5
            y=rep(rep(seq(miny,maxy),each=maxx-minx+1),3)-0.5
            hd<-(sin(cumsum(rnorm(mean=0, sd=0.3, n=length(x))))+1)/2*360
            pt@defaultXYSmoothing=0
            pt<-setPositrack(pt, pxX=x, pxY=y, hd=hd,
                             resSamplesPerWhlSample=400,samplingRateDat = 20000,pxPerCm = 1)

            st<-new("SpikeTrain")
            ## set the spike trains in the object
            spikeTimes<- (1:(maxx))*pt@resSamplesPerWhlSample
            st<-setSpikeTrain(st=st,res=spikeTimes,clu=rep(1,length(spikeTimes)),samplingRate=20000)
            st<-setIntervals(st,s=0,e=length(x)*pt@resSamplesPerWhlSample+pt@resSamplesPerWhlSample)

            sp<-new("SpatialProperties2d")
            sp@cmPerBin=1
            sp@smoothRateMapSd=0
            sp@smoothOccupancySd=0
            sp<-firingRateMap2d(sp,st,pt)
           # firingRateMapPlot(m=sp@maps[,,1])
            sp<-spikeTriggeredFiringRateMap2d(sp,st,pt,0,60)
            firingRateMapPlot(m=sp@maps[,,1])
            m<-sp@maps[,,1]
            expect_equal(length(which(m!=-1.0)),3) # within 60 ms time window, given whl sampling is 20 ms per sample,
                                                   # only 3 bins should have valid values.

            expect_equal(m[nrow(m)/2+1,ncol(m)/2+1],0) ## this is the middle of the map, and should be at 0 Hz because
                                                       ## distance between spike is one entire bin of the map

            rm(maxx,minx,maxy,miny,x,y,hd,st,sp,pt,spikeTimes,m)
})

test_that("grid score",
          {
            ### animal is at all location only once, from 1 to 50 in a 2d matrix
            pt<-new("Positrack")
            maxx=80
            minx=1
            maxy=80
            miny=1
            x=rep(rep(c(seq(minx,maxx),seq(maxx,minx)),(maxy-miny+1)/2),3)-0.5
            y=rep(rep(seq(miny,maxy),each=maxx-minx+1),3)-0.5
            hd<-(sin(cumsum(rnorm(mean=0, sd=0.3, n=length(x))))+1)/2*360
            pt@defaultXYSmoothing=0
            pt<-setPositrack(pt, pxX=x, pxY=y, hd=hd,
                             resSamplesPerWhlSample=400,samplingRateDat = 20000,pxPerCm = 1)

            ## get the coordinates of the fields of a grid cell
            delta<-seq(0,300,60)/180*pi ## angles of fields relative to center
            x_center<-(minx+(maxx-minx))/2
            y_center<-(miny+(maxy-miny))/2
            radius<-20
            xcoor<-c(x_center+radius*cos(delta),x_center)
            ycoor<-c(y_center+radius*sin(delta),y_center)
            coordinates<-matrix(data=c(xcoor,ycoor),ncol=2)
            t<-apply(coordinates,1,function(v,x,y){which.min(sqrt((v[1]-x)^2+(v[2]-y)^2))},x,y)

            st<-new("SpikeTrain")
            ## set the spike trains in the object
            spikeTimes<- t*pt@resSamplesPerWhlSample
            st<-setSpikeTrain(st=st,res=spikeTimes,clu=rep(1,length(spikeTimes)),samplingRate=20000)
            st<-setIntervals(st,s=0,e=length(x)*pt@resSamplesPerWhlSample+pt@resSamplesPerWhlSample)

            sp<-new("SpatialProperties2d")
            sp@cmPerBin=1
            sp@smoothRateMapSd=3.5
            sp@smoothOccupancySd=3.5
            sp<-firingRateMap2d(sp,st,pt)
            sp<-getMapStats(sp,st,pt)
            expect_gt(sp@gridScore,1.0) # perfect grid should have a high grid score

            ##########################
            ### test grid spacing ####
            expect_lt(abs(sp@gridSpacing-radius),expected=1) # radius is the spacing
            ### try a different radius
            delta<-seq(0,300,60)/180*pi ## angles of fields relative to center
            x_center<-(minx+(maxx-minx))/2
            y_center<-(miny+(maxy-miny))/2
            radius<-25
            xcoor<-c(x_center+radius*cos(delta),x_center)
            ycoor<-c(y_center+radius*sin(delta),y_center)
            coordinates<-matrix(data=c(xcoor,ycoor),ncol=2)
            t<-apply(coordinates,1,function(v,x,y){which.min(sqrt((v[1]-x)^2+(v[2]-y)^2))},x,y)
            ## set the spike trains in the object
            spikeTimes<- t*pt@resSamplesPerWhlSample
            st<-setSpikeTrain(st=st,res=spikeTimes,clu=rep(1,length(spikeTimes)),samplingRate=20000)
            st<-setIntervals(st,s=0,e=length(x)*pt@resSamplesPerWhlSample+pt@resSamplesPerWhlSample)
            sp<-getMapStats(sp,st,pt)
            expect_lt(abs(sp@gridSpacing-radius),expected=1) # radius is the spacing

            ### try a different radius
            delta<-seq(0,300,60)/180*pi ## angles of fields relative to center
            x_center<-(minx+(maxx-minx))/2
            y_center<-(miny+(maxy-miny))/2
            radius<-40
            xcoor<-c(x_center+radius*cos(delta),x_center)
            ycoor<-c(y_center+radius*sin(delta),y_center)
            coordinates<-matrix(data=c(xcoor,ycoor),ncol=2)
            t<-apply(coordinates,1,function(v,x,y){which.min(sqrt((v[1]-x)^2+(v[2]-y)^2))},x,y)
            ## set the spike trains in the object
            spikeTimes<- t*pt@resSamplesPerWhlSample
            st<-setSpikeTrain(st=st,res=spikeTimes,clu=rep(1,length(spikeTimes)),samplingRate=20000)
            st<-setIntervals(st,s=0,e=length(x)*pt@resSamplesPerWhlSample+pt@resSamplesPerWhlSample)
            sp<-getMapStats(sp,st,pt)
            expect_lt(abs(sp@gridSpacing-radius),expected=1) # radius is the spacing


            ### test grid orientation

            rm(maxx,minx,maxy,miny,x,y,hd,st,sp,pt,spikeTimes,
               delta,x_center,y_center,radius,xcoor,ycoor,coordinates,t)
          })
