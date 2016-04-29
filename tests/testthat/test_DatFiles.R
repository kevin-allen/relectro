library(relectro)
library(testthat)
context("DatFiles")
test_that("DatFiles",{
  
  nChannels<-5
  nSamples<-10
  f<-file("k01.dat","wb")
  writeBin(object=rep(c(1:nChannels),nSamples),con=f,size=2) # 5 channels, 6
  close(f)
  f<-file("k02.dat","wb")
  writeBin(object=rep(c((1:nChannels)),nSamples*2),con=f,size=2) # 5 channels, 6
  close(f)
  f1.size<-file.info("k01.dat")$size
  f2.size<-file.info("k02.dat")$size
  s<-c(f1.size,f2.size)
  
  
  
  ## check size of files and samples
  df<-new("DatFiles")
  df<-datFilesSet(df,fileNames=c("k01.dat","k02.dat"),path="",nChannels=5)
  expect_equal(sum(df@size),sum(s))
  expect_equal(sum(df@samples),nSamples+nSamples*2)
  
  ## check data read 
  ## channel 0 contains only 1s
  d<-datFilesGetOneChannel(df,channelNo=0,firstSample=0,lastSample=sum(df@samples)-1)
  expect_equal(sum(d),nSamples+2*nSamples)
  ## channel 4 contains only 5s
  d<-datFilesGetOneChannel(df,channelNo=4,firstSample=0,lastSample=sum(df@samples)-1)
  expect_equal(sum(d),(nSamples+2*nSamples)*5)
  
  f<-file("k01.dat","wb")
  writeBin(object=rep(c(1:5),each=nChannels),con=f,size=2) # 5 channels, 6
  close(f)
  f<-file("k02.dat","wb")
  writeBin(object=rep(c(6:11),each=nChannels),con=f,size=2) # 5 channels, 6
  close(f)
  
  ## check that it reads the correct samples
  df<-datFilesSet(df,fileNames=c("k01.dat","k02.dat"),path="",nChannels=5)
  d<-datFilesGetOneChannel(df,channelNo=0,firstSample=0,lastSample=sum(df@samples)-1)
  expect_equal(any(diff(d)!=1),FALSE)
  
  ## with different start index
  d<-datFilesGetOneChannel(df,channelNo=0,firstSample=3,lastSample=5)
  expect_equal(d,c(4,5,6))
  
  ## with several indices
  d<-datFilesGetOneChannel(df,channelNo=0,firstSample=c(3,2),lastSample=c(5,4))
  expect_equal(d,c(4,5,6,3,4,5))
  
  ## arrays out of bound, 11 samples, 10 is the last index
  expect_error(datFilesGetOneChannel(df,channelNo=0,firstSample=0,lastSample=11))
  
  file.remove("k01.dat","k02.dat")
  rm(nChannels,nSamples,f,f1.size,f2.size,s,df,d)
})