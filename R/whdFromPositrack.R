#' Create whd files for each trial (or .dat file) and for the recording session (all trials together). 
#' 
#' The tracking system called Positrack creates files with the extension .positrack 
#' in which the x, y position of the animal is recorded.
#' The head direction of the animal is also in the .positrack. file.
#' This function create files with .whd extension with x, y, hd is given at fixed time intervals
#' relative to the electrophysiological recording. 
#' During recording, positrack sends ttl pulses to the recording system and the ttl pulses
#' are used to synchronize the position and electrophysiological data. By default, the ttl
#' pulses are on the last channel in the .dat files.
#' This function runs on a single recording session
#' 
#' @param rs RecSession object
#' @param resSamplesPerWhdSample Number of samples in the electrophysiological recording for each whd sample.
#' @param ttlChannel Channel with the ttl signal in the .dat files. By default it is the last channel. Channel numbers start at 0.
#' If you want a different channel for each trial, give a numeric vector with the list of channel.
#' @param maxUpDiffRes Maximum difference between successive up for which position will be interpolated
#' @param overwrite Logical indicating whether to recreate a whd file if one already exists
#' @param checkUpIntegrity Logical indicating whether to check up integrity
#' @param checkPositrackIntegrity Logical indicating whether to check positrack integrity
#' @param minInterEventCor Minimum correlation between intervals of ttl pulses and frame capture
whdFromPositrack<-function(rs,
                           resSamplesPerWhdSample=400,
                           ttlChannel=NA,
                           maxUpDiffRes=4000,
                           overwrite=FALSE,
                           checkUpIntegrity=TRUE,
                           checkPositrackIntegrity=TRUE,
                           minInterEventCor=0.8)
  {
  if(rs@session=="")
    stop(paste("whdFromPositrack, rs@session == \"\""))
  if(length(rs@trialNames)==0)
    stop(paste("whdFromPositrack, rs@trialNames has a length of 0"))
  if(resSamplesPerWhdSample<=0)
    stop(paste("whdFromPositrack, resSamplesPerWhdSample <= 0", resSamplesPerWhdSample))
  if(maxUpDiffRes<=0)
    stop(paste("whdFromPositrack, maxUpDiffRes <= 0", maxUpDiffRes))
  if(rs@nChannels==0)
    stop(paste("whdFromPositrack, rs@nChannels equals 0"))
  if(is.na(rs@samplingRate))
    stop(paste("rs@samplingRate is NA"))
  if(file.exists(paste(rs@fileBase,"whd",sep="."))&overwrite==FALSE)
    return()
  
  df<-new("DatFiles")
  df<-datFilesSet(df,
                  fileNames=paste(rs@trialNames,"dat",sep="."),
                  path=rs@path,
                  nChannels=rs@nChannels)
  
  if(length(ttlChannel)!=length(rs@trialNames) & length(ttlChannel)!=1){
    stop(paste("length of ttlChannel should be 1 or the number of trial in the session"))
  }
  if(length(ttlChannel)==1 & any(is.na(ttlChannel))){
    ttlChannel=rep(rs@nChannels-1,length(rs@trialNames))
  }
  if(length(ttlChannel)==1 & any(!is.na(ttlChannel))){
    ttlChannel=rep(ttlChannel,length(rs@trialNames))
  }
  
  ext="whd"
  mainUp<-vector()
  mainPosi<-data.frame()
  
  # create the whd file for each .dat file
  for(tIndex in 1:length(rs@trialNames)){
    
    print(paste(tIndex, rs@trialNames[tIndex]))
    
    ################################
    ## get the data from dat file ##
    ################################
    print(paste("reading sycn channel",ttlChannel[tIndex],"from",rs@trialStartRes[tIndex],"to",rs@trialEndRes[tIndex]))
    x<-as.numeric(datFilesGetChannels(df,channels=ttlChannel[tIndex],
                                      firstSample = rs@trialStartRes[tIndex],
                                      lastSample = rs@trialEndRes[tIndex]))
    
    up<-detectUps(x) ## detect rising times of ttl pulses
    if(checkUpIntegrity){
    if(checkIntegrityUp(up,samplingRate=rs@samplingRate,datLengthSamples=rs@trialEndRes[tIndex]-rs@trialStartRes[tIndex])!=0)
      stop(paste("check of integrity of up failed"))
    }
    ######################################
    ## get the data from positrack file ##
    ######################################
    fn<-paste(paste(rs@path,rs@trialNames[tIndex],sep='/'),"positrack",sep='.')
    if(!file.exists(fn))
      stop(paste("whdFromPositrack, file missing:",fn))
    print(paste("reading",fn))
    posi<-read.table(fn,header=T)  ## now assumes that there is a header
    if(checkPositrackIntegrity){
    if(checkIntegrityPositrackData(posi)!=0)
      stop(paste("check of integrity of positrack file failed"))
    }
    #############################################
    ### compare the .dat and .positrack data  ###
    #############################################
    lup<-length(up)
    lposi<-length(posi$startProcTime)
    print(paste("Number of ttl pulses:",lup))
    print(paste("Number of frames in positrack file:",lposi))
    
    if(lup!=lposi)  
    { # try to align the frames using jitters
      print(paste("whdFromPositrack,",rs@trialNames[tIndex],"length of up (",lup,") and positrack (",lposi,") differs"))
      if(!is.list(x<-whdAlignedTtlPositrack(up,posi))){
        print(paste("alignment failed"))
        stop()
      }      
      up<-x$up
      posi<-x$posi
    }
    interEventCor<-cor(diff(up),diff(posi$startProcTime))
    print(paste("correlation between interUp and interPosi:",round(interEventCor,4)))
    if(interEventCor<minInterEventCor)
    {
      paste("The correlation between interUp and interPosi is below 0.8:",interEventCor)
      stop("Something is wrong with alignment")
    }
    
    ## the frame is capture before it is received by the computer
    ## up in .dat file is frame processing and not frame capture
    ## Therefore, we remove the capture-to-processing delay from the up values
    delay<-(posi$startProcTime-posi$capTime)*rs@samplingRate/1000
    ## remove the delay from the ttl pulse
    up<- up - delay
    
    ### create a main up and main posi for the entire recording session
    mainUp<-c(mainUp,up+rs@trialStartRes[tIndex])
    mainPosi<-rbind(mainPosi,posi[,c("startProcTime","no","x","y","hd")])
    
    
    ## call c function to make the whl file for this trial
    whd<- .Call("whd_file",
                posi$x,
                posi$y,
                posi$hd,
                as.integer(up),
                length(up),
                rs@trialEndRes[tIndex]-rs@trialStartRes[tIndex],
                resSamplesPerWhdSample,
                maxUpDiffRes)
    
    ## save a whd file in the session directory
    fn<-paste(paste(rs@path,rs@trialNames[tIndex],sep='/'),ext,sep='.')
    print(paste("writing",fn))
    write.table(x = whd,file=fn,quote = FALSE,row.names = FALSE, col.names = FALSE)
  }
  
  if(length(mainUp)!=length(mainPosi$x)){
    print(paste("whdFromPositrack,",rs@session,"length of mainUp (",length(mainUp),") and mainPosi (",length(mainPosi$x),") differs"))
    stop()      
  }
  
  ## make the main whd file
  whd<- .Call("whd_file",
              mainPosi$x,
              mainPosi$y,
              mainPosi$hd,
              as.integer(mainUp),
              length(mainUp),
              rs@trialEndRes[length(rs@trialEndRes)],
              resSamplesPerWhdSample,
              maxUpDiffRes)
  
  ## save a whd file in the session directory
  fn<-paste(paste(rs@path,rs@session,sep='/'),ext,sep='.')
  print(paste("writing",fn))
  write.table(x = whd,file=fn,quote = FALSE,row.names = FALSE, col.names = FALSE)
}



#' Create whd files for the entire session. Used when one positrack file covers several dat files. 
#' 
#' The tracking system called Positrack creates files with the extension .positrack 
#' in which the x, y position and head direction of the animal is recorded.
#' This function create files with .whd extension with x, y, hd is given at fixed time intervals
#' relative to the electrophysiological recording. 
#' During recording, positrack sends ttl pulses to the recording system and the ttl pulses
#' are used to synchronize the position and electrophysiological data. By default, the ttl
#' pulses are on the last channel in the .dat files.
#' This function runs on a single recording session
#' 
#' @param rs RecSession object
#' @param resSamplesPerWhdSample Number of samples in the electrophysiological recording for each whd sample.
#' @param ttlChannel Channel with the ttl signal in the .dat files. By default it is the last channel. Channel numbers start at 0.
#' If you want a different channel for each trial, give a numeric vector with the list of channel.
#' @param maxUpDiffRes Maximum difference between successive up for which position will be interpolated
#' @param overwrite Logical indicating whether to recreate a whd file if one already exists
#' @param checkUpIntegrity Logical indicating whether to check up integrity
#' @param checkPositrackIntegrity Logical indicating whether to check positrack integrity
#' @param minInterEventCor Minimum correlation between intervals of ttl pulses and frame capture
#' @param positrackDatMatchFile A file that describes how many dat file there is for each positrack file. Each line represent a positrack file, with consecutive indices (_01, _02, etc.). Each line contains a number indicating how many .dat files there is for the .positrack file. By default, the name is the session filebase.positrackDatMatch
whdFromLongPositrackFile<-function(rs,
                           resSamplesPerWhdSample=400,
                           ttlChannel=NA,
                           maxUpDiffRes=4000,
                           overwrite=FALSE,
                           checkUpIntegrity=TRUE,
                           checkPositrackIntegrity=TRUE,
                           minInterEventCor=0.8,
                           positrackDatMatchFile=NA)
{
  if(rs@session=="")
    stop(paste("whdFromPositrack, rs@session == \"\""))
  if(length(rs@trialNames)==0)
    stop(paste("whdFromPositrack, rs@trialNames has a length of 0"))
  if(resSamplesPerWhdSample<=0)
    stop(paste("whdFromPositrack, resSamplesPerWhdSample <= 0", resSamplesPerWhdSample))
  if(maxUpDiffRes<=0)
    stop(paste("whdFromPositrack, maxUpDiffRes <= 0", maxUpDiffRes))
  if(rs@nChannels==0)
    stop(paste("whdFromPositrack, rs@nChannels equals 0"))
  if(is.na(rs@samplingRate))
    stop(paste("rs@samplingRate is NA"))
  if(file.exists(paste(rs@fileBase,"whd",sep="."))&overwrite==FALSE)
    return()

  if(is.na(positrackDatMatchFile))
    positrackDatMatchFile<-paste(rs@fileBase,"positrackDatMatch",sep=".")
  if(!file.exists(positrackDatMatchFile))
    stop(paste(positrackDatMatchFile,"does not exist"))
  
  df<-new("DatFiles")
  df<-datFilesSet(df,
                  fileNames=paste(rs@trialNames,"dat",sep="."),
                  path=rs@path,
                  nChannels=rs@nChannels)
  
  if(length(ttlChannel)!=length(rs@trialNames) & length(ttlChannel)!=1){
    stop(paste("length of ttlChannel should be 1 or the number of trial in the session"))
  }
  if(length(ttlChannel)==1 & any(is.na(ttlChannel))){
    ttlChannel=rep(rs@nChannels-1,length(rs@trialNames))
  }
  if(length(ttlChannel)==1 & any(!is.na(ttlChannel))){
    ttlChannel=rep(ttlChannel,length(rs@trialNames))
  }
  
  ext="whd"
  mainUp<-vector()
  mainPosi<-data.frame()
  
  
  ## establish the matching of positrack and dat files
  ## should be done in a file and not here
  ## each line is a positrack file, enters the number of dat files covering the positrack file
  numDatFilesPerPositrack<-read.table(positrackDatMatchFile)$V1
  numDatFilesPerPositrack
  previousDatIndex=0
  # create the whd file for each .dat file
  for(tIndex in 1:length(numDatFilesPerPositrack)){
  
    positrackFileName<-paste(paste(rs@path,rs@trialNames[tIndex],sep="/"),"positrack",sep=".")
    firstDatIndex=previousDatIndex+1
    lastDatIndex=previousDatIndex+numDatFilesPerPositrack[tIndex]
    print(paste(positrackFileName, ", dat indices:",firstDatIndex,lastDatIndex))
    
    #################################
    ## get the data from dat files ##
    #################################
    print(paste("reading sycn channel",ttlChannel[tIndex],"from",rs@trialStartRes[firstDatIndex],"to",rs@trialEndRes[lastDatIndex]))
    x<-as.numeric(datFilesGetChannels(df,channels=ttlChannel[tIndex],
                                      firstSample = rs@trialStartRes[firstDatIndex],
                                      lastSample = rs@trialEndRes[lastDatIndex]))
    
    up<-detectUps(x) ## detect rising times of ttl pulses
    if(checkUpIntegrity){
      if(checkIntegrityUp(up,samplingRate=rs@samplingRate,datLengthSamples=rs@trialEndRes[lastDatIndex]-rs@trialStartRes[firstDatIndex])!=0)
        stop(paste("check of integrity of up failed"))
    }
    ######################################
    ## get the data from positrack file ##
    ######################################
    
    if(!file.exists(positrackFileName))
      stop(paste("whdFromPositrack, file missing:",positrackFileName))
    print(paste("reading",positrackFileName))
    posi<-read.table(positrackFileName,header=T)  ## now assumes that there is a header
    if(checkPositrackIntegrity){
      if(checkIntegrityPositrackData(posi)!=0)
        stop(paste("check of integrity of positrack file failed"))
    }
    #############################################
    ### compare the .dat and .positrack data  ###
    #############################################
    lup<-length(up)
    lposi<-length(posi$startProcTime)
    print(paste("Number of ttl pulses:",lup))
    print(paste("Number of frames in positrack file:",lposi))
    
    if(lup!=lposi)  
    { # try to align the frames using jitters
      print(paste("whdFromPositrack,",rs@trialNames[tIndex],"length of up (",lup,") and positrack (",lposi,") differs"))
      if(!is.list(x<-whdAlignedTtlPositrack(up,posi))){
        print(paste("alignment failed"))
        stop()
      }      
      up<-x$up
      posi<-x$posi
    }
    interEventCor<-cor(diff(up),diff(posi$startProcTime))
    print(paste("correlation between interUp and interPosi:",round(interEventCor,4)))
    if(interEventCor<minInterEventCor)
    {
      paste("The correlation between interUp and interPosi is below 0.8:",interEventCor)
      stop("Something is wrong with alignment")
    }
    
    ## the frame is capture before it is received by the computer
    ## up in .dat file is frame processing and not frame capture
    ## Therefore, we remove the capture-to-processing delay from the up values
    delay<-(posi$startProcTime-posi$capTime)*rs@samplingRate/1000
    ## remove the delay from the ttl pulse
    up<- up - delay
    
    ### create a main up and main posi for the entire recording session
    mainUp<-c(mainUp,up+rs@trialStartRes[tIndex])
    mainPosi<-rbind(mainPosi,posi[,c("startProcTime","no","x","y","hd")])
    
    previousDatIndex=lastDatIndex
  }
  
  if(length(mainUp)!=length(mainPosi$x)){
    print(paste("whdFromPositrack,",rs@session,"length of mainUp (",length(mainUp),") and mainPosi (",length(mainPosi$x),") differs"))
    stop()      
  }
  
  ## make the main whd file
  whd<- .Call("whd_file",
              mainPosi$x,
              mainPosi$y,
              mainPosi$hd,
              as.integer(mainUp),
              length(mainUp),
              rs@trialEndRes[length(rs@trialEndRes)],
              resSamplesPerWhdSample,
              maxUpDiffRes)
  
  ## save a whd file in the session directory
  fn<-paste(paste(rs@path,rs@session,sep='/'),ext,sep='.')
  print(paste("writing",fn))
  write.table(x = whd,file=fn,quote = FALSE,row.names = FALSE, col.names = FALSE)
}



#' Test alignment of one dat and positrack file. 
#' 
#' The tracking system called Positrack creates files with the extension .positrack 
#' in which the x, y position of the animal is recorded.
#' Positrack also sends ttl pulses to the electrophysiological system that are then saved into a .dat file 
#' This function test that the number of ttl pulses and the number of recorded video frames are the same
#' It can be used to test for synchronization after recording every trial
#' 
#' @param datFileName Name of the dat file
#' @param positrackFileName Name of the positrack file
#' @param numberChannelsDat Number of channels in the dat file
#' @param ttlChannel Channel with the ttl signal in the .dat files.
#' @param datSamplingRate Sampling rate for the .dat file.
positrackDatAlignmentCheck<-function(datFileName,
                           positrackFileName,
                           numberChannelsDat=NA,
                           ttlChannel=NA,
                           datSamplingRate=20000)
{

  # for testing  
  # positrackFileName="/home/kevin/test-16032017_01.positrack"
  # datFileName="/home/kevin/test-16032017_01.dat"
  # ttlChannel=48
  # numberChannelsDat=49
  # datSamplingRate=20000
  
  
  print(paste("dat file:",datFileName))
  print(paste("positrack file:",positrackFileName))
  print(paste("numberChannelsDat:",numberChannelsDat))
  print(paste("ttlChannel:",ttlChannel))
  print(paste("datSamplingRate:",datSamplingRate))
  
  if(!file.exists(datFileName)){
    stop(paste(datFileName, "is missing"))
  }
  if(!file.exists(positrackFileName)){
    stop(paste(positrackFileName, "is missing"))
  }
  
  if(length(unlist(strsplit(datFileName,split = "/")))>1){
   dfn<-tail(unlist(strsplit(datFileName,split = "/")),n=1)
   dfpath<-paste(unlist(strsplit(datFileName,split = "/"))[1:(length(unlist(strsplit(datFileName,split = "/")))-1)],collapse = "/")
  } else {
    dfn<-datFileName
    dfpath<-""
  }
  
  df<-new("DatFiles")
  df<-datFilesSet(df,fileNames=dfn,path=dfpath,nChannels=numberChannelsDat)
  
  print(paste("reading sycn channel",ttlChannel,"from",df@fileNames))
  x<-as.numeric(datFilesGetChannels(df,channels=ttlChannel,firstSample = 0,lastSample = df@samples-1))
  up<-detectUps(x) ## detect rising times of ttl pulses
  
  print("Check up integrity")
  if(checkIntegrityUp(up,samplingRate=datSamplingRate,datLengthSamples=df@samples-1)!=0)
    stop(paste("check of integrity of up failed"))
    
  ######################################
  ## get the data from positrack file ##
  ######################################
  if(!file.exists(positrackFileName))
      stop(paste("file missing:",positrackFileName))
  posi<-read.table(positrackFileName,header=T)  ## now assumes that there is a header
  print("Check positrack integrity")
  if(checkIntegrityPositrackData(posi)!=0)
      stop(paste("check of integrity of positrack file failed"))
    
  #############################################
  ### compare the .dat and .positrack data  ###
  #############################################
  lup<-length(up)
  lposi<-length(posi$startProcTime)
  print(paste("Number of ttl pulses:",lup))
  print(paste("Number of frames in positrack file:",lposi))

  if(lup==lposi){
  interEventCor<-cor(diff(up),diff(posi$startProcTime))
  print(paste("correlation between interUp and interPosi:",round(interEventCor,4)))
  plot(diff(up)/20,diff(posi$startProcTime))
  print("Difference between inter ttl and inter startProcTime")
  print(summary(diff(up)/20-diff(posi$startProcTime)))
  } else{
    interEventCor=0
  }
    
  if(lup!=lposi|interEventCor<0.70)  
  { # try to align the frames using jitters
    print(paste("length of up (",lup,") and positrack (",lposi,") differs"))
    if(!is.list(x<-whdAlignedTtlPositrack(up,posi))){
      print(paste("alignment failed"))
      stop()
    }      
  }
  print("Sychronization test passed")
}


#' Try to realign the up part of ttl pulses from .dat file and the positrack frames for a single trial. 
#' 
#' One assumption that is made is that there is variability in the intervals between frames in the positrack software
#' These time intervals can be measured from the positrack file (startProcTime) and the ttl pulses.
#' When there is no problem of alignment, the correlation between intervals of ttl and positrack is above .97
#' If there is an additional ttl detected, then the correlation goes down. 
#' 
#' The intervals between ups is calculated (dup)
#' The intervals between the posi$startProcTime is calculated (dposi)
#' We search for points where there is a large difference between dup and dposi.
#' The crosscorrelation function for the data points just after the large difference is calculated. 
#' If the peak is before or after 0, then an up or a posi is removed.
#' 
#' @param up Time of up phase of ttl pulses
#' @param posi Data from the positrack file
#' @param ccLength Lenght of the segment of data to use for crosscorrelation when a delay mismach is found
whdAlignedTtlPositrack<-function(up,posi,ccLength=20){
  dup<-diff(up)
  dposi<-diff(posi$startProcTime*20)
  minLength<-min(c(length(dup),length(dposi)))
  cor1=cor(head(dup,n=500),head(dposi,n=500))
  cor2=cor(dup[(minLength-500):minLength],dposi[(minLength-500):minLength])
  print(paste("correlation first 500:",round(cor1,3)))
  print(paste("correlation last 500:",round(cor2,3)))
  print(paste("length of up:",length(up)))
  print(paste("length of posi:",length(posi$no)))
  if(length(up)==length(posi$no)&cor1>0.95&cor2>0.95){
    print("Alignment appears ok")
    return(list(up=up,posi=posi))
  }
  if(abs(length(up)-length(posi$no))>10){
    print("The alignment problem is for more than 10 frames, no solution implemented for this yet")
    return(NA)
  }
  
  removedPosi=0
  removedUp=0
  latestRemoved=0
  attempted=0
  print(paste("Try to remove ttl or posi data points based on crosscorrelation"))
  while(TRUE){
    ## recalculate the differences
    dup<-diff(up/20)
    dposi<-diff(posi$startProcTime)
    print(paste("length up:", length(up),"length posi",length(posi$startProcTime)))
  
    ## compare the diff of two signal
    ## start at last index, to avoid retesting always the same data point
    lmin<-min(length(dup),length(dposi))
    dif<-head(dup,n=lmin)-head(dposi,n=lmin)
    index<-which(dif< -2.5 | dif > 2.5) ## if difference is larger than 2.5 ms
    index<-index[which(index>latestRemoved)][1+attempted]
    print(index)
    if(length(index)==0|is.na(index)){
      print("no more indices to test")
      break()
    }
    print(paste("large difference at index",index))
    
    plot(dup[(index-ccLength):(index+ccLength)],type='l')  
    lines(dposi[(index-ccLength):(index+ccLength)],type='l',col="red")  
    Sys.sleep(1)
    
    ## should we remove an up or posi line?
    ## if the crosscorrelation function of the 2 signals has a peak after 0, remove up[index]
    ## get a segment of data of approximately 250 data points
    if(index+ccLength>length(posi$no)){
      endIndex=length(posi$no) 
    }else{
      endIndex=index+ccLength
    }
    
    cc<-ccf(dup[(index+2):endIndex],dposi[(index+2):endIndex])
    print(paste("peak crosscorrelation between", index, "and", endIndex, "is" ,round(cc$acf[which.max(cc$acf)],3),"at lag",cc$lag[which.max(cc$acf)]))
    if(cc$lag[which.max(cc$acf)]>0 & cc$acf[which.max(cc$acf)] > 0.5){
      print(paste("Removing index in up",index))
      up<-up[0-index]
      removedUp=removedUp+1
      lastRemoved=index  
      attempted=0
    } else if(cc$lag[which.max(cc$acf)]<0 & cc$acf[which.max(cc$acf)] > 0.8){
      print(paste("Removing index in posi",index))
      posi<-posi[0-index,]
      removedPosi=removedPosi+1
      lastRevmoved=index
      attempted=0
    } else{
      print(paste("Index",index,"not removed"))
      attempted=attempted+1
    }
  }
    
  if(length(up)!=length(posi$no)){
    print(paste("alignment failed, number of up and posi still differ"))
    return(NA)
  }
  
  dup<-diff(up)
  dposi<-diff(posi$startProcTime*20)
  dif<-dup-dposi
  plot(dif,type='l')
  
  if(removedUp+removedPosi>10){
    print(paste("Removed up:",removedUp, "removed posi:",removedPosi))
    print(paste("More than 10 missalignments were detected. There is a problem with recording setup"))
    return(NA)
  }
  
  cor1<-cor(head(dup,n=2000),head(dposi,n=2000))
  cor2<-cor(dup[(length(dposi)-2000):length(dposi)],dposi[(length(dposi)-2000):length(dposi)])
  cor3<-cor(dup,dposi)
  print(paste("correlation first 2000:",round(cor1,4)))
  print(paste("correlation last 2000:",round(cor2,4)))
  print(paste("correlation all:",round(cor3,4)))
  print(paste("Removed up:",removedUp, "removed posi:",removedPosi))
  if(cor1<0.95|cor2<0.95| cor3<.95){
    print(paste("alignment failed because of correlation of jitter too low"))
    return(NA)
  }
  return(list(up=up,posi=posi))
}

#' Check the integrity of positrack data read from positrack file. 
#' 
#' @param posi Data frame containing the positrack file read via read.table(fn,header=T)
#' @param maxDelayCapProc Maximum allowed delay between frame capture and processing
#' @param maxInterCapDelay Maximum allowed delay between the time stamps of frame capture
#' @param maxInterProcDelay Maximum allowed delay between the time stamps of frame processing
#' @param maxProcDuration Maximum allowed processing of frame duration
#' @return Return 0 if all is ok and positive number if something is wrong
checkIntegrityPositrackData<-function(posi,maxDelayCapProc=1000,
                                      maxInterCapDelay=1000,
                                      maxInterProcDelay=1000,
                                      maxProcDuration=1000){
  print("Check integrity of positrack file")
  ## check for valid header
  validHeaderBeginning<-c("no","capTime","startProcTime")
  if(!all(validHeaderBeginning==names(posi)[1:3])){
    print(paste("Header of positrack is not starting with",validHeaderBeginning))
    print("Make sure you are using the latest version of positrack",validHeaderBeginning)
    return(1)
  }

  # there should not be any characters in the data
  if(any(is.character(posi))){
    print(paste("There are characters in the data"))
    return(2)
  }

  if(any(diff(posi$no)!=1)){
    print("According to posi$no, some frames are missing from the positrack file")
    return(2)
  }
  
  if(any(posi$x<0&posi$x!=-1.0)){
    print("Some x values are negative but not -1.0")
    return(2)
  }
  
  if(any(posi$y<0&posi$y!=-1.0)){
    print("Some y values are negative but not -1.0")
    return(2)
  }
  
  # check capture/processing delays
  delayCapProc<-posi$startProcTime-posi$capTime
  print(paste("Max delay between capture and processing of frame is",round(max(delayCapProc)),"ms at index",which.max(delayCapProc)))
  print(paste("Number of delay between capture and precessing of frame above 100 ms:",sum(delayCapProc>100)))
  if(max(delayCapProc)>maxDelayCapProc)
  {
    print(paste("There is a delay between frame capture and frame processing that is longer",maxDelayCapProc))
    return(2)
  }
  
  # Check inter-caputre delays
  # This would suggest that some frames might have been lost
  interCapDelay<-diff(posi$capTime)
  print(paste("Max inter caputre delay:",max(interCapDelay),"ms at index",which.max(interCapDelay)))
  print(paste("Number of inter capture time above 30 ms:",sum(interCapDelay>30)))
  if(max(interCapDelay)>maxInterCapDelay)
  {
    print(paste("There is a inter caputre delay that is larger than ",maxInterCapDelay))
    return(3)
  }
  
  ## Check inter-processing delays
  interProcDelay<-diff(posi$startProcTime)
  print(paste("Max inter processing delay:",round(max(interProcDelay)),"ms at index",which.max(interProcDelay)))
  print(paste("Number of inter processing time above 50 ms:",sum(interProcDelay>50)))
  if(max(interProcDelay)>maxInterProcDelay)
  {
    print(paste("There is a inter processing delay that is larger than ",maxInterProcDelay))
    return(3)
  }
  
  ## Check for processing time
  print(paste("Max processing duration:",max(posi$procDuration)))
  print(paste("Number of processing duration above 50 ms:",sum(posi$procDuration>50)))
  if(max(posi$procDuration)>maxProcDuration)
  {
    print(paste("There is a processing duration that is larger than ",maxProcDuration))
    return(4)
  }
  
  ## Check the percentage of invalid data point
  percentInvalid<-sum(posi$x==-1.0)/length(posi$x)*100
  print(paste("Percentage of invalid data points:",round(percentInvalid,2),"%"))
  if(percentInvalid>2.5)
  {
    print("********* The percentage of invalid data points is too high **********")
    print("Make sure your LED array and video camera are set properly")  
  }
  
  ## Check for jump in the x y position of the animal
  x<-posi$x[which(posi$x!=-1.0)]
  y<-posi$y[which(posi$x!=-1.0)]
  distance<-sqrt(diff(x)^2+diff(y)^2)
  thresholdJump<-30
  percentJump<-length(which(distance>thresholdJump))/length(distance)*100
  print(paste("Percentage of position jumps (diff >",thresholdJump,"px):",round(percentJump,2),"%"))
  
  ## Check for jump in the head direction of the animal
  hd<-posi$hd[which(posi$hd!=-1.0)]
  thresholdJump<-30
  diffHd<-abs(diff(hd))
  diffHd[which(diffHd>360-thresholdJump)]<-360-diffHd[which(diffHd>360-thresholdJump)]
  percentJump<-length(which(diffHd>thresholdJump))/length(diff(diffHd))*100
  print(paste("Percentage of head-direction jumps (diff >",thresholdJump,"deg):",round(percentJump,2),"%"))
  
  return(0)  
}

#' Check the integrity of the up. 
#' 
#' @param up Time stamps of the ttl pulses of tracking system
#' @param maxInterUpDelay Maximal delay between ups that will be allowed
#' @param samplingRate Sampling rate for the time point in up file
#' @param datLengthSamples Number of samples in the .dat file
#' @return Return 0 if all is ok and positive number if something is wrong
checkIntegrityUp<-function(up,maxInterUpDelay=1000,samplingRate=20000,datLengthSamples){
  d<-diff(up)
  print("Check integrity of up")
  print(paste("Minimum delay between up", min(d)/samplingRate*1000, "ms"))
  print(paste("Maximum delay between up", max(d)/samplingRate*1000, "ms"))
  if(any(diff(up)<1*samplingRate/1000)){
    print(paste("There are up intervals shorter than 1 ms"))
    print(paste("The first case occured at index",head(which(diff(up)<20),n=1)))
    print(paste("Assuming 50 hz, this is", head(which(diff(up)<20),n=1)/50,"sec in the trial"))
    return(1)
  }
  if(any(diff(up)>maxInterUpDelay*samplingRate/1000)){
    print(paste("There are up intervals larger than",maxInterUpDelay))
    print(paste("The first case occured at index",head(which(diff(up)>maxInterUpDelay*samplingRate/1000),n=1)))
    return(1)
  }
  fromStart=head(up,n=1)/samplingRate
  print(paste("First up is",fromStart, "seconds from the beginning of dat file"))
  if(fromStart<0.200){
    print(paste("The first up is less than 200 ms from the start of the dat file"))
    return(1)
  }
  fromEnd=(datLengthSamples-tail(up,n=1))/samplingRate
  print(paste("Last up is",fromEnd, "seconds from the end of the dat file"))
  if(fromEnd<0.200){
    print(paste("The last up",tail(up,n=1),"is less than 200 ms from the end of the dat file",datLengthSamples))
    return(1)
  }
  return(0)  
}
