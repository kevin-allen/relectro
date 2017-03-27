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
#' @param overwirte Logical indicating whether to recreate a whd file if one already exists
whdFromPositrack<-function(rs,
                           resSamplesPerWhdSample=400,
                           ttlChannel=NA,
                           maxUpDiffRes=4000,
                           overwrite=FALSE)
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
  {
    return()
  }
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
    if(checkIntegrityUp(up,samplingRate=rs@samplingRate)!=0)
      stop(paste("check of integrity of up failed"))
   
    
    ######################################
    ## get the data from positrack file ##
    ######################################
    fn<-paste(paste(rs@path,rs@trialNames[tIndex],sep='/'),"positrack",sep='.')
    if(!file.exists(fn))
      stop(paste("whdFromPositrack, file missing:",fn))
    print(paste("reading",fn))
    posi<-read.table(fn,header=T)  ## now assumes that there is a header
    if(checkIntegrityPositrackData(posi)!=0)
      stop(paste("check of integrity of positrack file failed"))
    
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
    if(interEventCor<0.6)
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
whdAlignedTtlPositrack<-function(up,posi){
  dup<-diff(up)
  dposi<-diff(posi$startProcTime*20)
  minLength<-min(c(length(dup),length(dposi)))
  print(paste("correlation first 500:",round(cor(head(dup,n=500),head(dposi,n=500)),3)))
  print(paste("correlation last 500:",round(cor(dup[(minLength-500):minLength],dposi[(minLength-500):minLength]),3)))
  print(paste("length of up:",length(up)))
  print(paste("length of posi:",length(posi$no)))
  
  if(length(up)==length(posi$no)){
    print("Alignment appears ok")
    return(list(up=up,posi=posi))
  }
  
  if(abs(length(up)-length(posi$no))>10){
    print("The alignment problem is for more than 10 frames, no solution implemented for this yet")
    return(NA)
  }
    
  ## visualize the situation
  dup<-diff(up)
  dposi<-diff(posi$startProcTime*20)
  lmin<-min(length(dup),length(dposi))
  dif<-head(dup,n=lmin)-head(dposi,n=lmin)
  plot(dif,type='l')
  removedPosi=0
  removedUp=0
  attempted=0
  print(paste("Try to remove ttl or posi data points based on crosscorrelation"))
  while(TRUE){
    ## recalculate the differences
    dup<-diff(up)
    dposi<-diff(posi$startProcTime*20)
    print(paste("length up:", length(up),"length posi",length(posi$startProcTime)))
    
    ## compare the diff of two signal
    ## start at last index, to avoid retesting always the same data point
    lmin<-min(length(dup),length(dposi))
    dif<-head(dup,n=lmin)-head(dposi,n=lmin)
    #plot(dif,type='l')
    if(attempted==0){
      index<-head(which(dif< -50),n=1)
    } else{
      index<-head(which(dif< -50)[-c(1:attempted)],n=1)
    }
    if(length(index)==0){
      print("no more indices to test")
      break()
    }
    print(paste("large difference at index",index))
    ## should we remove an up or posi line?
    ## if the crosscorrelation function of the 2 signals has a peak after 0, remove up[index]
    ## get a segment of data of approximately 250 data points
    if(index+250>length(posi$no)){
      endIndex=length(posi$no) 
    }else{
      endIndex=index+250
    }
    cc<-ccf(dup[index:endIndex],dposi[index:endIndex])
    print(paste("peak crosscorrelation between", index, "and", endIndex, "is" ,round(cc$acf[which.max(cc$acf)],3),"at lag",cc$lag[which.max(cc$acf)]))
    ## Warning: the peak in cc$acf might be below 0.90 if there are two data point to remove to remove
    if(cc$lag[which.max(cc$acf)]>0 & cc$acf[which.max(cc$acf)] > 0.7){
      print(paste("Removing index in up",index))
      up<-up[0-index]
      removedUp=removedUp+1
      attempted=0
    } else if(cc$lag[which.max(cc$acf)]<0 & cc$acf[which.max(cc$acf)] > 0.7){
      print(paste("Removing index in posi",index))
      posi<-posi[0-index,]
      removedPosi=removedPosi+1
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
  if(cor1<0.96|cor2<0.96| cor3<.98){
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

  # check capture/processing delays
  delayCapProc<-posi$startProcTime-posi$capTime
  print(paste("Max delay between capture and processing of frame is",max(delayCapProc),"at index",which.max(delayCapProc)))
  print(paste("Number of delay between capture and precessing of frame above 100 ms:",sum(delayCapProc>100)))
  if(max(delayCapProc)>maxDelayCapProc)
  {
    print(paste("There is a delay between frame capture and frame processing that is longer",maxDelayCapProc))
    return(2)
  }
  
  # Check inter-caputre delays
  # This would suggest that some frames might have been lost
  interCapDelay<-diff(posi$capTime)
  print(paste("Max inter caputre delay:",max(interCapDelay),"at index",which.max(interCapDelay)))
  print(paste("Number of inter capture time above 30 ms:",sum(interCapDelay>30)))
  if(max(interCapDelay)>maxInterCapDelay)
  {
    print(paste("There is a inter caputre delay that is larger than ",maxInterCapDelay))
    return(3)
  }
  
  ## Check inter-processing delays
  interProcDelay<-diff(posi$startProcTime)
  print(paste("Max inter processing delay:",max(interProcDelay),"at index",which.max(interProcDelay)))
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
  
  
  #plot((which.max(delayCapProc)-15):(which.max(delayCapProc)+100),
  #     delayCapProc[(which.max(delayCapProc)-15):(which.max(delayCapProc)+100)],ylim=c(0,max(delayCapProc)),
  #     ylab="Capture-Processing delays",
  #     xlab="Frame number")
  plot(delayCapProc,type='l',
       ylab="Capture-Processing delays",
       xlab="Frame number")
  lines(interCapDelay,col='red')
  lines(posi$procDuration,col='blue')
  return(0)  
}

#' Check the integrity of the up. 
#' 
#' @param up Time stamps of the ttl pulses of tracking system
#' @param maxInterUpDelay Maximal delay between ups that will be allowed
#' @param samplingRate Sampling rate for the time point in up file
#' @return Return 0 if all is ok and positive number if something is wrong
checkIntegrityUp<-function(up,maxInterUpDelay=1000,samplingRate=20000){
  
  d<-diff(up)
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
  
  
  return(0)  
}