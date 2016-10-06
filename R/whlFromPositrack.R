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
whdFromPositrack<-function(rs,
                           resSamplesPerWhdSample=400,
                           ttlChannel=NA,
                           maxUpDiffRes=4000)
  {
  
  if(rs@session==""){
    stop(paste("whdFromPositrack, rs@session == \"\""))
  }
  if(resSamplesPerWhdSample<=0)
    stop(paste("whdFromPositrack, resSamplesPerWhdSample <= 0", resSamplesPerWhdSample))
  if(maxUpDiffRes<=0)
    stop(paste("whdFromPositrack, maxUpDiffRes <= 0", maxUpDiffRes))
  
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
    ## get the data from dat file
    print(paste("reading sycn channel",ttlChannel[tIndex],"from",rs@trialStartRes[tIndex],"to",rs@trialEndRes[tIndex]))
    x<-as.numeric(datFilesGetChannels(df,channels=ttlChannel[tIndex],
                                      firstSample = rs@trialStartRes[tIndex],lastSample = rs@trialEndRes[tIndex]))
    
    up<-detectUps(x) ## detect rising times of ttl pulses
    fn<-paste(paste(rs@path,rs@trialNames[tIndex],sep='/'),"positrack",sep='.')
    if(!file.exists(fn))
      stop(paste("whdFromPositrack, file missing:",fn))
    
    con=file(fn,open="r")
    line=readLines(con,n=1) 
    close(con)
    
    if(substr(line, 1, 2)=="no"){ 
      # this is for the new positrack files with header
      print("positrack file with header")
      posi<-read.table(fn,header=T)  
      lup<-length(up)
      lposi<-length(posi$startProcTime)
      print(paste("Number of ttl pulses:",lup))
      print(paste("Number of frames in positrack file:",lposi))
      if(lup!=lposi)  
      { # try to align the frames using jitters
        print(paste("whdFromPositrack,",rs@trialNames[tIndex],"length of up (",lup,") and positrack (",lposi,") differs"))
        if(!is.list(x<-whdAlignedTtlPositrack(up,posi))){
          print(paste("aligmnet failed"))
          stop()
        }      
        up<-x$up
        posi<-x$posi
      }
      ## the frame is capture before it is received by the computer, this is the delay in .dat samples
      delay<-(posi$startProcTime-posi$capTime)*rs@samplingRate/1000
      ## remove the delay from the ttl pulse
      up<-up-delay
    } else{ 
      # this is for the old positrack files without header
      print("positrack file without header")
      posi<-read.table(fn)  
      colnames(posi)<-c("startProcTime","no","x","y","hd")
      lup<-length(up)
      lposi<-length(posi$startProcTime)
      print(paste("Number of ttl pulses:",lup))
      print(paste("Number of frames in positrack file:",lposi))
      if(lup!=lposi)  
      {
        print(paste("whdFromPositrack,",rs@trialNames[tIndex],"length of up (",lup,") and positrack (",lposi,") differs"))
        stop()      
      }
      ## because of the firewire camera buffer is one frame late ##
      ## we shift the spots forward by one and ignore the last up ##
      posi<-posi[2:length(posi$startProcTime),]
      up<-up[1:(length(up)-1)]  
    }
    
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
#' These time intervals can be measured from the positrack file (startProcTime) and the ttl pulse
#' When there is no problem of alignment, the correlation between intervals of ttl and positrack is above .97
#' If there is an additional ttl detected, then the correlation goes down. 
#' How low it is depends at which point the alignment brake
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
    
  if(length(up)<length(posi$no)){
    print("not implemented")  
    return(NA)
  }
  if(length(up)>length(posi$no))
  {
    print(paste("more ups (",length(up),") than lines in posi (",length(posi$no),")"))
    toRemove=length(up)-length(posi$no)
    print(paste("Try to remove",toRemove,"ttl pulses"))
    
    removed=0
    while(removed<toRemove){
      index<-head(which(head(dup,n=length(dposi))-dposi< -50),n=1)
    #  index
    #  dup[990:1003]
    #  dposi[990:1003]
      if(length(index)==1){
        print(paste("Removing index",index))
        up<-up[0-index]
        removed=removed+1
      }else{
        removed=toRemove
      }
    }
    
    if(length(up)!=length(posi$no)){
      print(paste("alignment failed, number of up and posi still differ"))
      return(NA)
    }
    
    dup<-diff(up)
    dposi<-diff(posi$startProcTime*20)
    
    cor1<-cor(head(dup,n=1000),head(dposi,n=1000))
    cor2<-cor(dup[(length(dposi)-1000):length(dposi)],dposi[(length(dposi)-1000):length(dposi)])
    print(paste("correlation first 1000:",round(cor1,4)))
    print(paste("correlation last 1000:",round(cor2,4)))
    if(cor1<0.98|cor2<0.98){
      print(paste("alignment failed because of correlation of jitter too low"))
      return(NA)
    }
    return(list(up=up,posi=posi))
  }
}
