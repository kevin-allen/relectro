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
  if(is.na(ttlChannel)){
    ttlChannel=rs@nChannels-1
  }
  
  ext="whd.r"
  mainUp<-vector()
  mainPosi<-data.frame()
  
  # create the whd file for each .dat file
  for(tIndex in 1:length(rs@trialNames)){
    print(paste(tIndex, rs@trialNames[tIndex]))
    ## get the data from dat file
    print(paste("reading sycn channel",ttlChannel,"from",rs@trialStartRes[tIndex],"to",rs@trialEndRes[tIndex]))
    x<-as.numeric(datFilesGetChannels(df,channels=ttlChannel,
                                     firstSample = rs@trialStartRes[tIndex],lastSample = rs@trialEndRes[tIndex]))
    up<-detectUps(x) ## detect rising times of ttl pulses
    fn<-paste(paste(rs@path,rs@trialNames[tIndex],sep='/'),"positrack",sep='.')
    if(!file.exists(fn))
      stop(paste("whdFromPositrack, file missing:",fn))
    posi<-read.table(fn)  
    colnames(posi)<-c("time","no","x","y","hd")
    lup<-length(up)
    lposi<-length(posi$time)
    
    if(lup!=lposi)  
    {
      print(paste("whdFromPositrack,",rs@trialNames[tIndex],"length of up (",lup,") and positrack (",lposi,") differs"))
      stop()      
    }
    
    ## because of the firewire camera buffer is one frame late ##
    ## we shift the spots forward by one and ignore the last up ##
    posi<-posi[2:length(posi$time),]
    up<-up[1:(length(up)-1)]  
  
    ### create a main up and main posi for the entire recording session
    mainUp<-c(mainUp,up+rs@trialStartRes[tIndex])
    mainPosi<-rbind(mainPosi,posi)
    
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
