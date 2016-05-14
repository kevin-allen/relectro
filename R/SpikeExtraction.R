#' Perform the spike extraction for a recording session
#' 
#' This function calls spikeExtractionTetrode for each tetrode of the recording session.
#' 
#' A RecSession object is used to get the information about the assignation of the recorded channels to tetrodes and the file names.
#' 
#' @param rs RecSession object

spikeExtractionSession<-function(rs){
  df<-new("DatFiles")
  df<-datFilesSet(df,
                  fileNames=paste(rs@trialNames,"dat",sep="."),
                  path=rs@path,
                  nChannels=rs@nChannels)
  for (i in 3:rs@nElectrodes)
  {
    spikeExtractionTetrode(rs,df,tetrodeNumber=i)
  }
}


#' Perform the spike extraction for a given tetrode
#' 
#' A RecSession object is used to get the information about the assignation of the recorded channels to tetrodes and the file names.
#' For each channel of the tetrode, 
#' the raw signal is band-pass filtered and the power is estimated with the root mean square in sliding windows.
#' Baseline variation in power are estimated by the standard deviation of power.
#' Time windows above a threshold contain a spike. 
#' The spike times are aligned to the most negative value within the adjacent windows with power above the threshold.
#' The spike times in sample number are saved in .res.tetrodeNumber.
#' The spike waveforms are extracted and saved in .spk.tetrodeNumber.
#' The features of spikes are obtained via principal component analysis and 3 features are kept for each channel.
#' The spike features are saved in .fet.tetrodeNumber.
#' 
#' @param rs RecSession object
#' @param df DatFile object
#' @param tetrodeNumber Tetrode number. Index starts at 1.
#' @param minPassHz Minimal frequency in Hz that is kept before calculating power
#' @param maxPassHz Maximal frequency in Hz that is kept
#' @param powerWindowSizeMs Size in ms of the sliding window to calculate power
#' @param powerWindowSlideMs Shift of the power window between calculation of power
#' @param SDThreshold Power threshold in standard deviation above the mean for detection of a spike
#' @param spikeDetectionRefractoryMs Period of refractory in the detection of spikes
#' @param waveformWindowSizeMs Window size used when extracting the spike waveforms.
#' @return Return 0 if sucessfull
spikeExtractionTetrode<-function(rs,df,tetrodeNumber,
                                 minPassHz=800,maxPassHz=5000,powerWindowSizeMs=0.4,powerWindowSlideMs=0.1,SDThreshold=2.5,
                                 simultaneousSpikeMaxJitterMs=0.2,                               
                                 waveformWindowSizeMs=1){
  if(minPassHz<0){
    stop(paste("spikeExtractionTetrode: minPassHz<0"))
  }
  if(maxPassHz<=minPassHz){
    stop(paste("spikeExtractionTetrode: maxPassHz<=minPassHz"))
  }
    
  if(powerWindowSlideMs*rs@samplingRate/1000<1){
    stop(paste("spikeExtractionTetrode: powerWindowSlideMs*rs@samplingRate/1000<1"))
  }
  
  ## get the channel of the electrode
  channels<-rs@channelsTetrode[tetrodeNumber,]
  window=32
  
  ## array to store the detected spike times
  spikes<-matrix(ncol=4)
  colnames(spikes)<-c("time","power","trough","channel")
  
  if(length(channels)==0)
    stop(paste("spikeExtractionTetrode, session:", rs@session, ", tetrode:",tetrodeNumber,", channels has length of 0"))
  
  print(paste("sampling rate:", rs@samplingRate,"Hz"))
  print(paste("minPassHz:", minPassHz,"Hz"))
  print(paste("maxPassHz:", maxPassHz,"Hz"))
  print(paste("powerWindowSizeMs:", powerWindowSizeMs,"ms"))
  print(paste("powerWindowSlideMs:", powerWindowSlideMs,"ms"))
  print(paste("SDThreshold:", SDThreshold))
  
  ## loop for each channel of the electrode
  for(chan in channels){
    print(paste("spike detection tetrode:",tetrodeNumber," channel:",chan))
    ## load all data from one channel
    data<-datFilesGetOneChannel(df,chan,
                              firstSample=0,
                              lastSample=df@samples[1])
    ## filter the data
    dataf<-bandPassFilter(data,rs@samplingRate,minPassHz=minPassHz,maxPassHz=maxPassHz)
    
    ## save a file with filtered data for each channel
    writeBin(dataf,con=paste(paste(rs@path,rs@session,sep="/"),"fil",chan,sep="."),size=4)
    
    sdetec<-detectSpikesFromTrace(data=dataf,
                          rs@samplingRate,
                          powerWindowSizeMs,
                          powerWindowSlideMs,
                          SDThreshold)
    
    sp<-cbind(sdetec$spikeTime,
              sdetec$spikePower,
              sdetec$spikeTrough,
              rep(chan,length(sdetec$spikeTime)))
    
    ## remove spikes at the very beginning and end of signal
    sp<-sp[which(sp[,1]<(length(dataf)-(window/2+1))),]
    sp<-sp[which(sp[,1]>(window/2+1)),]
    
    spikes<-rbind(spikes,sp)
    #   
    #   # plot data
#       m=47000
#       M=48000
#       plot(seq(m,M,1),dataf[m:M],type='l')
#       lines(rms.t[which(rms.t>m&rms.t<M)],
#             rms[which(rms.t>m&rms.t<M)],col='red')
#       lines(c(m,M),c(rms.threshold,rms.threshold),col="purple")
#       points(sp[which(sp[,"time"]>m&sp[,"time"]<M),"time"]+1,
#              sp[which(sp[,"time"]>m&sp[,"time"]<M),"trough"],
#              col="red")
#       
#     
    rm(sp,data,dataf,rms)
  
  }
  
  spikes<-spikes[order(spikes[,1]),] # sort spike matrix according to time
  spikes<-spikes[complete.cases(spikes),] # remove a row of NA that was there from creation of spike matrix
  
  ## join spikes that are within 0.2 ms of each other (4 samples)
  res<-mergeSimultaneousSpikes(spikes[,"time"],spikes[,"trough"],simultaneousSpikeMaxJitterMs*rs@samplingRate/1000)
  
  
  #########################################
  ## write the res file for this tetrode ##
  #########################################
  write(res,file=paste(paste(rs@path,rs@session,sep="/"),"res",tetrodeNumber,sep="."),ncolumns=1)
  write(rep(1,length(res)+1),file=paste(paste(rs@path,rs@session,sep="/"),"clu",tetrodeNumber,sep="."),ncolumns=1)
  
  ###################################
  ## visualize the detected spikes ##
  ###################################
  #plotSpikes(rs,channels,head(res,n=100))
  
  ###################################################
  # extract waveform of each spike, save .spk files #
  ###################################################
  createSpkFile(rs,channels,res,window,tetrodeNumber)
  
  ##########################################
  # PCA for each channel, save .fet files ##
  ##########################################
  window=20
  spikePCA(rs,channels,res,window,tetrodeNumber)
  
  ##############################
  ### delete the fil files  ####
  ##############################
  for(chan in channels){
    file.remove(paste(paste(rs@path,rs@session,sep="/"),"fil",chan,sep="."))
  }
}



#' Detect spikes in a signal
#' 
#' If there are low components in the signal, they should be filtered out.
#' 
#' @param data Vector containing the spikes and some noise
#' @param samplingRate Sampling rate of the trace
#' @param powerWindowSizeMs Window size when calculating power (root mean square)
#' @param powerWindowSlideMs Shift of the window in ms between estimation of power
#' @param SDThreshold Power threshold for spike detection.
#' @return list containing rms, rmsT, rmsSD, rmsMean, rmsThreshold, spikeTime and spikePower
detectSpikesFromTrace<-function(data,
                             samplingRate,
                             powerWindowSizeMs,
                             powerWindowSlideMs,
                             SDThreshold)
{
  rms<-powerRootMeanSquare(data=data,
                           windowSizeSamples=powerWindowSizeMs*samplingRate/1000,
                           windowSlide=powerWindowSlideMs*samplingRate/1000)
  # get the time point of power window
  rms.t<-seq(from=(powerWindowSizeMs*samplingRate/1000)/2, # middle of power window
             by=powerWindowSlideMs*samplingRate/1000,
             length.out=length(rms))
  # get the power baseline and threshold
  rms.sd<-sd(rms)
  rms.mean<-mean(rms)
  rms.threshold<- rms.mean+(rms.sd*SDThreshold)
  # identify spikes by detecting negative peaks in windows with power above threshold
  sp<-identifySpikeTimes(data,
                         power=rms,
                         powerWindowSize=powerWindowSizeMs*samplingRate/1000,
                         powerWindowSlide=powerWindowSlideMs*samplingRate/1000,
                         powerThreshold=rms.threshold)
  list(rms=rms,
       rmsT=rms.t,
       rmsSD=rms.sd,
       rmsMean=rms.mean,
       rmsThreshold=rms.threshold,
       spikeTime=sp[,1],
       spikePower=sp[,2],
       spikeTrough=sp[,2])
}




#'  Do PCA analysis on the spike waveforms
#'    
#'  Use the fil file for each channel. Work on each channel independently.
#'  Make a matrix with a spike per row.
#'  Do PCA analysis of spikes.
#'  Combine the first components of each wire.
#'  Save a fet file
#' 
#' @param rs RecSession
#' @param channels Numeric with the list of channels for the tetrode on which the spikes were detected
#' @param res Numeric with the time stamps of the spikes in sample value
#' @param window Number of ms that will be shown
#' @param tetrodeNumber Tetrode number used in the file name
spikePCA<-function(rs,channels,res,window=32,tetrodeNumber){
  if(length(channels)==0)
    stop(paste("spikePCA: length(channels)==0"))
  if(length(res)==0)
    stop(paste("spikePCA: length(res)==0"))
  
  if(window<=0)
    stop(paste("spikePCA: window<=0"))
  
  for(i in 1:length(channels)){
    if(!file.exists(paste(paste(rs@path,rs@session,sep="/"),"fil",channels[i],sep=".")))
      stop(paste("missing file:",paste(paste(rs@path,rs@session,sep="/"),"fil",channels[i],sep=".")))
  }
  
  if(length(channels)<4){
    featuresPerChannel=4
  } else {
    featuresPerChannel=3
  }
  
  fet<-matrix(ncol=featuresPerChannel*length(channels),nrow=length(res))
  ## for each channel
  for(i in 1:length(channels)){
    v<-readBin(what="numeric",con=paste(paste(rs@path,rs@session,sep="/"),"fil",channels[i],sep="."),size=4,n=1000000000)
    m<-getWaveformMatrix(as.integer(res),as.integer(v),window)
    spike.pca<-prcomp(m)
    fet[,((i-1)*featuresPerChannel+1):((i-1)*featuresPerChannel+featuresPerChannel)]<-predict(spike.pca,m)[,1:featuresPerChannel]
  }
  
  ## add time as last column
  fet<-cbind(fet,res)
  
  write(ncol(fet),ncolumns=1,
        file=paste(paste(rs@path,rs@session,sep="/"),"fet",tetrodeNumber,sep="."))
  write(t(fet), ncolumns = ncol(fet),
        file=paste(paste(rs@path,rs@session,sep="/"),"fet",tetrodeNumber,sep="."),
        append=T)
  return(fet)
}




#'  Get the waveform matrix, one spike per row
#' 
#' @param res Numeric with the time stamps of the spikes in sample value
#' @param v electrophysiological signal
#' @param window Number of data points for each spike
#' @return matrix
getWaveformMatrix<-function(res,v,window=20){
  if(length(res)==0)
    stop("getWaveformMatrix: length(res)==0")
  if(length(v)==0)
    stop("getWaveformMatrix: length(v)==0")
  if(window<=0)
    stop("getWaveformMatrix: window<=0")
  if(max(res)+window/2>length(v))
    stop(paste("getWaveformMatrix: max(res)+window/2>length(v)",max(res)+window/2,length(v)))
  
  results<- .Call("get_waveform_matrix",
                  v,
                  length(v),
                  as.integer(res),
                  length(res),
                  as.integer(window))
  return(results)
  
}
  



#'  Create the spk file for each tetrode
#'    
#'  Use the fil file for each channel
#' 
#' @param rs RecSession
#' @param channels Numeric with the list of channels for the tetrode on which the spikes were detected
#' @param res Numeric with the time stamps of the spikes in sample value
#' @param window Number of ms that will be shown
#' @param tetrodeNumber Tetrode number used in the file name
createSpkFile<-function(rs,channels,res,window=32,tetrodeNumber){
  
  if(length(channels)==0)
    stop(paste("createSpkFile: length(channels)==0"))
  
  if(length(res)==0)
    stop(paste("createSpkFile: length(res)==0"))

  if(window<=0)
    stop(paste("createSpkFile: window<=0"))
  
  for(i in 1:length(channels)){
    if(!file.exists(paste(paste(rs@path,rs@session,sep="/"),"fil",channels[i],sep=".")))
      stop(paste("missing file:",paste(paste(rs@path,rs@session,sep="/"),"fil",channels[i],sep=".")))
  }
  
  ## check the length of the fil files
  a<-readBin(what="numeric",con=paste(paste(rs@path,rs@session,sep="/"),"fil",channels[i],sep="."),size=4,n=1000000000)
  ## allocate a matrix to load the data
  m<-matrix(ncol=length(channels),nrow=length(a))
  ## fill the matrix
  for(i in 1:length(channels)){
    m[,i]<-readBin(what="numeric",con=paste(paste(rs@path,rs@session,sep="/"),"fil",channels[i],sep="."),size=4,n=1000000000)
  }
  
  results<- .Call("create_spk_file",
                  as.integer(m),nrow(m),ncol(m),
                  as.integer(res),
                  length(res),
                  as.integer(window),
                  paste(paste(rs@path,rs@session,sep="/"),"spk",tetrodeNumber,sep="."))

}


#'  Plot the waveform of spikes
#'  
#'  Use the fil file for each channel
#' 
#' @param rs RecSession
#' @param channels Numeric with the list of channels for the tetrode on which the spikes were detected
#' @param res Numeric with the time stamps of the spikes in sample value
#' @param windowMs Number of ms that will be shown
plotSpikes<-function(rs,channels,res,windowMs=2){
  
  if(length(channels)==0)
    stop(paste("plotSpikes: length(channels)==0"))
  
  for(i in 1:length(channels)){
    if(!file.exists(paste(paste(rs@path,rs@session,sep="/"),"fil",channels[i],sep=".")))
      stop(paste("missing file:",paste(paste(rs@path,rs@session,sep="/"),"fil",channels[i],sep=".")))
  }
  
  ## check the length of the fil files
  a<-readBin(what="numeric",con=paste(paste(rs@path,rs@session,sep="/"),"fil",channels[i],sep="."),size=4,n=1000000000)
  ## allocate a matrix to load the data
  m<-matrix(ncol=length(channels),nrow=length(a))
  ## fill the matrix
  for(i in 1:length(channels)){
    m[,i]<-readBin(what="numeric",con=paste(paste(rs@path,rs@session,sep="/"),"fil",channels[i],sep="."),size=4,n=1000000000)
  }
  window<-windowMs*rs@samplingRate/1000
  
  num.cols<-length(channels)
  num.rows<-1
  plot.per.page=num.cols*num.rows
  n<-matrix(c(rep(seq(0,1-(1/num.cols),1/num.cols),num.rows),
              rep(seq(1/num.cols,1,1/num.cols),num.rows),
              rep(seq(1-(1/num.rows),0,0-1/num.rows),each=num.cols),
              rep(seq(1,1/num.rows,0-1/num.rows),each=num.cols)),ncol=4)
  split.screen(n)  
  for(chanIndex in 1:length(channels)){
    screen(chanIndex)
    for(spikeIndex in 1:length(res)){
      wfi<-(res[spikeIndex]+1-window/2):(res[spikeIndex]+1+window/2)
      if(spikeIndex==1){
        par(mar=c(1,1,1,1), oma=c(2,1,0,0),cex.lab=0.6,cex.axis=0.6,mgp=c(1,0.3,0.2))
       plot((1:length(wfi))/rs@samplingRate*1000,m[wfi,chanIndex],type='l',ylim=c(-1500,1000),col=chanIndex,xlab="ms")
      } else{
        lines((1:length(wfi))/rs@samplingRate*1000,m[wfi,chanIndex],col=chanIndex)
      }
     }
    }
  close.screen(all.screens = TRUE)
}







#' Merge simulatenous spikes that were detected on different tetrodes
#' 
#' The spike time associated with the smallest trough is kept.
#' 
#' 
#' @param time Numeric containing the spike times
#' @param trough Numeric containing the trough of detected spikes
#' @param maxTimeDifference Maximal time difference to be considered simultaneous
#' @return Numeric containing the spike times
mergeSimultaneousSpikes<-function(time,
                                   trough,
                                   maxTimeDifference){
  if(length(time)==0)
    stop(paste("mergeSimultaneousSpike: length(time)==0"))
  if(length(trough)==0)
    stop(paste("mergeSimultaneousSpike: length(trough)==0"))
  if(length(time)!=length(trough))
    stop(paste("mergeSimultaneousSpike: length(trough)!-length(time)"))
  if(maxTimeDifference<=0)
    stop(paste("mergeSimultaneousSpike: maxTimeDifference<=0"))
  
  results<- .Call("merge_simultaneous_spikes",
                  as.integer(time),trough,length(time), as.integer(maxTimeDifference))
            
  return(results)
}


#' Detect spike time using filtered signal and root mean square arrays
#' 
#' 
#' @param dataf Filtered brain signal
#' @param power Power in sliding time windows
#' @param powerWindowSize Size of the power window
#' @param powerWindowSlide Amount of shift in the power window between subsequent power estimates
#' @param powerThreshold Threshold used to detect spikes
#' @return Matrix containing spike times and negative peak value of the spike. The spike times are in sample index starting at index 0.
identifySpikeTimes<-function(dataf,
                            power,
                            powerWindowSize,
                            powerWindowSlide,
                            powerThreshold){
  if(length(dataf)==0)
    stop(paste("identifySpikeTime: length(dataf)==0"))
  if(length(power)==0)
    stop(paste("identifySpikeTime: length(power)==0"))
  if(powerWindowSize<=0)
    stop(paste("identifySpikeTime: powerWindowSize<=0"))
  if(powerWindowSlide<=0)
    stop(paste("identifySpikeTime: powerWindowSlide<=0"))
  if(powerThreshold<0)
    stop(paste("identifySpikeTime: powerThreshold<0"))
  
  results<- .Call("identify_spike_times",
                  dataf, length(dataf),
                  power, length(power),
                  powerWindowSize,
                  powerWindowSlide,
                  powerThreshold)
  colnames(results)<-c("time","trough","power")
  
  return(results)
}

#' Generate a raw trace with some background gaussian noise and surimposed spikes
#' 
#' Several different waveforms can be used. They are all generated from the same generic spike by adding gaussian noise to it.
#' 
#' @param samplingRate Sampling rate of the trace
#' @param durationSec Total duration in second of the trace
#' @param noiseSD Standard deviation of gaussian noise
#' @param noiseMean Mean of noise
#' @param wavefromAmplitude Negative amplitude of the generic spike waveform
#' @param nClusters Number of different waveforms (neurons) in the trace
#' @param waveformDifferentiationSD Differentiation of the waveforms of different cluster (gaussian noise added in generic waveform)
#' @return list containing trace, spikeTimes and cluId
simulateRawTrace<-function(samplingRate=20000,
                           durationSec=1,
                           noiseSD=100,
                           noiseMean=0,
                           waveformAmplitude=700,
                           nClusters=3,
                           waveformDifferentiationSD=200,
                           maxSpikes=10000)
{
  # noise in signal gaussian
  noise<-rnorm(n=samplingRate*durationSec,mean=noiseMean,sd=noiseSD)
  signal<-rep(0,length(noise))
  
  # generic waveform
  genericWaveform<-c(0.0,0.2,0,-0.3,-1,-0.6,0,0.2,0.1,0.05,0)*waveformAmplitude
  # generate different spike patterns
  spikeWaveforms<-matrix(ncol=length(genericWaveform),nrow=nClusters)
  for(clu in 1:nClusters){
    spikeWaveforms[clu,] <- genericWaveform + rnorm(n=length(genericWaveform),mean=0,sd=waveformDifferentiationSD)
  }
  
  # get some spike times, minimum of 1 ms isi, one spikes every 10 ms, so 100 Hz on wire
  isi<-rpois(n=maxSpikes,lambda=10)*samplingRate/1000
  isi<-isi[which(isi>samplingRate/1000)]
  spikeTime<-(cumsum(isi))
  # remove spikes at the very end end
  spikeTime<-spikeTime[which(spikeTime<(length(noise)-length(spikeWaveforms)))]
  # remove spikes at the very beginning
  spikeTime<-spikeTime[which(spikeTime>length(spikeWaveforms))]
  # get clu id for each spike
  cluId<-sample(x=1:nClusters,size=length(spikeTime),replace=T)
  
  spikeLength<-length(genericWaveform)
  for(i in 1:length(spikeTime)){
    signal[spikeTime[i]:(spikeTime[i]+spikeLength-1)]<-spikeWaveforms[cluId[i],]
  }
  
  # set spikeTime to the trough of spikes
  spikeTime<-spikeTime+which.min(genericWaveform)-1
  trace<-signal+noise
  list(trace=trace,spikeTime=spikeTime,cluId=cluId)
}


#' Get the results of spike detection
#' 
#' @param sTimeD Spike times of detected spikes
#' @param sTimeT Spike times of simulated spikes
#' @param maxJitter Max jitter to considered simulated and detected spike as the same spikes
spikeDetectionAccuracy<-function(sTimeD,sTimeT,maxJitter=2)
{
## how many of the detected spikes were true spikes
detectedTrue<-sum(sapply(sTimeD,function(x,y){min(abs(x-y))},sTimeT)<=maxJitter)
detectedFalse<-length(sTimeD)-detectedTrue
## how many of the true spikes were detected
trueDetected<-sum(sapply(sTimeT,function(x,y){min(abs(x-y))},sTimeD)<=maxJitter)
## how many of the true spikes were not detected
trueNonDetected<-length(sTimeT)-trueDetected
list(detectedTrue=detectedTrue,
     detectedFalse=detectedFalse,
     trueDetected=trueDetected,
     trueNonDetected=trueNonDetected)
}





