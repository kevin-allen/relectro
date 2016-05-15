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
  for (i in 1:rs@nElectrodes)
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
  print(paste("spike extration for",rs@session,tetrodeNumber))
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
  
  if(length(channels)==0)
    stop(paste("spikeExtractionTetrode, session:", rs@session, ", tetrode:",tetrodeNumber,", channels has length of 0"))
  
  print(paste("sampling rate:", rs@samplingRate,"Hz"))
  print(paste("filters:", minPassHz,maxPassHz,"Hz"))
  print(paste("Power window:", powerWindowSizeMs,"ms, slide:", powerWindowSlideMs,"ms"))
  print(paste("Power SD threshold:", SDThreshold))
  traces<-matrix(ncol=length(channels),nrow=sum(df@samples))
  for(chan in 1:length(channels)){
    print(paste("reading channel",channels[chan]))
    ## load the raw trace from dat files
    traces[,chan]<-datFilesGetOneChannel(df,channels[chan],
                                      firstSample=0,
                                      lastSample=sum(df@samples)-1)
    ## filter
    traces[,chan]<-bandPassFilter(traces[,chan],rs@samplingRate,minPassHz=minPassHz,maxPassHz=maxPassHz)
  }
  res<-detectSpikesTetrodes(data=traces,
                            samplingRate=rs@samplingRate,
                            powerWindowSizeMs=powerWindowSizeMs,
                            powerWindowSlideMs=powerWindowSlideMs,
                            SDThreshold=SDThreshold,
                            simultaneousSpikeMaxJitterMs=simultaneousSpikeMaxJitterMs)
  
  #########################################
  ## write the res file for this tetrode ##
  #########################################
  write(res,file=paste(paste(rs@path,rs@session,sep="/"),"res",tetrodeNumber,sep="."),ncolumns=1)
  write(rep(1,length(res)+1),file=paste(paste(rs@path,rs@session,sep="/"),"clu",tetrodeNumber,sep="."),ncolumns=1)
  
  ###################################################
  # extract waveform of each spike, save .spk files #
  ###################################################
  createSpkFile(traces,rs,res,32,tetrodeNumber)
  
  #######################################
  # get spike waveforms in a 3D array ###
  #######################################
  swf<-spikeWaveformFromTraces(traces,res,20)
  
  rm(traces)
  ##########################################
  # PCA for each channel, save .fet files ##
  ##########################################
  fet<-spikePCA(swf)
  fet<-cbind(fet,res)
  write(ncol(fet),ncolumns=1,
        file=paste(paste(rs@path,rs@session,sep="/"),"fet",tetrodeNumber,sep="."))
  write(t(fet), ncolumns = ncol(fet),
        file=paste(paste(rs@path,rs@session,sep="/"),"fet",tetrodeNumber,sep="."),
        append=T)
  rm(res,fet,swf)
}


#'  Create the spike waveform array
#'    
#'  Returns a 3D array with the waveforms of each spike for each channel.
#'  array[spike,time,channel]
#' 
#' @param traces Matrix with the traces for each channel. Each column is a channel
#' @param res Numeric with the time stamps of the spikes in sample value
#' @param window Number of samples for each spike
#' @return array[spike,time,channel]
spikeWaveformFromTraces<-function(traces,res,window=20){
  if(class(traces)!="matrix")
    stop("spikeWaveformFromTraces: traces should be a matrix")
  if(length(res)==0)
    stop(paste("spikeWaveformFromTraces: length(res)==0"))
  if(window<=0)
    stop(paste("spikeWaveformFromTraces: window<=0"))
  if((length(traces[,1])-window/2)<max(res)){
    print(paste("length(traces[,1]):",length(traces[,1]),"max(res)",max(res)))
    stop(paste("spikeWaveformFromTraces:",
               "length(traces[,1]):",length(traces[,1]),"max(res)",max(res),
                "traces not long enough for some spike times"))
  }
  if(min(res)<window/2)
    stop(paste("spikeWaveformFromTraces: some res values start before the half window mark"))
  
  results<- .Call("spike_waveform_from_traces",
                  as.integer(traces),nrow(traces),ncol(traces),
                  as.integer(res),
                  length(res),
                  as.integer(window))
  return(array(data=results,dim=c(length(res),window,ncol(traces))))
}


#'  Create the spk file for each tetrode
#'    
#'  Use the fil file for each channel
#' 
#' @param traces Matrix with the electrophysiological traces
#' @param rs RecSession
#' @param res Numeric with the time stamps of the spikes in sample value
#' @param window Number of ms that will be shown
#' @param tetrodeNumber Tetrode number used in the file name
createSpkFile<-function(traces,rs,res,window=32,tetrodeNumber){
  
  if(class(traces)!="matrix")
    stop("createSpkFile: traces is not a matrix")
  if(length(res)==0)
    stop(paste("createSpkFile: length(res)==0"))
  if(window<=0)
    stop(paste("createSpkFile: window<=0"))
  
  results<- .Call("create_spk_file",
                  as.integer(traces),nrow(traces),ncol(traces),
                  as.integer(res),
                  length(res),
                  as.integer(window),
                  paste(paste(rs@path,rs@session,sep="/"),"spk",tetrodeNumber,sep="."))
  
}

#'  Do PCA analysis on the spike waveforms, treating each channel independently
#'  
#'  Do PCA analysis from spike waveform, keep the first 3 or 4 principal components
#'  Combine the first components of each wire.
#' 
#' @param swf Array with the spike waveforms [spike,time,channel]
#' @return Matrix with the different spike features, one spike per row

spikePCA<-function(swf){
  if(class(swf)!="array")
    stop(paste("spikePCA: swf is not an array but a",class(swf)))
  if(length(swf[,1,1])<2)
    stop("spikePCA: swf[,1,1]<2, too few spikes")
  if(length(dim(swf))!=3)
    stop("spikePCA: swf does not have 3 dimensions")
  if(dim(swf)[3]<4){
    featuresPerChannel=4
  } else {
    featuresPerChannel=3
  }
  
  fet<-matrix(ncol=featuresPerChannel*dim(swf)[3],nrow=length(swf[,1,1]))
  ## for each channel
  for(i in 1:dim(swf)[3]){
    spike.pca<-prcomp(swf[,,i])
    fet[,((i-1)*featuresPerChannel+1):((i-1)*featuresPerChannel+featuresPerChannel)]<-predict(spike.pca,swf[,,i])[,1:featuresPerChannel]
  }
  return(fet)
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

#' Detect spikes on the channels of a tetrode
#' 
#' If there are low components in the signal, they should be filtered out.
#' 
#' @param data Matrix containing the spikes and some noise, one channel per column
#' @param samplingRate Sampling rate of the trace
#' @param powerWindowSizeMs Window size when calculating power (root mean square)
#' @param powerWindowSlideMs Shift of the window in ms between estimation of power
#' @param SDThreshold Power threshold for spike detection.
#' @return list containing rms, rmsT, rmsSD, rmsMean, rmsThreshold, spikeTime and spikePower
detectSpikesTetrodes<-function(data,
                               samplingRate,
                               powerWindowSizeMs,
                               powerWindowSlideMs,
                               SDThreshold,
                               simultaneousSpikeMaxJitterMs,
                               noDetectionBeginEndMs=0.5)
{
  nChannels=ncol(data)
  ## array to store the detected spike times
  spikes<-matrix(ncol=4)
  colnames(spikes)<-c("time","power","trough","channel")
  results<-list()
  for(chan in 1:nChannels){
    sdetec<-detectSpikesFromTrace(data[,chan],
                          samplingRate,
                          powerWindowSizeMs,
                          powerWindowSlideMs,
                          SDThreshold)
    sp<-cbind(sdetec$spikeTime,
              sdetec$spikePower,
              sdetec$spikeTrough,
              rep(chan,length(sdetec$spikeTime)))
    
    ## remove spikes at the very beginning and end of signal
    sp<-sp[which(sp[,1]<(length(data[,chan])-noDetectionBeginEndMs*samplingRate/1000)),]
    sp<-sp[which(sp[,1]>(noDetectionBeginEndMs*samplingRate/1000)),]
    spikes<-rbind(spikes,sp)
  }
  spikes<-spikes[order(spikes[,1]),] # sort spike matrix according to time
  spikes<-spikes[complete.cases(spikes),] # remove a row of NA that was there from creation of spike matrix
  
  ## join spikes that are within 0.2 ms of each other (4 samples)
  res<-mergeSimultaneousSpikes(spikes[,"time"],spikes[,"trough"],simultaneousSpikeMaxJitterMs*samplingRate/1000)
  return(res)
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






#' Generate raw traces with some background gaussian noise and surimposed spikes
#' 
#' Several different waveforms can be used. They are all generated from the same generic spike by adding gaussian noise to it.
#' 
#' @param samplingRate Sampling rate of the trace
#' @param durationSec Total duration in second of the trace
#' @param noiseSD Standard deviation of gaussian noise
#' @param noiseMean Mean of noise
#' @param wavefromAmplitude Negative amplitude of the generic spike waveform
#' @param nClusters Number of different waveforms (neurons) in the trace
#' @param nChannels Number of channels, 4 in case of tetrodes
#' @param waveformDifferentiationSD Differentiation of the waveforms of different cluster (gaussian noise added in generic waveform)
#' @param minChannelScalling Minimal value (between 0 and 1) that can be used for creating the tetorde effect 
#' (scalling of waveform on different channels)
#' @return list containing trace, spikeTimes and cluId
simulateRawTrace<-function(samplingRate=20000,
                           durationSec=1,
                           noiseSD=100,
                           noiseMean=0,
                           waveformAmplitude=700,
                           nClusters=3,
                           nChannels=4,
                           waveformDifferentiationSD=200,
                           minChannelScalling=0.3,
                           maxSpikes=10000)
{
  # noise in signal gaussian
  noise<-matrix(data=rnorm(n=samplingRate*durationSec*nChannels,mean=noiseMean,sd=noiseSD),
                nrow=samplingRate*durationSec,ncol=nChannels)
  signal<-matrix(data=rep(0,ncol(noise)*nrow(noise)),
                 nrow=samplingRate*durationSec,ncol=nChannels)
  
  # generic waveform
  genericWaveform<-c(0.0,0.2,0,-0.3,-1,-0.6,0,0.2,0.1,0.05,0)*waveformAmplitude
  # generate different spike patterns, shared across channels
  spikeWaveforms<-matrix(ncol=length(genericWaveform),nrow=nClusters)
  for(clu in 1:nClusters){
    spikeWaveforms[clu,] <- genericWaveform + rnorm(n=length(genericWaveform),mean=0,sd=waveformDifferentiationSD)
  }
  # resize the spike so that amplitude vary across channels
  channelScaling<-matrix(data=runif(nChannels*nClusters, minChannelScalling, 1),nrow=nClusters,ncol=nChannels)
  # ensure that one random channel has an Scaling of 1 for each cluster
  for(i in nClusters){
    channelScaling[i,sample(1:4,1)]=1
  }
  
  # get some spike times, minimum of 1 ms isi, one spikes every 10 ms, so 100 Hz on wire
  isi<-rpois(n=maxSpikes,lambda=10)*samplingRate/1000
  isi<-isi[which(isi>samplingRate/1000)]
  spikeTime<-(cumsum(isi))
  # remove spikes at the very end end
  spikeTime<-spikeTime[which(spikeTime<(length(noise[,1])-length(spikeWaveforms)))]
  # remove spikes at the very beginning
  spikeTime<-spikeTime[which(spikeTime>length(spikeWaveforms))]
  # get clu id for each spike
  cluId<-sample(x=1:nClusters,size=length(spikeTime),replace=T)
  
  spikeLength<-length(genericWaveform)
  for(i in 1:length(spikeTime)){
    for(j in 1:nChannels){
      signal[spikeTime[i]:(spikeTime[i]+spikeLength-1),j]<-spikeWaveforms[cluId[i],]*channelScaling[cluId[i],j]
    }
  }
  
  # set spikeTime to the trough of spikes
  spikeTime<-spikeTime+which.min(genericWaveform)-1
  trace<-signal+noise
  list(trace=trace,spikeTime=spikeTime,cluId=cluId)
}