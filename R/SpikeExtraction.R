#' Perform the spike extraction for a recording session
#' 
#' This function calls spikeExtractionTetrode for each tetrode of the recording session.
#' 
#' A RecSession object is used to get the information about the assignation of the recorded channels to tetrodes and the file names.
#' 
#' @param rs RecSession object
#' @param minPassHz Minimal frequency in Hz that is kept before calculating power
#' @param maxPassHz Maximal frequency in Hz that is kept
#' @param powerWindowSizeMs Size in ms of the sliding window to calculate power
#' @param powerWindowSlideMs Shift in ms of the power window for power calculation
#' @param SDThreshold Power threshold in standard deviation above the mean for detection of a spike
#' @param spikeDetectionRefractoryMs Period of refractory in the detection of spikes
#' @param waveformWindowSizeMs Window size used when extracting the spike waveforms.
spikeExtractionSession<-function(rs,minPassHz=800,maxPassHz=5000,powerWindowSizeMs=0.5,powerWindowSlideMs=0.5,SDThreshold=5,
                                 spikeDetectionRefractoryMs=0.4,
                                 waveformWindowSizeMs=1){

  df<-new("DatFiles")
  df<-datFilesSet(df,
                  fileNames=paste(rs@trialNames,"dat",sep="."),
                  path=rs@path,
                  nChannels=rs@nChannels)
  for (i in 1:rs@nElectrodes)
  {
    spikeExtractionTetrode(rs,df,tetrodeNumber=i,
                           minPassHz=minPassHz,maxPassHz=maxPassHz,powerWindowSizeMs=powerWindowSizeMs,
                           powerWindowSlideMs=powerWindowSlideMs,
                           SDThreshold=SDThreshold,
                           spikeDetectionRefractoryMs=spikeDetectionRefractoryMs,
                           waveformWindowSizeMs=waveformWindowSizeMs)
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
                                 minPassHz=800,maxPassHz=5000,powerWindowSizeMs=0.5,powerWindowSlideMs=0.1,SDThreshold=3,
                                 spikeDetectionRefractoryMs=0.4,
                                 waveformWindowSizeMs=1){
  
  tetrodeNumber=1
  minPassHz=800
  maxPassHz=5000
  powerWindowSizeMs=0.5
  powerWindowSlideMs=0.10
  SDThreshold=3
  spikeDetectionRefractoryMs=0.4
  ## get the channel of the electrode
  channels<-rs@channelsTetrode[tetrodeNumber,]
  ##
  spikes<-matrix(ncol=4)
  colnames(spikes)<-c("time","trough","power","channel")
  
  if(length(channels)==0)
    stop(paste("spikeExtractionTetrode, session:", rs@session, ", tetrode:",tetrodeNumber,", channels has length of 0"))
  
  
  for(chan in channels){
  
    ## load all data from one channel
    system.time(data<-datFilesGetOneChannel(df,chan,
                              firstSample=0,
                              lastSample=sum(df@samples)-1))
    ## filter the data
    dataf<-bandPassFilter(data,rs@samplingRate,minPassHz=minPassHz,maxPassHz=maxPassHz)
    rm(data) # free memory
  
    # calculate power (root mean square)
    rms<-powerRootMeanSquare(data=dataf,
                             windowSizeSamples=powerWindowSizeMs*rs@samplingRate/1000,
                             windowSlide=powerWindowSlideMs*rs@samplingRate/1000)
  
    # get the time point of power window
    #rms.t<-seq(from=(powerWindowSizeMs*rs@samplingRate/1000)/2, # middle of power window
    #           by=powerWindowSlideMs*rs@samplingRate/1000,
    #           length.out=length(rms))
  
    # get the power baseline and threshold
    rms.sd<-sd(rms)
    rms.mean<-mean(rms)
    rms.threshold<- rms.mean+(rms.sd*SDThreshold)
  
    # identify spikes by detecting negative peaks in windows with power above threshold
    sp<-identifySpikeTimes(dataf,
                      power=rms,
                      powerWindowSize=powerWindowSizeMs*rs@samplingRate/1000,
                      powerWindowSlide=powerWindowSlideMs*rs@samplingRate/1000,
                      powerThreshold=rms.threshold)
    
    sp<-cbind(sp,rep(chan,length(sp[,1]))) # add channel number in the spike matrix
    spa<-rbind(spikes,sp) # 
  }
  
  
  
#   
#   # plot data
#   m=40000
#   M=42000
#   plot(seq(m,M,1),dataf[m:M],type='l')
#   lines(rms.t[which(rms.t>m&rms.t<M)],
#         rms[which(rms.t>m&rms.t<M)],col='red')
#   lines(c(m,M),c(rms.threshold,rms.threshold),col="purple")
#   points(spikes[which(spikes[,"time"]>m&spikes[,"time"]<M),"time"]+1,
#          spikes[which(spikes[,"time"]>m&spikes[,"time"]<M),"trough"],
#          col="red")
#   
#   
  
  
  
  
  
  # align spikes on the peak of waveform, save .res files
  # extract waveform of each spike, save .spk files
  # PCA for each channel, save .fet files
  
}



#' Detect spike time using filtered signal and root mean square arrays
#' 
#' 
#' @param dataf Filtered brain signal
#' @param power Power in sliding time windows
#' @param powerWindowSize Size of the power window
#' @param powerWindowSlide Amount of shift in the power window between subsequent power estimates
#' @param powerThreshold Threshold used to detect spikes
#' @return Matrix containing spike times and negative peak value of the spike
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
  colnames(results)<-c("time","trough")
  
  return(results)
}