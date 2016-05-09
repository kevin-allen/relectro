#' Apply a band-pass filter to a signal
#' 
#' Use a Butterworth filter.
#' 
#' This function was not tested extensively
#' 
#' @param data Numeric vector containing the data to filter
#' @param samplingRate Sampling rate of the data in Hz
#' @param minPassHz Minimal frequency in Hz that is kept
#' @param maxPassHz Maximal frequency in Hz that is kept
#' @return Filtered copy of the argument data
bandPassFilter<-function(data,samplingRate,minPassHz,maxPassHz){
  
  if(samplingRate<0|samplingRate>100000)
    stop(paste("bandPassFilter(), samplingRate should be between 0 and 100000 Hz but was",samplingRate))
  if(minPassHz<=0|minPassHz>samplingRate/2)
    stop(paste("bandPassFilter(), minPass should between 0 and Nyquist frequency but was",minPassHz))
  if(maxPassHz<=0|maxPassHz>samplingRate/2)
    stop(paste("bandPassFilter(), maxPass should between 0 and Nyquist frequency but was",maxPassHz))
  if(minPassHz>maxPassHz)
    stop(paste("minPass (",minPassHz,") should be smaller than maxPass(",maxPassHz,")"))
  
  results<- .Call("band_pass_filter_one_channel_fftw_cwrap",
                  data, length(data),
                  samplingRate,
                  minPassHz,maxPassHz)
  return(results)
}


