#' Apply a band-pass filter to a signal
#' 
#' Use a Butterworth filter.
#' 
#' NOT TESTED AT ALL
#' 
#' @param data Numeric vector containing the data to filter
#' @param samplingRate Sampling rate of the data in Hz
#' @param minPass Minimal frequency in Hz that is kept
#' @param maxPass Maximal frequency in Hz that is kept
#' @return Filtered version of the argument data
bandPassFilter<-function(data,samplingRate,minPass,maxPass){
  
  if(samplingRate<0|samplingRAte>100000)
    stop(paste("bandPassFilter(), samplingRate should be between 0 and 100000 Hz but was",samplingRate))
  if(minPass<=0|minPass>samplingRate/2)
    stop(paste("bandPassFilter(), minPass should between 0 and Nyquist frequency but was",minPass))
  if(maxPass<=0|maxPass>samplingRate/2)
    stop(paste("bandPassFilter(), maxPass should between 0 and Nyquist frequency but was",maxPass))
  if(minPass<=maxPass)
    stop(paste("minPass (",minPass,") should be smaller than maxPass(",maxPass,")"))
  
  results<- .Call("band_pass_filter_one_channel_fftw_cwrap",
                  data, length(data),
                  samplingRate,
                  minPass,maxPass)
}

