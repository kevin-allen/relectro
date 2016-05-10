#' Apply a band-pass filter to a signal
#' 
#' Use a Butterworth filter.
#' 
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

#' Calculate the root mean square within a sliding window
#' 
#' 
#' @param data Numeric vector containing the signal
#' @param windowSizeSamples Window size in number of samples within which the power will be calculated
#' @param windowSlide Shift of the window between calculation of root mean square
#' @return Numeric, root mean square of all windows
powerRootMeanSquare<-function(data,windowSizeSamples,windowSlide){
  
  if(length(data)==0)
    stop(paste("powerRootMeanSquare, length(data)==0"))
  if(windowSizeSamples<=0)
    stop(paste("powerRootMeanSquare, windowSizeSamples<=0"))
  if(windowSlide<=0)
    stop(paste("powerRootMeanSquare, windowSlide<=0"))
  
  results<- .Call("power_root_mean_square",
                  data, length(data),
                  windowSizeSamples,
                  windowSlide)
  return(results)
}