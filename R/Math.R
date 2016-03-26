#' Returns possible pairs from one or two vectors of values
#' 
#' If only on vector of values is given, all possible pairs within these values are returned.
#' If two vectors of values are given, only pairs across the two vectors are returned.
#' 
#' @param cl1 Numeric vector containing a vector of values
#' @param cl2 Optional argument containing a second vector of values
#' @return A data.frame containing the pairs
#' @examples MakePairs(cl1=1:10)
makePairs<-function(cl1="",cl2=NULL){
  if(is.null(cl2)){
    m<-combn(cl1,m=2)
    data.frame(Var1=m[1,],Var2=m[2,])
  }
  else{
    expand.grid(1:3,4:5)
  }
}

#' Smooth the values in a numeric vector using a Gaussian kernel
#' 
#' The values at -1.0 are consider invalid by defaults and are not used or changed. 
#' Set the value of the argument degrees to TRUE if the data are circular. 
#' For example, if the first and last data of the vector should be considered next to each other.
#' 
#' @param x Numeric vector
#' @param sd The standard deviation of the Gaussian kernel used to smooth
#' @param invalid Numeric indicating which value should be treated as NA, by default -1.0.
#' @param degrees Logical indicating if the smoothing should consider the vector as circular, i.e. the first and last values are adjacent.
#' @examples smoothGaussian(x=c(1:10,9:1),sd=2,invalid=-1.0,degrees=FALSE)
smoothGaussian<-function(x,sd=2,invalid=-1.0,degrees=FALSE)
{
  if(length(x)==0)
    return
  if(sd==0)
    return
  if(sd<0)
    stop(paste("sd is smaller than 0:",sd))
  if(class(x)=="integer"){
    x<-as.numeric(x)
  }
  if(degrees==FALSE){
    results<- .Call("smooth_double_gaussian_cwrap",
                    x, length(x), sd, invalid)
  }
  if(degrees==TRUE){
    if(any(x>360))
      stop(paste("x values larger than 360"))
    results<- .Call("smooth_double_gaussian_degrees_cwrap",
                    x, length(x), sd, invalid)  
  }
  return(results)
}

#' Shift values in a vector by a certain number of places and in a given direction
#' 
#' The vector is wrapped around so that values that would end up after the end of the vector are place at the beginning.
#' 
#' @param v A vector
#' @param places Number of places the values will be moved.
#' @param dir Direction of the shift, values should be right or left (or r or l).
#' @examples shift(v=1:10, place=2, dir="r")
shift <- function (v, places, dir = "right") 
{# shift a vector in place
  vnew <- NULL
  d <- substring(dir, 1, 1)
  if (d == "r" & places < 0) {
    d <- "l"
  }
  else {
    if (d == "l" & places < 0) {
      d <- "r"
    }
  }
  n <- length(v)
  p <- abs(places)
  if (p == 0) {
    vnew <- v
  }
  else {
    if (d == "r") {
      vnew <- c(v[(n - p + 1):n], v[1:(n - p)])
    }
    else {
      vnew <- c(v[(p + 1):n], v[1:p])
    }
  }
  return(vnew)
}

#' Shift a the values of a vector by a random amount that is at least as large as the argument minMvMs
#' 
#' @param x A vector
#' @param timePerSampleRes Time in sample values (from the .dat files) between the position sample
#' @param minMvMs Minimum shift in ms 
#' @param samplingRate Sampling rate of the .dat files.
#' @examples shift.position.vector(x=1:100, timePerSamplesRes=400, minMvMs = 1000, samplingRate=20000)
shiftPositionVector<-function(x,
                         timePerSampleRes,
                         minMvMs,
                         samplingRate){
  minMv<- minMvMs*(samplingRate/1000)/timePerSampleRes
  mv<-sample(minMv:length(x)-minMv,1)
  return(shift(x,mv))
}

#' Shift a the values of two vectors by a random amount that is at least as large as the argument minMvMs
#' 
#' The two vectors are shifted by the same amount.
#' 
#' @param x A vector
#' @param y A second vector
#' @param timePerSampleRes Time in sample values (from the .dat files) between the position sample
#' @param minMvMs Minimum shift in ms 
#' @param samplingRate Sampling rate of the .dat files.
#' @examples shift.position.vector(x=1:100,y=201:300, timePerSamplesRes=400, minMvMs = 1000, samplingRate=20000)
shiftPositionVectors<-function(x,y,
                         timePerSampleRes,
                         minMvMs,
                         samplingRate){
  minMv<- minMvMs*(samplingRate/1000)/timePerSampleRes
  mv<-sample(minMv:length(x)-minMv,1)
  x<-shift(x,mv)
  y<-shift(y,mv)
  return(list(x=x,y=y))
}
