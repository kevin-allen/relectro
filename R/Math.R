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

shuffle.vector<-function(x,
                         time.per.sample.res,
                         min.mv.ms,
                         sampling_rate){
  min.mv<- min.mv.ms*(sampling_rate/1000)/time.per.sample.res
  mv<-sample(min.mv:length(x)-min.mv,1)
  return(shift(x,mv))
}

shuffle.vectors<-function(x,y,
                         time.per.sample.res,
                         min.mv.ms,
                         sampling_rate){
  min.mv<- min.mv.ms*(sampling_rate/1000)/time.per.sample.res
  mv<-sample(min.mv:length(x)-min.mv,1)
  x<-shift(x,mv)
  y<-shift(y,mv)
  return(list(x=x,y=y))
}