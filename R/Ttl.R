#' Detect the positive deflection of a ttl pulse in a vector
#' 
#' @param x numeric vector containing the data
#' @param threshold amplitude of the positive deflections to be considered
#' @return indices at which the positive deflections were detected in \code{x}
#' @examples  detectUps(x=c(rep(0,10),rep(20000,10)))
detectUps<-function(x, threshold=10000)
{  # check if x[i+1]-x[i] > threshold
  
  if(class(x)=="integer") x<-as.numeric(x)
  .Call("detect_ttl_ups_cwrap",
          x,
          length(x),
          threshold)
}
#' Detect the negative deflection of a ttl pulse in a vector
#' 
#' @param x numeric vector containing the data
#' @param threshold amplitude of the negative deflections to be considered
#' @return indices at which the negative deflections were detected in \code{x}
#' @examples  detectDowns(x=c(rep(0,10),rep(20000,10)))
detectDowns<-function(x, threshold=10000)
{  # check if x[i+1]-x[i] < 0-threshold
  if(class(x)=="integer") x<-as.numeric(x)
  .Call("detect_ttl_downs_cwrap",
        x,
        length(x),
        threshold)
}
