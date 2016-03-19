#' Detect the positive deflection of a ttl pulse in a numiric vector
#' 
#' @param x numeric vector containing the data
#' @param threshold amplitude of the positive deflection to be considered
#' @return indices at which the positive deflections were detected in \code{x}
#' @examples  detectUps(x=c(rep(0,10),rep(20000,10)))
detectUps<-function(x, threshold=10000)
{  # check if x[i+1]-x[i] > threshold
  .Call("detect_ttl_ups_cwrap",
          x,
          length(x),
          threshold)
}
detectDowns<-function(x, threshold=10000)
{  # check if x[i+1]-x[i] < 0-threshold
  .Call("detect_ttl_downs_cwrap",
        x,
        length(x),
        threshold)
}
