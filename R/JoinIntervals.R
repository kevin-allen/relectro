#' Join two sets of intervals using an AND logic
#' 
#' To be included in the output interval, the time period needs to be included in both set of time intervals
#'  
#' @param s1 Numeric vector containing the start times of the first set of intervals
#' @param e1 Numeric vector containing the end times of the first set of intervals
#' @param s2 Numeric vector containing the start times of the second set of intervals
#' @param e2 Numeric vector containing the end times of the second set of intervals
#' @return matrix containing the start and end of the resulting intervals. First col is the start times. Second col is end times.
#' @examples 
#' s1<-c(0,20000)
#' e1<-c(500,25000)
#' s2<-c(250,24000)
#' e2<-c(600,30000)
#' joinIntervalsAND(s1,e1,s2,e2)
joinIntervalsAND<-function(s1,e1,s2,e2){
  ## check args
  if(length(s1)!=length(e1))
    stop("unequal length of first intervals")
  if(length(s2)!=length(e2))
    stop("unequal lenght of the second intervals")
  if(any(s1>e1))
    stop("problem with chronology of the first set of intervals")
  if(any(s2>e2))
    stop("problem with chronology of the first set of intervals")
  res<-.Call("joinIntervalsAND_cwrap",
        as.integer(s1),
        as.integer(e1),
        length(s1),
        as.integer(s2),
        as.integer(e2),
        length(s2))
  colnames(res)<-c("start","end")
  return(res)
}
