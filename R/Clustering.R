#' Check if clusters of a recording session are well isolated for other clusters
#' 
#' This uses spike-time autocorrelation, spike-time crosscorrelation and isolation Distance.
#'  
#' @param rs RecSession
#' @return A list including the clusters that need to be checked or removed
clusterIsolationCheck<-function(rs){
  if(class(rs)!="RecSession")
    stop("clusterIsolationCheck: rs is not a RecSession")
  
  myList<-getRecSessionObjects(rs)
  ## the object needed
  st<-myList$st
  cg<-myList$cg
  ##
  st<-meanFiringRate(st)
  st<-isolationDistance(st,cg)
  st<-refractoryRatio(st)
  st<-crossRefractoryRatio(st)
  st
  cg
  ## get crosscorrelation refractory

}

#' Read a fet file
#' 
#' Read a file with the spike features. Uses a c function for speed.
#'  
#' @param file File name
#' @return A matrix with the spike features, one spike per row
readFetFile<-function(fileName){
  if(!file.exists(fileName))
    stop("readFetFile: needs ",fileName)
  m<-.Call("read_fet_file_cwrap",
        fileName)
  return(m)
}



