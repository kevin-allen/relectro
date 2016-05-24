#' Check if clusters of a recording session are well isolated for other clusters
#' 
#' This function generate a list with the isolation distance and refractory ratio for all clusters.
#' The returned list with the values used to assess cluster isolation. 
#'  
#' @param rs RecSession
#' @return A list including isolation distance and refractory ratio of all clusters. The check element is a
#' logical indicating which clusters are likely to be poorly isolated.
clusterIsolationCheck<-function(rs){
  if(class(rs)!="RecSession")
    stop("clusterIsolationCheck: rs is not a RecSession")
  print(rs@session)
  myList<-getRecSessionObjects(rs)
  ## the object needed
  st<-myList$st
  cg<-myList$cg
  
  if(any(st@nSpikes<10)){
    print(paste("At least one cluster has fewer than 10 spikes.",rs@session))
    print("list of clusters with fewer than 10 spikes")
    print(paste(st@cellList[which(st@nSpikes<10)]))
    print(paste("tetrode:",cg@tetrode[which(st@nSpikes<10)],", clu:", cg@cluToTetrode[which(st@nSpikes<10)]))
    stop()
  }
    
  ##
  st<-meanFiringRate(st)
  st<-isolationDistance(st,cg)
  st<-refractoryRatio(st)
  st<-crossRefractoryRatio(st)
  check<-ifelse(st@refractoryRatio>0.15|st@crossRefractoryRatio<0.15,TRUE,FALSE)
  check[which(!is.na(st@isolationDistance)&st@isolationDistance<5)]<-TRUE
  
  ## plot the autocorrelation and crosscorrelation so we double check.
  
  cluster.check<-data.frame(cluId=cg@id,
       tetrode=cg@tetrode,
       cluToTetrode=cg@cluToTetrode,
       meanFiringRate=st@meanFiringRate,
       isolationDistance=st@isolationDistance,
       refractoryRatio=st@refractoryRatio,
       crossRefractoryRatio=st@crossRefractoryRatio,
       check=check)
  cluster.check
  return(list(cluster.check=cluster.check))
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



