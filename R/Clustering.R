#' Check if clusters of a recording session are well isolated from other clusters
#' 
#' This function generate a data.frame with the isolation distance and refractory period ratio for all clusters.
#' It is returned within a list to conform to runOnSessionList policy on returned object
#'  
#' @param rs RecSession
#' @return A list including isolation distance and refractory ratio of all clusters. The check element is a
#' logical indicating which clusters are likely to be poorly isolated.
clusterIsolationCheck<-function(rs){
  if(class(rs)!="RecSession")
    stop("clusterIsolationCheck: rs is not a RecSession")
  print(rs@session)
 
  st<-new("SpikeTrain",session=rs@session,path=rs@path)
  st<-loadSpikeTrain(st)
  cg<-new("CellGroup",session=rs@session,path=rs@path,nTetrodes=rs@nElectrodes)
  cg<-loadCellGroup(cg)
  
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
  check<-ifelse(st@refractoryRatio>0.125|st@crossRefractoryRatio<0.15,TRUE,FALSE)
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


#' Merge tetrode specific clu and res file into a main res and clu file
#' 
#' This is run after the clustering of all tetrodes is done. It create
#' the main .res and .clu files that are used when doing subsequent 
#' analysis. Cluster 1 is noise. 
#' The clusters 1 from different tetrodes are merged. The other
#' clusters are not merged. They will be numbered from 2 to N, where N is 
#' the total number of clusters in the recording session.
#' 
#' @param rs RecSession
#' 
mergeTetrodeSpecificResCluFiles<-function(rs){
  mclu<-vector()
  mres<-vector()
  preClu=0
  for(t in 1:rs@nElectrodes)
  {
    rfile=paste(paste(rs@path,rs@session,sep="/"),".res.",t,sep="")
    cfile=paste(paste(rs@path,rs@session,sep="/"),".clu.",t,sep="")
    if(!file.exists(rfile))
      stop(paste("mergeTetrodeSpecificResCluFiles:",rfile,"missing"))
    if(!file.exists(cfile))
      stop(paste("mergeTetrodeSpecificResCluFiles:",cfile,"missing"))
    r<-as.numeric(.Call("read_one_column_int_file_cwrap", rfile))
    c<-as.numeric(.Call("read_one_column_int_file_cwrap", cfile))
    nClu<-c[1]
    c<-c[-1]
    if(length(c)!=length(r))
      stop(paste("mergeTetrodeSpecificResCluFiles: problem with length of",rfile,"and",cfile))
    if(nClu>1)
      c[which(c>1)]<-c[which(c>1)]+preClu

    mres<-c(mres,r)
    mclu<-c(mclu,c)
    
    preClu=preClu+nClu-1        
  }
 
  ## reorder according to time
  o<-order(mres)
  mres<-mres[o]
  mclu<-mclu[o]
  
  ## add number of clu at the head of clu data
  mclu<-c(preClu+1,mclu)
  
  ## save main clu and res files
  write(x=as.integer(mres),
        file=paste(paste(rs@path,rs@session,sep="/"),".res",sep=""),
        ncolumns = 1,
        append=F)
  write(x=as.integer(mclu),
        file=paste(paste(rs@path,rs@session,sep="/"),".clu",sep=""),
        ncolumns = 1,
        append=F)
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



