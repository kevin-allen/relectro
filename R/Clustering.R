#' Check if clusters of a recording session are well isolated from one another
#' 
#' This function generate a data.frame with the isolation distance and refractory period ratio for all clusters.
#' It is returned within a list to conform to runOnSessionList policy on returned object.
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
  cluster.check<-data.frame(
      session=rs@session,
      clu=cg@clu,
      cluId=cg@id,
      tetrode=cg@tetrode,
      cluToTetrode=cg@cluToTetrode,
      meanFiringRate=st@meanFiringRate,
      isolationDistance=st@isolationDistance,
      refractoryRatio=st@refractoryRatio,
      crossRefractoryRatio=st@crossRefractoryRatio,
      check=check)
  return(list(cluster.check=cluster.check))
}


#' Delete a specific cluster in a recording session
#' 
#' This will delete the cluster in the tetrode specific clu file
#' The deleted cluster is assigned a value of 1 in the clu file.
#' Don't forget to recreate the new session wide clu and res file with
#' the function mergeTetrodeSpecificResCluFiles()
#' 
#' @param rs RecSession
#' @param clu Cluster number. This is the session wide cluster number.
#' Not the tetrode specific number
#' 
deleteCluster<-function(rs,clu)
{
  if(class(rs)!="RecSession")
    stop("deleteCluster, the class of rs should be RecSession")
  if(class(clu)!="numeric")
    stop("deleteCluster, the class of clu should be numeric")
  if(length(clu)!=1)
    stop("deleteCluster, length of clu should be 1")
  
  cg<-new("CellGroup",session=rs@session,path=rs@path,nTetrodes=rs@nElectrodes)
  cg<-loadCellGroup(cg)
  if(clu%in%cg@clu==FALSE)
    stop(paste("deleteCluster: cluster to delete is not in the list of clusters of the CellGroup object"))
  index<-which(cg@clu==clu)
  tetrode=cg@tetrode[index]
  cluToTet=cg@cluToTetrode[index]
  print(paste("Delete cluster",cluToTet, "from tetrode",tetrode))
  
  ## read clu file
  cfile=paste(paste(rs@path,rs@session,sep="/"),".clu.",tetrode,sep="")
  if(!file.exists(cfile))
    stop(paste("deleteCluster:",cfile,"missing"))
  c<-as.numeric(.Call("read_one_column_int_file_cwrap", cfile))
  nClu<-c[1]
  c<-c[-1]
  
  ## delete the cluster
  c[which(c==cluToTet)]<- 1
  
  ## add the new cluster number at the head of clu data
  c<-c(nClu-1,c)

  ## save the new clu file
  write(x=as.integer(c),
        file=cfile,
        ncolumns = 1,
        append=F)
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
  if(class(rs)!="RecSession")
    stop("mergeTetrodeSpecificResCluFiles, the class of rs should be RecSession")
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
#' @param fileName File name
#' @return A matrix with the spike features, one spike per row
readFetFile<-function(fileName){
  if(!file.exists(fileName))
    stop("readFetFile: needs ",fileName)
  m<-.Call("read_fet_file_cwrap",
        fileName)
  return(m)
}
