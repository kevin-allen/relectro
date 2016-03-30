############################################
### definition of ElectroProject Class   ###
############################################
#' A S4 class to represent an electrophysiological project containing several recording sessions.
#' 
#' An ElectroProject object can be used to run programs on several recording sessions.
#' It contains a list of recording sessions. You can call lapply to run a function on each recording session.
#' If you want to speed up things, you can try parLapply of the snow package.
#' It can also be used to give you an overview of the progress during the data acquisition period.
#' 
#' @slot directory A path to a directory tree where the recording sessions can be found
#' @slot nSessions Number of recording sessions in the project
#' @slot sessionNameList Name of the recording sessions
#' @slot sessionList List of RecSession objects.
ElectroProject <- setClass(
  "ElectroProject", ## name of the class
  slots=c(directory="character",
          nSessions="numeric",
          sessionNameList="character",
          sessionList="list"),
  prototype = list(directory=""))

### setSessionList ###
setGeneric(name="setSessionList",
           def=function(ep)
           {standardGeneric("setSessionList")}
)
setMethod(f="setSessionList",
          signature="ElectroProject",
          definition=function(ep)
          {
            if(ep@directory=="")
              stop("ep@directory not set")
            
            ## list all directories in the project path
            dirs<-list.dirs(path=ep@directory)
            
            ## only keep the directory with a hyphen in the name
            dirs<-dirs[grepl(pattern="-",dirs)]
            ep@nSessions<-length(dirs)
            
            dirDepth<-length(strsplit(dirs,split="/")[[1]])
            ep@sessionNameList<-unlist(strsplit(dirs,split="/"))[seq(from=dirDepth,to=dirDepth*ep@nSessions,by=dirDepth)]
              
            for (i in 1:length(dirs))
            {
              ep@sessionList[i]<-new("RecSession",session=ep@sessionNameList[i],path=dirs[i])
            }
          ep@sessionList<-lapply(ep@sessionList,loadRecSession)  
          return(ep)
          })


### show ###
setMethod("show", "ElectroProject",
          function(object){
            print(paste("directory:",object@directory))
            print(paste("nSessions:",object@nSessions))
            print(paste("sessionNameList:"))
            print(object@sessionNameList)
            print(paste("Clustered sessions:",sum(sapply(ep@sessionList,getIsClustered))))
            print(object@sessionNameList[sapply(ep@sessionList,getIsClustered)])
            print(paste("Not clustered sessions:",sum(!sapply(ep@sessionList,getIsClustered))))
            print(object@sessionNameList[!sapply(ep@sessionList,getIsClustered)])
            print(paste("Early processed sessions:", sum(sapply(ep@sessionList,getIsEarlyProcessed))))
            print(paste("Not early processed sessions:", sum(!sapply(ep@sessionList,getIsEarlyProcessed))))
            print(object@sessionNameList[!sapply(ep@sessionList,getIsEarlyProcessed)])
          })

setGeneric(name="getClusteredSessionList",
           def=function(ep)
           {standardGeneric("getClusteredSessionList")}
)
setMethod(f="getClusteredSessionList",
          signature="ElectroProject",
          definition=function(ep)
          {
            if(ep@directory=="")
              stop("ep@directory not set")
            return(ep@sessionList[sapply(ep@sessionList,getIsClustered)])
          })
# 
# setGeneric(name="runOnSessionList",
#            def=function(ep,...)
#            {standardGeneric("runOnSessionList")}
# )
# setMethod(f="runOnSessionList",
#           signature="ElectroProject",
#           definition=function(ep,)
#           {
#             
#           })

## should be a method of ElectroProject
runOnSessionList<-function(sessionList,fnct=function(x){NA},overwrite=T,parallel=F,cluster=""){
  if(!is.list(sessionList))
    stop("runOnSessionList needs a list as first argument")
  if(length(sessionList)==0)
    stop("runOnSessionList, sessionList has size 0")
  if(!is.function(fnct))
    stop("runOnSessionList, fnct needs to be a function")
  if(parallel==T&cluster="")
    stop("runOnSessionList, give a valid snow cluster if you want to run the function on several threads")
  if(parallel==T){
    list.res<-parLapply(cluster,sessionList,fnct)   
  } else {
    list.res<-lapply(sessionList,fnct)
  }
  
  ## check that list.res is a list of list
  if(!is.list(list.res))
    stop("runOnSessionList, list.res is not a list")
  if(!is.list(list.res[[1]]))
    stop(paste("runOnSessionList, first item of the list list.res is not a list, fnct should return a list"))
  
  ## list of objects to merge and save
  objectNames<-names(list.res[[1]])
   
  if(overwrite==T){
    for(n in objectNames){
      #merge and assig
      #save in results directory
    }
  } else{
    for(n in objectNames){
      # merge and assign to a new name
      # load the existing object
      # remove data from same cells
      # merge new and existing object
      # save object
    }
  }
}
