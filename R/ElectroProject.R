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
          resultsDirectory="character",
          nSessions="numeric",
          sessionNameList="character",
          sessionPathList="character",
          sessionList="list"),
  prototype = list(directory=""))




#' Create a list of RecSession objects for the ElectroProject
#'
#' Will list directories in the project directories
#' The one with a hyphen will be considered recSession directories.
#'
#' @param ep ElectroProject object
#' @param loadSessions logical, whether the RecSession object should be created
#' @return ElectroProject object with a list of session names and directories
#' 
#' @docType methods
#' @rdname setSessionList-methods
setGeneric(name="setSessionList",
           def=function(ep,loadSessions=TRUE)
           {standardGeneric("setSessionList")}
)
#' @rdname setSessionList-methods
#' @aliases setSessionList,ANY,ANY-method
setMethod(f="setSessionList",
          signature="ElectroProject",
          definition=function(ep,loadSessions=TRUE)
          {
            if(ep@directory=="")
              stop("ep@directory not set")
            
            ep@resultsDirectory=paste(ep@directory,"results",sep="/")
            
            ## list all directories in the project path
            dirs<-list.dirs(path=ep@directory)
            
            ## only keep the directory with a hyphen in the name
            ep@sessionPathList<-dirs<-dirs[grepl(pattern="-",dirs)]
            ep@sessionPathList
            ep@nSessions<-length(dirs)
            
            ## make sure all directories have the same depth
            x<-sapply(ep@sessionPathList,function(x){length(unlist(strsplit(x,split="/")))})
            if(any(x!=x[1])){
              stop("the depth of some session directories differs")
            }
            
            dirDepth<-as.numeric(x[1])
            ep@sessionNameList<-unlist(strsplit(dirs,split="/"))[seq(from=dirDepth,to=dirDepth*ep@nSessions,by=dirDepth)]
            if(loadSessions){
              loadSessionsInList(ep)
            }
            return(ep)
          })


#' Create a List of RecSessions objects from the SessionNameList
#' 
#' The RecSession objects will be initialized
#'
#' @param ep ElectroProject object
#' @return ElectroProject object with RecSession loaded
#' 
#' @docType methods
#' @rdname setSessionList-methods
setGeneric(name="loadSessionsInList",
           def=function(ep)
           {standardGeneric("loadSessionsInList")}
)
#' @rdname loadSessionsInList-methods
#' @aliases loadSessionsInList,ANY,ANY-method
setMethod(f="loadSessionsInList",
          signature="ElectroProject",
          definition=function(ep)
          {
            if(ep@directory=="")
              stop("ep@directory not set")
            
            if(length(ep@sessionNameList)==0)
              stop("ep@sessionNameList has length of 0")
            if(length(ep@sessionPathList)==0)
              stop("ep@sessionPathList has length of 0")
            
            for (i in 1:length(ep@sessionNameList))
            {
              ep@sessionList[i]<-new("RecSession",session=ep@sessionNameList[i],path=ep@sessionPathList[i])
            }
            ep@sessionList<-lapply(ep@sessionList,loadRecSession)     
            
            return(ep)
          })



#' Return a list of clustered RecSession objects
#'
#'
#' @param ep ElectroProject object
#' @return list of clustered RecSession objects
#' 
#' @docType methods
#' @rdname getClusteredSessionList-methods
setGeneric(name="getClusteredSessionList",
           def=function(ep)
           {standardGeneric("getClusteredSessionList")}
)
#' @rdname getClusteredSessionList-methods
#' @aliases getClusteredSessionList,ANY,ANY-method
setMethod(f="getClusteredSessionList",
          signature="ElectroProject",
          definition=function(ep)
          {
            if(ep@directory=="")
              stop("ep@directory not set")
            return(ep@sessionList[sapply(ep@sessionList,getIsClustered)])
          })




#' Return a list of RecSession objects that have some common properties
#'
#'
#' @param ep ElectroProject object
#' @param clustered logical indicating whether the session should be clustered or not
#' @param region Set to a given brain region to select only sessions with tetrodes in this brain region
#' @param env Set to a given environment code to select only sessions during which this environment was presented
#' @return list of RecSession objects
#' 
#' @docType methods
#' @rdname getSessionList-methods
setGeneric(name="getSessionList",
           def=function(ep,clustered="",region="",env="")
           {standardGeneric("getSessionList")}
)
#' @rdname getSessionList-methods
#' @aliases getSessionList,ANY,ANY-method
setMethod(f="getSessionList",
          signature="ElectroProject",
          definition=function(ep,clustered="",region="",env="")
          {
            if(ep@directory=="")
              stop("ep@directory not set")
            if(length(ep@sessionList)==0)
              stop("sp@sessionList has a length of 0")
            myList<-ep@sessionList
            length(myList)
            if(clustered==T)
              myList<-myList[sapply(myList,getIsClustered)]
            length(myList)
            if(region!=""){
              myList<-myList[sapply(myList,containsElectrodeLocation,location=region)]
            }
            if(env!=""){
              myList<-myList[sapply(myList,containsEnvironment,environment=env)]
            }
            return(myList)
          })

















 
#' Running a function on a set of recording sessions
#'
#' This applies a function to a list of RecSession objects.
#' If save is set to TRUE, the results returned by the function will be saved. 
#' Not that the function should return the results in a list. 
#' The results are saved in the resultsDirectory of the ElectroProject object.
#' The names of the files saved will be the name of the elements in the list returned by the function
#' The data from each recording session will be concatenated.
#' You can use parLapply instead of lapply by setting parallel to TRUE and passing a valid cluster to the function
#' If save is set to FALSE, the data returned by the function will not be saved.
#'
#' @param ep ElectroProject object
#' @param sessionList List of RecSession objects on which the function will be applied
#' @param fnct A function to run on each RecSession
#' @param save Whether you want to save the data returned by the function
#' @param overwrite Whether you want to overwrite the previous data when saving the results
#' @param parallel Whether you want to run the function in parallel
#' @param cluster A cluster generated from the makeCluster function of the snow package
#'
#' 
#' @docType methods
#' @rdname runOnSessionList-methods
setGeneric(name="runOnSessionList",
            def=function(ep,sessionList,fnct=function(x){NA},save=T,overwrite=T,parallel=F,cluster="")
            {standardGeneric("runOnSessionList")}
          )

#' @rdname runOnSessionList-methods
#' @aliases runOnSessionList,ANY,ANY-method
setMethod(f="runOnSessionList",
         signature="ElectroProject",
         definition=function(ep,sessionList,fnct=function(x){NA},save=T,overwrite=T,parallel=F,cluster="")
         {
           if(class(ep)!="ElectroProject")
             stop("runOnSessionList, ep should be a ElectroProject object")
           if(!is.list(sessionList))
             stop("runOnSessionList, sessionList has to be a list")
           if(length(sessionList)==0)
             stop("runOnSessionList, sessionList has size 0")
           if(class(sessionList[[1]])!="RecSession")
             stop("runOnSessionList, sessionList should contain RecSession objects")
           if(!is.function(fnct))
             stop("runOnSessionList, fnct needs to be a function")
           if(parallel==T)
             if(class(cluster)[2]!="cluster")
               stop("runOnSessionList, give a valid snow cluster if you want to run the function on several threads")
           if(parallel==T){
             list.res<-parLapply(cluster,sessionList,fnct)   
           } else {
             list.res<-lapply(sessionList,fnct)
           }
           
           if(save==T){
             ## check that list.res is a list of list
             if(!is.list(list.res))
               stop("runOnSessionList, list.res is not a list")
             if(!is.list(list.res[[1]]))
               stop(paste("runOnSessionList, first item of the list list.res is not a list, fnct should return a list"))
             
             ## list of objects to merge and save
             objectNames<-names(list.res[[1]])
             for(n in objectNames){
               
               if(overwrite==T){assign(n,do.call("rbind", sapply(list.res,function(x){x[n]})))
               }else{## concatonate to existing data
                 assign(paste(n,"new",sep="."),do.call("rbind", sapply(list.res,function(x){x[n]})))
                 load(file=paste(ep@resultsDirectory,n,sep="/"))
                 assign(n,rbind(get(n),get(paste(n,"new",sep="."))))
               }
               print(paste("saving",paste(ep@resultsDirectory,n,sep="/")))
               save(list=n,file=paste(ep@resultsDirectory,n,sep="/"))
             }
           }
         })

### show ###
setMethod("show", "ElectroProject",
          function(object){
            print(paste("directory:",object@directory))
            print(paste("nSessions:",object@nSessions))
            if(length(object@sessionNameList)!=0){
              print(paste("sessionNameList:"))
              print(object@sessionNameList)
              if(length(object@sessionList)!=0){
                m<-matrix(c(sapply(object@sessionList,getIsClustered),sapply(object@sessionList,getIsEarlyProcessed)),ncol=2)
                print(paste("Clustered sessions:",sum(m[,1])))
                print(object@sessionNameList[m[,1]])
                print(paste("Not clustered, but early processed sessions:",length(which(m[,1]==F&m[,2]==T))))
                print(object@sessionNameList[which(m[,1]==F&m[,2]==T)])
                print(paste("Not early processed sessions:", length(which(m[,1]==F&m[,2]==F))))
                print(object@sessionNameList[which(m[,1]==F&m[,2]==F)])
              }
            }
          })
