#' A S4 class to represent a recording session.
#' 
#' A RecSession object contains a description of the recording session.
#' 
#' Sessions are usually divided in several trials with one .dat file for each trial.
#' 
#' 
#' @slot session Name of the recording session
#' @slot path Directory where the recording session data are located
#' @slot fileBase Filebase of the session
#' @slot animalName Name of the animal
#' @slot samplingRate Sampling rate of the electrophysiological data
#' @slot resofs Number of samples in each trial
#' @slot env List of environment for each trial
#' @slot electrodeLocation List of electrode location, one per electrode
#' @slot trialStartRes Sample at which a trial starts. Index starts at 0
#' @slot trialEndRes Sample at which a trial ends. Index starts at 0
#' @slot trialNames Filebase of the individual trials
#' @slot trialDurationSec Length of the trials in sec
#' @slot sessionDurationSec Total length of the session in sec
#' @slot nElectrodes Number of electrodes
#' @slot nChannels Number of channels in the electrophysiological recordings
#' @slot nTrials Number of trials in the session
#' @slot channelsTetrode Matrix containing the channel numbers associated with each tetrode
#' @slot clustered Logical indicating if the spikes are clustered
#' @slot earlyProcessed Logical indicating if spike extraction has been done
RecSession <- setClass(
  "RecSession", ## name of the class
  slots=c(session="character",
          path="character",
          fileBase="character", # path+session
          animalName="character",
          samplingRate="numeric",
          resofs="numeric",
          env="character",
          electrodeLocation="character",
          trialStartRes="numeric",
          trialEndRes="numeric",
          trialNames="character",
          trialDurationSec="numeric",
          sessionDurationSec="numeric",
          nElectrodes="numeric",
          nChannels="numeric",
          nTrials="numeric",
          channelsTetrode="matrix",
          clustered="logical",
          earlyProcessed="logical"),  # cell list to limit the analysis to these cells
  prototype = list(session="",path=""))

#' Load the data regarding a recording session
#'
#' You need to set the session and  path slots of your object before calling this function
#' This will read the .par, .desen, .desel, .sampling_rate_dat, .resofs files
#' to fill up the RecSession object
#'
#' @param rs A RecSession object
#' @return RecSession
#' 
#' @docType methods
#' @rdname loadRecSession-methods
setGeneric(name="loadRecSession",
           def=function(rs)
           {standardGeneric("loadRecSession")}
)
#' @rdname loadRecSession-methods
#' @aliases loadRecSession,ANY,ANY-method
setMethod(f="loadRecSession",
          signature="RecSession",
          definition=function(rs)
          {
            if(rs@session=="")
              stop("rs@session is empty")
            if(rs@path=="")
              rs@path=getwd()
            
            rs@clustered=FALSE
            rs@earlyProcessed<-FALSE
            rs@fileBase<-paste(rs@path,rs@session,sep="/")
            rs@animalName<-unlist(strsplit(rs@session,"-"))[1]
            
            if(!file.exists(paste(rs@fileBase,"par",sep=".")))
              stop("needs ",paste(rs@fileBase,"par",sep="."))
            
            ## read the par file line per line## shitty format
            conn <- file(paste(rs@fileBase,"par",sep="."),open="r")
            par<-readLines(conn)
            close(conn)
            rs@nChannels<-as.numeric(unlist(strsplit(par[1], split=" "))[1])
            rs@nElectrodes  <-as.numeric(unlist(strsplit(par[3], split=" "))[1])
            rs@nTrials<-as.numeric(par[rs@nElectrodes+4])
            rs@trialNames<-par[(rs@nElectrodes+5):(rs@nElectrodes+5+rs@nTrials-1)]
            
            ## map of channel and tetrodes
            chan<-strsplit(par[4:(4+rs@nElectrodes-1)], split=" ")
            max.channelsTetrode<-max(unlist(lapply(chan,length))-1)
            rs@channelsTetrode<-matrix(nrow=rs@nElectrodes,ncol=max.channelsTetrode)
            
            if(rs@nElectrodes>0) { 
            for(i in 1:rs@nElectrodes) {
              l1<-length((rs@channelsTetrode[i,]))
              l2<-length((as.numeric(chan[[i]][-1])))
              if(l2>l1)
                stop(paste("loadRecSession, problem with number of channels on a tetrode",rs@session))
              rs@channelsTetrode[i,1:l2]<-as.numeric(chan[[i]][-1])
            }
            } else {
              print("No channelsTetrode")
              rs@channelsTetrode<-matrix(NA)
            }
            
            if(file.exists(paste(rs@fileBase,"desen",sep="."))){
              rs@env<-as.character(read.table(paste(rs@fileBase,"desen",sep="."))$V1)
              if(length(rs@env)!=length(rs@trialNames))
                stop(paste("loadRecSession, problem with length of par and desen files",rs@session))
            }
            if(file.exists(paste(rs@fileBase,"desel",sep="."))){
                rs@electrodeLocation<-as.character(read.table(paste(rs@fileBase,"desel",sep="."))$V1)
              if(rs@nElectrodes!=length(rs@electrodeLocation))
                stop(paste("loadRecSession, problem with length of par and desel files",rs@session))
            }
            
            if(rs@nTrials!=length(rs@trialNames))
              stop(paste("loadRecSession, problem with number of trials in par file",rs@session))

            if(file.exists(paste(rs@fileBase,"sampling_rate_dat",sep="."))){
              rs@samplingRate<-read.table(paste(rs@fileBase,"sampling_rate_dat",sep="."))$V1
              if(rs@samplingRate<1 | rs@samplingRate > 100000)
                stop(paste("loadRecSession, samplingRate is out of range:",rs@samplingRate,rs@session))
            }
            
            ## if early process was run on this one, get more informaiton
            if(file.exists(paste(rs@fileBase,"resofs",sep=".")))
            {
              ## read the resofs file  
              rs@resofs<-read.table(paste(rs@fileBase,"resofs",sep="."))$V1
              if(length(rs@resofs)!=rs@nTrials)
                stop(paste("loadRecSession, problem with length of resofs",rs@session))
              ## trial times
              rs@trialStartRes<-c(0,rs@resofs[-length(rs@resofs)])
              rs@trialEndRes<-rs@resofs-1 # resofs is the number of samples in one file, index is -1 that
              rs@trialDurationSec<-(rs@trialEndRes-rs@trialStartRes)/rs@samplingRate
              rs@sessionDurationSec<-sum(rs@trialDurationSec)
            }
            if(file.exists(paste(rs@fileBase,"clu",sep="."))&
               file.exists(paste(rs@fileBase,"res",sep="."))) rs@clustered=T
            
            if(file.exists(paste(rs@fileBase,"resofs",sep=".")))rs@earlyProcessed=T
            return(rs)
          }
)

#' Is the recording session clustered?
#'
#' @param rs A RecSession object
#' @return TRUE or FALSE
#' 
#' @docType methods
#' @rdname getIsClustered-methods
setGeneric(name="getIsClustered",
           def=function(rs)
           {standardGeneric("getIsClustered")}
)
#' @rdname getIsClustered-methods
#' @aliases getIsClustered,ANY,ANY-method
setMethod(f="getIsClustered",
          signature="RecSession",
          definition=function(rs)
          {
            return(rs@clustered)
          })


#' Has the spike extraction been run on the recording session?
#'
#' @param rs A RecSession object
#' @return TRUE or FALSE
#' 
#' @docType methods
#' @rdname getIsEarlyProcessed-methods
setGeneric(name="getIsEarlyProcessed",
           def=function(rs)
           {standardGeneric("getIsEarlyProcessed")}
)
#' @rdname getIsEarlyProcessed-methods
#' @aliases getIsEarlyProcessed,ANY,ANY-method
setMethod(f="getIsEarlyProcessed",
          signature="RecSession",
          definition=function(rs)
          {
            return(rs@earlyProcessed)
          })


#' Check if the session had an electrode in a particular brain area
#'
#' This will check whether the value of location is in the electrodeLocation vector.
#' 
#' @param rs A RecSession object
#' @param location A brain area of interest
#' @return TRUE or FALSE
#' 
#' @docType methods
#' @rdname containsElectrodeLocation-methods
setGeneric(name="containsElectrodeLocation",
           def=function(rs,location="")
           {standardGeneric("containsElectrodeLocation")}
)
#' @rdname containsElectrodeLocation-methods
#' @aliases containsElectrodeLocation,ANY,ANY-method
setMethod(f="containsElectrodeLocation",
          signature="RecSession",
          definition=function(rs,location="")
          {
            return(any(rs@electrodeLocation==location))
          })

#' Check if the session had a trial in a given environment
#'
#' This will check whether the value of environment is in the env vector.
#' 
#' @param rs A RecSession object
#' @param environment The name of an environment
#' @return TRUE or FALSE
#' 
#' @docType methods
#' @rdname containsEnvironment-methods
setGeneric(name="containsEnvironment",
           def=function(rs,environment="")
           {standardGeneric("containsEnvironment")}
)
#' @rdname containsEnvironment-methods
#' @aliases containsEnvironment,ANY,ANY-method
setMethod(f="containsEnvironment",
          signature="RecSession",
          definition=function(rs,environment="")
          {
            return(any(rs@env==environment))
          })

#' Get the recording date of a recSession, taken from session name
#'
#' @param rs A RecSession object
#' @return Date
#' 
#' @docType methods
#' @rdname recordingDate-methods
setGeneric(name="recordingDate",
           def=function(rs)
           {standardGeneric("recordingDate")}
)
#' @rdname recordingDate-methods
#' @aliases recordingDate,ANY,ANY-method
setMethod(f="recordingDate",
          signature="RecSession",
          definition=function(rs)
          {
            d<-unlist(strsplit(x=rs@session,split="-"))[2]
            if(nchar(d)==8)
            {
              rDate<-as.Date(d, "%d%m%Y")
            } else{
              rDate<-NA
            }
            return(rDate)
})


#' Get the time intervals in sample values for trials in a given environment
#'
#' @param rs A RecSession object
#' @param environment The name of an environment
#' @return matrix with 2 columns containing the start and end of each trial in the environment
#' 
#' @docType methods
#' @rdname getIntervalsEnvironment-methods
setGeneric(name="getIntervalsEnvironment",
           def=function(rs,environment="lt")
           {standardGeneric("getIntervalsEnvironment")}
)
#' @rdname getIntervalsEnvironment-methods
#' @aliases getIntervalsEnvironment,ANY,ANY-method
setMethod(f="getIntervalsEnvironment",
          signature="RecSession",
          definition=function(rs,environment="lt")
          {
            if(length(rs@trialStartRes)==0){
              print("trialStartRes is not set")
              return()
            }
            if(!environment%in%rs@env){
              print("environment not used in the session")
              return()
            }
          return(matrix(data=c(rs@trialStartRes[which(rs@env==environment)],rs@trialEndRes[which(rs@env==environment)]),ncol=2,
                 dimnames=list(rep(environment,length(which(rs@env==environment))),c("start","end"))))
          })

#' Load a set of objects that are session specific
#'
#' This is used to get the objects that are most commonly needed when doing analysis. 
#' Instead of creating the object each at a time and having to set some values manually,
#' just call this function and a list of objects are returned.
#' 
#' If you just want to use one or two objects, the code might run faster if you load what you need manually.
#'
#' @param rs A RecSession object
#' @return a list of objects containing spike trains, position data, data files, cell groups, 
#' spatial properties, 1D spatial properties, head direction data.
#' 
#' @docType methods
#' @rdname getRecSessionObjects-methods
setGeneric(name="getRecSessionObjects",
           def=function(rs)
           {standardGeneric("getRecSessionObjects")}
)
#' @rdname getRecSessionObjects-methods
#' @aliases getRecSessionObjects,ANY,ANY-method
setMethod(f="getRecSessionObjects",
          signature="RecSession",
          definition=function(rs)
          {
            if(rs@session=="")
              stop("rs@session is empty")
            if(rs@path=="")
              stop("rs@path not set")
            if(rs@clustered==FALSE)
              stop("rs is not clustered")
              
            st<-new("SpikeTrain",session=rs@session,path=rs@path)
            st<-loadSpikeTrain(st)
            pt<-new("Positrack",session=rs@session,path=rs@path)
            pt<-loadPositrack(pt)
            df<-new("DatFiles",fileNames=paste(rs@trialNames,"dat",sep="."),path=rs@path,nChannels=rs@nChannels)
            cg<-new("CellGroup",session=rs@session,path=rs@path,nTetrodes=rs@nElectrodes)
            cg<-loadCellGroup(cg)
            sp<-new("SpatialProperties2d",session=rs@session)
            sp1<-new("SpatialProperties1d",session=rs@session)
            hd<-new("HeadDirection",session=rs@session)
            return(list(st=st,pt=pt,df=df,cg=cg,sp=sp,sp1=sp1,hd=hd))
          })


#' Make a copy of the files of the recording session in another directory
#'
#' This is used to export the data of an experiment. 
#'
#' @param rs A RecSession object
#' @param destination A directory in which to do the backup. 
#' If /data is given for destination, the data go in the directory /data/animalName/sessionName
#' @param sessionSpecificExtensions List of file extensions. These are in the format session.extension
#' @param tetrodeSpecificExtensions List of file extensions. These are in the format session.extension.tetrodeNo
#' 
#' @docType methods
#' @rdname copyRecSessionFiles-methods
setGeneric(name="copyRecSessionFiles",
           def=function(rs,destination,sessionSpecificExtensions,tetrodeSpecificExtensions)
           {standardGeneric("copyRecSessionFiles")}
)
#' @rdname copyRecSessionFiles-methods
#' @aliases copyRecSessionFiles,ANY,ANY-method
setMethod(f="copyRecSessionFiles",
          signature="RecSession",
          definition=function(rs,destination,sessionSpecificExtensions,tetrodeSpecificExtensions)
          {
            if(rs@session=="")
              stop("rs@session is empty")
            if(rs@path=="")
              stop("rs@path not set")
           
            print(paste("copyRecSessionFiles",rs@session))
            ## create mouse directory
            mouseDestination=paste(destination,rs@animalName,sep="/")
            if(!dir.exists(mouseDestination)){
              print(paste("Creating",mouseDestination))
              dir.create(mouseDestination)
            }
            ## create session directory
            sessionDestination=paste(destination,rs@animalName,rs@session,sep="/")
            if(!dir.exists(sessionDestination)){
              print(paste("Creating",sessionDestination))
              dir.create(sessionDestination)
            }  
            ## check that session specific files exists
            fileNames<-paste(rs@fileBase,sessionSpecificExtensions,sep=".")
            if(any(!file.exists(fileNames))){
              print(paste("file missing:",fileNames[!file.exists(fileNames)]))
              stop("copyRecSessionFiles, missing files")
            }
            ## check that tetrode specific files exists
            fileNames<-paste(rs@fileBase,tetrodeSpecificExtensions,sep=".")
            tet<-rep(1:rs@nElectrodes,each=length(fileNames))
            fileNames<-paste(fileNames,tet,sep='.')
            if(any(!file.exists(fileNames))){
              print(paste("file missing:",fileNames(!file.exists(fileNames))))
              stop("copyRecSessionFiles, missing files")
            }
            ## copy session specific files
            fileNames<-paste(rs@fileBase,sessionSpecificExtensions,sep=".")
            newFileNames<-paste(paste(destination,rs@animalName,rs@session,rs@session,sep="/"),sessionSpecificExtensions,sep=".")
            if(any(!file.copy(from = fileNames,to = newFileNames,overwrite = T,recursive = F)))
              stop(paste("copy of files failed in", rs@session))
            ## copy tetrode specific files
            fileNames<-paste(rs@fileBase,tetrodeSpecificExtensions,sep=".")
            tet<-rep(1:rs@nElectrodes,each=length(fileNames))
            fileNames<-paste(fileNames,tet,sep='.')
            newFileNames<-paste(paste(destination,rs@animalName,rs@session,rs@session,sep="/"),tetrodeSpecificExtensions,sep=".")
            tet<-rep(1:rs@nElectrodes,each=length(newFileNames))
            newFileNames<-paste(newFileNames,tet,sep='.')
            if(any(!file.copy(from = fileNames,to = newFileNames,overwrite = T,recursive = F)))
              stop(paste("copy of files failed in", rs@session))
})



### show ###
setMethod("show", "RecSession",
          function(object){
            print(paste("session:",object@session))
            print(paste("path:",object@path))
            print(paste("samplingRate:",object@samplingRate,"Hz"))
            print(paste("nChannels:",object@nChannels))
            print(paste("nTrials:",object@nTrials))
            print(paste("env:"))
            print(paste(object@env))
            print(paste("nElectrodes:",object@nElectrodes))
            print(paste("electrodeLocation:"))
            print(paste(object@electrodeLocation))
            print(paste("trialNames:"))
            print(paste(object@trialNames))
            print(paste("trialDurationSec:"))
            print(paste(object@trialDurationSec,"sec"))
            print(paste("session duration:",object@sessionDurationSec,"sec"))
            print(paste("trialStartRes and trialEndRes:"))
            print(paste(object@trialStartRes,object@trialEndRes))
            print(paste("Map of channels per tetrode (channelsTetrode)"))
            print(object@channelsTetrode)
            print(paste("clustered:",object@clustered))
            print(paste("earlyProcessed:",object@earlyProcessed))
          })



#' Get animal name from session name
#' 
#' Assumes the session name is in the format name-date-rest
#' 
#' @param sessionName Character vector with the session name
#' @return Character vector with animal name
animalNameFromSessionName<-function(sessionName=NULL){
    return(unlist(lapply(strsplit(as.character(sessionName),split="-"),function(x){return(x[[1]])})))
}

#' Get session name from cluId
#' 
#' Assumes the session name is in the format name-date-rest and cluId is name-data-rest_cluId
#' 
#' @param cluId Character vector with the session name
#' @return Character vector with session name
sessionNameFromCluId<-function(cluId=NULL){
  return(unlist(lapply(strsplit(as.character(cluId),split="_"),function(x){return(x[[1]])})))
}

