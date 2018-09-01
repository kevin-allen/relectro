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
#' @slot environment List of environments for each trial
#' @slot stimulation List of stimulation types
#' @slot setup List of recording setups used in the recording session.
#' @slot environmentFamiliarity Indicate whether the environment was familiar (fam) or novel (nov) for the animal
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
#' @slot pxPerCm Numeric representing the number of pixels per cm in the position data
#' @slot resSamplesPerWhdSample Number of electrophysiological sample between two samples in the position file (whd)
RecSession <- setClass(
  "RecSession", ## name of the class
  slots=c(session="character",
          path="character",
          fileBase="character", # path+session
          animalName="character",
          samplingRate="numeric",
          resofs="numeric",
          environment="character",
          stimulation="character",
          setup="character",
          environmentFamiliarity="character",
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
          earlyProcessed="logical",
          pxPerCm="numeric",
          resSamplesPerWhdSample="numeric"),  # cell list to limit the analysis to these cells
  prototype = list(session="",path=""))

#' Load the data regarding a recording session
#'
#' You need to set the session and  path slots of your object before calling this function
#' This will read the .par, .desen, .desel, .sampling_rate_dat, .resofs files.
#' If the resofs is missing will try to get trial times from a datFiles object.
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
            
            ## add tests in case of weird .par file
            if(length(rs@nChannels)!=1)
              stop(paste("rs@nChannels is not set correctly for",rs@session))
            if(length(rs@nElectrodes)!=1)
              stop(paste("rs@nElectrodes is not set correctly for",rs@session))
            if(is.na(rs@nChannels))
              stop(paste("rs@nChannels is na for",rs@session))
            if(is.na(rs@nElectrodes))
              stop(paste("rs@nChannels is na for",rs@session))
            if(length(rs@nTrials)!=1)
              stop(paste("rs@nTrials is not set correctly for",rs@session))
            if(is.na(rs@nTrials))
              stop(paste("rs@nTrials is na for",rs@session))

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
              try(
                rs@environment<-as.character(read.table(paste(rs@fileBase,"desen",sep="."))$V1),
                silent=F)
              if(length(rs@environment)!=length(rs@trialNames))
                stop(paste("loadRecSession, problem with length of par and desen files",rs@session))
            }
            
            if(file.exists(paste(rs@fileBase,"stimulation",sep="."))){
              try(
                rs@stimulation<-as.character(read.table(paste(rs@fileBase,"stimulation",sep="."))$V1),
                silent=F)
              if(length(rs@stimulation)!=length(rs@trialNames))
                stop(paste("loadRecSession, problem with length of par and stimulation files",rs@session))
            }

            if(file.exists(paste(rs@fileBase,"setup",sep="."))){
              try(
                rs@setup<-as.character(read.table(paste(rs@fileBase,"setup",sep="."))$V1),
                silent=F)
              if(length(rs@setup)!=length(rs@trialNames))
                stop(paste("loadRecSession, problem with length of par and setup files",rs@session))
            }
            if(file.exists(paste(rs@fileBase,"environmentFamiliarity",sep="."))){
              try(
                rs@environmentFamiliarity<-as.character(read.table(paste(rs@fileBase,"environmentFamiliarity",sep="."))$V1),
                silent=F)
              if(length(rs@environmentFamiliarity)!=length(rs@trialNames))
                stop(paste("loadRecSession, problem with length of par and environmentFamiliarity files",rs@session))
            }
            
            if(file.exists(paste(rs@fileBase,"desel",sep="."))){
              try(
                rs@electrodeLocation<-as.character(read.table(paste(rs@fileBase,"desel",sep="."))$V1),
                silent=F)
              if(rs@nElectrodes!=length(rs@electrodeLocation))
                stop(paste("loadRecSession, problem with length of par and desel files",rs@session))
            }


            if(rs@nTrials!=length(rs@trialNames))
              stop(paste("loadRecSession, problem with number of trials in par file",rs@session))

            if(file.exists(paste(rs@fileBase,"sampling_rate_dat",sep="."))){
              try(rs@samplingRate<-read.table(paste(rs@fileBase,"sampling_rate_dat",sep="."))$V1,
                  silent=F)
              if(length(rs@samplingRate)>1)
                stop(paste("loadRecSession, samplingRate has a length > 1, check",paste(rs@fileBase,"sampling_rate_dat",sep=".")))
              if(rs@samplingRate<1 | rs@samplingRate > 100000)
                stop(paste("loadRecSession, samplingRate is out of range:",rs@samplingRate,rs@session))
            }

            ## get the number of pixels per cm in tracking data
            if(file.exists(paste(rs@fileBase,"px_per_cm",sep="."))){
              rs@pxPerCm<-read.table(paste(rs@fileBase,"px_per_cm",sep="."))$V1
            }
            
            ## get the number of res samples per whd sample in the tracking data
            if(file.exists(paste(rs@fileBase,"res_samples_per_whd_sample",sep="."))){
              rs@resSamplesPerWhdSample<-read.table(paste(rs@fileBase,"res_samples_per_whd_sample",sep="."))$V1
            }
            
            
            ## if early process was run on this one, get more informaiton from resofs file
            if(file.exists(paste(rs@fileBase,"resofs",sep=".")))
            {
              ## read the resofs file
              try(
              rs@resofs<-read.table(paste(rs@fileBase,"resofs",sep="."))$V1,
              silent=F)
              if(length(rs@resofs)!=rs@nTrials)
                stop(paste("loadRecSession, problem with length of resofs",rs@session))
              ## trial times
              rs@trialStartRes<-c(0,rs@resofs[-length(rs@resofs)])
              rs@trialEndRes<-rs@resofs-1 # resofs is the number of samples in one file, index is -1 that
              if(length(rs@samplingRate)!=0){
                rs@trialDurationSec<-(rs@trialEndRes-rs@trialStartRes)/rs@samplingRate
                rs@sessionDurationSec<-sum(rs@trialDurationSec)
              }
            }else
            {
              # try to get the info from DatFiles object, but don't return error if not there
              df<-new("DatFiles")
              try(
                df<-datFilesSet(df,
                              fileNames=paste(rs@trialNames,"dat",sep="."),
                              path=rs@path,
                              nChannels=rs@nChannels),silent=TRUE)
              if(df@path!=""){
                rs@trialStartRes<-head(c(0,cumsum(df@samples)),rs@nTrials)
                rs@trialEndRes<-head(cumsum(df@samples),rs@nTrials)-1  # -1 because we want the index
                if(length(rs@samplingRate)!=0){
                  rs@trialDurationSec<-(rs@trialEndRes-rs@trialStartRes)/rs@samplingRate
                  rs@sessionDurationSec<-sum(rs@trialDurationSec)
                }
              }
            }

            if(file.exists(paste(rs@fileBase,"clu",sep="."))&
               file.exists(paste(rs@fileBase,"res",sep="."))) rs@clustered=T

            if(file.exists(paste(rs@fileBase,"resofs",sep=".")))rs@earlyProcessed=T
            return(rs)
          }
)


#' Create the session configuration files from a RecSession object
#'
#' This function create files with the following extension:
#' .par .desen .desel .px_per_cm .sampling_rate_dat .stimulation
#' 
#' @param rs A RecSession object
#' @param kiloSortConfig logical indicating whether to create configuration files for KiloSort
#'
#' @docType methods
#' @rdname saveRecSessionParameterFiles-methods
setGeneric(name="saveRecSessionParameterFiles",
           def=function(rs,kiloSortConfig=FALSE)
           {standardGeneric("saveRecSessionParameterFiles")}
)
#' @rdname saveRecSessionParameterFiles-methods
#' @aliases saveRecSessionParameterFiles,ANY,ANY-method
setMethod(f="saveRecSessionParameterFiles",
          signature="RecSession",
          definition=function(rs,kiloSortConfig=FALSE)
          {
            print(paste("saveRecSessionParameterFiles",rs@session))
            
            if(rs@session=="")
              stop("rs@session is empty")
            if(class(kiloSortConfig)!="logical")
              stop("kiloSortConfig is not a logical")
            if(!dir.exists(rs@path))
              stop(paste("saveRecSessionParameterFiles:",rs@path,"does not exist"))
            if(rs@fileBase=="")
              stop(paste("saveRecSessionParameterFiles: rs@fileBase is empty"))
            # write a par file
            print(paste("create",paste(rs@fileBase,"par",sep=".")))
            write(x=c(rs@nChannels, 16),file=paste(rs@fileBase,"par",sep="."),append = F,ncolumns = 2)
            write(x=c(1000000/rs@samplingRate, 800),file=paste(rs@fileBase,"par",sep="."),append = T,ncolumns = 2)
            write(x=c(rs@nElectrodes,0),file=paste(rs@fileBase,"par",sep="."),append = T,ncolumns = 2)
            for(t in 1:rs@nElectrodes){
              write(x =c(length(rs@channelsTetrode[t,which(!is.na(rs@channelsTetrode[t,]))]),
                         rs@channelsTetrode[t,which(!is.na(rs@channelsTetrode[t,]))]),
                    file=paste(rs@fileBase,"par",sep="."),append = T,ncolumns =length(rs@channelsTetrode[t,which(!is.na(rs@channelsTetrode[t,]))])+1 )
            }
            write(x=rs@nTrials,file=paste(rs@fileBase,"par",sep="."),append = T,ncolumns = 1)
            write(x=rs@trialNames,file=paste(rs@fileBase,"par",sep="."),append = T,ncolumns = 1)
            # write .desen file
            
            if(length(rs@environment)!=0){
              print(paste("create",paste(rs@fileBase,"desen",sep=".")))
              write(x=rs@environment,
                    file=paste(rs@fileBase,"desen",sep="."),
                    ncolumns = 1)
            }
            # write .desel file
            if(length(rs@electrodeLocation)!=0){
              print(paste("create",paste(rs@fileBase,"desel",sep=".")))
                write(x=rs@electrodeLocation,
                    file=paste(rs@fileBase,"desel",sep="."),
                    ncolumns = 1)
            }
            # write .sampling_rate_dat file
            print(paste("create",paste(rs@fileBase,"sampling_rate_dat",sep=".")))
            write(x=rs@samplingRate,
                  file=paste(rs@fileBase,"sampling_rate_dat",sep="."),
                  ncolumns = 1)
            # write .px_per_cm file
            print(paste("create",paste(rs@fileBase,"px_per_cm",sep=".")))
            write(x=rs@pxPerCm,
                  file=paste(rs@fileBase,"px_per_cm",sep="."),
                  ncolumns = 1)
            # write .res_samples_per_whd_sample
            print(paste("create",paste(rs@fileBase,"res_samples_per_whd_sample",sep=".")))
            write(x=rs@resSamplesPerWhdSample,
                  file=paste(rs@fileBase,"res_samples_per_whd_sample",sep="."),
                  ncolumns = 1)
            # write .stimulation file if needed
            if(length(rs@stimulation)!=0){
              print(paste("create",paste(rs@fileBase,"stimulation",sep=".")))
            write(x=rs@stimulation,
                  file=paste(rs@fileBase,"stimulation",sep="."),
                  ncolumns = 1)
            }
            # write .setup file if needed
            if(length(rs@setup)!=0){
              print(paste("create",paste(rs@fileBase,"setup",sep=".")))
              write(x=rs@setup,
                    file=paste(rs@fileBase,"setup",sep="."),
                    ncolumns = 1)
            }
            # write .environmentFamiliarity file if needed
            if(length(rs@environmentFamiliarity)!=0){
              print(paste("create",paste(rs@fileBase,"environmentFamiliarity",sep=".")))
              write(x=rs@environmentFamiliarity,
                    file=paste(rs@fileBase,"environmentFamiliarity",sep="."),
                    ncolumns = 1)
            }
            #############################################
            ## create configuration files for KiloSort ##
            #############################################
            if(kiloSortConfig==TRUE)
            { # function located in KiloSort.R
              writeKiloSortConfigurationFiles(rs)
            }
        }
)


#' Set a RecSession object with data passed as arguments
#'
#'
#' @param rs A RecSession object
#' @param session Session name
#' @param path Session directory path
#' @param samplingRate Sampling rate in Hz
#' @param nChannels Number of channels in the dat files.
#' @param nTrials Number of trials
#' @param nElectrodes Number of electrodes
#' @param trialNames Names of each trials
#' @param channelsTetrode Matrix containing the map of channel number for each tetrode, has 4 columns
#' @param environment Recording environment names during each trial
#' @param stimulation Codes for stimulation during each trial
#' @param setup Recording setups during each trial
#' @param environmentFamiliarity Familiarity of the recording environment during each trial
#' @param electrodeLocation Brain region for each electrode
#' @param pxPerCm Pixels per cm in position data
#' @param resSamplesPerWhdSample Number of res samples between whd sample
#' @return RecSession
#'
#' @docType methods
#' @rdname setRecSession-methods
setGeneric(name="setRecSession",
           def=function(rs,session,path,samplingRate,nChannels,nTrials,nElectrodes,
                        trialNames,channelsTetrode,environment,stimulation,setup,environmentFamiliarity,electrodeLocation,
                        pxPerCm,resSamplesPerWhdSample)
           {standardGeneric("setRecSession")}
)
#' @rdname setRecSession-methods
#' @aliases setRecSession,ANY,ANY-method
setMethod(f="setRecSession",
          signature="RecSession",
          definition=function(rs,session,path,samplingRate,nChannels,nTrials,
                              nElectrodes,trialNames,channelsTetrode,environment,stimulation,setup,environmentFamiliarity,electrodeLocation,
                              pxPerCm,resSamplesPerWhdSample)
          {
            if(session=="")
              stop("session is empty, you need to set a session name with session argument")
            rs@session<-session
            rs@path<-path
            rs@fileBase<-paste(rs@path,rs@session,sep="/")
            rs@animalName<-unlist(strsplit(rs@session,"-"))[1]
            if(samplingRate!="")
            {
              if(samplingRate<1|samplingRate>48000)
              {stop(paste("samplingRate is out of range:",samplingRate))}
              rs@samplingRate<-samplingRate  
            }
            if(nChannels!=""){
              if(nChannels<1){
                stop(paste("nChannels is out of range:",nChannels))
              }
              rs@nChannels<-nChannels
            }
            if(nTrials!=""){
              if(nTrials<1){
                stop(paste("nTrials is out of range:",nTrials))
              }
              rs@nTrials<-nTrials
            }
            if(nElectrodes!=""){
              if(nElectrodes<1){
                stop(paste("nElectrodes is out of range:",nElectrodes))
              }
              rs@nElectrodes<-nElectrodes
            }
            if(length(trialNames)!=0){
              if(length(trialNames)!=rs@nTrials)
                stop(paste("length of trialNames is not equal to rs@nTrials"))
              rs@trialNames<-trialNames
            }
            if(ncol(channelsTetrode)!=0){
              if(class(channelsTetrode)!="matrix")
                stop(paste("channelsTetrode should be a matrix but is a",class(channelsTetrode)))
              if(ncol(channelsTetrode)!=4)
                stop(paste("ncol(channelsTetrode) should be 4 but is",ncol(channelsTetrode)))
              if(nrow(channelsTetrode)!=rs@nElectrodes)
                stop(paste("nrow(channelsTetrode) should be rs@nElectrodes (",rs@nElectrodes,") but is", 
                           nrow(channelsTetrode)))
              rs@channelsTetrode<-channelsTetrode
            }              
            if(length(environment)!=0)
            {
              if(length(environment)!=rs@nTrials)
                stop(paste("length(env) should be rs@nTrials (",rs@nTrials,") but is",length(environment)))
              rs@environment<-environment
            }
            if(length(stimulation)!=0)
            {
              if(length(stimulation)!=rs@nTrials)
                stop(paste("length(stimulation) should be rs@nTrials (",rs@nTrials,") but is",length(stimulation)))
              rs@stimulation<-stimulation
            }
            if(length(setup)!=0)
            {
              if(length(setup)!=rs@nTrials)
                stop(paste("length(setup) should be rs@nTrials (",rs@nTrials,") but is",length(setup)))
              rs@setup<-setup
            }
            if(length(environmentFamiliarity)!=0)
            {
              if(length(environmentFamiliarity)!=rs@nTrials)
                stop(paste("length(environmentFamiliarity) should be rs@nTrials (",rs@nTrials,") but is",length(environmentFamiliarity)))
              rs@environmentFamiliarity<-environmentFamiliarity
            }  
            
            if(length(electrodeLocation)!=0)
            {
              if(length(electrodeLocation)!=rs@nElectrodes)
                stop(paste("length(electrodeLocation) should be rs@nElectrodes (",rs@nElectrodes,") but is",length(electrodeLocation)))
              rs@electrodeLocation<-electrodeLocation
            }
            
            if(pxPerCm!=""){
              if(pxPerCm<1){
                stop(paste("pxPerCm is out of range:",pxPerCm))
              }
              rs@pxPerCm<-pxPerCm
            }
          
            if(resSamplesPerWhdSample!=""){
              if(resSamplesPerWhdSample<1){
                stop(paste("resSamplesPerWhdSample is out of range:",resSamplesPerWhdSample))
              }
              rs@resSamplesPerWhdSample<-resSamplesPerWhdSample
            }
              
            rs@clustered=FALSE
            rs@earlyProcessed<-FALSE
            rs@fileBase<-paste(rs@path,rs@session,sep="/")
            rs@animalName<-unlist(strsplit(rs@session,"-"))[1]
            
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
            return(any(rs@environment==environment))
          })

#' Check if the session had a trial in a given stimulation type
#'
#' This will check whether the value of stimulation is in the stimulation vector.
#' 
#' @param rs A RecSession object
#' @param stimulation The name of a stimulation
#' @return TRUE or FALSE
#' 
#' @docType methods
#' @rdname containsStimulation-methods
setGeneric(name="containsStimulation",
           def=function(rs,stimulation="")
           {standardGeneric("containsStimulation")}
)
#' @rdname containsStimulation-methods
#' @aliases containsStimulation,ANY,ANY-method
setMethod(f="containsStimulation",
          signature="RecSession",
          definition=function(rs,stimulation="")
          {
            return(any(rs@stimulation==stimulation))
          })



#' Check if the session directory contains a file ending with the value of the argument extension
#'
#' By default test whether paste(rs@fileBase,extension,sep=".") exists
#'
#' @param rs A RecSession object
#' @param extension The extension of the file you are looking for.
#' @return TRUE or FALSE
#'
#' @docType methods
#' @rdname fileExists-methods
setGeneric(name="fileExists",
           def=function(rs,extension="")
           {standardGeneric("fileExists")}
)
#' @rdname fileExists-methods
#' @aliases fileExists,ANY,ANY-method
setMethod(f="fileExists",
          signature="RecSession",
          definition=function(rs,extension="")
          {
            return(file.exists(paste(rs@fileBase,extension,sep=".")))
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
            if(!environment%in%rs@environment){
              print("environment not used in the session")
              return()
            }
          return(matrix(data=c(rs@trialStartRes[which(rs@environment==environment)],rs@trialEndRes[which(rs@environment==environment)]),ncol=2,
                 dimnames=list(rep(environment,length(which(rs@environment==environment))),c("start","end"))))
          })

#' Get the time intervals in sample values for trials for a given stimulation
#'
#' @param rs A RecSession object
#' @param stimulation The name of a stimulation
#' @return matrix with 2 columns containing the start and end of each trial in the stimulation
#' 
#' @docType methods
#' @rdname getIntervalsStimulation-methods
setGeneric(name="getIntervalsStimulation",
           def=function(rs,stimulation="lt")
           {standardGeneric("getIntervalsStimulation")}
)
#' @rdname getIntervalsStimulation-methods
#' @aliases getIntervalsStimulation,ANY,ANY-method
setMethod(f="getIntervalsStimulation",
          signature="RecSession",
          definition=function(rs,stimulation="lt")
          {
            if(length(rs@trialStartRes)==0){
              print("trialStartRes is not set")
              return()
            }
            if(!stimulation%in%rs@stimulation){
              print("stimulation not used in the session")
              return()
            }
          return(matrix(data=c(rs@trialStartRes[which(rs@stimulation==stimulation)],rs@trialEndRes[which(rs@stimulation==stimulation)]),ncol=2,
                 dimnames=list(rep(stimulation,length(which(rs@stimulation==stimulation))),c("start","end"))))
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
            if(file.exists(paste(df@path,df@fileNames[1],sep="/")))
            { # make it possible to use this function with recSession that have no .dat files
              df<-datFilesSet(df,
                              fileNames=paste(rs@trialNames,"dat",sep="."),
                              path=rs@path,
                              nChannels=rs@nChannels)
            }

            cg<-new("CellGroup",session=rs@session,path=rs@path,nTetrodes=rs@nElectrodes)
            cg<-loadCellGroup(cg)
            sp<-new("SpatialProperties2d",session=rs@session)
            sp1<-new("SpatialProperties1d",session=rs@session)
            hd<-new("HeadDirection",session=rs@session)
            sw<-new("SpikeWaveform",session=rs@session)

            if(st@nCells!=cg@nCells){
              print(paste("st@nCells is not equal to cg@nCells for",rs@session))
              print("There is probably a cluster with no spike that was not removed at clustering time")
              print("cg object")
              print(cg)
              print("st object")
              print(st)
              stop()
            }

            return(list(st=st,pt=pt,df=df,cg=cg,sp=sp,sp1=sp1,hd=hd,sw=sw))
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
            if(length(object@environment)!=0){
              print(paste("environment:"))
              print(paste(object@environment))
            }
            if(length(object@stimulation)!=0){
              print(paste("stimulation:"))
              print(paste(object@stimulation))
            }
            if(length(object@setup)!=0){
              print(paste("setup:"))
              print(paste(object@setup))
            }
            if(length(object@environmentFamiliarity)!=0){
              print(paste("environmentFamiliarity:"))
              print(paste(object@environmentFamiliarity))
            }
            
            print(paste("nElectrodes:",object@nElectrodes))
            if(length(object@electrodeLocation)!=0){
              print(paste("electrodeLocation:"))
              print(paste(object@electrodeLocation))
            }
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
            print(paste("pxPerCm:",object@pxPerCm))
            print(paste("resSamplesPerWhdSample:",object@resSamplesPerWhdSample))
          })



#' Get the tetrode number of a channel
#'
#' @param rs A RecSession object
#' @param channelNo Channel no of interest, can be a vector of length > 1
#' @return integer representing the tetrode no for each channel
#' 
#' @docType methods
#' @rdname getTetrodeNumberOfChannel-methods
setGeneric(name="getTetrodeNumberOfChannel",
           def=function(rs,channelNo)
           {standardGeneric("getTetrodeNumberOfChannel")}
)
#' @rdname getTetrodeNumberOfChannel-methods
#' @aliases getTetrodeNumberOfChannel,ANY,ANY-method
setMethod(f="getTetrodeNumberOfChannel",
          signature="RecSession",
          definition=function(rs,channelNo)
          {
            if(sum(dim(rs@channelsTetrode))==0)
              stop("rs@channelsTetrodes has dimension of 0 0")
            if(any(channelNo<0) | any(channelNo>=rs@nChannels))
              stop("channelNo is out of range")
            return(as.integer(apply(apply(rs@channelsTetrode,1,function(x,y) y %in% x ,channelNo),1,function(x) head(which(x),n = 1))))
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
#' Assumes the session name is in the format name-date-rest and cluId is name-data-rest_cluNo
#'
#' @param cluId Character vector with the session name
#' @return Character vector with session name
sessionNameFromCluId<-function(cluId=NULL){
  return(unlist(lapply(strsplit(as.character(cluId),split="_"),function(x){return(x[[1]])})))
}

#' Get cluNo from cluId
#'
#' Assumes the cluId is in the format name-date-rest_cluNo
#'
#' @param cluId Character vector with the cluIds
#' @return Numeric vectors with cluNo
cluNoFromCluId<-function(cluId=NULL){
  return(as.numeric(unlist(lapply(strsplit(as.character(cluId),split="_"),function(x){return(x[[2]])}))))
}

