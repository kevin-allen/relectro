############################################
#### definition of RecSession Class      ###
############################################
RecSession <- setClass(
  "RecSession", ## name of the class
  slots=c(session="character",
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
          channelsTetrode="matrix"),  # cell list to limit the analysis to these cells
  prototype = list(session=""))


### loadRecSession ###
setGeneric(name="loadRecSession",
           def=function(rs)
           {standardGeneric("loadRecSession")}
)
setMethod(f="loadRecSession",
          signature="RecSession",
          definition=function(rs)
          {
            #setwd("~/repo/r_packages/data")
            #session="jp4298-15022016-0106"
            if(rs@session=="")
              stop("rs@session is empty")
            if(!file.exists(paste(rs@session,"par",sep=".")))
              stop("need",paste(rs@session,"par",sep="."))
            if(!file.exists(paste(rs@session,"desen",sep=".")))
              stop("need",paste(rs@session,"desen",sep="."))
            if(!file.exists(paste(rs@session,"desel",sep=".")))
              stop("need",paste(rs@session,"desel",sep="."))
            if(!file.exists(paste(rs@session,"resofs",sep=".")))
              stop("need",paste(rs@session,"resofs",sep="."))
            if(!file.exists(paste(rs@session,"sampling_rate_dat",sep=".")))
              stop("need",paste(rs@session,"sampling_rate_dat",sep="."))
            
            ## get sampling rate
            rs@samplingRate<-read.table(paste(rs@session,"sampling_rate_dat",sep="."))$V1
            ## read the par file line per line## shitty format
            conn <- file(paste(rs@session,"par",sep="."),open="r")
            par<-readLines(conn)
            close(conn)
            rs@nChannels<-as.numeric(unlist(strsplit(par[1], split=" "))[1])
            rs@nElectrodes  <-as.numeric(unlist(strsplit(par[3], split=" "))[1])
            rs@nTrials<-as.numeric(par[rs@nElectrodes+4])
            rs@trialNames<-par[(rs@nElectrodes+5):(rs@nElectrodes+5+rs@nTrials-1)]
            ## read the desen file
            rs@env<-as.character(read.table(paste(rs@session,"desen",sep="."))$V1)
            ## read the desel file
            rs@electrodeLocation<-as.character(read.table(paste(rs@session,"desel",sep="."))$V1)
            ## read the resofs file  
            rs@resofs<-read.table(paste(rs@session,"resofs",sep="."))$V1
            
            ## check that things add up
            if(length(rs@env)!=length(rs@trialNames))
              stop("Problem with length of par and desen files")
            if(rs@nElectrodes!=length(rs@electrodeLocation))
              stop("Problem with length of par and desel files")
            if(rs@nTrials!=length(rs@trialNames))
              stop("Problem with number of trials in par file")
            if(length(rs@resofs)!=rs@nTrials)
              stop("Problem with length of resofs")
            if(rs@samplingRate<1 | rs@samplingRate > 100000)
              stop(paste("samplingRate is out of range:",rs@samplingRate))
            
            ## trial times
            rs@trialStartRes<-c(0,rs@resofs[-length(rs@resofs)])
            rs@trialEndRes<-rs@resofs
            rs@trialDurationSec<-(rs@trialEndRes-rs@trialStartRes)/rs@samplingRate
            rs@sessionDurationSec<-sum(rs@trialDurationSec)
            
            ## map of channel and tetrodes
            chan<-strsplit(par[4:(4+rs@nElectrodes-1)], split=" ")
            max.channelsTetrode<-max(unlist(lapply(chan,length))-1)
            rs@channelsTetrode<-matrix(nrow=rs@nElectrodes,ncol=max.channelsTetrode)
            
            for(i in 1:rs@nElectrodes)
              rs@channelsTetrode[i,]<-as.numeric(chan[[i]][-1])
            
            return(rs)
          }
)

### getTrialIntervalsWithEnv ###
setGeneric(name="getIntervalsEnvironment",
           def=function(rs,env="lt")
           {standardGeneric("getIntervalsEnvironment")}
)
setMethod(f="getIntervalsEnvironment",
          signature="RecSession",
          definition=function(rs,env="lt")
          {
            if(!env%in%rs@env){
              print("environment not used in the session")
              return()
            }
          return(matrix(data=c(rs@trialStartRes[which(rs@env==env)],rs@trialEndRes[which(rs@env==env)]),ncol=2,
                 dimnames=list(rep(env,length(which(rs@env==env))),c("start","end"))))
          })
            

### show ###
setMethod("show", "RecSession",
          function(object){
            print(paste("session:",object@session))
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
          
          })
