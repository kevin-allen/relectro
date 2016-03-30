############################################
#### definition of RecSession Class      ###
############################################
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


### loadRecSession ###
setGeneric(name="loadRecSession",
           def=function(rs)
           {standardGeneric("loadRecSession")}
)
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
            
            for(i in 1:rs@nElectrodes)
              rs@channelsTetrode[i,]<-as.numeric(chan[[i]][-1])
            
            if(file.exists(paste(rs@fileBase,"desen",sep="."))){
              rs@env<-as.character(read.table(paste(rs@fileBase,"desen",sep="."))$V1)
              if(length(rs@env)!=length(rs@trialNames))
                stop("Problem with length of par and desen files")
            }
            if(file.exists(paste(rs@fileBase,"desel",sep="."))){
                rs@electrodeLocation<-as.character(read.table(paste(rs@fileBase,"desel",sep="."))$V1)
              if(rs@nElectrodes!=length(rs@electrodeLocation))
                stop("Problem with length of par and desel files")
            }
            
            if(rs@nTrials!=length(rs@trialNames))
              stop("Problem with number of trials in par file")

            if(file.exists(paste(rs@fileBase,"sampling_rate_dat",sep="."))){
              rs@samplingRate<-read.table(paste(rs@fileBase,"sampling_rate_dat",sep="."))$V1
              if(rs@samplingRate<1 | rs@samplingRate > 100000)
                stop(paste("samplingRate is out of range:",rs@samplingRate))
            }
            
            ## if early process was run on this one, get more informaiton
            if(file.exists(paste(rs@fileBase,"resofs",sep=".")))
            {
              ## read the resofs file  
              rs@resofs<-read.table(paste(rs@fileBase,"resofs",sep="."))$V1
              if(length(rs@resofs)!=rs@nTrials)
                stop("Problem with length of resofs")
              ## trial times
              rs@trialStartRes<-c(0,rs@resofs[-length(rs@resofs)])
              rs@trialEndRes<-rs@resofs
              rs@trialDurationSec<-(rs@trialEndRes-rs@trialStartRes)/rs@samplingRate
              rs@sessionDurationSec<-sum(rs@trialDurationSec)
            }
            if(file.exists(paste(rs@fileBase,"clu",sep="."))) rs@clustered=T
            
            if(file.exists(paste(rs@fileBase,"resofs",sep=".")))rs@earlyProcessed=T
            return(rs)
          }
)


### getIsClustered ###
setGeneric(name="getIsClustered",
           def=function(rs)
           {standardGeneric("getIsClustered")}
)
setMethod(f="getIsClustered",
          signature="RecSession",
          definition=function(rs)
          {
            return(rs@clustered)
          })


### getIsEarlyProcessed ###
setGeneric(name="getIsEarlyProcessed",
           def=function(rs)
           {standardGeneric("getIsEarlyProcessed")}
)
setMethod(f="getIsEarlyProcessed",
          signature="RecSession",
          definition=function(rs)
          {
            return(rs@earlyProcessed)
          })


### containsElectrodeLocation ###
setGeneric(name="containsElectrodeLocation",
           def=function(rs,location="")
           {standardGeneric("containsElectrodeLocation")}
)
setMethod(f="containsElectrodeLocation",
          signature="RecSession",
          definition=function(rs,location="")
          {
            return(any(rs@electrodeLocation==location))
          })

### containsEnvironment ###
setGeneric(name="containsEnvironment",
           def=function(rs,environment="")
           {standardGeneric("containsEnvironment")}
)
setMethod(f="containsEnvironment",
          signature="RecSession",
          definition=function(rs,environment="")
          {
            return(any(rs@env==environment))
          })




### getTrialIntervalsWithEnv ###
setGeneric(name="getIntervalsEnvironment",
           def=function(rs,env="lt")
           {standardGeneric("getIntervalsEnvironment")}
)
setMethod(f="getIntervalsEnvironment",
          signature="RecSession",
          definition=function(rs,env="lt")
          {
            
            
            if(length(rs@trialStartRes)==0){
              print("trialStartRes is not set")
              return()
            }
            if(!env%in%rs@env){
              print("environment not used in the session")
              return()
            }
          return(matrix(data=c(rs@trialStartRes[which(rs@env==env)],rs@trialEndRes[which(rs@env==env)]),ncol=2,
                 dimnames=list(rep(env,length(which(rs@env==env))),c("start","end"))))
          })

### getRecSessionObjects ###
setGeneric(name="getRecSessionObjects",
           def=function(rs)
           {standardGeneric("getRecSessionObjects")}
)
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
