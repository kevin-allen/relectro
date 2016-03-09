############################################
#### definition of RecSession Class      ###
############################################
RecSession <- setClass(
  "RecSession", ## name of the class
  slots=c(session="character",
          samplingRate="numeric",
          resofs="numeric",
          env="character",
          electrode.location="character",
          trial.start.res="numeric",
          trial.end.res="numeric",
          trial.names="character",
          trial.duration.sec="numeric",
          session.duration.sec="numeric",
          n.electrodes="numeric",
          n.channels="numeric",
          n.trials="numeric",
          channels.tetrode="matrix"),  # cell list to limit the analysis to these cells
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
            rs@n.channels<-as.numeric(unlist(strsplit(par[1], split=" "))[1])
            rs@n.electrodes  <-as.numeric(unlist(strsplit(par[3], split=" "))[1])
            rs@n.trials<-as.numeric(par[rs@n.electrodes+4])
            rs@trial.names<-par[(rs@n.electrodes+5):(rs@n.electrodes+5+rs@n.trials-1)]
            ## read the desen file
            rs@env<-as.character(read.table(paste(rs@session,"desen",sep="."))$V1)
            ## read the desel file
            rs@electrode.location<-as.character(read.table(paste(rs@session,"desel",sep="."))$V1)
            ## read the resofs file  
            rs@resofs<-read.table(paste(rs@session,"resofs",sep="."))$V1
            
            ## check that things add up
            if(length(rs@env)!=length(rs@trial.names))
              stop("Problem with length of par and desen files")
            if(rs@n.electrodes!=length(rs@electrode.location))
              stop("Problem with length of par and desel files")
            if(rs@n.trials!=length(rs@trial.names))
              stop("Problem with number of trials in par file")
            if(length(rs@resofs)!=rs@n.trials)
              stop("Problem with length of resofs")
            if(rs@samplingRate<1 | rs@samplingRate > 100000)
              stop(paste("samplingRate is out of range:",rs@samplingRate))
            
            
            ## trial times
            rs@trial.start.res<-c(0,rs@resofs[-length(rs@resofs)])
            rs@trial.end.res<-rs@resofs
            rs@trial.duration.sec<-(rs@trial.end.res-rs@trial.start.res)/rs@samplingRate
            rs@session.duration.sec<-sum(rs@trial.duration.sec)
            
            ## map of channel and tetrodes
            chan<-strsplit(par[4:(4+rs@n.electrodes-1)], split=" ")
            max.channels.tetrode<-max(unlist(lapply(chan,length))-1)
            rs@channels.tetrode<-matrix(nrow=rs@n.electrodes,ncol=max.channels.tetrode)
            
            for(i in 1:rs@n.electrodes)
              rs@channels.tetrode[i,]<-as.numeric(chan[[i]][-1])
            
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
          return(matrix(data=c(rs@trial.start.res[which(rs@env==env)],rs@trial.end.res[which(rs@env==env)]),ncol=2,
                 dimnames=list(rep(env,length(which(rs@env==env))),c("start","end"))))
          })
            









### show ###
setMethod("show", "RecSession",
          function(object){
            print(paste("session:",object@session))
            print(paste("samplingRate:",object@samplingRate,"Hz"))
            print(paste("n.channels:",object@n.channels))
            print(paste("n.trials:",object@n.trials))
            print(paste("env:"))
            print(paste(object@env))
            print(paste("n.electrodes:",object@n.electrodes))
            print(paste("electrode.location:"))
            print(paste(object@electrode.location))
            print(paste("trial.names:"))
            print(paste(object@trial.names))
            print(paste("trial.duration.sec:"))
            print(paste(object@trial.duration.sec,"sec"))
            print(paste("session duration:",object@session.duration.sec,"sec"))
            print(paste("trial.start.res and trial.end.res:"))
            print(paste(object@trial.start.res,object@trial.end.res))
            print(paste("Map of channels per tetrode (channels.tetrode)"))
            print(object@channels.tetrode)
          
          })
