###############################################
#### definition of SpikeTrain class         ###
###############################################
SpikeTrain <- setClass(
  "SpikeTrain", ## name of the class
  slots=c(session="character",samplingRate="numeric",res="numeric",time="numeric",clu="numeric",
          nCells="numeric",nSpikes="numeric",nSpikesPerCell="numeric",
          startInterval="numeric",endInterval="numeric", # to limit analysis to these intervals
          startResIndexc="numeric",endResIndexc="numeric",
          events="numeric",
          cellList="numeric",cellPairList="data.frame"),  # cell list to limit the analysis to these cells
  prototype = list(session=""))

### spikeTimeAutocorrelation ###
# function using .Call() to call c wrapper
# It calls the the c function with the cwrap suffix
setGeneric(name="spikeTimeAutocorrelation",
           def=function(st,...)
             {standardGeneric("spikeTimeAutocorrelation")})
setMethod(f="spikeTimeAutocorrelation",
          signature = "SpikeTrain",
          definition=function(st,bin.size.ms=1,window.size.ms=200,probability=FALSE,...)
            {
            n.bins=(window.size.ms*2)/bin.size.ms
            window.size=2*window.size.ms*st@samplingRate/1000 # window size in res value from - to + extrems
            # call cwrapper function
            
            results<- .Call("autocorrelation_cwrap",
                        st@cellList,
                        length(st@cellList),
                        st@clu,
                        st@res,
                        st@nSpikes,
                        n.bins,
                        window.size,
                        st@startResIndexc,
                        st@endResIndexc,
                        length(st@startResIndexc),
                        probability)
            
            # return a data fram with clu time count
            if(probability==F)
              data.frame(clu=rep(st@cellList,each=n.bins),
                         time=seq(-window.size.ms+bin.size.ms,window.size.ms,bin.size.ms)-(bin.size.ms/2),
                         count=results)
            else{
              data.frame(clu=rep(st@cellList,each=n.bins),
                         time=seq(-window.size.ms+bin.size.ms,window.size.ms,bin.size.ms)-(bin.size.ms/2),
                         prob=results)
            }
          }
          )




### spikeTimeCrosscorrelationEvents ###
# function using .Call() to call c wrapper
# It calls the the c function with the cwrap suffix
setGeneric(name="spikeTimeCrosscorrelationEvents",
           def=function(st,bin.size.ms=1,window.size.ms=200,probability=FALSE)
           {standardGeneric("spikeTimeCrosscorrelationEvents")})

setMethod(f="spikeTimeCrosscorrelationEvents",
          signature = "SpikeTrain",
          definition=function(st,bin.size.ms=1,window.size.ms=200,probability=FALSE)
          {
            
            if(length(st@events)==0)
              stop("events is empty")
            
            
            n.bins=(window.size.ms*2)/bin.size.ms
            window.size=2*window.size.ms*st@samplingRate/1000 # window size in res value from - to + extrems
            # call cwrapper function
            
            dyn.load("~/repo/r_packages/relectro/src/relectro.so")
            results<- .Call("crosscorrelationEvents_cwrap",
                            st@cellList,
                            length(st@cellList),
                            st@clu,
                            st@res,
                            st@nSpikes,
                            n.bins,
                            window.size,
                            st@startResIndexc,
                            st@endResIndexc,
                            length(st@startResIndexc),
                            probability,
                            st@events,
                            length(st@events))
            
            
            # return a data fram with clu time count
            if(probability==F){
              data.frame(clu=rep(st@cellList,each=n.bins),
                         time=seq(-window.size.ms+bin.size.ms,window.size.ms,bin.size.ms)-(bin.size.ms/2),
                         count=results)
            }
            else{
              data.frame(clu=rep(st@cellList,each=n.bins),
                         time=seq(-window.size.ms+bin.size.ms,window.size.ms,bin.size.ms)-(bin.size.ms/2),
                         prob=results)
            }
          
        })

### spikeTimeCrosscorrelation ###
# function using .Call() to call c wrapper
# It calls the the c function with the cwrap suffix
setGeneric(name="spikeTimeCrosscorrelation",
           def=function(st,bin.size.ms=1,window.size.ms=200,probability=FALSE)
           {standardGeneric("spikeTimeCrosscorrelation")})

setMethod(f="spikeTimeCrosscorrelation",
          signature = "SpikeTrain",
          definition=function(st,bin.size.ms=1,window.size.ms=200,probability=FALSE)
          {
            
            if(length(st@cellPairList[,1])==0)
              stop("cellPairList is empty")
            
           # bin.size.ms=1
          #  window.size.ms=200
          #  probability=FALSE
            
            n.bins=(window.size.ms*2)/bin.size.ms
            window.size=2*window.size.ms*st@samplingRate/1000 # window size in res value from - to + extrems
            # call cwrapper function
            window.size
            
            results<- .Call("crosscorrelation_cwrap",
                            st@cellPairList[,1],
                            st@cellPairList[,2],
                            length(st@cellPairList[,1]),
                            st@clu,
                            st@res,
                            st@nSpikes,
                            n.bins,
                            window.size,
                            st@startResIndexc,
                            st@endResIndexc,
                            length(st@startResIndexc),
                            probability)
            
            
            # return a data fram with clu time count
            if(probability==F)
              data.frame(clu1=rep(st@cellPairList[,1],each=n.bins),
                         clu2=rep(st@cellPairList[,2],each=n.bins),
                         time=seq(-window.size.ms+bin.size.ms,window.size.ms,bin.size.ms)-(bin.size.ms/2),
                         count=results)
            else
              data.frame(clu1=rep(st@cellPairList[,1],each=n.bins),
                         clu2=rep(st@cellPairList[,2],each=n.bins),
                         time=seq(-window.size.ms+bin.size.ms,window.size.ms,bin.size.ms)-(bin.size.ms/2),
                         prob=results)
            
          }
          )


### loadSpikeTrain ###
setGeneric(name="loadSpikeTrain",
           def=function(st)
           {standardGeneric("loadSpikeTrain")}
)
setMethod(f="loadSpikeTrain",
          signature="SpikeTrain",
          definition=function(st)
          {
            if(!file.exists(paste(st@session,"res",sep=".")))
              stop("need",paste(st@session,"res",sep="."))
            if(!file.exists(paste(st@session,"clu",sep=".")))
              stop("need",paste(st@session,"clu",sep="."))
            if(!file.exists(paste(st@session,"sampling_rate_dat",sep=".")))
              stop("need",paste(st@session,"sampling_rate_dat",sep="."))
            
            st@res<-read.table(paste(st@session,"res",sep="."))$V1
            st@clu<-read.table(paste(st@session,"clu",sep="."))$V1
            st@samplingRate<-read.table(paste(st@session,"sampling_rate_dat",sep="."))$V1
            if(st@samplingRate<1 | st@samplingRate > 100000)
              stop(paste("samplingRate is out of range:",st@samplingRate))
            st@clu<-st@clu[-1] ## remove first number
            if(length(st@res)!=length(st@clu))
              stop(paste("length of res and clu files not equals:",length(st@res),length(st@clu)))
            # remove noise cluster 1 or 0
            df<-data.frame(res=st@res,clu=st@clu)
            df<-df[which(df$clu>1),]
            st@res<-df$res
            st@clu<-df$clu
            # get time in sec
            st@time<-st@res/st@samplingRate
            st@nSpikes<-length(st@res)
            st@nCells<-length(unique(st@clu))
            st@nSpikesPerCell<-as.numeric(table(st@clu))
            # by default analysis on all cells and all recording period
            st@cellList<-sort(unique(st@clu))
            # by default get all the possible pairs
            if(length(st@cellList>1))
              st@cellPairList<-make.pairs(st@cellList)
            st@startInterval<-0
            st@endInterval<-max(st@res)
            st@startResIndexc<-0
            st@endResIndexc<-length(st@res)-1
            return(st)
          
          }
)

### setSpikeTrain ###
setGeneric(name="setSpikeTrain",
           def=function(st,...)
           {standardGeneric("setSpikeTrain")}
)


setMethod(f="setSpikeTrain",
          signature="SpikeTrain",
          definition=function(st,res,clu,sampling.rate)
          {
            
            #res=df$res
            #clu=df$clu
            #sampling.rate=20000
            
            st@res<-res
            st@clu<-clu
            st@samplingRate<-sampling.rate
          
            if(st@samplingRate<1 | st@samplingRate > 100000)
              stop(paste("samplingRate is out of range:",st@samplingRate))
          
            if(length(st@res)!=length(st@clu))
              stop(paste("length of res and clu files not equals:",length(st@res),length(st@clu)))
          
            # get time in sec
            st@time<-st@res/st@samplingRate
            st@nSpikes<-length(st@res)
            st@nCells<-length(unique(st@clu))
            st@nSpikesPerCell<-as.numeric(table(st@clu))
            # by default analysis on all cells and all recording period
            st@cellList<-sort(unique(st@clu))
            
            if(length(st@cellList)>1)
              st@cellPairList<-make.pairs(st@cellList)
            st@startInterval<-0
            st@endInterval<-max(st@res)
            st@startResIndexc<-0
            st@endResIndexc<-length(st@res)-1
            return(st)
          }
)

### meanFiringRate ###
setGeneric(name="meanFiringRate",
           def=function(st)
           {standardGeneric("meanFiringRate")})

setMethod(f="meanFiringRate",
          signature = "SpikeTrain",
          definition=function(st)
          {
            # call cwrapper function
            results<- .Call("meanFiringRate_cwrap",
                            st@cellList,
                            length(st@cellList),
                            st@clu,
                            st@res,
                            st@nSpikes,
                            st@startInterval,
                            st@endInterval,
                            st@startResIndexc,
                            st@endResIndexc,
                            length(st@startResIndexc),
                            st@samplingRate)
            return(results)
          }
)


### setIntervals ###
setGeneric(name="setIntervals",
           def=function(st,s,e)
           {standardGeneric("setIntervals")})

setMethod(f="setIntervals",
          signature = "SpikeTrain",
          definition=function(st,s,e)
          {
            
            ## if s is a matrix, then e is ignored
            if(class(s)=="matrix"){
              start.intervals<-as.numeric(s[,1])
              end.intervals<-as.numeric(s[,2])
              
            }else{
              start.intervals<-s
              end.intervals<-e
            }
            
            if(length(start.intervals)!=length(end.intervals))
              stop("problem with length of start.intervals and end.intervals")
            if(any(start.intervals>end.intervals))
              stop("start.intervals>end.intervals")
            if(any(diff(start.intervals)<0))
              stop("problem with chonology within start.intervals")
            if(any(diff(end.intervals)<0))
              stop("problem with chonology within end.intervals")
            if(any(start.intervals[-1]-end.intervals[-length(end.intervals)]<0))
              stop("problem with chronology between intervals, from end to next start")
            
            st@startInterval<-start.intervals
            st@endInterval<-end.intervals
            
            st@startInterval
            st@endInterval
            # call cwrapper function
            results<- .Call("resIndexForIntervals_cwrap",
                            length(st@startInterval),
                            st@startInterval,
                            st@endInterval,
                            st@nSpikes,
                            st@res)
            
            st@startResIndexc<-results[1,]
            st@endResIndexc<-results[2,]
            st@startInterval<-st@startInterval[1:length(st@startResIndexc)]
            st@endInterval<-st@endInterval[1:length(st@endResIndexc)]
            return(st)
          }
)




### setEvents ###
setGeneric(name="setEvents",
           def=function(st,events)
           {standardGeneric("setEvents")})

setMethod(f="setEvents",
          signature = "SpikeTrain",
          definition=function(st,events=NULL)
          {
            if(is.null(events))
              stop("events is empty")
            if(any(events<0))
              stop("negative values as events")
            if(any(diff(events)<0))
              stop("problem with the chronology of the events")
            
            st@events<-events
            return(st)
          }
)

### show ###
setMethod("show", "SpikeTrain",
          function(object){
            print(paste("session:",object@session))
            print(paste("samplingRate:",object@samplingRate))
            print(paste("nCells:",object@nCells))
            print(paste("nSpikes:",object@nSpikes))
            print(paste("nSpikesPerCell:"))
            print(object@nSpikesPerCell)
            print(paste("cellList:"))
            print(object@cellList)
            print(paste("n cellPairList:",length(object@cellPairList[,1])))
            print(paste("nIntervals:",length(object@startInterval))) 
            print(paste("Interval time:", sum(object@endInterval-object@startInterval)/object@samplingRate,"sec"))
            print(paste(object@startInterval,object@endInterval))
            print(paste("nIntervalsc:",length(object@startResIndexc)))
            print(paste(object@startResIndexc,object@endResIndexc))
            print(paste("nEvents:",length(object@events)))
            
                  })