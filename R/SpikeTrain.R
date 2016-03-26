###############################################
#### definition of SpikeTrain class         ###
###############################################
SpikeTrain <- setClass(
  "SpikeTrain", ## name of the class
  slots=c(session="character",
          path="character",
          samplingRate="numeric",
          res="numeric",
          time="numeric",
          clu="numeric",
          nCells="numeric",
          nSpikes="numeric",
          nSpikesPerCell="numeric",
          startInterval="numeric",endInterval="numeric", # to limit analysis to these intervals
          startResIndexc="numeric",endResIndexc="numeric",
          events="numeric",
          cellList="numeric",cellPairList="data.frame", # cell list to limit the analysis to these cells
          # ifr
          ifrKernelSdMs="numeric",
          ifrWindowSizeMs="numeric",
          ifrSpikeBinMs="numeric",
          ifr="matrix",
          ifrTime="numeric",
          #
          meanFiringRate="numeric"
          ),  
  prototype = list(session="",ifrKernelSdMs=50,ifrWindowSizeMs=50,ifrSpikeBinMs=1))


### ifr ###
setGeneric(name="ifr",
           def=function(st)
           {standardGeneric("ifr")})

setMethod(f="ifr",
          signature = "SpikeTrain",
          definition=function(st)
          {
            if(st@nSpikes==0)
              stop("ifr(): there are no spike in the SpikeTrain object")

            results<-.Call("ifr_from_spike_density",
                  as.integer(st@res),as.integer(st@clu), as.integer(st@nSpikes), 
                  st@ifrWindowSizeMs, st@ifrKernelSdMs, st@ifrSpikeBinMs,
                  as.integer(st@cellList), length(st@cellList),
                  as.integer(st@startInterval),as.integer(st@endInterval),length(st@startInterval),
                  st@samplingRate)
            
            st@ifrTime<-results[1,]
            st@ifr<-matrix(results[-1,],nrow=st@nCells)
            
            return(st)
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
            
            if(st@session=="")
              stop("st@session is not set")
            if(st@path=="") ## path is given or is getwd()
              st@path=getwd()
            pathSession=paste(st@path,st@session,sep="/")
            
            
            if(!file.exists(paste(pathSession,"res",sep=".")))
              stop("need ",paste(pathSession,"res",sep="."))
            if(!file.exists(paste(pathSession,"clu",sep=".")))
              stop("need ",paste(pathSession,"clu",sep="."))
            if(!file.exists(paste(pathSession,"sampling_rate_dat",sep=".")))
              stop("need ",paste(pathSession,"sampling_rate_dat",sep="."))
            
            st@res<-as.numeric(.Call("read_one_column_int_file_cwrap", paste(pathSession,"res",sep=".")))
            st@clu<-as.numeric(.Call("read_one_column_int_file_cwrap", paste(pathSession,"clu",sep=".")))
            
            st@samplingRate<-as.numeric(readLines(paste(pathSession,"sampling_rate_dat",sep=".")))
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
              st@cellPairList<-makePairs(st@cellList)
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
          definition=function(st,res,clu,samplingRate)
          {
            
            st@res<-res
            st@clu<-clu
            st@samplingRate<-samplingRate
            
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
              st@cellPairList<-makePairs(st@cellList)
            st@startInterval<-0
            st@endInterval<-max(st@res)
            st@startResIndexc<-0
            st@endResIndexc<-length(st@res)-1
            return(st)
          }
)





### spikeTimeAutocorrelation ###
# function using .Call() to call c wrapper
# It calls the the c function with the cwrap suffix
setGeneric(name="spikeTimeAutocorrelation",
           def=function(st,...)
             {standardGeneric("spikeTimeAutocorrelation")})
setMethod(f="spikeTimeAutocorrelation",
          signature = "SpikeTrain",
          definition=function(st,binSizeMs=1,windowSizeMs=200,probability=FALSE,...)
            {
            nBins=(windowSizeMs*2)/binSizeMs
            windowSize=2*windowSizeMs*st@samplingRate/1000 # window size in res value from - to + extrems
            # call cwrapper function
            
            results<- .Call("autocorrelation_cwrap",
                        st@cellList,
                        length(st@cellList),
                        st@clu,
                        st@res,
                        st@nSpikes,
                        nBins,
                        windowSize,
                        st@startResIndexc,
                        st@endResIndexc,
                        length(st@startResIndexc),
                        probability)
            
            # return a data fram with clu time count
            if(probability==F)
              data.frame(clu=rep(st@cellList,each=nBins),
                         time=seq(-windowSizeMs+binSizeMs,windowSizeMs,binSizeMs)-(binSizeMs/2),
                         count=results)
            else{
              data.frame(clu=rep(st@cellList,each=nBins),
                         time=seq(-windowSizeMs+binSizeMs,windowSizeMs,binSizeMs)-(binSizeMs/2),
                         prob=results)
            }
          }
          )




### spikeTimeCrosscorrelationEvents ###
# function using .Call() to call c wrapper
# It calls the the c function with the cwrap suffix
setGeneric(name="spikeTimeCrosscorrelationEvents",
           def=function(st,binSizeMs=1,windowSizeMs=200,probability=FALSE)
           {standardGeneric("spikeTimeCrosscorrelationEvents")})

setMethod(f="spikeTimeCrosscorrelationEvents",
          signature = "SpikeTrain",
          definition=function(st,binSizeMs=1,windowSizeMs=200,probability=FALSE)
          {
            
            if(length(st@events)==0)
              stop("events is empty")
            
            
            nBins=(windowSizeMs*2)/binSizeMs
            windowSize=2*windowSizeMs*st@samplingRate/1000 # window size in res value from - to + extrems
            # call cwrapper function
            results<- .Call("crosscorrelationEvents_cwrap",
                            st@cellList,
                            length(st@cellList),
                            st@clu,
                            st@res,
                            st@nSpikes,
                            nBins,
                            windowSize,
                            st@startResIndexc,
                            st@endResIndexc,
                            length(st@startResIndexc),
                            probability,
                            st@events,
                            length(st@events))
            
            
            # return a data fram with clu time count
            if(probability==F){
              data.frame(clu=rep(st@cellList,each=nBins),
                         time=seq(-windowSizeMs+binSizeMs,windowSizeMs,binSizeMs)-(binSizeMs/2),
                         count=results)
            }
            else{
              data.frame(clu=rep(st@cellList,each=nBins),
                         time=seq(-windowSizeMs+binSizeMs,windowSizeMs,binSizeMs)-(binSizeMs/2),
                         prob=results)
            }
        })

### spikeTimeCrosscorrelation ###
# function using .Call() to call c wrapper
# It calls the the c function with the cwrap suffix
setGeneric(name="spikeTimeCrosscorrelation",
           def=function(st,binSizeMs=1,windowSizeMs=200,probability=FALSE)
           {standardGeneric("spikeTimeCrosscorrelation")})

setMethod(f="spikeTimeCrosscorrelation",
          signature = "SpikeTrain",
          definition=function(st,binSizeMs=1,windowSizeMs=200,probability=FALSE)
          {
            
            if(length(st@cellPairList[,1])==0)
              stop("cellPairList is empty")
          
            nBins=(windowSizeMs*2)/binSizeMs
            window.size=2*windowSizeMs*st@samplingRate/1000 # window size in res value from - to + extrems
            # call cwrapper function
            window.size
            
            results<- .Call("crosscorrelation_cwrap",
                            st@cellPairList[,1],
                            st@cellPairList[,2],
                            length(st@cellPairList[,1]),
                            st@clu,
                            st@res,
                            st@nSpikes,
                            nBins,
                            window.size,
                            st@startResIndexc,
                            st@endResIndexc,
                            length(st@startResIndexc),
                            probability)
            
            
            # return a data fram with clu time count
            if(probability==F)
              data.frame(clu1=rep(st@cellPairList[,1],each=nBins),
                         clu2=rep(st@cellPairList[,2],each=nBins),
                         time=seq(-windowSizeMs+binSizeMs,windowSizeMs,binSizeMs)-(binSizeMs/2),
                         count=results)
            else
              data.frame(clu1=rep(st@cellPairList[,1],each=nBins),
                         clu2=rep(st@cellPairList[,2],each=nBins),
                         time=seq(-windowSizeMs+binSizeMs,windowSizeMs,binSizeMs)-(binSizeMs/2),
                         prob=results)
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
            st@meanFiringRate<- .Call("meanFiringRate_cwrap",
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
            return(st)
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
              startIntervals<-as.numeric(s[,1])
              endIntervals<-as.numeric(s[,2])
              
            }else{
              startIntervals<-s
              endIntervals<-e
            }
            
            if(length(startIntervals)!=length(endIntervals))
              stop("problem with length of startIntervals and endIntervals")
            if(any(startIntervals>endIntervals))
              stop("startIntervals>endIntervals")
            if(any(diff(startIntervals)<0))
              stop("problem with chonology within startIntervals")
            if(any(diff(endIntervals)<0))
              stop("problem with chonology within endIntervals")
            if(any(startIntervals[-1]-endIntervals[-length(endIntervals)]<0))
              stop("problem with chronology between intervals, from end to next start")
            
            st@startInterval<-startIntervals
            st@endInterval<-endIntervals
            
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
            if(length(object@cellList)!=0){
              print(paste("n cellPairList:",length(object@cellPairList[,1])))
            }
            if(length(object@meanFiringRate)!=0){
              print(paste(object@meanFiringRate))
            }
            
            if(length(object@startInterval)!=0){
              print(paste("nIntervals:",length(object@startInterval))) 
              print(paste("Interval time:", sum(object@endInterval-object@startInterval)/object@samplingRate,"sec"))
              print(paste(object@startInterval,object@endInterval))
              print(paste("nIntervalsc:",length(object@startResIndexc)))
              print(paste(object@startResIndexc,object@endResIndexc))
            }
            if(length(object@events)!=0)
              print(paste("nEvents:",length(object@events)))
                  })