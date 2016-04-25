#' A S4 class to represent spike trains of one recording session
#' 
#' Set intervals in this object to limit the analysis to some parts of the 
#' recording session.
#' 
#' @slot session Name of the recording session
#' @slot path Directory where the recording session is located
#' @slot samplingRate Sampling rate of the electrophysiological data
#' @slot res Spike time in sample value
#' @slot time Spike time in seconds
#' @slot clu Cluster id for each spike
#' @slot nCells Number of cells in the recording session
#' @slot nSpikes Number of spikes in the recording session
#' @slot nSpikesPerCell Number of spikes for each cell
#' @slot startInterval Sample values of the beginning of intervals. To limit the analysis to some intervals
#' @slot endInterval Sample values of the end of intervals
#' @slot startResIndexc Index in the spike arrays for the start of intervals. Index is for a C array with 0 indexing.
#' @slot endResIndexc Index in the spike arrays for the start of intervals. Index is for a C array with 0 indexing.
#' @slot events Time in sample number for some events
#' @slot cellList Cell list
#' @slot auto Matrix holding spike-time autocorrelation
#' @slot autoMsPerBin Ms per bin in spike-time autocorrelation
#' @slot autoProbability Logical, whether the spike-time autocorrelation contains probability values or spike count
#' @slot cellPairList Data frame containing pairs of cells
#' @slot ifrKernelSdMs Standard deviation of a gaussian kernel used to calculate the instantaneous firing rate
#' @slot ifrWindowSizeMs Window size of the instantaneous firing rate array
#' @slot ifrSpikeBinMs Size of the bins in which the number of spikes are counted
#' @slot ifr Matrix containing the instantaneous firing rate of the neurons
#' @slot ifrTime Time associated with each window of the instantaneous firing rate
#' @slot meanFiringRate Mean firing rate of the neurons
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
          auto="matrix",
          autoMsPerBin="numeric",
          autoProbability="logical",
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



#' Calculate the instantaneous firing rate from the spike trains
#'
#'
#' @param st SpikeTrain object
#' @return SpikeTrain object with the instantaneous firing rate
#' 
#' @docType methods
#' @rdname ifr-methods
setGeneric(name="ifr",
           def=function(st)
           {standardGeneric("ifr")})

#' @rdname ifr-methods
#' @aliases ifr,ANY,ANY-method
setMethod(f="ifr",
          signature = "SpikeTrain",
          definition=function(st)
          {
            if(st@nSpikes==0)
              stop(paste("ifr(): there are no spike in the SpikeTrain object",st@session))

            if(st@ifrWindowSizeMs<st@ifrSpikeBinMs)
              stop(paste("ifr(): ifrWindowSizeMs should not be smaller than ifrSpikeBinMs",
                         st@ifrWindowSizeMs,st@ifrSpikeBinMs))
            
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


#' Load the spike train from the .clu and .res files
#'
#'
#' @param st SpikeTrain object
#' @return SpikeTrain object with the instantaneous firing rate
#' 
#' @docType methods
#' @rdname loadSpikeTrain-methods
setGeneric(name="loadSpikeTrain",
           def=function(st)
           {standardGeneric("loadSpikeTrain")}
)
#' @rdname loadSpikeTrain-methods
#' @aliases loadSpikeTrain,ANY,ANY-method
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
              stop(paste("samplingRate is out of range:",st@samplingRate,st@session))
            st@clu<-st@clu[-1] ## remove first number
            if(length(st@res)!=length(st@clu))
              stop(paste("length of res and clu files not equals:",length(st@res),length(st@clu),st@session))
            # remove noise cluster 1 or 0
            df<-data.frame(res=st@res,clu=st@clu)
            df<-df[which(df$clu>1),]
            st@res<-df$res
            st@clu<-df$clu
            if(length(st@res)==0)
              stop(paste("No valid spike in session",st@session))
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

#' Set new spike trains
#' 
#' This function is used to set spike trains 
#'
#'
#' @param st SpikeTrain object
#' @param res Spike times in sample number
#' @param clu Cluster id for each spike
#' @param samplingRate The number of spamples per second (Hz)
#' @return SpikeTrain object with the instantaneous firing rate
#' 
#' @docType methods
#' @rdname setSpikeTrain-methods
setGeneric(name="setSpikeTrain",
           def=function(st,res,clu,samplingRate)
           {standardGeneric("setSpikeTrain")}
)

#' @rdname setSpikeTrain-methods
#' @aliases setSpikeTrain,ANY,ANY-method
setMethod(f="setSpikeTrain",
          signature="SpikeTrain",
          definition=function(st,res,clu,samplingRate)
          {
            st@res<-res
            st@clu<-clu
            st@samplingRate<-samplingRate
            
            if(st@samplingRate<1 | st@samplingRate > 100000)
              stop(paste("samplingRate is out of range:",st@samplingRate,st@session))
            
            if(length(st@res)!=length(st@clu))
              stop(paste("length of res and clu files not equals:",length(st@res),length(st@clu),st@session))
            
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



#' Calculate the spike-time autocorrelation
#' 
#' Each spike is treated in turn as a reference spike.
#' The number of spikes or probability to observe a spike around the reference spike is calculated.
#' You can set the bins size in ms and and the time window for which you want to do the analysis on.
#'
#' @param st SpikeTrain object
#' @param binSizeMs Default is 1
#' @param windowSizeMs Default is 200. This means that autocorrelation ranges from -windowSizeMs to windowSizeMs
#' @param probability If TRUE, will calculate the probability of a spike in a given bin instead of the spike count
#' @return SpikeTrain object with autocorrelation in slot auto
#' 
#' @docType methods
#' @rdname spikeTimeAutocorrelation-methods
setGeneric(name="spikeTimeAutocorrelation",
           def=function(st,binSizeMs,windowSizeMs,probability)
             {standardGeneric("spikeTimeAutocorrelation")})
#' @rdname spikeTimeAutocorrelation-methods
#' @aliases spikeTimeAutocorrelation,ANY,ANY-method
setMethod(f="spikeTimeAutocorrelation",
          signature = "SpikeTrain",
          definition=function(st, binSizeMs=1,windowSizeMs=200,probability=FALSE)
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
            
            st@auto<-matrix(results,nrow=nBins,ncol=length(st@cellList))
            st@autoMsPerBin=binSizeMs
            st@autoProbability=probability
            return(st)
            }
)



#' Get spike-time autocorrelation as data.frame
#' 
#'
#' @param st SpikeTrain object
#' @return data.frame with spike-time autocorrelation
#' 
#' @docType methods
#' @rdname spikeTimeAutocorrelationAsDataFrame-methods
setGeneric(name="spikeTimeAutocorrelationAsDataFrame",
           def=function(st)
           {standardGeneric("spikeTimeAutocorrelationAsDataFrame")})
#' @rdname spikeTimeAutocorrelationAsDataFrame-methods
#' @aliases spikeTimeAutocorrelationAsDataFrame,ANY,ANY-method
setMethod(f="spikeTimeAutocorrelationAsDataFrame",
          signature = "SpikeTrain",
          definition=function(st)
          {
            nBins=dim(st@auto)[1]
            # return a data fram with clu time count
            if(st@autoProbability==F)
              data.frame(clu=rep(st@cellList,each=nBins),
                         time=rep(seq(-st@autoMsPerBin*nBins/2+st@autoMsPerBin/2,
                              st@autoMsPerBin*nBins/2,
                              st@autoMsPerBin),length(st@cellList)),
                         count=as.numeric(st@auto))
            else{
              data.frame(clu=rep(st@cellList,each=nBins),
                         time=rep(seq(-st@autoMsPerBin*nBins/2+st@autoMsPerBin/2,
                                  st@autoMsPerBin*nBins/2,
                                  st@autoMsPerBin),length(st@cellList)),
                         prob=as.numeric(st@auto))
            }
          }
)


#' Calculate the spike-time crosscorrelation between the spike trains and a list of events
#' 
#' Each event is treated in turn as a reference event.
#' The number of spikes or probability to observe a spike around the reference event is calculated.
#' You can set the bins size in ms and and the time window for which you want to do the analysis on.
#'
#'
#' @param st SpikeTrain object
#' @param binSizeMs Default is 1
#' @param windowSizeMs Default is 200
#' @param probability If TRUE, will calculate the probability of a spike in a given bin instead of the spike count
#' @return a data.frame with the spike-time crosscorrelation
#' 
#' @docType methods
#' @rdname spikeTimeCrosscorrelationEvents-methods
setGeneric(name="spikeTimeCrosscorrelationEvents",
           def=function(st,binSizeMs=1,windowSizeMs=200,probability=FALSE)
           {standardGeneric("spikeTimeCrosscorrelationEvents")})
#' @rdname spikeTimeCrosscorrelationEvents-methods
#' @aliases spikeTimeCrosscorrelationEvents,ANY,ANY-method
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


#' Calculate the spike-time crosscorrelation between the spike trains of cell pairs
#' 
#' Each spike of cell 1 are used in turn as a reference spike.
#' The number of spikes or probability to observe a spike of cell 2 around the reference spikes is calculated.
#' You can set the bins size in ms and and the time window for which you want to do the analysis on.
#'
#'
#' @param st SpikeTrain object
#' @param binSizeMs Default is 1
#' @param windowSizeMs Default is 200
#' @param probability If TRUE, will calculate the probability of a spike in a given bin instead of the spike count
#' @return a data.frame with the spike-time crosscorrelation of cell pairs
#' 
#' @docType methods
#' @rdname spikeTimeCrosscorrelation-methods
setGeneric(name="spikeTimeCrosscorrelation",
           def=function(st,binSizeMs=1,windowSizeMs=200,probability=FALSE)
           {standardGeneric("spikeTimeCrosscorrelation")})
#' @rdname spikeTimeCrosscorrelation-methods
#' @aliases spikeTimeCrosscorrelation,ANY,ANY-method
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


#' Calculate the mean firing rate (Hz) of each neuron in a SpikeTrain object
#' 
#' This is simply the number of spikes divided by the time within the intervals
#' set in the SpikeTrain object.
#'
#'
#' @param st SpikeTrain object
#' @return a SpikeTrain object with the mean firing rate in slot meanFiringRate
#' 
#' @docType methods
#' @rdname meanFiringRate-methods
setGeneric(name="meanFiringRate",
           def=function(st)
           {standardGeneric("meanFiringRate")})

#' @rdname meanFiringRate-methods
#' @aliases meanFiringRate,ANY,ANY-method
setMethod(f="meanFiringRate",
          signature = "SpikeTrain",
          definition=function(st)
          {
            
            # call cwrapper function
            st@meanFiringRate<- .Call("meanFiringRate_cwrap",
                                      as.integer(st@cellList),
                                      length(st@cellList),
                                      as.integer(st@clu),
                                      as.integer(st@res),
                                      as.integer(st@nSpikes),
                                      as.integer(st@startInterval),
                                      as.integer(st@endInterval),
                                      as.integer(st@startResIndexc),
                                      as.integer(st@endResIndexc),
                                      length(st@startResIndexc),
                                      as.integer(st@samplingRate))
            return(st)
          }
)


#' Set time intervals to limit the period used in the analysis
#' 
#' Only the data within the intervals are used for analysis. 
#' For example with interval 0-20000, a spike at time 0 or 20000 is not included.
#' Spikes between time 0 and 20000 are used.
#' These intervals are used in SpikeTrain methods and also in 
#' methods of other classes (SaptialProperties2d, SpatialProperties1d, HeadDirection, etc.).
#' By default, the intervals are set from 0 to time point of the last recorded spike.
#'
#' @param st SpikeTrain object
#' @param s Vector or matrix. If a matrix, should have 2 columns (beginning and end of intervals).
#' If a vector, beginning of intervals.
#' @param e Vector, end of the intervals. If s is a matrix, e is not used.
#' @return a SpikeTrain object with the intervals set.
#' 
#' @docType methods
#' @rdname setIntervals-methods
setGeneric(name="setIntervals",
           def=function(st,s,e)
           {standardGeneric("setIntervals")})

#' @rdname setIntervals-methods
#' @aliases setIntervals,ANY,ANY-method
setMethod(f="setIntervals",
          signature = "SpikeTrain",
          definition=function(st,s,e)
          {
            
            ## if s is a matrix, then e is ignored
            if(class(s)=="matrix"){
              if(dim(s)[2]!=2){
                stop("matrix should have 2 columns")
              }
              startIntervals<-as.numeric(s[,1])
              endIntervals<-as.numeric(s[,2])
            }else{
              startIntervals<-s
              endIntervals<-e
            }
            if(length(startIntervals)!=length(endIntervals))
              stop(paste("problem with length of startIntervals and endIntervals in set intervals",st@session))
            if(any(startIntervals>endIntervals))
              stop(paste("startIntervals>endIntervals",st@session))
            if(any(diff(startIntervals)<0))
              stop(paste("problem with chonology within startIntervals in set intervals",st@session))
            if(any(diff(endIntervals)<0))
              stop(paste("problem with chonology within endIntervals in set intervals",st@session))
            if(any(startIntervals[-1]-endIntervals[-length(endIntervals)]<0))
              stop(paste("problem with chronology between intervals, from end to next start in set intervals",st@session))
            
            st@startInterval<-startIntervals
            st@endInterval<-endIntervals
            
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


#' Set some events for spike-time crosscorrelation to the events. 
#' 
#' These events could be laser stimulation or some behavioural events
#'
#' @param st SpikeTrain object
#' @param events Time in sample number of the events
#' @return a SpikeTrain object with the events set
#' 
#' @docType methods
#' @rdname setEvents-methods
setGeneric(name="setEvents",
           def=function(st,events)
           {standardGeneric("setEvents")})
#' @rdname setEvents-methods
#' @aliases setEvents,ANY,ANY-method
setMethod(f="setEvents",
          signature = "SpikeTrain",
          definition=function(st,events=NULL)
          {
            if(is.null(events))
              stop(paste("events is empty",st@session))
            if(any(events<0))
              stop(paste("negative values as events",st@session))
            if(any(diff(events)<0))
              stop(paste("problem with the chronology of the events",st@session))
            
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
            if(length(object@cellList)>1){
              print(paste("n cellPairList:",length(object@cellPairList[,1])))
            }
            if(length(object@meanFiringRate)!=0){
              print(paste("firing rate (Hz):"))
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