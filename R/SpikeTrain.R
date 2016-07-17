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
#' @slot startInterval Sample values of the beginning of intervals. To limit the analysis to some intervals.
#' Only data within the intervals are considered
#' @slot endInterval Sample values of the end of intervals. Only data within the intervals are considered.
#' @slot startResIndexc Index in the spike arrays for the start of intervals. Index is for a C array with 0 indexing.
#' The spike at the index is to be considered for analysis
#' @slot endResIndexc Index in the spike arrays for the start of intervals. Index is for a C array with 0 indexing.
#' The spike at the intdex is to be considered for analysis
#' @slot events Time in sample number for some events
#' @slot cellList Cell list
#' @slot auto Matrix holding spike-time autocorrelation
#' @slot autoMsPerBin Ms per bin in spike-time autocorrelation
#' @slot autoProbability Logical, whether the spike-time autocorrelation contains probability values or spike count
#' @slot cross Matrix holding spike-time crosscorrelation between spikes and events
#' @slot crossMsPerBin Ms per bin in spike-time crosscorrelation with events
#' @slot crossProbability Logical, whether the spike-time crosscorrelation to events contains probability values or spike count
#' @slot crossEvents Matrix holding spike-time crosscorrelation between spikes and events
#' @slot crossEventsMsPerBin Ms per bin in spike-time crosscorrelation with events
#' @slot crossEventsProbability Logical, whether the spike-time crosscorrelation to events contains probability values or spike count
#' @slot cellPairList Data frame containing pairs of cells
#' @slot ifrKernelSdMs Standard deviation of a gaussian kernel used to calculate the instantaneous firing rate
#' @slot ifrWindowSizeMs Window size of the instantaneous firing rate array
#' @slot ifrSpikeBinMs Size of the bins in which the number of spikes are counted
#' @slot ifr Matrix containing the instantaneous firing rate of the neurons
#' @slot ifrTime Time associated with each window of the instantaneous firing rate
#' @slot meanFiringRate Mean firing rate of the neurons
#' @slot isolationDistance Isolation distance of each cluster
#' @slot refractoryRatio Refractory ratio of each cluster
#' @slot crossRefractoryRatio Refractory ratio found in the spike-time crosscorrelation between a cluster 
#' and all other clusters. The smallest ratio is used.

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
          cross="matrix",
          crossMsPerBin="numeric",
          crossProbability="logical",
          crossEvents="matrix",
          crossEventsMsPerBin="numeric",
          crossEventsProbability="logical",
          # ifr
          ifrKernelSdMs="numeric",
          ifrWindowSizeMs="numeric",
          ifrSpikeBinMs="numeric",
          ifr="matrix",
          ifrTime="numeric",
          #
          meanFiringRate="numeric",
          isolationDistance="numeric",
          refractoryRatio="numeric",
          crossRefractoryRatio="numeric"
          ),  
  prototype = list(session=""))



#' Calculate the instantaneous firing rate from the spike trains.
#'
#' The ifr for the periods within the intervals set within the SpikeTrain object will be calculated.
#' A vector with the spike count in each bin is calculated
#' Then a gaussian kernel is applied to the spike count vector
#' Finally, the firing probability is integrated over a set window size and transform to a firing rate.
#' 
#'
#' @param st SpikeTrain object
#' @param windowSizeMs The window size for each ifr value
#' @param spikeBinMs The bin size for the spike count array
#' @param kernelSdMs Standard deviation of the gaussian kernel used to smooth the spike count vector
#' @return SpikeTrain object with the instantaneous firing rate
#' 
#' @docType methods
#' @rdname ifr-methods
setGeneric(name="ifr",
           def=function(st,windowSizeMs=100,spikeBinMs=1,kernelSdMs=100)
           {standardGeneric("ifr")})

#' @rdname ifr-methods
#' @aliases ifr,ANY,ANY-method
setMethod(f="ifr",
          signature = "SpikeTrain",
          definition=function(st,windowSizeMs=100,spikeBinMs=1,kernelSdMs=100)
          {
            if(st@nSpikes==0)
              stop(paste("ifr(): there are no spike in the SpikeTrain object",st@session))

            if(windowSizeMs<spikeBinMs)
              stop(paste("ifr(): windowSizeMs should not be smaller than spikeBinMs",
                         windowSizeMs,spikeBinMs))
            
            st@ifrWindowSizeMs<-windowSizeMs
            st@ifrSpikeBinMs<-spikeBinMs
            st@ifrKernelSdMs<-kernelSdMs
            
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


#' Calculate the correlation coefficient between the instantaneous firing rate of the neurons.
#'
#' Call ifr(st) method before calling this function.
#' Only the ifr falling withing the intervals set in the SpikeTime object will be considered.
#'
#' @param st SpikeTrain object with ifr already calculated
#' @return data.frame with the clu.ids of pairs and correlation coefficient between the instantaneous firing rate.
#' 
#' @docType methods
#' @rdname ifrAssociation-methods
setGeneric(name="ifrAssociation",
           def=function(st)
           {standardGeneric("ifrAssociation")})

#' @rdname ifrAssociation-methods
#' @aliases ifrAssociation,ANY,ANY-method
setMethod(f="ifrAssociation",
          signature = "SpikeTrain",
          definition=function(st)
          {
            if(st@nSpikes==0)
              stop(paste("ifrAssociation(): there are no spike in the SpikeTrain object",st@session))
            if(length(st@ifr)==0)
              stop(paste("ifrAssociation(): st@ifr has a size of 0",st@session))
            
            ## get the ifr and ifrTime inside st@intervals
            resTime<-st@ifrTime*st@samplingRate
            index<-as.logical(.Call("resWithinIntervals",
                                    length(st@startInterval),
                                    as.integer(st@startInterval),
                                    as.integer(st@endInterval),
                                    length(resTime),
                                    as.integer(resTime)))
            
            ifrSel<-matrix(st@ifr[,index],nrow=length(st@cellList))
            m<-cor(t(ifrSel))
            r<-m[which(lower.tri(m,diag=FALSE))]
            if(length(r)!=length(st@cellPairList[,1]))
              stop(paste("ifrAssociation(): length of r is not equal to the number of cell pair list"))
            return(data.frame(clu.id1=st@cellPairList[,1],
                       clu.id2=st@cellPairList[,2],
                       r=r))
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
            st@endInterval<-max(st@res)+1 ## add one so that the last spike is within the intervals by default
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
            st@endInterval<-max(st@res)+1 # add one so that the last spike is considered
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
#' @param windowSizeMs Default is 200, meaning that it ranges from + and - 200
#' @param probability If TRUE, will calculate the probability of a spike in a given bin instead of the spike count
#' @return st SpikeTrain object with the slot crossEvents filled
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
            
            st@crossEvents=matrix(results,nrow=nBins,ncol=length(st@cellList))
            st@crossEventsMsPerBin=binSizeMs
            st@crossEventsProbability=probability
            return(st) 
          }
          )
            
            
          


#' Get spike-time crosscorelation to events as data.frame
#' 
#'
#' @param st SpikeTrain object
#' @return data.frame with spike-time crosscorrelation to events
#' 
#' @docType methods
#' @rdname spikeTimeCrosscorrelationEventsAsDataFrame-methods
setGeneric(name="spikeTimeCrosscorrelationEventsAsDataFrame",
           def=function(st)
           {standardGeneric("spikeTimeCrosscorrelationEventsAsDataFrame")})
#' @rdname spikeTimeCrosscorrelationEventsAsDataFrame-methods
#' @aliases spikeTimeCrosscorrelationEventsAsDataFrame,ANY,ANY-method
setMethod(f="spikeTimeCrosscorrelationEventsAsDataFrame",
          signature = "SpikeTrain",
          definition=function(st)
          {
            nBins=dim(st@crossEvents)[1]
            # return a data fram with clu time count
            if(st@crossEventsProbability==F)
              data.frame(clu=rep(st@cellList,each=nBins),
                         time=rep(seq(-st@crossEventsMsPerBin*nBins/2+st@crossEventsMsPerBin/2,
                                      st@crossEventsMsPerBin*nBins/2,
                                      st@crossEventsMsPerBin),length(st@cellList)),
                         count=as.numeric(st@crossEvents))
            else{
              data.frame(clu=rep(st@cellList,each=nBins),
                         time=rep(seq(-st@crossEventsMsPerBin*nBins/2+st@crossEventsMsPerBin/2,
                                      st@crossEventsMsPerBin*nBins/2,
                                      st@crossEventsMsPerBin),length(st@cellList)),
                         prob=as.numeric(st@crossEvents))
            }
          }
)
 

#' Calculate the spike-time crosscorrelation between the spike trains of cell pairs
#' 
#' The spikes of cell 1 are used in turn as a reference spike.
#' The number of spikes or probability to observe a spike of cell 2 around the reference spikes is calculated.
#' You can set the bins size in ms and and the time window for which you want to do the analysis on.
#'
#'
#' @param st SpikeTrain object
#' @param binSizeMs Default is 1
#' @param windowSizeMs Default is 200. Will span from -windowSizeMs to windowSizeMs
#' @param probability If TRUE, will calculate the probability of a spike in a given bin instead of the spike count
#' @return SpikeTrain object with spike-time crosscorrelation of cell pairs in slot cross
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
            
            results<- .Call("crosscorrelation_cwrap",
                            as.integer(st@cellPairList[,1]),
                            as.integer(st@cellPairList[,2]),
                            length(st@cellPairList[,1]),
                            as.integer(st@clu),
                            as.integer(st@res),
                            st@nSpikes,
                            as.integer(nBins),
                            as.integer(window.size),
                            as.integer(st@startResIndexc),
                            as.integer(st@endResIndexc),
                            length(st@startResIndexc),
                            probability)
            
            st@cross=matrix(results,nrow=nBins,ncol=length(st@cellPairList[,1]))
            st@crossMsPerBin=binSizeMs
            st@crossProbability=probability
            return(st)
          })




#' Get spike-time crosscorelation between cells as data.frame
#' 
#'
#' @param st SpikeTrain object
#' @return data.frame with spike-time crosscorrelation between neurons
#' 
#' @docType methods
#' @rdname spikeTimeCrosscorrelationAsDataFrame-methods
setGeneric(name="spikeTimeCrosscorrelationAsDataFrame",
           def=function(st)
           {standardGeneric("spikeTimeCrosscorrelationAsDataFrame")})
#' @rdname spikeTimeCrosscorrelationAsDataFrame-methods
#' @aliases spikeTimeCrosscorrelationAsDataFrame,ANY,ANY-method
setMethod(f="spikeTimeCrosscorrelationAsDataFrame",
          signature = "SpikeTrain",
          definition=function(st)
          {
            nBins=dim(st@cross)[1]
            # return a data fram with clu time count
            if(st@crossProbability==F)
              data.frame(clu1=rep(st@cellPairList[,1],each=nBins),
                         clu2=rep(st@cellPairList[,2],each=nBins),
                         time=rep(seq(-st@crossMsPerBin*nBins/2+st@crossMsPerBin/2,
                                      st@crossMsPerBin*nBins/2,
                                      st@crossMsPerBin),length(st@cellPairList[,1])),
                         count=as.numeric(st@cross))
            else{
              data.frame(clu1=rep(st@cellPairList[,1],each=nBins),
                         clu2=rep(st@cellPairList[,2],each=nBins),
                         time=rep(seq(-st@crossMsPerBin*nBins/2+st@crossMsPerBin/2,
                                      st@crossMsPerBin*nBins/2,
                                      st@crossMsPerBin),length(st@cellPairList[1])),
                         prob=as.numeric(st@cross))
            }
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
                            as.integer(st@startInterval),
                            as.integer(st@endInterval),
                            as.integer(st@nSpikes),
                            as.integer(st@res),
                            as.integer(0))# will keep the intervals after the last spikes
            
            st@startResIndexc<-results[1,]
            st@endResIndexc<-results[2,]
            st@startInterval<-st@startInterval[1:length(st@startResIndexc)]
            st@endInterval<-st@endInterval[1:length(st@endResIndexc)]
            
            return(st)
          }
)


#' Set the list of cells to limit the analysis to these cells
#' 
#' Only these cells will be considered for analysis.
#'
#' @param st SpikeTrain object
#' @param cellList Numiric vector containing the clu id of the neurons
#' @return a SpikeTrain object with a new cellList.
#' 
#' @docType methods
#' @rdname setCellList-methods
setGeneric(name="setCellList",
           def=function(st,cellList)
           {standardGeneric("setCellList")})

#' @rdname setCellList-methods
#' @aliases setCellList,ANY,ANY-method
setMethod(f="setCellList",
          signature = "SpikeTrain",
          definition=function(st,cellList)
          {
            if(class(cellList)!="numeric")
              stop("setCellList: cellList is not a numeric")
            if(length(cellList)==0)
              stop("setCellList: length(cellList==0)")
            st@cellList<-cellList
            st@nCells<-length(cellList)
            
            ## spikes per cell
            st@nSpikesPerCell<-as.numeric(table(st@clu)[as.character(cellList)])
            
            if(any(is.na(st@nSpikesPerCell)))
              stop("setCellList: a cell has no spike")
            
            # all the possible pairs
            if(length(st@cellList>1))
              st@cellPairList<-makePairs(st@cellList)
            
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

#' Get the isolation distance of each cluster
#' 
#' For each cluster, the mahalanobis distance of all spikes relative to the cluster
#' is calculated. The isolation distance is the minimal distance from the cluster center
#' at which there are as many spikes from other cluster than from the cluster of interest.
#' Distance is not defined for cases in which the number of cluster spikes is 
#' greater than the number of noise spikes.
#' 
#' From Schmitzer-Torbert et al. 2005
#'
#' @param st SpikeTrain object
#' @param cg CellGroup object
#' @return a SpikeTrain object with the isolationDistance set.
#' 
#' @docType methods
#' @rdname isolationDistance-methods
setGeneric(name="isolationDistance",
           def=function(st,cg)
           {standardGeneric("isolationDistance")})
#' @rdname isolationDistance-methods
#' @aliases isolationDistance,ANY,ANY-method
setMethod(f="isolationDistance",
          signature = "SpikeTrain",
          definition=function(st,cg)
          {
            if(class(cg)!="CellGroup")
              stop("isolationDistance: cg should be an object of class CellGroup")
            if(st@session=="")
              stop("st@session is not set")
            if(st@path=="") ## path is given or is getwd()
              st@path=getwd()
            pathSession=paste(st@path,st@session,sep="/")
            index=1
            for(tetrode in unique(cg@tetrode)){
            
              clus<-cg@cluToTetrode[which(cg@tetrode==tetrode)]
            
              ## check if the files are there
              if(!file.exists(paste(pathSession,"fet",tetrode,sep=".")))
                stop("isolationDistance needs",paste(pathSession,"fet",tetrode,sep="."))
              if(!file.exists(paste(pathSession,"clu",tetrode,sep=".")))
                stop("isolationDistance needs",paste(pathSession,"clu",tetrode,sep="."))
              
              ## load the tetrode clu file
              cluTet<-as.numeric(.Call("read_one_column_int_file_cwrap", paste(pathSession,"clu",tetrode,sep=".")))
              cluTet<-cluTet[-1]# remove first line
              
              ## load the tetrode fet file
              fet<-readFetFile(paste(pathSession,"fet",tetrode,sep="."))
              fet<-fet[,1:(ncol(fet)-1)] # remove time, last column
              
              if(length(cluTet)!=nrow(fet))
                stop("isolationDistance, length of cluTet and fet not equal")
              
              for(clu in clus){
                # get the mahalanobis distance of each spikes relative to clu
                fetClu<-fet[which(cluTet==clu),]
                maha<-mahalanobis(x=fet,colMeans(fetClu),cov=cov(fetClu))
                cm<-cbind(cluTet,maha) # cluTet - mahalanobis distance
                cm<-cm[order(cm[,2]),] ## order spikes according to distance
                propSpikeFromClu<-cumsum(cm[,1]==clu)/1:length(cm[,1]) ## probability that spikes are from clu
                if(any(propSpikeFromClu<.5)){
                  st@isolationDistance[index]=min(cm[which(propSpikeFromClu<.5),2])  
                } else{
                  st@isolationDistance[index]=NA
                }
                index<-index+1
              }
            }
            return(st)
          }
)

#' Get refractory ratio of each cluster from its spike-time autocorrelation
#' 
#' This is the ratio between the mean number of spikes falling in the bins of the refractory period
#' compared to the max number of spikes falling in one bin of the control period outside the refractory period.
#' Note that the spike-time autocorrelations in the SpikeTrain object will be modified.
#'
#' @param st SpikeTrain object
#' @param refractoryMs Length of the refractory period in ms
#' @param binSizeMs Size of the bins in the spike-time autocorrelation
#' @param windowSizeMs Size of the window used to construct the spike-time autocorrelation
#' @param minControlWindowMs Minimum time of the control window
#' @param maxControlWindowMs Maximal time of the control window
#' @return a SpikeTrain object with the refractoryRation set.
#' 
#' @docType methods
#' @rdname refractoryRatio-methods
setGeneric(name="refractoryRatio",
           def=function(st,refractoryMs=1.5,binSizeMs=0.5,windowSizeMs=25,
                        minControlWindowMs=5.0,maxControlWindowMs=25)
           {standardGeneric("refractoryRatio")})
#' @rdname refractoryRatio-methods
#' @aliases refractoryRatio,ANY,ANY-method
setMethod(f="refractoryRatio",
          signature = "SpikeTrain",
          definition=function(st,refractoryMs=1.5,binSizeMs=0.5,windowSizeMs=25,
                              minControlWindowMs=5.0,maxControlWindowMs=25)
          {
            
            if(st@session=="")
              stop("st@session is not set")
            if(st@path=="") ## path is given or is getwd()
              st@path=getwd()
            if(refractoryMs>windowSizeMs)
              stop("refractoryRatio, refractoryMs>windowSizeMs")
            if(minControlWindowMs>windowSizeMs)
              stop("refractoryRatio, minControlWindowMs>windowSizeMs")
            if(maxControlWindowMs>windowSizeMs)
              stop("refractoryRatio, maxControlWindowMs>windowSizeMs")
            if(minControlWindowMs>=maxControlWindowMs)
              stop("refractoryRatio, minControlWindowMs>=maxControlWindowMs")
            st<-spikeTimeAutocorrelation(st,binSizeMs=binSizeMs,
                                      windowSizeMs=windowSizeMs,probability=F)
            st@refractoryRatio<-apply(st@auto,2,function(x,autoMsPerBin,refractoryMs,minControlWindowMs,maxControlWindowMs)
              {
                time<-seq(-st@autoMsPerBin*length(x)/2+st@autoMsPerBin/2,
                             st@autoMsPerBin*length(x)/2,st@autoMsPerBin)
                ref<-mean(x[which(time>=0&time<=refractoryMs)])
                con<-max(x[which(time>=minControlWindowMs&time<=maxControlWindowMs)])
                if(con==0)
                {
                  return(NA)
                } else{
                  return(ref/con)
                  
                }
            },
            st@autoMsPerBin,
            refractoryMs,
            minControlWindowMs,
            maxControlWindowMs)
            return(st)
          }
)

#' Get the cross refractory ratio of each cluster from its spike-time crosscorrelation
#' 
#' This tells you whether a cluster has a common refractory period with other cluster. If this is the 
#' case, perhaps the two clusters are from the same neuron.
#' 
#' This is the ratio between the mean number of spikes falling in the bins of the refractory period
#' compared to the max number of spikes falling in one bin of the control period outside the refractory period.
#' Note that the spike-time crosscorrelations in the SpikeTrain object will be modified.
#'
#' @param st SpikeTrain object
#' @param refractoryMs Length of the refractory period in ms
#' @param binSizeMs Size of the bins in the spike-time autocorrelation
#' @param windowSizeMs Size of the window used to construct the spike-time autocorrelation
#' @param minControlWindowMs Minimum time of the control window
#' @param maxControlWindowMs Maximal time of the control window
#' @return a SpikeTrain object with the refractoryRation set.
#' 
#' @docType methods
#' @rdname crossRefractoryRatio-methods
setGeneric(name="crossRefractoryRatio",
           def=function(st,refractoryMs,binSizeMs,windowSizeMs,
                        minControlWindowMs,maxControlWindowMs)
           {standardGeneric("crossRefractoryRatio")})
#' @rdname crossRefractoryRatio-methods
#' @aliases crossRefractoryRatio,ANY,ANY-method
setMethod(f="crossRefractoryRatio",
          signature = "SpikeTrain",
          definition=function(st,refractoryMs=1.5,binSizeMs=0.5,windowSizeMs=25,
                              minControlWindowMs=5.0,maxControlWindowMs=25)
          {
            if(st@session=="")
              stop("st@session is not set")
            if(st@path=="") ## path is given or is getwd()
              st@path=getwd()
            if(refractoryMs>windowSizeMs)
              stop("refractoryRatio, refractoryMs>windowSizeMs")
            if(minControlWindowMs>windowSizeMs)
              stop("refractoryRatio, minControlWindowMs>windowSizeMs")
            if(maxControlWindowMs>windowSizeMs)
              stop("refractoryRatio, maxControlWindowMs>windowSizeMs")
            if(minControlWindowMs>=maxControlWindowMs)
              stop("refractoryRatio, minControlWindowMs>=maxControlWindowMs")
            if(st@nCells==1){
              st@crossRefractoryRatio[1]<-NA
              return(st)
            }
            st<-spikeTimeCrosscorrelation(st,binSizeMs=binSizeMs,
                                         windowSizeMs=windowSizeMs,probability=F)
            
            refractoryRatio<-apply(st@cross,2,function(x,crossMsPerBin,refractoryMs,minControlWindowMs,maxControlWindowMs)
            {
              time<-seq(-st@crossMsPerBin*length(x)/2+st@crossMsPerBin/2,
                        st@crossMsPerBin*length(x)/2,st@crossMsPerBin)
              ref<-mean(x[which(time>=0&time<=refractoryMs)])
              con<-max(x[which(time>=minControlWindowMs&time<=maxControlWindowMs)])
              if(con==0)
              {
                return(NA)
              } else{
                return(ref/con)
                
              }
            },
            st@crossMsPerBin,
            refractoryMs,
            minControlWindowMs,
            maxControlWindowMs)
            
            # get the sum of spikes in the positive side of the crosscorrelation
            if(ncol(st@cross)>1){
              totalSpikes<-colSums(st@cross[ (nrow(st@cross)/2):(nrow(st@cross))  , ])
            } else {
              totalSpikes<-sum(st@cross[(length(st@cross)/2):(length(st@cross))])
            }
            
            # don't consider the cc with very few spikes
            refractoryRatio[which(totalSpikes<200)]<-NA
            
            for(clu in 1:st@nCells){
              if(all(is.na(refractoryRatio[which(st@cellPairList[,1]==st@cellList[clu]|
                                                 st@cellPairList[,2]==st@cellList[clu])]))){
                st@crossRefractoryRatio[clu]<-NA
              }else{
                st@crossRefractoryRatio[clu]<-min(refractoryRatio[which(st@cellPairList[,1]==st@cellList[clu]|
                                                                      st@cellPairList[,2]==st@cellList[clu])],na.rm=T)
              }
            }
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
            if(length(object@isolationDistance)!=0){
              print(paste("isolation distance:"))
              print(paste(object@isolationDistance))
            }
            if(length(object@refractoryRatio)!=0){
              print(paste("refractory ratio:"))
              print(paste(object@refractoryRatio))
            }
            if(length(object@crossRefractoryRatio)!=0){
              print(paste("crosscorrelation refractory ratio:"))
              print(paste(object@crossRefractoryRatio))
            }
            
            
            if(length(object@startInterval)!=0){
              print(paste("nIntervals:",length(object@startInterval))) 
              print(paste("Interval time:", sum(object@endInterval-object@startInterval)/object@samplingRate,"sec"))
              if(length(object@startInterval<500)){
                print(paste(object@startInterval,object@endInterval))
              }else{
                print(head(paste(object@startInterval,object@endInterval),n=200))
              }
              print(paste("nIntervalsc:",length(object@startResIndexc)))
              if(length(object@startInterval<500)){
                print(paste(object@startResIndexc,object@endResIndexc))
              }else{
                print(head(paste(object@startResIndexc,object@endResIndexc),n=200))
              }
            }
            if(length(object@events)!=0)
              print(paste("nEvents:",length(object@events)))
                  })