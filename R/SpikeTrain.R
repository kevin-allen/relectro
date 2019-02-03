#' A S4 class to represent the spike trains of one recording session
#'
#' This class is used to analyse the spike trains of neurons. You can calculate mean firing rates,
#' instantaneous firing rates, spike-time autocorrelations, spike-time crosscorrelations,
#' crosscorrelation between spike trains and events, etc.
#' You can limit the analysis to some time interval by using the method setInterval.
#' Similarly, you can limit the analysis to some cells with the method setCellList.
#'
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
#' @slot autoTimePoints Time points for data points in the spike-time autocorrelation
#' @slot autoProbability Logical, whether the spike-time autocorrelation contains probability values or spike count
#' @slot autoCOM Center of mass of the positive half of each spike-time autocorrelation
#' @slot cross Matrix holding spike-time crosscorrelation between spikes and events
#' @slot crossMsPerBin Ms per bin in spike-time crosscorrelation with events
#' @slot crossTimePoints Time points for data points in the spike-time crosscorrelation
#' @slot crossProbability Logical, whether the spike-time crosscorrelation to events contains probability values or spike count
#' @slot crossEvents Matrix holding spike-time crosscorrelation between spikes and events
#' @slot crossEventsMsPerBin Ms per bin in spike-time crosscorrelation with events
#' @slot crossEventsTimePoints Time points for data points in the spike-time crosscorrelation to events
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
          autoTimePoints="numeric",
          autoProbability="logical",
          autoCOM="numeric",
          cross="matrix",
          crossMsPerBin="numeric",
          crossTimePoints="numeric",
          crossProbability="logical",
          crossEvents="matrix",
          crossEventsMsPerBin="numeric",
          crossEventsTimePoints="numeric",
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
#' A vector with the spike count in each bin is calculated.
#' Then a gaussian kernel is applied to the spike count vector.
#' Finally, the firing probability is integrated over a set window size and transform to a firing rate.
#'
#' If you want to return the ifr matrix in a function called by runOnSessions, you need to transpose the ifr array so that the 
#' data from all recording sessions can be joined together. The matrices of different sessions also need to have the same number of rows.
#' This could be accomplished by setting the same intervals in the SpikeTrain object for all recording sessions
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


#' Calculate the power spectrum of instantaneous firing rate from the spike trains.
#'
#' This function can use different function to get the power spectrum
#' of the instantaneous firing rate. One is the pwelch function of the oce R package.
#' 
#'
#' @param st SpikeTrain object
#' @param nfft The size of the fast fourier transform
#' @param method Method used to calculate the power spectrum. Currently supported: pwelch
#' @return list with the frequency and the power for each cell
#'
#' @docType methods
#' @rdname ifrPowerSpectrum-methods
setGeneric(name="ifrPowerSpectrum",
           def=function(st,nfft=1024,method="pwelch")
           {standardGeneric("ifrPowerSpectrum")})

#' @rdname ifrPowerSpectrum-methods
#' @aliases ifrPowerSpectrum,ANY,ANY-method
setMethod(f="ifrPowerSpectrum",
          signature = "SpikeTrain",
          definition=function(st,nfft=1024,method="pwelch")
          {
            if(st@nSpikes==0)
              stop(paste("ifr(): there are no spike in the SpikeTrain object",st@session))
            if(dim(st@ifr)[1]<1)
              stop(paste("st@ifr has less then 1 column, run ifr() before calling this function"))
            if(st@ifrWindowSizeMs<=0)
              stop(paste("st@ifrWindowSizeMs is smaller or equal to  0"))
            ## number of observation by unit of time (second)
            Fs=1000/st@ifrWindowSizeMs
            
            if(method=="pwelch"){
              ## calculate the power spectrum of each time serie
              ## we could speed up things with the snow package here 
              out<-apply(st@ifr,1,oce::pwelch,nfft=nfft,plot=FALSE,log='no',fs=Fs)
              # get the frequencies
              freq<-out[[1]]$freq
              # get the power values of each list
              ps<-matrix(unlist(lapply(out,function(x) x$spec)),ncol=nrow(st@ifr))
            }
            if(method==""){
              ## alternative method to confirm the results obtained by pwelch
            }
            
            return(list(freq=freq,ps=ps))
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

            ## replace ~ by the home directory as c code does not know ~
            if(grepl(pattern = "~",pathSession))
            {
              pathSession<-gsub("~",replacement = Sys.getenv("HOME"),x = pathSession)
            }

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
            cluNoFromCluFile=st@clu[1]
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
            if(st@nCells!=cluNoFromCluFile-1){
              stop(paste(st@session,"
                         There are fewer unique clusters in the spike train than
                         it should based on the first line of the clu file.
                         Could affect how tetrode_id are assigned to cluster if using cellGroup object.
                         Consider reclustering and deleting cluster without spikes."))
            }

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
#' The results are saved in the slot auto and autoTimePoints.
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
            st@autoTimePoints=seq(-st@autoMsPerBin*nBins/2+st@autoMsPerBin/2,
                                  st@autoMsPerBin*nBins/2,
                                  st@autoMsPerBin)
            st@autoProbability=probability
            rownames(st@auto)<-st@autoTimePoints
            return(st)
            }
)

#' Calculate spike-time autocorrelation center of mass
#'
#' Only the positive half of the spike-time autocorrelation is used.
#' Only positive values should be returned by this function.
#'
#' @param st SpikeTrain object
#'
#' @docType methods
#' @rdname spikeTimeAutocorrelationCenterOfMass-methods
setGeneric(name="spikeTimeAutocorrelationCenterOfMass",
           def=function(st)
           {standardGeneric("spikeTimeAutocorrelationCenterOfMass")})
#' @rdname spikeTimeAutocorrelationCenterOfMass-methods
#' @aliases spikeTimeAutocorrelationCenterOfMass,ANY,ANY-method
setMethod(f="spikeTimeAutocorrelationCenterOfMass",
          signature = "SpikeTrain",
          definition=function(st)
          {
            if(length(st@auto)==0)
              stop(paste("length(st@auto) == 0, call spikeTimeAutocorrelation() before spikeTimeAutocorrelationCenterOfMass()"))
            
            ## get center of mass of positive part of st auto
            autoPositive<-st@auto[which(st@autoTimePoints>=0),]
            if(class(autoPositive)=="numeric"|class(autoPositive)=="integer")
              {
              #only one cell
              com<-centerOfMass(autoPositive)
            }else
              {
              com<-apply(autoPositive,2,centerOfMass) # the values are in indices from 1 to length of positive auto
              }
            
            ## transform the com from indices to ms
            m<-min(st@autoTimePoints[which(st@autoTimePoints>=0)])
            M<-max(st@autoTimePoints[which(st@autoTimePoints>=0)])
            propRange<-(com-1)/(length(st@autoTimePoints[which(st@autoTimePoints>=0)])-1) ## proportion of the range
            st@autoCOM<-m+propRange*(M-m)
                        return(st)
          })





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
            st@crossEventsTimePoints=seq(-st@crossEventsMsPerBin*nBins/2+st@crossEventsMsPerBin/2,
                                   st@crossEventsMsPerBin*nBins/2,
                                   st@crossEventsMsPerBin)
            st@crossEventsProbability=probability
            rownames(st@crossEvents)<-st@crossEventsTimePoints
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
            st@crossTimePoints=seq(-st@crossMsPerBin*nBins/2+st@crossMsPerBin/2,
                                   st@crossMsPerBin*nBins/2,
                                   st@crossMsPerBin)
            rownames(st@cross)<-st@crossTimePoints
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



#' Set the list of cell pairs to limit the analysis to these cell pairs
#'
#' Only these cell pairs will be considered for analysis.
#'
#' @param st SpikeTrain object
#' @param cellPairList Data.frame with 2 columns containing the clu id of the neurons, one pair per row
#' @return a SpikeTrain object with a new cellPairList.
#'
#' @docType methods
#' @rdname setCellPairList-methods
setGeneric(name="setCellPairList",
           def=function(st,cellPairList)
           {standardGeneric("setCellPairList")})

#' @rdname setCellPairList-methods
#' @aliases setCellPairList,ANY,ANY-method
setMethod(f="setCellPairList",
          signature = "SpikeTrain",
          definition=function(st,cellPairList)
          {
            if(class(cellPairList)!="data.frame")
              stop("setCellPairList: cellPairList is not a data.frame")
            if(length(cellPairList[,1])==0)
              stop("setCellPairList: length(cellPairList[,1]==0)")
            if(dim(cellPairList)[2]!=2)
              stop("setCellPairList: dim(cellPairList)[2] != 2")
            st@cellPairList<-cellPairList
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
           def=function(st,refractoryMs=1.5,binSizeMs=0.5,windowSizeMs=25,
                        minControlWindowMs=5.0,maxControlWindowMs=25)
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
            print(paste("path:",object@path))
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
              if(length(object@startInterval<100)){
                print(paste(object@startInterval,object@endInterval))
              }else{
                print(head(paste(object@startInterval,object@endInterval),n=100))
                print("only the first 100 are printed")
              }
              print(paste("nIntervals for c function:",length(object@startResIndexc)))
              if(length(object@startInterval<100)){
                print(paste(object@startResIndexc,object@endResIndexc))
              }else{
                print(head(paste(object@startResIndexc,object@endResIndexc),n=100))
                print("only the first 100 are printed")
              }
            }
            if(length(object@events)!=0)
              print(paste("nEvents:",length(object@events)))
                  })


#' Plot the ifr over time of a SpikeTrain object
#'
#' You first need to run ifr(st) before calling this function.
#'
#' @param st SpikeTrain object
#' @param cluId cluId of the cell you want to plot. cluId here is just a number
#' @param name name for the plot
#' @param axis.y.pos position of the y axis
#' @param axis.x.pos position of the x axis
#' @param axis.y.las orientation of numbers on y axis, 1 or 2
#' @param mgp.x mgp of x axis
#' @param mgp.y mgp of y axis
#' @param xlab label of x axis
#' @param ylab label of y axis
#' @param plotxlim limits of x axis
#' @param plotylim limits of y axis
#' @param outma outma for graph
#' @param margin margin of graph
#' @param xaxis.at x axis tic location
#' @param yaxis.at y axis tics location
#' @param add.text text you can add on the graph
#' @param add.text.pos position of the text you can add on the graph
#' 
#' @docType methods
#' @rdname ifrPlot-methods
setGeneric(name="ifrPlot",
           def=function(st,cluId,
                        name="",
                        axis.y.pos=NA,
                        axis.x.pos=NA,
                        axis.y.las=2,
                        mgp.x=c(0.5,0.1,0.1),
                        mgp.y=c(1.1,0.2,0.1),
                        xlab="Time (sec)",
                        ylab="Firing rate (Hz)",
                        plotxlim=NA,
                        plotylim=NA,
                        outma=c(0.5,0.5,0.5,0.5),
                        margin=c(1.5,1.7,1,0.3),
                        xaxis.at=NA,
                        yaxis.at=NA,
                        add.text="",
                        add.text.pos=c(0,0.5))
           {standardGeneric("ifrPlot")})
#' @rdname ifrPlot-methods
#' @aliases ifrPlot,ANY,ANY-method
setMethod(f="ifrPlot",
          signature = "SpikeTrain",
          definition=function(st,cluId,
                              name="",
                              axis.y.pos=NA,
                              axis.x.pos=NA,
                              axis.y.las=2,
                              mgp.x=c(0.5,0.1,0.1),
                              mgp.y=c(1.1,0.2,0.1),
                              xlab="Time (sec)",
                              ylab="Firing rate (Hz)",
                              plotxlim=NA,
                              plotylim=NA,
                              outma=c(0.5,0.5,0.5,0.5),
                              margin=c(1.5,1.7,1,0.3),
                              xaxis.at=NA,
                              yaxis.at=NA,
                              add.text="",
                              add.text.pos=c(0,0.5))
          {
            if(dim(st@ifr)[1]==0)
            {
              stop("you need to call ifr() before calling ifrPlot()")
            }
            
            if(cluId%in%st@cellList==FALSE)
            {
              stop(paste("cluId",cluId, "is not in the SpikeTrain object"))
            }
            ## get the ifr for this CluId
            m<-st@ifr[which(st@cellList==cluId),]
            
            par(mar=margin, oma=outma,cex.lab=0.6,cex.axis=0.6)
            if(any(is.na(plotxlim)))
              plotxlim=c(min(0),max(st@ifrTime))
            if(any(is.na(plotylim)))
              plotylim=c(0,max(m))
            if(any(is.na(axis.y.pos)))
              axis.y.pos<-min(st@ifrTime)
            if(any(is.na(axis.x.pos)))
              axis.x.pos<-0
            graphics::plot(x=plotxlim,y=plotylim,type='n', axes=FALSE, pch=20,lwd=1,xlab="",ylab="")
            
            #keep only data within xlim
            
            m<-m[which(st@ifrTime>=plotxlim[1]&st@ifrTime<=plotxlim[2])]
            t<-st@ifrTime[which(st@ifrTime>=plotxlim[1]&st@ifrTime<=plotxlim[2])]
            
            lines(t,m)
            par(mgp=mgp.x)
            if(is.na(xaxis.at)){
              graphics::axis(side = 1, pos=axis.x.pos, tck=-0.05, cex.axis=0.6)
            }else{
              graphics::axis(side = 1, pos=axis.x.pos, at=xaxis.at, tck=-0.05, cex.axis=0.6)
            }
            par(mgp=mgp.y)
            if(is.na(yaxis.at)){
              graphics::axis(side = 2, las=axis.y.las, pos=axis.y.pos,tck=-0.05,cex.axis=0.6)
            } else{
              graphics::axis(side = 2, at=yaxis.at, las=axis.y.las, pos=axis.y.pos,tck=-0.05,cex.axis=0.6)
            }
            graphics::title(xlab=xlab,mgp=mgp.x)
            graphics::title(ylab=ylab,mgp=mgp.y)
            if(name!=""){
              graphics::title(main=name,cex.main=0.5)
            }
            if(add.text!=""){
              graphics::text(labels=add.text,x=add.text.pos[1],y=add.text.pos[2],cex=0.6)
            }
            return()
          }
)



#' Plot the ifr vector over time
#'
#'
#' @param ifr ifr vector
#' @param timePoints vector with the time points for ifr values
#' @param name name for the plot
#' @param axis.y.pos position of the y axis
#' @param axis.x.pos position of the x axis
#' @param axis.y.las orientation of numbers on y axis, 1 or 2
#' @param mgp.x mgp of x axis
#' @param mgp.y mgp of y axis
#' @param xlab label of x axis
#' @param ylab label of y axis
#' @param plotxlim limits of x axis
#' @param plotylim limits of y axis
#' @param outma outma for graph
#' @param margin margin of graph
#' @param xaxis.at x axis tic location
#' @param yaxis.at y axis tics location
#' @param add.text text you can add on the graph
#' @param add.text.pos position of the text you can add on the graph
ifrPlotVec<-function(ifr,timePoints,
                  name="",
                  axis.y.pos=NA,
                  axis.x.pos=NA,
                  axis.y.las=2,
                  mgp.x=c(0.5,0.1,0.1),
                  mgp.y=c(1.1,0.2,0.1),
                  xlab="Time (sec)",
                  ylab="Firing rate (Hz)",
                  plotxlim=NA,
                  plotylim=NA,
                  outma=c(0.5,0.5,0.5,0.5),
                  margin=c(1.5,1.7,1,0.3),
                  xaxis.at=NA,
                  yaxis.at=NA,
                  add.text="",
                  add.text.pos=c(0,0.5))
{
  if(length(ifr)==0)
  {
    stop("ifr vector is empty")
  }
  if(length(timePoints)==0)
  {
    stop("timePoints vector is empty")
  }
  if(length(timePoints)!=length(ifr))
  {
    stop("length of timePoints differs from ifr")
  }
  
  m<-ifr
  
  par(mar=margin, oma=outma,cex.lab=0.6,cex.axis=0.6)
  if(any(is.na(plotxlim)))
    plotxlim=c(min(0),max(timePoints))
  if(any(is.na(plotylim)))
    plotylim=c(0,max(m))
  if(any(is.na(axis.y.pos)))
    axis.y.pos<-min(timePoints)
  if(any(is.na(axis.x.pos)))
    axis.x.pos<-0
  graphics::plot(x=plotxlim,y=plotylim,type='n', axes=FALSE, pch=20,lwd=1,xlab="",ylab="")
  
  #keep only data within xlim
  
  m<-m[which(timePoints>=plotxlim[1]&timePoints<=plotxlim[2])]
  t<-timePoints[which(timePoints>=plotxlim[1]&timePoints<=plotxlim[2])]
  
  lines(t,m)
  par(mgp=mgp.x)
  if(is.na(xaxis.at)){
    graphics::axis(side = 1, pos=axis.x.pos, tck=-0.05, cex.axis=0.6)
  }else{
    graphics::axis(side = 1, pos=axis.x.pos, at=xaxis.at, tck=-0.05, cex.axis=0.6)
  }
  par(mgp=mgp.y)
  if(is.na(yaxis.at)){
    graphics::axis(side = 2, las=axis.y.las, pos=axis.y.pos,tck=-0.05,cex.axis=0.6)
  } else{
    graphics::axis(side = 2, at=yaxis.at, las=axis.y.las, pos=axis.y.pos,tck=-0.05,cex.axis=0.6)
  }
  graphics::title(xlab=xlab,mgp=mgp.x)
  graphics::title(ylab=ylab,mgp=mgp.y)
  if(name!=""){
    graphics::title(main=name,cex.main=0.5)
  }
  if(add.text!=""){
    graphics::text(labels=add.text,x=add.text.pos[1],y=add.text.pos[2],cex=0.6)
  }
  return()
}

