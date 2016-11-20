#' An S4 class used to get head direction histograms and circular statistics for each neuron.
#'
#' The user creates positrack and spike train objects first.
#' Then the user calls the different methods of this class with the spike train and positrack objects as arguments.
#' It can get the vector length of the head direction histograms and do shuffling to know the threshold for significance.
#' 
#' @slot session A charactor object with the name of the recording session
#' @slot degPerBin Degrees per bin in the head direction histograms
#' @slot smoothOccupancySd The standard deviation of the gaussian kernel used to smooth the occupancy histogram
#' @slot smoothRateHistoSd The standard deviation of the gaussian kernel used to smooth the rate histograms
#' @slot nBinHisto Number of bins in the histograms
#' @slot hdSpikes Head direction of all valid spikes
#' @slot histo An array holding the histograms of all neurons
#' @slot occupancy A numeric holding the occupancy histogram
#' @slot cellList A list of cluster id of the neurons
#' @slot cellPairList Data frame containing pairs of cells
#' @slot histoRepetitions Number of repetition of the 0-360 range in the histogram. Default is 0 repetition
#' @slot peakRates Peak firing rate in Hz in each firing rate histogram
#' @slot vectorLength Mean vector length of each firing rate histogram
#' @slot meanDirection Mean direction in each firing rate histogram
#' @slot nShufflings Number of shufflings to get a distribution of vector length that would be obtained by chance
#' @slot minShiftMs Minimum time shift of the head direction data used during the shuffling procedure
#' @slot peakRatesShuffle Numeric to hold the peak firing rates obtained during the shuffling procedure
#' @slot vectorLengthShuffle Numeric to hold the vector length obtained during the shuffling procedure
#' @slot crossHisto Array containing crosscorrelation rate histogram
HeadDirection<- setClass(
  "HeadDirection", ## name of the class
  slots=c(session="character",
          degPerBin="numeric",
          smoothOccupancySd="numeric", ## in deg
          smoothRateHistoSd="numeric", ## in deg
          nBinHisto="integer",
          hdSpikes="numeric",
          histo="array",# 2d
          occupancy="numeric",
          cellList="numeric",
          cellPairList="data.frame",       
          histoRepetitions="numeric", # 0 is only once
          peakRates="numeric",
          vectorLength="numeric",
          meanDirection="numeric",
          nShufflings="numeric",
          minShiftMs="numeric",
          peakRatesShuffle="numeric",
          vectorLengthShuffle="numeric",
          crossHisto="array"
  ),
  prototype = list(session="",degPerBin=10,nBinHisto=as.integer(36),smoothOccupancySd=10,smoothRateHistoSd=10,nShufflings=100,
                   minShiftMs=20000,histoRepetitions=0))


### show ###
setMethod("show", "HeadDirection",
          function(object){
            print(paste("session:",object@session))
            print(paste("degPerBin:",object@degPerBin))
            print(paste("smoothOccupancySd:",object@smoothOccupancySd))
            print(paste("smoothRateHistoSd:",object@smoothRateHistoSd))
            print(paste("nBinHisto:",object@nBinHisto))
            if(length(object@cellList)!=0){
              print(paste("cellList:"))
              print(object@cellList)
            }
            print(paste("nShufflings:",object@nShufflings))
            print(paste("shuffled values:",length(object@vectorLengthShuffle)))
            if(length(object@meanDirection)!=0){
              print("meanDirction:")
              print(paste(round(object@meanDirection,1)))
              print("vectorLength:")
              print(paste(round(object@vectorLength,3)))
              print("peakRates:")
              print(paste(round(object@peakRates,3)))
            }
          })


#' Calculate the head direction rate histograms of neurons using a HeadDirection, SpikeTrain and Positrack objects
#'
#' @param hd HeadDirection object
#' @param st SpikeTrain object
#' @param pt Positrack
#' @return HeadDirection object with the firing rate histograms in the slot histo
#' 
#' @docType methods
#' @rdname headDirectionHisto-methods
setGeneric(name="headDirectionHisto",
           def=function(hd,st,pt)
           {standardGeneric("headDirectionHisto")}
)
#' @rdname headDirectionHisto-methods
#' @aliases headDirectionHisto,ANY,ANY-method
setMethod(f="headDirectionHisto",
          signature="HeadDirection",
          definition=function(hd,st,pt)
          {
            if(length(pt@hd)==0)
              stop("pt@hd is empty")
            if(length(pt@x)==0)
              stop("pt@x has length of 0")
            if(st@nSpikes==0)
              stop("st@nSpikes==0")
            
            hd@cellList<-st@cellList
            
            ## use -1 as invalid values in c functions
            hdir<-pt@hd
            hdir[is.na(hdir)]<- -1.0
            hd@nBinHisto<-as.integer(ceiling(360/hd@degPerBin))
            
            ## get spike head direction
            hd@hdSpikes<-.Call("spike_head_direction_cwrap",
                           hdir,
                           length(hdir),
                           as.integer(st@res),
                           as.integer(st@nSpikes),
                           as.integer(pt@resSamplesPerWhlSample),
                           as.integer(st@startInterval),
                           as.integer(st@endInterval),
                           length(st@startInterval))
            
            ## make the occupancy map
            hd@occupancy<-.Call("occupancy_histogram_cwrap",
                                hd@nBinHisto,
                                hd@degPerBin,
                                hdir,
                                length(hdir),
                                pt@resSamplesPerWhlSample/pt@samplingRateDat*1000, ## ms per whl samples
                                as.integer(st@startInterval),
                                as.integer(st@endInterval),
                                length(st@startInterval),
                                as.integer(pt@resSamplesPerWhlSample),
                                hd@histoRepetitions+1)
            
            hd@occupancy<- .Call("smooth_double_gaussian_circular_cwrap",
                                 as.numeric(hd@occupancy),
                                 hd@nBinHisto,
                                 hd@smoothOccupancySd/hd@degPerBin,
                                 -1.0)
            
            ## make the 2d maps
            results<- .Call("firing_rate_histo_cwrap",
                            hd@nBinHisto,
                            hd@degPerBin,
                            hd@hdSpikes,
                            as.integer(st@clu),
                            as.integer(st@nSpikes),
                            as.integer(hd@cellList),
                            length(hd@cellList),
                            as.numeric(hd@occupancy),
                            hd@smoothRateHistoSd/hd@degPerBin,
                            hd@histoRepetitions+1,
                            1)
            hd@histo<-array(data=results,dim=(c(hd@nBinHisto,length(hd@cellList))))
            
            return(hd)
          }
)

#' Calculate the spike-triggered head-direction histogram using a HeadDirection, SpikeTrain and Positrack objects
#'
#' Each spike is treated as a reference spike in turn, and set to phase 180 degree. 
#' The histogram is constructed from the data following the reference spikes 
#' by shifting the head direction data so that the head direction of the 
#' agent at the time of the reference spike is 180 degree.
#' 
#' The occupancy histogram and the firing rate histogram are smoothed with a Gaussian kernel.
#' The amount of smoothing is determined by slots smoothOccupancySd and smoothRateHistoSd of the HeadDirection object.
#' 
#' You can set the temporal limit for the data used to construct the histogram with minIsiMs and maxIsiMs
#' 
#' 
#' @param hd HeadDirection object
#' @param st SpikeTrain object
#' @param pt Positrack object
#' @param minIsiMs Minimal interspike interval to consider in ms
#' @param maxIsiMs Maximal interspike interval to consider in ms
#' @return HeadDirection object with the spike-triggered rate histograms in the hist slot
#' 
#' @docType methods
#' @rdname spikeTriggeredHeadDirectionHisto-methods
setGeneric(name="spikeTriggeredHeadDirectionHisto",
           def=function(hd,st,pt,minIsiMs,maxIsiMs)
           {standardGeneric("spikeTriggeredHeadDirectionHisto")}
)
#' @rdname spikeTriggeredHeadDirectionHisto-methods
#' @aliases spikeTriggeredHeadDirectionHisto,ANY,ANY-method
setMethod(f="spikeTriggeredHeadDirectionHisto",
          signature="HeadDirection",
          definition=function(hd,st,pt,minIsiMs,maxIsiMs)
          {
            if(length(pt@x)==0)
              stop(paste("pt@x has length of 0 in firingRateMap2d",st@session))
            if(st@nSpikes==0)
              stop(paste("st@nSpikes==0 in firingRateMap2d",st@session))
            if(minIsiMs<0)
              stop(paste("minIsiMs should be 0 or larger than 0"))
            if(maxIsiMs<0)
              stop(paste("maxIsiMs should larger than 0"))
            if(maxIsiMs<=minIsiMs)
              stop(paste("maxIsiMs should be larger than minIsiMs"))
            hd@cellList<-st@cellList
            ## use -1 as invalid values in c functions
            hdir<-pt@hd
            hdir[is.na(hdir)]<- -1.0
            hd@nBinHisto<-as.integer(ceiling(360/hd@degPerBin))
            results<- .Call("spike_triggered_head_direction_histo_cwrap",
                            as.integer(hd@nBinHisto),
                            hd@degPerBin,
                            as.integer(hd@cellList),
                            length(hd@cellList),
                            hdir,
                            length(pt@hd),
                            as.integer(st@res),
                            as.integer(st@clu),
                            st@nSpikes,
                            as.integer(st@startInterval),
                            as.integer(st@endInterval),
                            length(st@startInterval),
                            pt@resSamplesPerWhlSample/pt@samplingRateDat*1000,
                            as.integer(pt@resSamplesPerWhlSample),
                            hd@smoothOccupancySd,
                            hd@smoothRateHistoSd,
                            minIsiMs,
                            maxIsiMs,
                            as.integer(st@samplingRate))
            hd@histo<-array(data=results,dim=(c(hd@nBinHisto,length(hd@cellList))))
            return(hd)
          }
)

#' Calculate the head direction statistics for the head direction rate histograms with a HeadDirection, SpikeTrain and Positrack objects
#'
#' @param hd HeadDirection object
#' @param st SpikeTrain object
#' @param pt Positrack
#' @param newHisto logical indicating if new rate histograms are calculated when calling the function.
#' Default value is TRUE
#' @return HeadDirection object with the statistics in the slots peakRates, vectorLength and meanDirection
#' 
#' @docType methods
#' @rdname headDirectionStats-methods
setGeneric(name="headDirectionStats",
           def=function(hd,st,pt,newHisto=TRUE)
           {standardGeneric("headDirectionStats")}
)
#' @rdname headDirectionStats-methods
#' @aliases headDirectionStats,ANY,ANY-method
setMethod(f="headDirectionStats",
          signature="HeadDirection",
          definition=function(hd,st,pt,newHisto=TRUE)
          {
            ### create the histo
            if(newHisto==TRUE)
              hd<-headDirectionHisto(hd,st,pt)
            
            ### get peak rates
            hd@peakRates<-apply(hd@histo,2,max)
            ### get vector length and mean direction
            results<- .Call("circular_stats_rate_histogram_cwrap",
                                 length(hd@cellList),
                                 as.numeric(hd@histo),
                                 as.integer(hd@nBinHisto))
            hd@vectorLength<-results[1,]
            hd@meanDirection<-results[2,]
            return(hd)
          }
)

#' Calculate the random head direction statistics for the head direction rate histograms with 
#' a HeadDirection, SpikeTrain and Positrack objects
#'
#' @param hd HeadDirection object
#' @param st SpikeTrain object
#' @param pt Positrack
#' @return HeadDirection object with the random statistics in the slots 
#' peakRatesShuffle, vectorLengthShuffle and meanDirectionShuffle
#' 
#' @docType methods
#' @rdname headDirectionStatsShuffle-methods
setGeneric(name="headDirectionStatsShuffle",
           def=function(hd,st,pt)
           {standardGeneric("headDirectionStatsShuffle")}
)
#' @rdname headDirectionStatsShuffle-methods
#' @aliases headDirectionStatsShuffle,ANY,ANY-method
setMethod(f="headDirectionStatsShuffle",
          signature="HeadDirection",
          definition=function(hd,st,pt){
            
            if(hd@nShufflings==0)            
              stop("sp@nShufflings==0")
            
            if(hd@session=="")
              stop("pt@session is empty")
            if(length(pt@x)==0)
              stop("pt@x has length of 0")
            if(st@nSpikes==0)
              stop("st@nSpikes==0")
            if(length(hd@vectorLengthShuffle)!=0){
              hd@vectorLengthShuffle=vector("numeric")
              hd@peakRatesShuffle=vector("numeric")
            }
            for(i in 1:hd@nShufflings){

              pts<-shiftHdRandom(pt)
              
              ### create the histo
              hd<-headDirectionHisto(hd,st,pts)
              
              ### get peak rates
              hd@peakRatesShuffle<-c(hd@peakRatesShuffle,apply(hd@histo,2,max))
              
              ### get vector length and mean direction
              results<- .Call("circular_stats_rate_histogram_cwrap",
                              length(hd@cellList),
                              as.numeric(hd@histo),
                              as.integer(hd@nBinHisto))
              hd@vectorLengthShuffle <-c(hd@vectorLengthShuffle,results[1,])
            }       
            return(hd)
          })

#' Return the hd histos as a data.frame
#' 
#' 
#' @param hd HeadDirection object
#' @return data.frame with the rate-head direction histograms. Column names are clu.id, deg, rate
#' 
#' @docType methods
#' @rdname hdHistoAsDataFrame-methods
setGeneric(name="hdHistoAsDataFrame",
           def=function(hd)
           {standardGeneric("hdHistoAsDataFrame")}
)
#' @rdname hdHistoAsDataFrame-methods
#' @aliases hdHistoAsDataFrame,ANY,ANY-method
setMethod(f="hdHistoAsDataFrame",
          signature="HeadDirection",
          definition=function(hd)
          {
            if(length(hd@histo)==0)
              stop("Need to call headDirectionHisto first to run hdHistoAsDataFrame()")
         
            data.frame(clu.id=rep(paste(hd@session,hd@cellList,sep="_"),each=hd@nBinHisto),
                       deg=rep(seq(hd@degPerBin/2,360*(hd@histoRepetitions+1),hd@degPerBin),length(hd@cellList)),
                       rate=as.numeric(hd@histo))
          }
)



#' Calculate the spike-triggered cross head-direction histogram using a HeadDirection, SpikeTrain and Positrack objects
#'
#' This works on pairs of cells listed in slot cellPairList of the st object.
#' Each spike of cell A is treated as a reference spike in turn. 
#' The histogram is constructed from the data following the reference spikes 
#' by shifting the head direction data so that the head direction of the 
#' agent at the time of the reference spike is 180 degree.
#' The spikes of cell B are used to calculate the sampling rates
#' The results are stored in the crossHisto slot
#' 
#' The occupancy histogram and the firing rate histogram are smoothed with a Gaussian kernel.
#' The amount of smoothing is determined by slots smoothOccupancySd and smoothRateHistoSd of the HeadDirection object.
#' 
#' You can set the temporal limit for the data used to construct the histogram with minIsiMs and maxIsiMs
#' 
#' @param hd HeadDirection object
#' @param st SpikeTrain object
#' @param pt Positrack object
#' @param minIsiMs Minimal interspike interval to consider in ms
#' @param maxIsiMs Maximal interspike interval to consider in ms
#' @return HeadDirection object with the spike-triggered cross rate histograms
#' 
#' @docType methods
#' @rdname spikeTriggeredCrossHeadDirectionHisto-methods
setGeneric(name="spikeTriggeredCrossHeadDirectionHisto",
           def=function(hd,st,pt,minIsiMs,maxIsiMs)
           {standardGeneric("spikeTriggeredCrossHeadDirectionHisto")}
)
#' @rdname spikeTriggeredCrossHeadDirectionHisto-methods
#' @aliases spikeTriggeredCrossHeadDirectionHisto,ANY,ANY-method
setMethod(f="spikeTriggeredCrossHeadDirectionHisto",
          signature="HeadDirection",
          definition=function(hd,st,pt,minIsiMs,maxIsiMs)
          {
            if(length(pt@x)==0)
              stop(paste("pt@x has length of 0 in firingRateMap2d",st@session))
            if(st@nSpikes==0)
              stop(paste("st@nSpikes==0 in firingRateMap2d",st@session))
            if(minIsiMs<0)
              stop(paste("minIsiMs should be 0 or larger than 0"))
            if(maxIsiMs<0)
              stop(paste("maxIsiMs should larger than 0"))
            if(maxIsiMs<=minIsiMs)
              stop(paste("maxIsiMs should be larger than minIsiMs"))
            if(sum(dim(st@cellPairList))==0)
              stop(paste("st@cellPairList has 0 dimension. set the cellPairList in the st object"))
            hd@cellPairList<-st@cellPairList
            ## use -1 as invalid values in c functions
            hdir<-pt@hd
            hdir[is.na(hdir)]<- -1.0
            hd@nBinHisto<-as.integer(ceiling(360/hd@degPerBin))
            results<- .Call("spike_triggered_head_direction_histo_cwrap",
                            as.integer(hd@nBinHisto),
                            hd@degPerBin,
                            as.integer(hd@cellList),
                            length(hd@cellList),
                            hdir,
                            length(pt@hd),
                            as.integer(st@res),
                            as.integer(st@clu),
                            st@nSpikes,
                            as.integer(st@startInterval),
                            as.integer(st@endInterval),
                            length(st@startInterval),
                            pt@resSamplesPerWhlSample/pt@samplingRateDat*1000,
                            as.integer(pt@resSamplesPerWhlSample),
                            hd@smoothOccupancySd,
                            hd@smoothRateHistoSd,
                            minIsiMs,
                            maxIsiMs,
                            as.integer(st@samplingRate))
            hd@crossHisto<-array(data=results,dim=(c(hd@nBinHisto,length(hd@cellList))))
            return(hd)
          }
)
