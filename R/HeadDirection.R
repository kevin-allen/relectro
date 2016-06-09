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
#' @slot histoRepetitions Number of repetition of the 0-360 range in the histogram. Default is 0 repetition
#' @slot peakRates Peak firing rate in Hz in each firing rate histogram
#' @slot vectorLength Mean vector length of each firing rate histogram
#' @slot meanDirection Mean direction in each firing rate histogram
#' @slot nShufflings Number of shufflings to get a distribution of vector length that would be obtained by chance
#' @slot minShiftMs Minimum time shift of the head direction data used during the shuffling procedure
#' @slot peakRatesShuffle Numeric to hold the peak firing rates obtained during the shuffling procedure
#' @slot vectorLengthShuffle Numeric to hold the vector length obtained during the shuffling procedure
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
          histoRepetitions="numeric", # 0 is only once
          peakRates="numeric",
          vectorLength="numeric",
          meanDirection="numeric",
          nShufflings="numeric",
          minShiftMs="numeric",
          peakRatesShuffle="numeric",
          vectorLengthShuffle="numeric"
  ),
  prototype = list(session="",degPerBin=10,smoothOccupancySd=10,smoothRateHistoSd=10,nShufflings=100,
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
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@x)==0)
              stop("pt@x has length of 0")
            if(st@session=="")
              stop("st@session is empty")
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

#' Calculate the head direction statistics for the head direction rate histograms with a HeadDirection, SpikeTrain and Positrack objects
#'
#' @param hd HeadDirection object
#' @param st SpikeTrain object
#' @param pt Positrack
#' @return HeadDirection object with the statistics in the slots peakRates, vectorLength and meanDirection
#' 
#' @docType methods
#' @rdname headDirectionStats-methods
setGeneric(name="headDirectionStats",
           def=function(hd,st,pt)
           {standardGeneric("headDirectionStats")}
)
#' @rdname headDirectionStats-methods
#' @aliases headDirectionStats,ANY,ANY-method
setMethod(f="headDirectionStats",
          signature="HeadDirection",
          definition=function(hd,st,pt)
          {
            ### create the histo
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
