#################################################
#### definition of HeadDirection Class  ###
#################################################
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
          ## stats
          peakRates="numeric",
          vectorLength="numeric",
          meanDirection="numeric",
          ##
          nShufflings="numeric",
          minShiftMs="numeric",
          peakRatesShuffle="numeric",
          vectorLengthShuffle="numeric"
  ),
  prototype = list(session="",degPerBin=10,smoothOccupancySd=10,smoothRateHistoSd=10,nShufflings=100,minShiftMs=20000,histoRepetitions=0))


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
            }
          })

### make firing rate histo
setGeneric(name="headDirectionHisto",
           def=function(hd,st,pt)
           {standardGeneric("headDirectionHisto")}
)
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
            
            ## smooth the occupancy map
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
                            hd@histoRepetitions+1)
            hd@histo<-array(data=results,dim=(c(hd@nBinHisto,length(hd@cellList))))
            
            return(hd)
          }
)


#### getHistoStats
setGeneric(name="getHistoStats",
           def=function(hd,st,pt)
           {standardGeneric("getHistoStats")}
)
setMethod(f="getHistoStats",
          signature="HeadDirection",
          definition=function(hd,st,pt)
          {
            ### create the histo
            hd<-headDirectionHisto(hd,st,ptsqr70)
            ### get peak rates
            hd@peakRates<-apply(hd@histo,2,max)
            ### get vector length and mean direction
            results<- .Call("circular_stats_rate_histogram_cwrap",
                                 as.integer(hd@cellList),
                                 length(hd@cellList),
                                 as.numeric(hd@histo),
                                 as.integer(hd@nBinHisto))
            hd@vectorLength<-results[1,]
            hd@meanDirection<-results[2,]
            return(hd)
          }
)


#### getHistoStatsShuffle 
setGeneric(name="getHistoStatsShuffle",
           def=function(hd,st,pt)
           {standardGeneric("getHistoStatsShuffle")}
)
setMethod(f="getHistoStatsShuffle",
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
          
            
            for(i in 1:sp@nShufflings){
              print(paste(i,"of",sp@nShufflings))
              
              pts<-shiftHdRandom(pt)
              
              ### create the histo
              hd<-headDirectionHisto(hd,st,pts)
              
              ### get peak rates
              hd@peakRatesShuffle<-c(hd@peakRatesShuffle,apply(hd@histo,2,max))
              
              ### get vector length and mean direction
              results<- .Call("circular_stats_rate_histogram_cwrap",
                              as.integer(hd@cellList),
                              length(hd@cellList),
                              as.numeric(hd@histo),
                              as.integer(hd@nBinHisto))
              hd@vectorLengthShuffle <-c(hd@vectorLengthShuffle,results[1,])
            }       
            return(hd)
          })



