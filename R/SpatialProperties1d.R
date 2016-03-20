#################################################
#### definition of SpatialProperties2d Class  ###
#################################################
SpatialProperties1d<- setClass(
  "SpatialProperties1d", ## name of the class
  slots=c(session="character",
          cmPerBin="numeric",
          smoothOccupancySd="numeric", ## in cm
          smoothRateHistoSd="numeric", ## in cm
          nBinRateHisto="integer",
          xSpikes="numeric",
          rateHisto="array",
          occupancy="numeric",
          cellList="numeric",
          reduceSize="logical",
          ##
          peakRate="numeric",
          infoScore="numeric",
          sparsity="numeric",
          ##
          nShufflings="numeric",
          minShiftMs="numeric",
          peakRateShuffle="numeric",
          infoScoreShuffle="numeric",
          sparsityShuffle="numeric"
          ),
  
  prototype = list(session="",cmPerBin=2,smoothOccupancySd=3,smoothRateHistoSd=3,reduceSize=T,
                   nShufflings=100,minShiftMs=20000))



### show ###
setMethod("show", "SpatialProperties1d",
          function(object){
            print(paste("session:",object@session))
            print(paste("cmPerBin:",object@cmPerBin))
            print(paste("smoothOccupancySd:",object@smoothOccupancySd))
            print(paste("smoothRateHistoSd:",object@smoothRateHistoSd))
            print(paste("reduceSize:",object@reduceSize))
            print(paste("nBinRateHisto:",object@nBinRateHisto))
            if(length(object@cellList)!=0){
              print(paste("cellList:"))
              print(object@cellList)
            }
            if(length(object@peakRate)!=0){
              print("peakRate:")
              print(paste(object@peakRate))
              print("infoScore:")
              print(paste(object@infoScore))
              print("sparsity:")
              print(paste(object@sparsity))
            }
            print(paste("nShufflings:",object@nShufflings))
            print(paste("shuffled values:",length(object@infoScoreShuffle)))
          })


### make firing rate histo
setGeneric(name="firingRateHisto",
           def=function(sp1,st,pt)
           {standardGeneric("firingRateHisto")}
)
setMethod(f="firingRateHisto",
          signature="SpatialProperties1d",
          definition=function(sp1,st,pt)
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@x)==0)
              stop("pt@x has length of 0")
            if(st@session=="")
              stop("st@session is empty")
            if(st@nSpikes==0)
              stop("st@nSpikes==0")
            if(length(pt@lin)==0)
              stop("pt@lin has length of 0")

            ## get list of cells from SpikeTrain object
            sp1@cellList<-st@cellList
            ## reduce the size of histo
            if(sp1@reduceSize==T){
              x<-pt@lin-min(pt@lin,na.rm=T)+sp1@cmPerBin
            }else{
              x<-pt@x
            }
            sp1@nBinRateHisto<-as.integer(ceiling(max(x,na.rm=T)/sp1@cmPerBin))
            ## use -1 as invalid values in c functions
            x[is.na(x)]<- -1.0
            
            ## get spike head direction
            sp1@xSpikes<-.Call("spike_position_1d_cwrap",
                               x,
                               length(x),
                               as.integer(st@res),
                               as.integer(st@nSpikes),
                               as.integer(pt@resSamplesPerWhlSample),
                               as.integer(st@startInterval),
                               as.integer(st@endInterval),
                               length(st@startInterval))
            ## make the occupancy map
            sp1@occupancy<-.Call("occupancy_histogram_cwrap",
                                sp1@nBinRateHisto,
                                sp1@cmPerBin,
                                x,
                                length(x),
                                pt@resSamplesPerWhlSample/pt@samplingRateDat*1000, ## ms per whl samples
                                as.integer(st@startInterval),
                                as.integer(st@endInterval),
                                length(st@startInterval),
                                as.integer(pt@resSamplesPerWhlSample),
                                1)
            ## smooth the occupancy map
            sp1@occupancy<- .Call("smooth_double_gaussian_circular_cwrap",
                                 as.numeric(sp1@occupancy),
                                 sp1@nBinRateHisto,
                                 sp1@smoothOccupancySd/sp1@cmPerBin,
                                 -1.0)
            ## make the 1d maps
            results<- .Call("firing_rate_histo_cwrap",
                            sp1@nBinRateHisto,
                            sp1@cmPerBin,
                            sp1@xSpikes,
                            as.integer(st@clu),
                            as.integer(st@nSpikes),
                            as.integer(sp1@cellList),
                            length(sp1@cellList),
                            as.numeric(sp1@occupancy),
                            sp1@smoothRateHistoSd/sp1@cmPerBin,
                            1)
            sp1@rateHisto<-array(data=results,dim=(c(sp1@nBinRateHisto,length(sp1@cellList))))
            return(sp1)
          }
)


