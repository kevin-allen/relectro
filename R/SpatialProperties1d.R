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
          rateHisto="numeric",
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
