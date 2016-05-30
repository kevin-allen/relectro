#' A S4 class to analyze spatial properties in 2D environments
#' 
#' Use to get firing rate maps of neurons in an open field.
#' It can get information scores, sparsity, etc.
#' You can use it to analyze the firing properties of grid cells or border cells.
#' Can do some shuffling to get the spatial properties that you would get by chance.
#' 
#' @slot session Name of the recording session
#' @slot cmPerBin Number of cm in each bin of a firing rate map, 2 by default
#' @slot smoothOccupancySd Standard deviation in cm of the gaussian kernel used to smooth the occupancy map, 3 by default
#' @slot smoothRateMapSd Standard deviation in cm of the gaussian kernel used to smooth the firing rate histograms, 3 by default
#' @slot nColMap Number of columns in the firing rate maps
#' @slot nRowMap Number of rows in the firing rate maps
#' @slot nColAuto Number of columns in the spatial autocorrelation
#' @slot nRowAuto Number of rows in the spatial autocorrelation
#' @slot xSpikes x position of the animal for each spike time
#' @slot ySpikes y position of the animal for each spike time
#' @slot maps Array containing the firing rate maps of the neurons
#' @slot occupancy Matrix containing the occupancy map
#' @slot autos Array containing the spatial autocorrelation of the firing rate maps
#' @slot autosDetect Array with the spatial autocorrelation once the firing fields have been removed
#' @slot autosDoughnut Array with the spatial autocorrelation containing only the ring of the 6 closest fields (excluding central field)
#' @slot autosDoughnutRotate Array with the rotated spatial autocorrelations
#' @slot cellList List of cells
#' @slot reduceSize Logical indicating if the size of map should be minimized by setting smallest x and y to 0
#' @slot minValidBinsAuto Number of valid bins to calculate a r value in the spatial autocorrelation maps
#' @slot peakRate Peak firing rates in the maps
#' @slot infoScore Information score of the firing rate maps
#' @slot sparsity Sparsity of the firing rate maps
#' @slot borderScore Border score for each firing rate maps
#' @slot borderCM Value of CM when calculating the border score
#' @slot borderDM Value of DM when calculating the border score
#' @slot borderNumFieldsDetected Number of detected fields when calculating the border score
#' @slot borderPercentageThresholdField Threshold for a bin being part of a field, expressed as percentage of peak
#' @slot borderMinBinsInField Minimum number of bins to be considered a field
#' @slot mapPolarity How polarized is the firing rate of the cell in the map. Each bin of the map is given an angle and firing rate.
#' mapPolarity is the resultant mean vector length
#' @slot gridScore Grid score calculated form the spatial autocorrelation
#' @slot gridOrientation Grid orientation
#' @slot gridSpacing Grid spacing in cm
#' @slot gridScoreNumberFieldsToDetect Maximal number of fields that will be detected
#' @slot gridScoreMinNumBinsPerField Minimum number of bins to be considered a field
#' @slot gridScoreFieldThreshold Minimum value of firing fields in the spatial autocorrelation map
#' @slot nAutoRotations Number of rotations when rotating the spatial autocorrelation map
#' @slot AutoRotationDegree Angle of rotation of the spatial autocorrelation map
#' @slot speedScore Speed score of the neurons
#' @slot speedRateSlope Slope of the linear regression line between speed and rate
#' @slot speedRateIntercept Rate intercept of the linear regression line between speed and rate
#' @slot nShufflings Number of shufflings to get chance levels, by default 100
#' @slot minShiftMs Minimum time shift of the position data when doing shuffling
#' @slot peakRateShuffle peak firing rate in the shuffling analysis
#' @slot infoScoreShuffle Information score from the shuffling analysis
#' @slot sparsityShuffle Sparsity from the shuffling analysis
#' @slot borderScoreShuffle Border score form the shuffling analysis
#' @slot borderCMShuffle Border CM from the shuffling analysis
#' @slot borderDMShuffle Border DM from the shuffling analysis
#' @slot mapPolarityShuffle Map polarity from the shuffling analysis
#' @slot gridScoreShuffle Grid score from the shuffling analysis
#' @slot speedScoreShuffle Speed score from the shuffling analysis
SpatialProperties2d<- setClass(
  "SpatialProperties2d", ## name of the class
    slots=c(session="character",
            cmPerBin="numeric",
            smoothOccupancySd="numeric", ## in cm
            smoothRateMapSd="numeric", ## in cm
            nColMap="integer",
            nRowMap="integer",
            nColAuto="integer",
            nRowAuto="integer",
            xSpikes="numeric",
            ySpikes="numeric",
            maps="array",
            occupancy="matrix",
            autos="array",
            autosDetect="array",
            autosDoughnut="array",
            autosDoughnutRotate="array",
            cellList="numeric",
            reduceSize="logical",
            minValidBinsAuto="numeric",
            ##
            peakRate="numeric",
            infoScore="numeric",
            sparsity="numeric",
            ##
            borderScore="numeric",
            borderCM="numeric",
            borderDM="numeric",
            borderNumFieldsDetected="numeric",
            borderPercentageThresholdField="numeric",
            borderMinBinsInField="numeric",
            mapPolarity="numeric",
            ##
            gridScore="numeric",
            gridOrientation="numeric",
            gridSpacing="numeric",
            gridScoreNumberFieldsToDetect="numeric",
            gridScoreMinNumBinsPerField="numeric",
            gridScoreFieldThreshold="numeric",
            nAutoRotations="numeric",
            AutoRotationDegree="numeric",
            ##
            speedScore="numeric",
            speedRateSlope="numeric",
            speedRateIntercept="numeric",
            ##
            nShufflings="numeric",
            minShiftMs="numeric",
            peakRateShuffle="numeric",
            infoScoreShuffle="numeric",
            sparsityShuffle="numeric",
            borderScoreShuffle="numeric",
            borderCMShuffle="numeric",
            borderDMShuffle="numeric",
            mapPolarityShuffle="numeric",
            gridScoreShuffle="numeric",
            speedScoreShuffle="numeric"
            ),
  prototype = list(session="",cmPerBin=2,smoothOccupancySd=3,smoothRateMapSd=3,minValidBinsAuto=20,reduceSize=T,
                   gridScoreNumberFieldsToDetect=40,gridScoreMinNumBinsPerField=50,gridScoreFieldThreshold=0.1,
                   borderPercentageThresholdField=20,borderMinBinsInField=10,nShufflings=100,minShiftMs=20000,
                   nAutoRotations=5,AutoRotationDegree=30))


### show ###
setMethod("show", "SpatialProperties2d",
          function(object){
            print(paste("session:",object@session))
            print(paste("cmPerBin:",object@cmPerBin))
            print(paste("smoothOccupancySd:",object@smoothOccupancySd))
            print(paste("smoothRateMapSd:",object@smoothRateMapSd))
            print(paste("reduceSize:",object@reduceSize))
            print(paste("nColMap:",object@nColMap))
            print(paste("nRowMap:",object@nRowMap))
            if(length(object@cellList)!=0){
              print(paste("nCells:",length(object@cellList)))
              print(paste("cellList:"))
              print(object@cellList)
            }
            if(length(object@peakRate)!=0){
              print("peakRate:")
              print(paste(object@peakRate))
              print("infoScore:")
              print(paste(object@infoScore))
              print("borderScore:")
              print(paste(object@borderScore))
              print("DM:")
              print(paste(object@borderDM))
              print("CM:")
              print(paste(object@borderCM))
              print("gridScore:")
              print(paste(object@gridScore))
            }
            if(length(object@speedScore)!=0){
              print("speedScore:")
              print(paste(object@speedScore))
            }
            if(length(object@speedRateSlope)!=0){
              print("speedRateSlope:")
              print(paste(object@speedRateSlope))
              print("speedRateIntercept:")
              print(paste(object@speedRateIntercept))
            }
            print(paste("nShufflings:",object@nShufflings))
            print(paste("shuffled values:",length(object@infoScoreShuffle)))
              
          })

#' Calculate the firing rate maps of neurons using a SpikeTrain and Positrack objects
#'
#' The occupancy map and the firing rate maps are smoothed with a Gaussian kernel
#' The amount of smoothing is determined by slots smoothOccupancySd and smoothRateMapSd of sp
#' 
#' 
#' @param sp SpatialProperties1d object
#' @param st SpikeTrain object
#' @param pt Positrack object
#' @return SpatialProperties2d object with the firing rate maps
#' 
#' @docType methods
#' @rdname firingRateMap2d-methods
setGeneric(name="firingRateMap2d",
           def=function(sp,st,pt)
           {standardGeneric("firingRateMap2d")}
)
#' @rdname firingRateMap2d-methods
#' @aliases firingRateMap2d,ANY,ANY-method
setMethod(f="firingRateMap2d",
          signature="SpatialProperties2d",
          definition=function(sp,st,pt)
          {
            if(length(pt@x)==0)
              stop(paste("pt@x has length of 0 in firingRateMap2d",st@session))
            if(st@nSpikes==0)
              stop(paste("st@nSpikes==0 in firingRateMap2d",st@session))
            
            sp@cellList<-st@cellList
            
            ## reduce the size of maps and map autocorrelation
            if(sp@reduceSize==T){
              x<-pt@x-min(pt@x,na.rm=T)+sp@cmPerBin
              y<-pt@y-min(pt@y,na.rm=T)+sp@cmPerBin
            }else{
              x<-pt@x
              y<-pt@y
            }
            
            ## use -1 as invalid values in c functions
            x[is.na(x)]<- -1.0
            y[is.na(y)]<- -1.0
            
            ## get the dimensions of the map
            sp@nRowMap=as.integer(((max(x)+1)/sp@cmPerBin)+1) # x in R is a row
            sp@nColMap=as.integer(((max(y)+1)/sp@cmPerBin)+1) # y in R is a col
            
            ## get spike position
            results<-.Call("spike_position_cwrap",
                  x,
                  y,
                  length(x),
                  as.integer(st@res),
                  as.integer(st@nSpikes),
                  as.integer(pt@resSamplesPerWhlSample),
                  as.integer(st@startInterval),
                  as.integer(st@endInterval),
                  length(st@startInterval))
            sp@xSpikes<-results[1,]
            sp@ySpikes<-results[2,]
            
            #plot(head(sp@xSpikes[which(sp@xSpikes!=-1.0)],n=20000),head(sp@ySpikes[which(sp@xSpikes!=-1.0)],n=20000))
            
            ## make the occupancy map
            sp@occupancy<-.Call("occupancy_map_cwrap",
                           sp@nRowMap,
                            sp@nColMap,
                            sp@cmPerBin,
                            sp@cmPerBin,
                            x,
                            y,
                            length(x),
                            pt@resSamplesPerWhlSample/pt@samplingRateDat*1000, ## ms per whl samples
                            as.integer(st@startInterval),
                            as.integer(st@endInterval),
                            length(st@startInterval),
                            as.integer(pt@resSamplesPerWhlSample))
            
            sp@occupancy
            
            #image((sp@occupancy),zlim=c(0,max(sp@occupancy,na.rm=T)))
            
            ## smooth the occupancy map
            sp@occupancy<- .Call("smooth_double_gaussian_2d_cwrap",
                  as.numeric(sp@occupancy),
                  sp@nColMap, # because C has a different way to order matrix as my c code
                  sp@nRowMap, #
                  sp@smoothOccupancySd/sp@cmPerBin,
                  -1.0)
            #image((sp@occupancy),zlim=c(0,max(sp@occupancy,na.rm=T)))
            
            ## make the 2d maps
            results<- .Call("firing_rate_map_2d_cwrap",
                            sp@nRowMap,
                            sp@nColMap,
                            sp@cmPerBin,
                            sp@cmPerBin,
                            sp@xSpikes,
                            sp@ySpikes,
                            as.integer(st@clu),
                            as.integer(st@nSpikes),
                            as.integer(sp@cellList),
                            length(sp@cellList),
                            as.numeric(sp@occupancy),
                            sp@smoothRateMapSd/sp@cmPerBin)
            sp@maps<-array(data=results,dim=(c(sp@nRowMap,sp@nColMap,length(sp@cellList))))
            return(sp)
          }
)

#' Calculate spatial statistics of the firing rate maps of neurons
#'
#' @param sp SpatialProperties2d object
#' @param st SpikeTrain object
#' @param pt Positrack object
#' @param border Set how the border of the environment will be detected. Value can be set to "rectangular" or "circular. 
#' Default value is "rectangular"
#' @return SpatialProperties2d object with the stats in the following slots: peakRate, infoScore, sparsity, borderScore and gridScore
#' 
#' @docType methods
#' @rdname getMapStats-methods
setGeneric(name="getMapStats",
           def=function(sp,st,pt,...)
           {standardGeneric("getMapStats")}
)
#' @rdname getMapStats-methods
#' @aliases getMapStats,ANY,ANY-method
setMethod(f="getMapStats",
          signature="SpatialProperties2d",
          definition=function(sp,st,pt,border="rectangular")
          {
            if(border!="rectangular" & border!= "circular")
              stop(paste("getMapstats, value of border can be \"rectangular\" or \"circular\" but is", border))
            
            ## make the maps
            sp<-firingRateMap2d(sp,st,pt)

            ### get peak rates
            sp@peakRate<-apply(sp@maps,3,max)

            ### get info scores
            sp@infoScore<- .Call("information_score_cwrap",
                length(sp@cellList),
                as.numeric(sp@maps),
                as.numeric(sp@occupancy),
                as.integer(sp@nColMap*sp@nRowMap))
            
            ### get sparsity scores
            sp@sparsity<- .Call("sparsity_score_cwrap",
                                length(sp@cellList),
                                as.numeric(sp@maps),
                                as.numeric(sp@occupancy),
                                as.integer(sp@nColMap*sp@nRowMap))
            
            if(border=="rectangular"){
              ### get border score
              results<-.Call("border_score_rectangular_environment_cwrap", ## if no field is detected, you get NaN values
                             as.integer(sp@cellList),
                             length(sp@cellList),
                             sp@nRowMap,
                             sp@nColMap,
                             sp@occupancy,
                             sp@maps,
                             sp@borderPercentageThresholdField,
                             as.integer(sp@borderMinBinsInField))
              sp@borderScore<-results[1,]
              sp@borderCM<-results[2,]
              sp@borderDM<-results[3,]
              sp@borderNumFieldsDetected<-results[4,]
            }
            
            if(border=="circular")
            {
              ### get border score
              results<-.Call("border_score_circular_environment_cwrap", ## if no field is detected, you get NaN values
                             as.integer(sp@cellList),
                             length(sp@cellList),
                             sp@nRowMap,
                             sp@nColMap,
                             sp@occupancy,
                             sp@maps,
                             sp@borderPercentageThresholdField,
                             as.integer(sp@borderMinBinsInField))
              sp@borderScore<-results[1,]
              sp@borderCM<-results[2,]
              sp@borderDM<-results[3,]
              sp@borderNumFieldsDetected<-results[4,]
              sp@mapPolarity<-results[5,]
            }

            # make spatial autocorrelations
            sp<-mapSpatialAutocorrelation(sp)
            sp@gridScore<-.Call("grid_score_cwrap",
                                length(sp@cellList),
                                sp@autos,
                                sp@nRowAuto,
                                sp@nColAuto,
                                sp@cmPerBin,
                                as.integer(sp@gridScoreNumberFieldsToDetect),
                                as.integer(sp@gridScoreMinNumBinsPerField),
                                sp@gridScoreFieldThreshold,
                                -2.0) #invalid
          return(sp)
          }
)



#' Calculate random spatial statistics of the firing rate maps of neurons
#'
#' The random values are obtained by shifting the position data.
#' The number of shufflings is nShufflings of the SpatialProperties2d object.
#'
#' @param sp SpatialProperties2d object
#' @param st SpikeTrain object
#' @param pt Positrack object
#' @param border Set how the border of the environment will be detected. Value can be set to "rectangular" or "circular.
#' @return SpatialProperties2d object with the random stats in the following slots: peakRate, infoScore, sparsity, borderScore and gridScore
#' 
#' @docType methods
#' @rdname getMapStatsShuffle-methods
setGeneric(name="getMapStatsShuffle",
           def=function(sp,st,pt,...)
           {standardGeneric("getMapStatsShuffle")}
)
#' @rdname getMapStatsShuffle-methods
#' @aliases getMapStatsShuffle,ANY,ANY-method
setMethod(f="getMapStatsShuffle",
          signature="SpatialProperties2d",
          definition=function(sp,st,pt,border="rectangular"){
        
            if(border!="rectangular" & border!= "circular")
              stop(paste("getMapstatsShuffle, value of border can be \"rectangular\" or \"circular\" but is", border))
            
            if(sp@nShufflings==0)            
              stop("sp@nShufflings==0")

            
            sp@peakRateShuffle=numeric()
            sp@infoScoreShuffle=numeric()
            sp@sparsityShuffle=numeric()
            sp@borderScoreShuffle=numeric()
            sp@borderCMShuffle=numeric()
            sp@borderDMShuffle=numeric()
            sp@gridScoreShuffle=numeric()
            
            print(paste("SpatialProperties2d shuffle",sp@nShufflings))
            for(i in 1:sp@nShufflings){
              pts<-shiftPositionRandom(pt)
              sp<-firingRateMap2d(sp,st,pts)
              sp@peakRateShuffle <- c(sp@peakRateShuffle,apply(sp@maps,3,max))
              sp@infoScoreShuffle <- c(sp@infoScoreShuffle,
                                       .Call("information_score_cwrap",
                                             length(sp@cellList),
                                             as.numeric(sp@maps),
                                             as.numeric(sp@occupancy),
                                             as.integer(sp@nColMap*sp@nRowMap)))
              sp@sparsityShuffle <- c(sp@sparsityShuffle, .Call("sparsity_score_cwrap",
                                                                length(sp@cellList),
                                                                as.numeric(sp@maps),
                                                                as.numeric(sp@occupancy),
                                                                as.integer(sp@nColMap*sp@nRowMap)))
              
              if(border=="rectangular"){
                ### get border score
                results<-.Call("border_score_rectangular_environment_cwrap", # if no field is detected, you get NaN values
                               as.integer(sp@cellList),
                               length(sp@cellList),
                               sp@nRowMap,
                               sp@nColMap,
                               sp@occupancy,
                               sp@maps,
                               sp@borderPercentageThresholdField,
                               as.integer(sp@borderMinBinsInField))
                sp@borderScoreShuffle <-c(sp@borderScoreShuffle,results[1,])
                sp@borderCMShuffle<-c(sp@borderCMShuffle,results[2,])
                sp@borderDMShuffle<-c(sp@borderDMShuffle,results[3,])
              }
              
              if(border=="circular")
              {
                ### get border score
                results<-.Call("border_score_circular_environment_cwrap", ## if no field is detected, you get NaN values
                               as.integer(sp@cellList),
                               length(sp@cellList),
                               sp@nRowMap,
                               sp@nColMap,
                               sp@occupancy,
                               sp@maps,
                               sp@borderPercentageThresholdField,
                               as.integer(sp@borderMinBinsInField))
                sp@borderScoreShuffle<-results[1,]
                sp@borderCMShuffle<-results[2,]
                sp@borderDMShuffle<-results[3,]
                sp@mapPolarityShuffle<-results[5,]
              }
              
              
              # make spatial autocorrelations
              sp<-mapSpatialAutocorrelation(sp)
              sp@gridScoreShuffle<-c(sp@gridScoreShuffle,
                                     .Call("grid_score_cwrap",
                                  length(sp@cellList),
                                  sp@autos,
                                  sp@nRowAuto,
                                  sp@nColAuto,
                                  sp@cmPerBin,
                                  as.integer(sp@gridScoreNumberFieldsToDetect),
                                  as.integer(sp@gridScoreMinNumBinsPerField),
                                  sp@gridScoreFieldThreshold,
                                  -2.0))
            }
            return(sp)
})


#' Calculate the spatial autocorrelation of the firing rate maps of neurons
#'
#' The autocorrelation is performed on the firing rate maps of the SpatialProperties2d object
#' 
#' @param sp SpatialProperties2d object
#' @return SpatialProperties2d object with the spatial autocorrelation in slot auto
#' 
#' @docType methods
#' @rdname mapSpatialAutocorrelation-methods
setGeneric(name="mapSpatialAutocorrelation",
           def=function(sp)
           {standardGeneric("mapSpatialAutocorrelation")}
)
#' @rdname mapSpatialAutocorrelation-methods
#' @aliases mapSpatialAutocorrelation,ANY,ANY-method
setMethod(f="mapSpatialAutocorrelation",
          signature="SpatialProperties2d",
          definition=function(sp)
          {
            if(length(sp@maps)==0)
              stop("Need to call firingRateMap2d first to run getMapStats")
            if(length(sp@occupancy)==0)
              stop("sp@occupancy length ==0")
            sp@nColAuto = as.integer((sp@nColMap*2)+1)
            sp@nRowAuto = as.integer((sp@nRowMap*2)+1)
            results<- .Call("map_autocorrelation_cwrap",
                  length(sp@cellList),
                  as.numeric(sp@maps),
                  sp@nRowMap,
                  sp@nColMap,
                  sp@nRowAuto,
                  sp@nColAuto,
                  as.integer(sp@minValidBinsAuto))
            sp@autos<-array(data=results,dim=(c(sp@nRowAuto,sp@nColAuto,length(sp@cellList))))
            return(sp)
          }
)

#' Calculate the spatial autocorrelation of the firing rate maps of neurons and removed the detected fields
#'
#' The autocorrelation is performed on the firing rate maps of the SpatialProperties2d object.
#' Fields are detected as when one wants to calculate grid scores
#' Usefull to test the different steps of the calculation of the grid score.
#' 
#' @param sp SpatialProperties2d object
#' @return SpatialProperties2d object with the spatial autocorrelation in slot autosDetect
#' 
#' @docType methods
#' @rdname autocorrelationNoFields-methods
setGeneric(name="autocorrelationNoFields",
           def=function(sp)
           {standardGeneric("autocorrelationNoFields")}
)
#' @rdname autocorrelationNoFields-methods
#' @aliases autocorrelationNoFields,ANY,ANY-method
setMethod(f="autocorrelationNoFields",
          signature="SpatialProperties2d",
          definition=function(sp)
          {
            if(length(sp@autos)==0)
              stop("Need to call mapSpatialAutocorrelation first to run gridScore()")
            results<-.Call("detect_and_remove_field_cwrap",
                           as.integer(sp@cellList),
                           length(sp@cellList),
                           sp@autos,
                           sp@nRowAuto,
                           sp@nColAuto,
                           as.integer(sp@gridScoreNumberFieldsToDetect),
                           as.integer(sp@gridScoreMinNumBinsPerField),
                           sp@gridScoreFieldThreshold,
                           -2.0) #invalid
            sp@autosDetect<-array(data=results,dim=(c(sp@nRowAuto,sp@nColAuto,length(sp@cellList))))
            return(sp)
          }
)


#' Calculate the spatial autocorrelation of the firing rate maps of neurons and get the region use for grid score calculation
#' Autocorrelation is performed on the firing rate maps of the SpatialProperties2d object.
#' Fields are detected as when one wants to calculate grid scores. 
#' Then the region with the 6 fields surrounding the center is kept
#' Usefull to test the different steps of the calculation of the grid score.
#' 
#' @param sp SpatialProperties2d object
#' @return SpatialProperties2d object with the spatial autocorrelation in slot autosDoughnut
#' 
#' @docType methods
#' @rdname autocorrelationDoughnut-methods
setGeneric(name="autocorrelationDoughnut",
           def=function(sp)
           {standardGeneric("autocorrelationDoughnut")}
)
#' @rdname autocorrelationDoughnut-methods
#' @aliases autocorrelationDoughnut,ANY,ANY-method
setMethod(f="autocorrelationDoughnut",
          signature="SpatialProperties2d",
          definition=function(sp)
          {
            if(length(sp@autos)==0)
              stop("Need to call mapSpatialAutocorrelation first to run gridScore()")
            results<-.Call("autocorrelation_doughnut_cwrap",
                           as.integer(sp@cellList),
                           length(sp@cellList),
                           sp@autos,
                           sp@nRowAuto,
                           sp@nColAuto,
                           as.integer(sp@gridScoreNumberFieldsToDetect),
                           as.integer(sp@gridScoreMinNumBinsPerField),
                           sp@gridScoreFieldThreshold,
                           sp@cmPerBin,
                           -2.0) #invalid
            sp@autosDoughnut<-array(data=results,dim=(c(sp@nRowAuto,sp@nColAuto,length(sp@cellList))))
            return(sp)
          }
)


#' Get the rotated region of the spatial autocorrelation of the firing rate maps when calculating grid scores
#' Autocorrelation is performed on the firing rate maps of the SpatialProperties2d object.
#' Fields are detected as when one wants to calculate grid scores. 
#' Then the region with the 6 fields surrounding the center is kept
#' The rotated copy of the region are stored in slot autosDoughnutRotate
#' Usefull to test the different steps of the calculation of the grid score.
#' 
#' @param sp SpatialProperties2d object
#' @return SpatialProperties2d object with the spatial autocorrelation in slot autosDoughnutRotate
#' 
#' @docType methods
#' @rdname autocorrelationDoughnutRotate-methods
setGeneric(name="autocorrelationDoughnutRotate",
           def=function(sp)
           {standardGeneric("autocorrelationDoughnutRotate")}
)
#' @rdname autocorrelationDoughnutRotate-methods
#' @aliases autocorrelationDoughnutRotate,ANY,ANY-method
setMethod(f="autocorrelationDoughnutRotate",
          signature="SpatialProperties2d",
          definition=function(sp)
          {
            if(length(sp@autos)==0)
              stop("Need to call mapSpatialAutocorrelation first to run gridScore()")
            results<-.Call("autocorrelation_doughnut_rotate_cwrap",
                           as.integer(sp@cellList),
                           length(sp@cellList),
                           sp@autos,
                           sp@nRowAuto,
                           sp@nColAuto,
                           as.integer(sp@gridScoreNumberFieldsToDetect),
                           as.integer(sp@gridScoreMinNumBinsPerField),
                           sp@gridScoreFieldThreshold,
                           sp@cmPerBin,
                           sp@nAutoRotations,
                           sp@AutoRotationDegree,
                           -2.0) #invalid
            sp@autosDoughnutRotate<-array(data=results,dim=(c(sp@nRowAuto,sp@nColAuto,length(sp@cellList)*sp@nAutoRotations)))
            return(sp)
          }
)

#' Calculate the grid scores from the spatial autocorrelations
#' 
#' 
#' @param sp SpatialProperties2d object
#' @return SpatialProperties2d object with the grid scores in the slot gridScore
#' 
#' @docType methods
#' @rdname gridScore-methods
setGeneric(name="gridScore",
           def=function(sp)
           {standardGeneric("gridScore")}
)
#' @rdname gridScore-methods
#' @aliases gridScore,ANY,ANY-method
setMethod(f="gridScore",
          signature="SpatialProperties2d",
          definition=function(sp)
          {
            if(length(sp@autos)==0)
              stop("Need to call mapSpatialAutocorrelation first to run gridScore()")
            sp@gridScore<-.Call("grid_score_cwrap",
                  length(sp@cellList),
                  sp@autos,
                  sp@nRowAuto,
                  sp@nColAuto,
                  sp@cmPerBin,
                  as.integer(sp@gridScoreNumberFieldsToDetect),
                  as.integer(sp@gridScoreMinNumBinsPerField),
                  sp@gridScoreFieldThreshold,
                  -2.0) #invalid
            return(sp)
          }
)

# setGeneric(name="gridOrientation",
#            def=function(sp)
#            {standardGeneric("gridOrientation")}
# )
# setMethod(f="gridOrientation",
#           signature="SpatialProperties2d",
#           definition=function(sp)
#           {
#             if(length(sp@autos)==0)
#               stop("Need to call mapSpatialAutocorrelation first to run gridOrientation()")
#             
# # 
# #             sp@gridOrientation<- .Call("grid_orientation_cwrap",
# #                                  as.integer(sp@cellList),
# #                                  length(sp@cellList),
# #                                  sp@autos,
# #                                  sp@nColAuto,
# #                                  sp@nRowAuto,
# #                                  sp@cmPerBin,
# #                                  as.integer(sp@gridScoreNumberFieldsToDetect),
# #                                  as.integer(sp@gridScoreMinNumBinsPerField),
# #                                  sp@gridScoreFieldThreshold,
# #                                  -2.0)
#             return(sp)
#           }
# )
# setGeneric(name="gridSpacing",
#            def=function(sp)
#            {standardGeneric("gridSpacing")}
# )
# setMethod(f="gridSpacing",
#           signature="SpatialProperties2d",
#           definition=function(sp)
#           {
#             if(length(sp@autos)==0)
#               stop("Need to call mapSpatialAutocorrelation first to run gridSpacing()")
#             
# #             sp@gridSpacing<- .Call("grid_spacing_cwrap",
# #                                         as.integer(sp@cellList),
# #                                         length(sp@cellList),
# #                                         sp@autos,
# #                                         sp@nColAuto,
# #                                         sp@nRowAuto,
# #                                         sp@cmPerBin,
# #                                         as.integer(sp@gridScoreNumberFieldsToDetect),
# #                                         as.integer(sp@gridScoreMinNumBinsPerField),
# #                                         sp@gridScoreFieldThreshold,
# #                                         -1.0)
#             return(sp)
#           }
# )

#' Calculate border scores from the firing maps 
#' 
#' 
#' @param sp SpatialProperties2d object
#' @return SpatialProperties2d object with the border scores in the slot borderScore
#' 
#' @docType methods
#' @rdname borderScore-methods
setGeneric(name="borderScore",
           def=function(sp,...)
           {standardGeneric("borderScore")}
)
#' @rdname borderScore-methods
#' @aliases borderScore,ANY,ANY-method
setMethod(f="borderScore",
          signature="SpatialProperties2d",
          definition=function(sp,border="rectangular")
          {
            
            
            if(length(sp@maps)==0)
              stop("Need to call firingRateMap2d first to run borderScore()")
            
            results<-.Call("border_score_rectangular_environment_cwrap",
                    as.integer(sp@cellList),
                    length(sp@cellList),
                    sp@nRowMap,
                    sp@nColMap,
                    sp@occupancy,
                    sp@maps,
                    sp@borderPercentageThresholdField,
                    as.integer(sp@borderMinBinsInField))
            
            sp@borderScore<-results[1,]
            sp@borderCM<-results[2,]
            sp@borderDM<-results[3,]
            sp@borderNumFieldsDetected<-results[4,]
            return(sp)
          }
)

#' Detect the border pixels in an environment and return
#' a matrix with border pixels identified.
#' 
#' Call firingRateMap2d() to construct the map first before calling this function
#' 
#' @param sp SpatialProperties2d object
#' @param border Set to rectangular or circular depending of the type of environment
#' Default value is rectangular 
#' @return A matrix with the border pixels set to non negative positive values
#'  
#'
#' @docType methods
#' @rdname borderDetection-methods
setGeneric(name="borderDetection",
           def=function(sp,...)
           {standardGeneric("borderDetection")}
)
#' @rdname borderDetection-methods
#' @aliases borderDetection,ANY,ANY-method
setMethod(f="borderDetection",
          signature="SpatialProperties2d",
          definition=function(sp,border="rectangular")
          {
            if(border!="rectangular" & border!= "circular")
              stop(paste("borderDetection, value of border can be \"rectangular\" or \"circular\" but is", border))
            
            if(length(sp@maps)==0)
              stop("Need to call firingRateMap2d first to run borderDetection()")
            if(border=="rectangular"){
              results<-.Call("border_detection_rectangular_environment_cwrap",
                             sp@nRowMap,
                             sp@nColMap,
                             sp@occupancy)
            }
            if(border=="circular"){
              results<-.Call("border_detection_circular_environment_cwrap",
                             sp@nRowMap,
                             sp@nColMap,
                             sp@occupancy)
            }
            return(results)
          })

#' Speed scores from spike train and positrack
#' 
#' 
#' @param sp SpatialProperties2d object
#' @param st SpikeTrain object with ifr
#' @param pt Positrack object
#' @param minSpeed Minimal speed to be considered
#' @param maxSpeed Maximal speed to be considered
#' @param runLm Logical, if TRUE a linear model will be build and slope and intercept calculated
#' @return SpatialProperties2d object with the speed scores in the slots SpeedScore
#' 
#' @docType methods
#' @rdname speedScore-methods
setGeneric(name="speedScore",
           def=function(sp,st,pt,minSpeed,maxSpeed,runLm)
           {standardGeneric("speedScore")}
)
#' @rdname speedScore-methods
#' @aliases speedScore,ANY,ANY-method
setMethod(f="speedScore",
          signature="SpatialProperties2d",
          definition=function(sp,st,pt,minSpeed=3,maxSpeed=100,runLm=F)
          {
            if(length(pt@speed)==0)
              stop("pt@speed has length of 0")
            if(dim(st@ifr)[1]!=length(sp@cellList))
              stop(paste("ifr dim",dim(st@ifr)[1], "is not equal to number of cells", length(sp@cellList)))
            
            ## get the ifr and ifrTime inside st@interval
            resTime<-st@ifrTime*st@samplingRate
            index<-as.logical(.Call("resWithinIntervals",
                                    length(st@startInterval),
                                    as.integer(st@startInterval),
                                    as.integer(st@endInterval),
                                    length(resTime),
                                    as.integer(resTime)))
            
            ifrSel<-matrix(st@ifr[,index],nrow=length(st@cellList))
            resTime<-resTime[index]
            
            ## get the speed for the res values
            speed<-getSpeedAtResValues(pt,resTime)
            
            ## speed filter
            index<-which(speed>minSpeed&speed<maxSpeed)
            ifrSel<-matrix(ifrSel[,index],nrow=length(st@cellList))
            speed<-speed[index]
            
            ## do the correlation
            sp@speedScore<-apply(ifrSel,1,cor,speed)
            if(runLm){
              c<-apply(ifrSel,1,function(x,y){lm(x~y)$coefficients},speed)
              sp@speedRateSlope=c[2,]
              sp@speedRateIntercept=c[1,]
            }
            return(sp)
          }
)


#' Random speed scores from spike train and positrack
#' 
#' 
#' @param sp SpatialProperties2d object
#' @param st SpikeTrain object with ifr
#' @param pt Positrack object
#' @param minSpeed Minimal speed to be considered
#' @param maxSpeed Maximal speed to be considered
#' @return SpatialProperties2d object with the random speed scores in the slots SpeedScoreShuffle
#' 
#' @docType methods
#' @rdname speedScoreShuffle-methods
setGeneric(name="speedScoreShuffle",
           def=function(sp,st,pt,minSpeed,maxSpeed)
           {standardGeneric("speedScoreShuffle")}
)
#' @rdname speedScoreShuffle-methods
#' @aliases speedScoreShuffle,ANY,ANY-method
setMethod(f="speedScoreShuffle",
          signature="SpatialProperties2d",
          definition=function(sp,st,pt,minSpeed=3,maxSpeed=100)
          {
            if(length(pt@speed)==0)
              stop("pt@speed has length of 0")
            if(dim(st@ifr)[1]!=length(sp@cellList))
              stop(paste("ifr dim",dim(st@ifr)[1], "is not equal to number of cells", length(sp@cellList)))
            if(sp@nShufflings==0)            
              stop("sp@nShufflings==0")
            
            ## get the ifr and ifrTime inside st@interval
            resTime<-st@ifrTime*st@samplingRate
            index<-as.logical(.Call("resWithinIntervals",
                                    length(st@startInterval),
                                    as.integer(st@startInterval),
                                    as.integer(st@endInterval),
                                    length(resTime),
                                    as.integer(resTime)))
            ifrSel<-matrix(st@ifr[,index],nrow=length(st@cellList))
            resTime<-resTime[index]
            #clear from previous data
            sp@speedScoreShuffle<-numeric()
            print(paste("speed shuffle",sp@nShufflings))
            for(i in 1:sp@nShufflings){
            
              ifrSels<-ifrSel
              resTimes<-resTime
              # shuffle speed
              pts<-shiftSpeedRandom(pt)
              ## get the speed for the res values
              speed<-getSpeedAtResValues(pts,resTimes)
              ## speed filter
              index<-which(speed>minSpeed&speed<maxSpeed)
              ifrSels<-matrix(ifrSels[,index],nrow=length(st@cellList))
              speed<-speed[index]
              ## do the correlation
              sp@speedScoreShuffle<-c(sp@speedScoreShuffle,
                                     apply(ifrSels,1,cor,speed))
            }
            return(sp)
          }
)


#' Calculate the center of mass of the firing rate maps in a SpatialProperties2d object
#' 
#' 
#' @param sp1 SpatialProperties2d object
#' @return matrix containing the center of mass of the firing rate maps
#' 
#' @docType methods
#' @rdname firingRateMapCenterOfMass-methods
setGeneric(name="firingRateMapCenterOfMass",
           def=function(sp)
           {standardGeneric("firingRateMapCenterOfMass")}
)
#' @rdname firingRateMapCenterOfMass-methods
#' @aliases firingRateMapCenterOfMass,ANY,ANY-method
setMethod(f="firingRateMapCenterOfMass",
          signature="SpatialProperties2d",
          definition=function(sp)
          {
            if(length(sp@maps)==0)
              stop("firingRateMapCenterOfMass: length(sp@maps)==0")
            if(length(dim(sp@maps))!=3)
              stop("firingRateMapCenterOfMass: length(dim(sp@maps))!=3")
            
            sp@maps[which(sp@maps==-1.0)]<-NA
            
            cm<-t(apply(sp@maps,3,centerOfMass))
            colnames(cm)<-c("row","col")
            return(cm)
          }
)




#' Calculate a Pearson correlation coefficients between the firing rate maps of two SpatialProperties2d objects
#' 
#' Used to estimate map stability between the firing rate maps in two conditions.
#' 
#' 
#' @param sp1 SpatialProperties2d object
#' @param sp2 SpatialProperties2d object
#' @return Numeric containing the Pearson correlation coefficients
#' 
#' @docType methods
#' @rdname firingRateMapCorrelation-methods
setGeneric(name="firingRateMapCorrelation",
           def=function(sp1,sp2)
           {standardGeneric("firingRateMapCorrelation")}
)
#' @rdname firingRateMapCorrelation-methods
#' @aliases firingRateMapCorrelation,ANY,ANY-method
setMethod(f="firingRateMapCorrelation",
          signature="SpatialProperties2d",
          definition=function(sp1,sp2)
          {
            if(class(sp1)!="SpatialProperties2d")
              stop("firingRateMapCorrelation: sp1 is not a SpatialProperties2d object")
            if(class(sp2)!="SpatialProperties2d")
              stop("firingRateMapCorrelation: sp2 is not a SpatialProperties2d object")
            if(length(sp1@maps)==0)
              stop("firingRateMapCorrelation: length(sp1@maps)==0")
            if(length(dim(sp1@maps))!=3)
              stop("firingRateMapCorrelation: length(dim(sp1@maps))!=3")
            if(!identical(dim(sp1@maps),dim(sp2@maps)))
              stop("firingRateMapCorrelation: the dimensions of the maps arrays in sp1 and sp2 are not equal")
            
            sp1@maps[which(sp1@maps==-1.0)]<-NA
            sp2@maps[which(sp2@maps==-1.0)]<-NA
  
            return(sapply(1:dim(sp1@maps)[3], 
                   function(i) cor(as.numeric(sp1@maps[,,i]), as.numeric(sp2@maps[,,i]),use="pairwise.complete.obs")))
          }
)




#' Return the firing rate maps in a data.frame
#' 
#' 
#' @param sp SpatialProperties2d object
#' @return data.frame with the firing rate maps. Column names are clu.id, x, y, rate
#' 
#' @docType methods
#' @rdname mapsAsDataFrame-methods
setGeneric(name="mapsAsDataFrame",
           def=function(sp)
           {standardGeneric("mapsAsDataFrame")}
)
#' @rdname mapsAsDataFrame-methods
#' @aliases mapsAsDataFrame,ANY,ANY-method
setMethod(f="mapsAsDataFrame",
          signature="SpatialProperties2d",
          definition=function(sp)
          {
            if(length(sp@maps)==0)
              stop("Need to call firingRateMap2d first to run mapsAsDataFrame()")
            
            data.frame(clu.id=rep(paste(sp@session,sp@cellList,sep="_"),each=sp@nRowMap*sp@nColMap),
                       x=rep(1:sp@nRowMap,sp@nColMap), ## warning the x and y were swapped
                       y=rep(1:sp@nColMap,each=sp@nRowMap), ## warning the x and y were swapped
                       rate=as.numeric(sp@maps))
          }
)

#' Return firing rate map statisitc as a data.frame
#' 
#' 
#' @param sp SpatialProperties2d object
#' @param shuffle Logical, if TRUE returns the stats of the shuffling analysis
#' @return data.frame with the firing rate map stats.
#' 
#' @docType methods
#' @rdname statsAsDataFrame-methods
setGeneric(name="statsAsDataFrame",
           def=function(sp,shuffle)
           {standardGeneric("statsAsDataFrame")}
)
#' @rdname statsAsDataFrame-methods
#' @aliases statsAsDataFrame,ANY,ANY-method
setMethod(f="statsAsDataFrame",
          signature="SpatialProperties2d",
          definition=function(sp,shuffle=FALSE)
          {
            
            if(shuffle==FALSE){
              if(length(sp@maps)==0)
                stop("Need to call getMapStats() before statsAsDataFrame")
              if(length(sp@infoScore)==0)
                stop("Need to call getMapStats() before statsAsDataFrame")
              df<-data.frame(clu.id=paste(sp@session,sp@cellList,sep="_"),
                         peakRate=sp@peakRate,
                         infoScore=sp@infoScore,
                         sparsity=sp@sparsity,
                         borderScore=sp@borderScore,
                         borderCM=sp@borderCM,
                         borderDM=sp@borderDM,
                         gridScore=sp@gridScore)
            }
            if(shuffle==TRUE){
              if(length(sp@maps)==0)
                stop("Need to call getMapStatsShuffle() before statsAsDataFrame with shuffle=T")
              if(length(sp@infoScoreShuffl)==0)
                stop("Need to call getMapStatsShuffle() before statsAsDataFrame with shuffle=T")
              df<-data.frame(clu.id=paste(sp@session,sp@cellList,sep="_"),
                         peakRate=sp@peakRateShuffle,
                         infoScore=sp@infoScoreShuffle,
                         sparsity=sp@sparsityShuffle,
                         borderScore=sp@borderScoreShuffle,
                         borderCM=sp@borderCMShuffle,
                         borderDM=sp@borderDMShuffle,
                         gridScore=sp@gridScoreShuffle)
            }
            return(df)
          }
)


