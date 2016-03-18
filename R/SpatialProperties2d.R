#################################################
#### definition of SpatialProperties2d Class  ###
#################################################
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
            spatialAutocorrelation="numeric",
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
            ##
            gridScore="numeric",
            gridOrientation="numeric",
            gridSpacing="numeric",
            gridScoreNumberFieldsToDetect="numeric",
            gridScoreMinNumBinsPerField="numeric",
            gridScoreFieldThreshold="numeric",
            ##
            nShufflings="numeric",
            minShiftMs="numeric",
            peakRateShuffle="numeric",
            infoScoreShuffle="numeric",
            sparsityShuffle="numeric",
            borderScoreShuffle="numeric",
            borderCMShuffle="numeric",
            borderDMShuffle="numeric",
            gridScoreShuffle="numeric"
            ),
      
  prototype = list(session="",cmPerBin=2,smoothOccupancySd=3,smoothRateMapSd=3,minValidBinsAuto=20,reduceSize=T,
                   gridScoreNumberFieldsToDetect=40,gridScoreMinNumBinsPerField=50,gridScoreFieldThreshold=0.1,
                   borderPercentageThresholdField=20,borderMinBinsInField=10,nShufflings=100,minShiftMs=20000))


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
              print("gridScore:")
              print(paste(object@gridScore))
            }
              
            print(paste("nShufflings:",object@nShufflings))
            print(paste("shuffled values:",length(object@infoScoreShuffle)))
              
          })

#### make firing rate maps 
setGeneric(name="firingRateMap2d",
           def=function(sp,st,pt)
           {standardGeneric("firingRateMap2d")}
)
setMethod(f="firingRateMap2d",
          signature="SpatialProperties2d",
          definition=function(sp,st,pt)
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@x)==0)
              stop("pt@x has length of 0")
            if(st@session=="")
              stop("st@session is empty")
            if(st@nSpikes==0)
              stop("st@nSpikes==0")
            
            sp@cellList<-st@cellList
            
            ## reduce the size of maps and map autocorrelation
            if(sp@reduceSize==T){
              x<-pt@x-min(pt@x,na.rm=T)+sp@cmPerBin
              y<-pt@y-min(pt@y,na.rm=T)+sp@cmPerBin
            }else{
              x<-pt@x
              y<-pt@y
            }
            
            #plot(x,y)
            ## use -1 as invalid values in c functions
            x[is.na(x)]<- -1.0
            y[is.na(y)]<- -1.0
            
            ## get the dimensions of the map
            sp@nColMap=as.integer(((max(x)+1)/sp@cmPerBin)+1) # x 
            sp@nRowMap=as.integer(((max(y)+1)/sp@cmPerBin)+1) # y
            
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
                           sp@nColMap,
                            sp@nRowMap,
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
            #image(t(sp@occupancy),zlim=c(0,max(sp@occupancy,na.rm=T)))
            
            ## smooth the occupancy map
            sp@occupancy<- .Call("smooth_double_gaussian_2d_cwrap",
                  as.numeric(sp@occupancy),
                  sp@nColMap,
                  sp@nRowMap,
                  sp@smoothOccupancySd/sp@cmPerBin,
                  -1.0)
            #image(t(sp@occupancy),zlim=c(0,max(sp@occupancy,na.rm=T)))
            
            ## make the 2d maps
            results<- .Call("firing_rate_map_2d_cwrap",
                            sp@nColMap,
                            sp@nRowMap,
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



#### getMapStats
setGeneric(name="getMapStats",
           def=function(sp)
           {standardGeneric("getMapStats")}
)
setMethod(f="getMapStats",
          signature="SpatialProperties2d",
          definition=function(sp)
          {
            
            if(length(sp@maps)==0)
              stop("Need to call firingRateMap2d first to run getMapStats")
            if(length(sp@occupancy)==0)
              stop("sp@occupancy length ==0")
          
            ### get peak rates
            sp@peakRate<-apply(sp@maps,3,max)

            ### get info scores
            sp@infoScore<- .Call("information_score_cwrap",
                as.integer(sp@cellList),
                length(sp@cellList),
                as.numeric(sp@maps),
                as.numeric(sp@occupancy),
                as.integer(sp@nColMap*sp@nRowMap))
            
            ### get sparsity scores
            sp@sparsity<- .Call("sparsity_score_cwrap",
                                as.integer(sp@cellList),
                                length(sp@cellList),
                                as.numeric(sp@maps),
                                as.numeric(sp@occupancy),
                                as.integer(sp@nColMap*sp@nRowMap))
            ### get border score
            results<-.Call("border_score_rectangular_environment_cwrap",
                           as.integer(sp@cellList),
                           length(sp@cellList),
                           sp@nColMap,
                           sp@nRowMap,
                           sp@occupancy,
                           sp@maps,
                           sp@borderPercentageThresholdField,
                           as.integer(sp@borderMinBinsInField))
            sp@borderScore<-results[1,]
            sp@borderCM<-results[2,]
            sp@borderDM<-results[3,]
            sp@borderNumFieldsDetected<-results[4,]
            
            sp@nColAuto = as.integer((sp@nColMap*2)+1)
            sp@nRowAuto = as.integer((sp@nRowMap*2)+1)
            results<- .Call("map_autocorrelation_cwrap",
                            as.integer(sp@cellList),
                            length(sp@cellList),
                            as.numeric(sp@maps),
                            sp@nColMap,
                            sp@nRowMap,
                            sp@nColAuto,
                            sp@nRowAuto,
                            as.integer(sp@minValidBinsAuto))
            sp@autos<-array(data=results,dim=(c(sp@nRowAuto,sp@nColAuto,length(sp@cellList))))
            
            sp@gridScore<-.Call("grid_score_cwrap",
                                as.integer(sp@cellList),
                                length(sp@cellList),
                                sp@autos,
                                sp@nColAuto,
                                sp@nRowAuto,
                                sp@cmPerBin,
                                as.integer(sp@gridScoreNumberFieldsToDetect),
                                as.integer(sp@gridScoreMinNumBinsPerField),
                                sp@gridScoreFieldThreshold,
                                -1.0) #invalid
          return(sp)
          }
)


#### getMapStatsShuffle 
setGeneric(name="getMapStatsShuffle",
           def=function(sp,st,pt)
           {standardGeneric("getMapStatsShuffle")}
)
setMethod(f="getMapStatsShuffle",
          signature="SpatialProperties2d",
          definition=function(sp,st,pt){
            
            if(sp@nShufflings==0)            
              stop("sp@nShufflings==0")
            
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@x)==0)
              stop("pt@x has length of 0")
            if(st@session=="")
              stop("st@session is empty")
            if(st@nSpikes==0)
              stop("st@nSpikes==0")
            
            sp@cellList<-st@cellList
            
            ## reduce the size of maps and map autocorrelation
            if(sp@reduceSize==T){
              x<-pt@x-min(pt@x,na.rm=T)+sp@cmPerBin
              y<-pt@y-min(pt@y,na.rm=T)+sp@cmPerBin
            }else{
              x<-pt@x
              y<-pt@y
            }
            
            #plot(x,y)
            ## use -1 as invalid values in c functions
            x[is.na(x)]<- -1.0
            y[is.na(y)]<- -1.0
            
            ## get the dimensions of the map
            sp@nColMap=as.integer(((max(x)+1)/sp@cmPerBin)+1) # x 
            sp@nRowMap=as.integer(((max(y)+1)/sp@cmPerBin)+1) # y
            ## dimensions of autos
            sp@nColAuto = as.integer((sp@nColMap*2)+1)
            sp@nRowAuto = as.integer((sp@nRowMap*2)+1)
            
            for(i in 1:sp@nShufflings){
              print(paste(i,"of",sp@nShufflings))
              
              ################################
              # shuffle the x and y position #
              # need to be changed           #
              # valid position values should #
              # not change environments      #
              ################################
              
              results<-shuffle.vectors(x=x,y=y,
                             time.per.sample.res=pt@resSamplesPerWhlSample,
                             sp@minShiftMs,
                             rs@samplingRate)
              x<-results[[1]]
              y<-results[[2]]
              
              
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
              
              ## make the occupancy map
              sp@occupancy<-.Call("occupancy_map_cwrap",
                                  sp@nColMap,
                                  sp@nRowMap,
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
              
              ## smooth the occupancy map
              sp@occupancy<- .Call("smooth_double_gaussian_2d_cwrap",
                                   as.numeric(sp@occupancy),
                                   sp@nColMap,
                                   sp@nRowMap,
                                   sp@smoothOccupancySd/sp@cmPerBin,
                                   -1.0)
              ## make the 2d maps
              results<- .Call("firing_rate_map_2d_cwrap",
                              sp@nColMap,
                              sp@nRowMap,
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
              
              
              ### get peak rates
              sp@peakRateShuffle <- c(sp@peakRateShuffle,apply(sp@maps,3,max))
              
              ### get info scores
              sp@infoScoreShuffle <- c(sp@infoScoreShuffle,
                                       .Call("information_score_cwrap",
                                             as.integer(sp@cellList),
                                             length(sp@cellList),
                                             as.numeric(sp@maps),
                                             as.numeric(sp@occupancy),
                                             as.integer(sp@nColMap*sp@nRowMap)))
              
              ### get sparsity scores
              sp@sparsityShuffle <- c(sp@sparsityShuffle, .Call("sparsity_score_cwrap",
                                                                as.integer(sp@cellList),
                                                                length(sp@cellList),
                                                                as.numeric(sp@maps),
                                                                as.numeric(sp@occupancy),
                                                                as.integer(sp@nColMap*sp@nRowMap)))
              ### get the border scores
              results<-.Call("border_score_rectangular_environment_cwrap",
                             as.integer(sp@cellList),
                             length(sp@cellList),
                             sp@nColMap,
                             sp@nRowMap,
                             sp@occupancy,
                             sp@maps,
                             sp@borderPercentageThresholdField,
                             as.integer(sp@borderMinBinsInField))
              
              sp@borderScoreShuffle <-c(sp@borderScoreShuffle,results[1,])
              sp@borderCMShuffle<-c(sp@borderCMShuffle,results[2,])
              sp@borderDMShuffle<-c(sp@borderDMShuffle,results[3,])
              
              ### get auto
              results<- .Call("map_autocorrelation_cwrap",
                              as.integer(sp@cellList),
                              length(sp@cellList),
                              as.numeric(sp@maps),
                              sp@nColMap,
                              sp@nRowMap,
                              sp@nColAuto,
                              sp@nRowAuto,
                              as.integer(sp@minValidBinsAuto))
              sp@autos<-array(data=results,dim=(c(sp@nRowAuto,sp@nColAuto,length(sp@cellList))))
              
              sp@gridScoreShuffle <-.Call("grid_score_cwrap",
                                  as.integer(sp@cellList),
                                  length(sp@cellList),
                                  sp@autos,
                                  sp@nColAuto,
                                  sp@nRowAuto,
                                  sp@cmPerBin,
                                  as.integer(sp@gridScoreNumberFieldsToDetect),
                                  as.integer(sp@gridScoreMinNumBinsPerField),
                                  sp@gridScoreFieldThreshold,
                                  -1.0) #invalid
            }
            return(sp)
})







setGeneric(name="mapSpatialAutocorrelation",
           def=function(sp)
           {standardGeneric("mapSpatialAutocorrelation")}
)
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
                  as.integer(sp@cellList),
                  length(sp@cellList),
                  as.numeric(sp@maps),
                  sp@nColMap,
                  sp@nRowMap,
                  sp@nColAuto,
                  sp@nRowAuto,
                  as.integer(sp@minValidBinsAuto))
            sp@autos<-array(data=results,dim=(c(sp@nRowAuto,sp@nColAuto,length(sp@cellList))))
            return(sp)
          }
)



setGeneric(name="gridScore",
           def=function(sp)
           {standardGeneric("gridScore")}
)
setMethod(f="gridScore",
          signature="SpatialProperties2d",
          definition=function(sp)
          {
            if(length(sp@autos)==0)
              stop("Need to call mapSpatialAutocorrelation first to run gridScore()")
            
            sp@gridScore<-.Call("grid_score_cwrap",
                  as.integer(sp@cellList),
                  length(sp@cellList),
                  sp@autos,
                  sp@nColAuto,
                  sp@nRowAuto,
                  sp@cmPerBin,
                  as.integer(sp@gridScoreNumberFieldsToDetect),
                  as.integer(sp@gridScoreMinNumBinsPerField),
                  sp@gridScoreFieldThreshold,
                  -1.0) #invalid
            return(sp)
          }
)
setGeneric(name="gridOrientation",
           def=function(sp)
           {standardGeneric("gridOrientation")}
)
setMethod(f="gridOrientation",
          signature="SpatialProperties2d",
          definition=function(sp)
          {
            if(length(sp@autos)==0)
              stop("Need to call mapSpatialAutocorrelation first to run gridOrientation()")
            

            sp@gridOrientation<- .Call("grid_orientation_cwrap",
                                 as.integer(sp@cellList),
                                 length(sp@cellList),
                                 sp@autos,
                                 sp@nColAuto,
                                 sp@nRowAuto,
                                 sp@cmPerBin,
                                 as.integer(sp@gridScoreNumberFieldsToDetect),
                                 as.integer(sp@gridScoreMinNumBinsPerField),
                                 sp@gridScoreFieldThreshold,
                                 -1.0)
            return(sp)
          }
)
setGeneric(name="gridSpacing",
           def=function(sp)
           {standardGeneric("gridSpacing")}
)
setMethod(f="gridSpacing",
          signature="SpatialProperties2d",
          definition=function(sp)
          {
            if(length(sp@autos)==0)
              stop("Need to call mapSpatialAutocorrelation first to run gridSpacing()")
            
            sp@gridSpacing<- .Call("grid_spacing_cwrap",
                                        as.integer(sp@cellList),
                                        length(sp@cellList),
                                        sp@autos,
                                        sp@nColAuto,
                                        sp@nRowAuto,
                                        sp@cmPerBin,
                                        as.integer(sp@gridScoreNumberFieldsToDetect),
                                        as.integer(sp@gridScoreMinNumBinsPerField),
                                        sp@gridScoreFieldThreshold,
                                        -1.0)
            return(sp)
          }
)


setGeneric(name="borderScore",
           def=function(sp)
           {standardGeneric("borderScore")}
)
setMethod(f="borderScore",
          signature="SpatialProperties2d",
          definition=function(sp)
          {
            if(length(sp@maps)==0)
              stop("Need to call firingRateMap2d first to run borderScore()")
            
            results<-.Call("border_score_rectangular_environment_cwrap",
                    as.integer(sp@cellList),
                    length(sp@cellList),
                    sp@nColMap,
                    sp@nRowMap,
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





