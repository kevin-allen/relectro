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
            autosDetect="array",
            autosDoughnut="array",
            autosDoughnutRotate="array",
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



#### getMapStats
setGeneric(name="getMapStats",
           def=function(sp,st,pt)
           {standardGeneric("getMapStats")}
)
setMethod(f="getMapStats",
          signature="SpatialProperties2d",
          definition=function(sp,st,pt)
          {
            ## make the maps
            sp<-firingRateMap2d(sp,st,pt)

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
            
            # make spatial autocorrelations
            sp<-mapSpatialAutocorrelation(sp)
            sp@gridScore<-.Call("grid_score_cwrap",
                                as.integer(sp@cellList),
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
                                             as.integer(sp@cellList),
                                             length(sp@cellList),
                                             as.numeric(sp@maps),
                                             as.numeric(sp@occupancy),
                                             as.integer(sp@nColMap*sp@nRowMap)))
              sp@sparsityShuffle <- c(sp@sparsityShuffle, .Call("sparsity_score_cwrap",
                                                                as.integer(sp@cellList),
                                                                length(sp@cellList),
                                                                as.numeric(sp@maps),
                                                                as.numeric(sp@occupancy),
                                                                as.integer(sp@nColMap*sp@nRowMap)))
              ### get border score
              results<-.Call("border_score_rectangular_environment_cwrap",
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
              
              # make spatial autocorrelations
              sp<-mapSpatialAutocorrelation(sp)
              sp@gridScoreShuffle<-c(sp@gridScoreShuffle,
                                     .Call("grid_score_cwrap",
                                    as.integer(sp@cellList),
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
                  sp@nRowMap,
                  sp@nColMap,
                  sp@nRowAuto,
                  sp@nColAuto,
                  as.integer(sp@minValidBinsAuto))
            sp@autos<-array(data=results,dim=(c(sp@nRowAuto,sp@nColAuto,length(sp@cellList))))
            return(sp)
          }
)

setGeneric(name="autocorrelationNoFields",
           def=function(sp)
           {standardGeneric("autocorrelationNoFields")}
)
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
setGeneric(name="autocorrelationDoughnut",
           def=function(sp)
           {standardGeneric("autocorrelationDoughnut")}
)
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

setGeneric(name="autocorrelationDoughnutRotate",
           def=function(sp)
           {standardGeneric("autocorrelationDoughnutRotate")}
)
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
            
# 
#             sp@gridOrientation<- .Call("grid_orientation_cwrap",
#                                  as.integer(sp@cellList),
#                                  length(sp@cellList),
#                                  sp@autos,
#                                  sp@nColAuto,
#                                  sp@nRowAuto,
#                                  sp@cmPerBin,
#                                  as.integer(sp@gridScoreNumberFieldsToDetect),
#                                  as.integer(sp@gridScoreMinNumBinsPerField),
#                                  sp@gridScoreFieldThreshold,
#                                  -2.0)
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
            
#             sp@gridSpacing<- .Call("grid_spacing_cwrap",
#                                         as.integer(sp@cellList),
#                                         length(sp@cellList),
#                                         sp@autos,
#                                         sp@nColAuto,
#                                         sp@nRowAuto,
#                                         sp@cmPerBin,
#                                         as.integer(sp@gridScoreNumberFieldsToDetect),
#                                         as.integer(sp@gridScoreMinNumBinsPerField),
#                                         sp@gridScoreFieldThreshold,
#                                         -1.0)
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


setGeneric(name="speedScore",
           def=function(sp,st,pt,minSpeed,maxSpeed,runLm)
           {standardGeneric("speedScore")}
)
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



setGeneric(name="speedScoreShuffle",
           def=function(sp,st,pt,minSpeed,maxSpeed)
           {standardGeneric("speedScoreShuffle")}
)
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



setGeneric(name="mapsAsDataFrame",
           def=function(sp)
           {standardGeneric("mapsAsDataFrame")}
)
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


setGeneric(name="statsAsDataFrame",
           def=function(sp,shuffle)
           {standardGeneric("statsAsDataFrame")}
)
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




