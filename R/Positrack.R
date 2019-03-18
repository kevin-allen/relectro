#' An S4 class representing the path of an animal during a recording session
#' 
#' This class is used to manipulate and represent the position data of the animal.
#' The first data point is at time resSamplesPerWhlSample and not 0.
#' The data are usually loaded from .whd, .res_samples_per_whl_sample, .sampling_rate_dat and .px_per_cm
#' 
#' The position data likely come from a tracking system (e.g. positrack) in which the y-axis has its origin
#' at the top-left of the screen. In contrast, most plotting functions in R have the y-axis origin at the bottom-left
#' of the screen. This is important to know if you want to find out where landmarks relative to the animal.
#' 
#' @slot session A character vector containing the names of the recording session.
#' @slot path The directory in which the files of the session are located.
#' @slot pxPerCm Pixels per centimeter in the x and y position data.
#' @slot samplingRateDat Sampling rate of the .dat files in this recording session
#' @slot resSamplesPerWhlSample Number of samples in the .dat files between each position values.
#' @slot xWhl A numeric vector containing x data loaded from the .whl file.
#' @slot yWhl A numeric vector containing y data loaded from the .whl file.
#' @slot hdWhd A numeric vector containing the head direction of the animal, loaded from .whd file.
#' @slot x x position of the animal in cm
#' @slot y y position of the animal in cm
#' @slot hd Head direction of the animal
#' @slot lin Linearized position of the animal. Used for linear track data.
#' @slot dir Direction in the linearized position data (0 or 1)
#' @slot speed Linear speed in cm per sec
#' @slot angularSpeed Angular speed in degrees per sec
#' @slot res Sample number associated with each position value
#' @slot defaultXYSmoothing Default smoothing apply to the x and y position.
#' @slot minShiftMs Minimum shift of the position vector during a shuffling procedure
#' @slot percentInvalidPosition Percentage of data points with invalid position
#' 
Positrack <- setClass(
  "Positrack", ## name of the class
  slots=c(session="character",
          path="character",
          pxPerCm="numeric",
          samplingRateDat="numeric",
          resSamplesPerWhlSample="numeric",
          xWhl="numeric",
          yWhl="numeric",
          hdWhd="numeric",
          x="numeric", # in cm
          y="numeric",
          hd="numeric",
          lin="numeric", # linear position
          dir="numeric", # direction in linear position
          speed="numeric",
          angularSpeed="numeric",
          res="numeric",
          defaultXYSmoothing="numeric",
          minShiftMs="numeric",
          percentInvalidPosition="numeric"
          ),  # cell list to limit the analysis to these cells
  prototype = list(session="",path="",defaultXYSmoothing=2,minShiftMs=20000))


#' Load position data from the .whl and .whd files
#'
#' This function also use the .res_samples_per_whl_sample, .sampling_rate_dat and .px_per_cm files
#' Some smoothing is done on the x and y data.
#' Slots x and y are in cm.
#' 
#' @param pt Positrack object
#' @return Positrack object
#' 
#' @docType methods
#' @rdname loadPositrack-methods
setGeneric(name="loadPositrack",
           def=function(pt)
           {standardGeneric("loadPositrack")}
)

#' @rdname loadPositrack-methods
#' @aliases loadPositrack,ANY,ANY-method
setMethod(f="loadPositrack",
          signature="Positrack",
          definition=function(pt)
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(pt@path==""){
              pt@path=getwd()
            }
            pathSession=paste(pt@path,pt@session,sep="/")
            if(!file.exists(paste(pathSession,"whd",sep=".")))
              stop("needs ",paste(pathSession,"whd",sep="."))
            if(!file.exists(paste(pathSession,"res_samples_per_whd_sample",sep="."))&
               !file.exists(paste(pathSession,"res_samples_per_whl_sample",sep=".")))
              stop("needs ",paste(pathSession,"res_samples_per_whd_sample",sep="."))
            if(!file.exists(paste(pathSession,"sampling_rate_dat",sep=".")))
              stop("needs ",paste(pathSession,"sampling_rate_dat",sep="."))
            if(!file.exists(paste(pathSession,"px_per_cm",sep=".")))
              stop("needs ",paste(pathSession,"px_per_cm",sep="."))
            ## get sampling rate
            whd<-read.table(paste(pathSession,"whd",sep="."))
            pt@xWhl<-whd$V1
            pt@yWhl<-whd$V2
            pt@hdWhd<-whd$V3
            
            if(file.exists(paste(pathSession,"res_samples_per_whd_sample",sep=".")))
            {
              pt@resSamplesPerWhlSample<-read.table(paste(pathSession,"res_samples_per_whd_sample",sep="."))$V1  
            } else {
              pt@resSamplesPerWhlSample<-read.table(paste(pathSession,"res_samples_per_whl_sample",sep="."))$V1  
            }
            
            pt@pxPerCm<-read.table(paste(pathSession,"px_per_cm",sep="."))$V1
            pt@samplingRateDat<-read.table(paste(pathSession,"sampling_rate_dat",sep="."))$V1
            
            ## check that things add up
            if(length(pt@xWhl)!=length(pt@hdWhd)){
              print("length whl:",paste(length(pt@xWhl),"whd:",length(pt@hdWhd)))
              stop("Problem with length of whl and whd files")
            }
            if(pt@samplingRateDat<1|pt@samplingRateDat>100000)
              stop(paste("pt@sampingRateDat is out of bound:",pt@samplingRateDat))
            if(pt@pxPerCm<1|pt@pxPerCm>10000)
              stop(paste("pt@pxPerCm is out of bound:", pt@pxPerCm))
            if(max(pt@xWhl)>5000|min(pt@xWhl)< -1.0){
              print(paste("min x:",min(pt@xWhl),"max x:",max(pt@xWhl)))
              stop(paste("values of pt@xWhl are out of bound"))
            }
            if(max(pt@yWhl)>5000|min(pt@yWhl)< -1.0){
              print(paste("min y:",min(pt@yWhl),"max y:",max(pt@yWhl)))
              stop(paste("values of pt@xWhl are out of bound"))
            }
            if(max(pt@hdWhd>360))
              stop(paste("max value of pt@hdWhd > 360:",max(pt@hdWhd)))
            if(min(pt@hdWhd< -1.0))
              stop(paste("min value of pt@hdWhd < -1.0:",min(pt@hdWhd)))
           
            
             ## pt@xWhl will never be changed, 
            pt@x<-pt@xWhl
            pt@y<-pt@yWhl
            pt@hd<-pt@hdWhd
            
            pt@x[which(pt@x!=-1.0)]<-pt@x[which(pt@x!=-1.0)]/pt@pxPerCm
            pt@y[which(pt@y!=-1.0)]<-pt@y[which(pt@y!=-1.0)]/pt@pxPerCm
            
            ## set the res value associated with eack whl sample
            pt@res<-seq(from=pt@resSamplesPerWhlSample,by=pt@resSamplesPerWhlSample,length.out=length(pt@x))
            
            ## apply a bit a smoothing to the x and y
            pt@x<-smoothGaussian(x=pt@x,sd=pt@defaultXYSmoothing,invalid=-1.0)
            pt@y<-smoothGaussian(x=pt@y,sd=pt@defaultXYSmoothing,invalid=-1.0)
         
            ## apply no smoothing to head direction
               
            ## get the speed from position
            pt@speed<- .Call("speed_from_whl_cwrap",
                            as.numeric(pt@x),
                            as.numeric(pt@y),
                            length(pt@x),
                            1.0, # already in cm
                            pt@samplingRateDat, 
                            pt@resSamplesPerWhlSample)
            
            
            ## get the angular speed from position
            pt@angularSpeed<- .Call("angular_speed_from_hd_cwrap",
                                     as.numeric(pt@hd),
                                     length(pt@hd),
                                     4,
                                     4,
                                     pt@samplingRateDat,
                                     pt@resSamplesPerWhlSample)
            
            # set -1.0 to NA for calculation in R       
            pt@x[which(pt@x==-1.0)]<-NA
            pt@y[which(pt@y==-1.0)]<-NA
            pt@hd[which(pt@hd==-1.0)]<-NA
            pt@speed[which(pt@speed==-1.0)]<-NA
            pt@angularSpeed[which(pt@angularSpeed==-1.0)]<-NA
            pt@percentInvalidPosition<-sum(is.na(pt@x))/length(pt@x)*100
            
            return(pt)
          }
)





#' set position values from arguments
#'
#' This function also use the .res_samples_per_whl_sample, .sampling_rate_dat and .px_per_cm files
#' 
#' @param pt Positrack object
#' @param pxX Numeric with the x position of animal. x is in pixels. invalid is -1.0
#' @param pxY Numeric with the y position of animal. y is in pixels. invalid is -1.0
#' @param hd Numeric with the head direction of animal in degrees. invalid is -1.0
#' @param resSamplesPerWhlSample Number of electrophysiological samples per position sample
#' @param samplingRateDat Sampling rate of electrophysiological data
#' @param pxPerCm Number of pixels per cm
#' @return Positrack object
#' 
#' @docType methods
#' @rdname setPositrack-methods
setGeneric(name="setPositrack",
           def=function(pt,pxX,pxY,hd,resSamplesPerWhlSample,samplingRateDat,pxPerCm)
           {standardGeneric("setPositrack")}
)
#' @rdname setPositrack-methods
#' @aliases setPositrack,ANY,ANY-method
setMethod(f="setPositrack",
          signature="Positrack",
          definition=function(pt,pxX,pxY,hd,resSamplesPerWhlSample,samplingRateDat,pxPerCm)
          {
          
            ## get sampling rate
            pt@xWhl<-pxX
            pt@yWhl<-pxY
            pt@hdWhd<-hd
            pt@resSamplesPerWhlSample<-resSamplesPerWhlSample
            pt@samplingRateDat<-samplingRateDat
            pt@pxPerCm<-pxPerCm
            
            ## check that things add up
            if(length(pt@xWhl)!=length(pt@hdWhd))
              stop("Problem with length of xWhl and hdWhd")
            if(pt@samplingRateDat<1|pt@samplingRateDat>100000)
              stop(paste("pt@sampingRateDat is out of bound:",pt@samplingRateDat))
            if(pt@pxPerCm<1|pt@pxPerCm>10000)
              stop(paste("pt@pxPerCm is out of bound:", pt@pxPerCm))
            if(max(pt@xWhl)>2000|min(pt@xWhl)< -1.0){
              print(paste("min x:",min(pt@xWhl),"max x:",max(pt@xWhl)))
              stop(paste("values of pt@xWhl are out of bound"))
            }
            if(max(pt@yWhl)>2000|min(pt@yWhl)< -1.0){
              print(paste("min y:",min(pt@yWhl),"max y:",max(pt@yWhl)))
              stop(paste("values of pt@xWhl are out of bound"))
            }
            if(max(pt@hdWhd>360))
              stop(paste("max value of pt@hdWhd > 360:",max(pt@hdWhd)))
            if(min(pt@hdWhd< -1.0))
              stop(paste("min value of pt@hdWhd < -1.0:",min(pt@hdWhd)))
            ## pt@xWhl will never be changed, 
            
            pt@x<-pt@xWhl
            pt@y<-pt@yWhl
            pt@hd<-pt@hdWhd
            
            pt@x[which(pt@x!=-1.0)]<-pt@x[which(pt@x!=-1.0)]/pt@pxPerCm
            pt@y[which(pt@y!=-1.0)]<-pt@y[which(pt@y!=-1.0)]/pt@pxPerCm
            
            ## set the res value associated with eack whl sample
            pt@res<-seq(from=pt@resSamplesPerWhlSample,by=pt@resSamplesPerWhlSample,length.out=length(pt@x))
            
            ## apply a bit a smoothing to the x and y
            pt@x<-smoothGaussian(x=pt@x,sd=pt@defaultXYSmoothing,invalid=-1.0)
            pt@y<-smoothGaussian(x=pt@y,sd=pt@defaultXYSmoothing,invalid=-1.0)
            
            ## apply no smoothing to head direction
            
            ## get the speed from position
            pt@speed<- .Call("speed_from_whl_cwrap",
                             pt@x,
                             pt@y,
                             length(pt@x),
                             1.0, # speed is already in cm
                             pt@samplingRateDat, 
                             pt@resSamplesPerWhlSample)
            
            ## get the angular speed from position
            pt@angularSpeed<- .Call("angular_speed_from_hd_cwrap",
                                    pt@hd,
                                    length(pt@hd),
                                    4,
                                    4,
                                    pt@samplingRateDat,
                                    pt@resSamplesPerWhlSample)
            
            # set -1.0 to NA for calculation in R       
            pt@x[which(pt@x==-1.0)]<-NA
            pt@y[which(pt@y==-1.0)]<-NA
            pt@hd[which(pt@hd==-1.0)]<-NA
            pt@speed[which(pt@speed==-1.0)]<-NA
            pt@angularSpeed[which(pt@angularSpeed==-1.0)]<-NA
            pt@percentInvalidPosition<-sum(is.na(pt@x))/length(pt@x)*100
            return(pt)
          }
)


#' Apply some smoothing to the x and y position data
#'
#' A Gaussian kernel is used for smoothing
#' 
#' @param pt Positrack object
#' @param sd Standard deviation of the smoothing kernel. The unit is position samples.
#' @return Positrack object with slots x and y smoothed.
#' 
#' @docType methods
#' @rdname smoothxy-methods
setGeneric(name="smoothxy",
           def=function(pt,sd)
           {standardGeneric("smoothxy")}
)

#' @rdname smoothxy-methods
#' @aliases smoothxy,ANY,ANY-method
setMethod(f="smoothxy",
          signature="Positrack",
          definition=function(pt,sd=2)
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@x)==0)
              stop("pt@x has length of 0")
            ## smooth x and y
            pt@x[which(is.na(pt@x))]<- -1.0
            pt@y[which(is.na(pt@y))]<- -1.0
            pt@x<-smoothGaussian(pt@x,sd=sd,invalid=-1.0)
            pt@y<-smoothGaussian(pt@y,sd=sd,invalid=-1.0)
            pt@x[which(pt@x==-1.0)]<- NA
            pt@y[which(pt@y==-1.0)]<- NA
            return(pt)
          }
)

#' Apply some smoothing to the hd array
#'
#' A Gaussian kernel is used for smoothing
#' This function is not working at the moment. A warning message will appear
#' The c function needs to be written
#' 
#' @param pt Positrack object
#' @param sd Standard deviation of the smoothing kernel. The unit is position samples.
#' @return Positrack object with slot hd smoothed.
#' 
#' @docType methods
#' @rdname smoothhd-methods
setGeneric(name="smoothhd",
           def=function(pt,sd)
           {standardGeneric("smoothhd")}
)
#' @rdname smoothhd-methods
#' @aliases smoothhd,ANY,ANY-method
setMethod(f="smoothhd",
          signature="Positrack",
          definition=function(pt,sd=2)
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@hd)==0)
              stop("pt@hd has length of 0")
            ## smooth hd
            pt@hd[which(is.na(pt@hd))]<- -1.0
            
            stop("Smoothing of hd data with an function that does not take into account that data are circular")
            pt@hd<-smoothGaussian(pt@hd,sd=sd,invalid=-1.0)
            pt@hd[which(pt@hd==-1.0)]<- NA
            return(pt)
          }
)


#' Set position data to NA for which the animal speed was not between a given range
#'
#' The slots x, y, hd, speed are affected.
#' Please note that the speed array is affected. 
#' 
#' @param pt Positrack object
#' @param minSpeed Minimal running speed
#' @param maxSpeed Maximal running speed
#' @return Positrack object with invalid values outside the speed range.
#' 
#' @docType methods
#' @rdname speedFilter-methods
setGeneric(name="speedFilter",
           def=function(pt,minSpeed,maxSpeed)
           {standardGeneric("speedFilter")}
)

#' @rdname speedFilter-methods
#' @aliases speedFilter,ANY,ANY-method
setMethod(f="speedFilter",
          signature="Positrack",
          definition=function(pt,minSpeed=3,maxSpeed=100)
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@x)==0)
              stop("pt@x has length of 0")
            
            ## set time when outside of speed to NA
            index<-which(pt@speed<minSpeed|pt@speed>maxSpeed)
            pt@x[index]<- NA
            pt@y[index]<- NA
            pt@hd[index]<- NA
            pt@speed[index]<- NA
            return(pt)
          }
)



#' Set position data to NA if animal's direction was not the one given in arguments
#'
#' The slots x, y, hd, speed and lin are affected.
#' Please note that this is not head direction but rather the direction of movement in 1D environment 
#' 
#' @param pt Positrack object
#' @param direction Value of 1 or 0 refereing to direction
#' @return Positrack object with only valid value in the specific direction.
#' 
#' @docType methods
#' @rdname directionFilter-methods
setGeneric(name="directionFilter",
           def=function(pt,direction)
           {standardGeneric("directionFilter")}
)
#' @rdname directionFilter-methods
#' @aliases directionFilter,ANY,ANY-method
setMethod(f="directionFilter",
          signature="Positrack",
          definition=function(pt,direction=0)
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@x)==0)
              stop("pt@x has length of 0")
            if(length(pt@lin)==0)
              stop(paste("pt@lin has length of 0 in direction filter. Need to linearized the position data first  with linearzeLinearTrack(ptlt)"))
            if(direction!=0&direction!=1)
              stop(paste("direction should have value of 0 or 1"))
            
            ## set time when outside of speed to NA
            index<-which(pt@dir!=direction)
            pt@x[index]<- NA
            pt@y[index]<- NA
            pt@hd[index]<- NA
            pt@speed[index]<- NA
            pt@lin[index]<-NA
            return(pt)
          }
)


#' Set position data outside time intervals to NA
#'
#' The slots x, y, hd, speed, lin and dir are affected.
#' The intervals are in sample values of the electrophysiological data 
#' 
#' @param pt Positrack object
#' @param s Vector or matrix, If matrix, then should have 2 columns, for beginning and end of intervals. 
#' If vector, then it is the beginning of intervals.
#' @param e Vector, end of the intervals. Not used if s is a matrix
#' @return Positrack object with only valid value within the intervals.
#' 
#' @docType methods
#' @rdname setInvalidOutsideInterval-methods
setGeneric(name="setInvalidOutsideInterval",
           def=function(pt,s,e)
           {standardGeneric("setInvalidOutsideInterval")}
)
#' @rdname setInvalidOutsideInterval-methods
#' @aliases setInvalidOutsideInterval,ANY,ANY-method
setMethod(f="setInvalidOutsideInterval",
          signature="Positrack",
          definition=function(pt,s,e="")
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@x)==0)
              stop("pt@x has length of 0")

            ## if s is a matrix, then e is ignored
            if(class(s)=="matrix"){
              startInterval<-as.numeric(s[,1])
              endInterval<-as.numeric(s[,2])
            }else{
              startInterval<-s
              endInterval<-e
            }
            if(length(startInterval)!=length(endInterval))
              stop("problem with length of startInterval and endInterval")
            if(any(startInterval>endInterval))
              stop("startInterval>endInterval")
            if(any(diff(startInterval)<0))
              stop("problem with chonology within startInterval")
            if(any(diff(endInterval)<0))
              stop("problem with chonology within endInterval")
            if(any(startInterval[-1]-endInterval[-length(endInterval)]<0))
              stop("problem with chronology between intervals, from end to next start")

            ## slow, 1 sec for a list of 2 intervals
            #isin<-sapply(pt@res,function(x)any(x>=startInterval&x<=endInterval)) # takes seconds
            isin<-as.logical(.Call("resWithinIntervals", ## about 1000 times faster
                  length(startInterval),
                  as.integer(startInterval),
                  as.integer(endInterval),
                  length(pt@res),
                  as.integer(pt@res)))
            pt@x[which(!isin)]<-NA
            pt@y[which(!isin)]<-NA
            pt@hd[which(!isin)]<-NA
            pt@speed[which(!isin)]<-NA
            pt@angularSpeed[which(!isin)]<-NA
            if(length(pt@lin)!=0) pt@lin[which(!isin)]<-NA
            if(length(pt@dir)!=0) pt@dir[which(!isin)]<-NA
            return(pt)
          }
)


#' Get time intervals at which the speed of the animal was within a specified ranged
#'
#' The intervals are in sample values of the electrophysiological data 
#' 
#' @param pt Positrack object
#' @param minSpeed Minimal speed of the animal. 
#' @param maxSpeed Maximal speed of the animal.
#' @return matrix with the time intervals
#' 
#' @docType methods
#' @rdname getIntervalsAtSpeed-methods
setGeneric(name="getIntervalsAtSpeed",
           def=function(pt,minSpeed,maxSpeed)
           {standardGeneric("getIntervalsAtSpeed")}
)
#' @rdname getIntervalsAtSpeed-methods
#' @aliases getIntervalsAtSpeed,ANY,ANY-method
setMethod(f="getIntervalsAtSpeed",
          signature="Positrack",
          definition=function(pt,minSpeed,maxSpeed)
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@x)==0)
              stop("pt@x has length of 0")
            if(minSpeed>maxSpeed)
              stop("minSpeed>maxSpeed")
            results<-.Call("speed_intervals_cwrap", 
                  pt@speed,
                  length(pt@speed),
                  pt@resSamplesPerWhlSample,
                  minSpeed,
                  maxSpeed)
            colnames(results)<-c("start","end")
            return(results)
          }
)


#' Get time intervals at which the animal ran in a given direction (0 or 1)
#'
#' The intervals are in sample values of the electrophysiological data.
#' The direction is based on the direction of movement in 1D environment.
#' 
#' @param pt Positrack object
#' @param direction 0 or 1 indicating the direction of movment
#' @return matrix with the time intervals
#' 
#' @docType methods
#' @rdname getIntervalsAtDirection-methods
setGeneric(name="getIntervalsAtDirection",
           def=function(pt,direction)
           {standardGeneric("getIntervalsAtDirection")}
)
#' @rdname getIntervalsAtDirection-methods
#' @aliases getIntervalsAtDirection,ANY,ANY-method
setMethod(f="getIntervalsAtDirection",
          signature="Positrack",
          definition=function(pt,direction=0)
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@x)==0)
              stop("pt@x has length of 0")
            if(length(pt@dir)==0)
              stop("pt@dir has length of 0")
            if(direction!=0& direction!=1)
              stop("direction should have a value of 0 or 1")
            x<-pt@dir
            x[is.na(x)]<- -1.0
            results<-.Call("direction_intervals_cwrap", 
                           as.integer(x),
                           length(x),
                           as.integer(pt@resSamplesPerWhlSample),
                           as.integer(direction))
            colnames(results)<-c("start","end")
            return(results)
          }
)


#' Get time intervals at which the animal's head direction is in a given range
#' 
#' There is a tric to deal with the circularity of the head direction data
#' If hdMin is larger than hdMax, then it is assumed that you want the range from hdMin-360 and 0-hdMax
#' This wrap around the 360-0.
#' 
#' @param pt Positrack object
#' @param hdMin Minimal head direction that will be considered (in degree, from 0 to 360)
#' @param hdMax Maximal head direction that will be considered (in degree, from 0 to 360)
#' @return matrix with the time intervals
#'
#' @docType methods
#' @rdname getIntervalsAtHeadDirection-methods
setGeneric(name="getIntervalsAtHeadDirection",
           def=function(pt,hdMin,hdMax)
           {standardGeneric("getIntervalsAtHeadDirection")}
)
#' @rdname getIntervalsAtHeadDirection-methods
#' @aliases getIntervalsAtHeadDirection,ANY,ANY-method
setMethod(f="getIntervalsAtHeadDirection",
          signature="Positrack",
          definition=function(pt,hdMin,hdMax)
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@x)==0)
              stop("pt@x has a length of 0")
            if(length(pt@hd)==0)
              stop("pt@hd has a length of 0")
            if(hdMin<0)
              stop("hdMin < 0")
            if(hdMin>360)
              stop("hdMin > 360")
            if(hdMax<0)
              stop("hdMax < 0")
            if(hdMax>360)
              stop("hdMax > 360")

            x<-pt@hd
            x[is.na(x)]<- -1.0
            results<-.Call("head_direction_intervals_cwrap",
                           x,
                           length(x),
                           as.integer(pt@resSamplesPerWhlSample),
                           hdMin,
                           hdMax)
            colnames(results)<-c("start","end")
            return(results)
          }
)


#' Get the speed of the animal at given time points
#' 
#' The time points are in sample values of the electrophysiological data
#' 
#' @param pt Positrack object
#' @param res Time values in electrophysiological samples.
#' @param angular Logical, if true will use the angular velocity instead of the running speed
#' @return Vector with the speed at the different time points
#' 
#' @docType methods
#' @rdname getSpeedAtResValues-methods
setGeneric(name="getSpeedAtResValues",
           def=function(pt,res,angular=FALSE)
           {standardGeneric("getSpeedAtResValues")}
)
#' @rdname getSpeedAtResValues-methods
#' @aliases getSpeedAtResValues,ANY,ANY-method
setMethod(f="getSpeedAtResValues",
          signature="Positrack",
          definition=function(pt,res,angular=FALSE)
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@x)==0)
              stop("pt@x has length of 0")
            
            if(any(res<0))
              stop("some res<0")
            
            if(angular==FALSE){
              results<-.Call("speed_at_res_values_cwrap", 
                             pt@speed,
                             length(pt@speed),
                             as.integer(res),
                             length(res),
                             as.integer(pt@resSamplesPerWhlSample))
            }
            else{
              print("angular")
              results<-.Call("speed_at_res_values_cwrap", 
                             pt@angularSpeed,
                             length(pt@angularSpeed),
                             as.integer(res),
                             length(res),
                             as.integer(pt@resSamplesPerWhlSample))
            }
            return(results)
          }
)


#' Get the x, y and hd of the animal at given time points
#' 
#' The time points are in sample values of the electrophysiological data
#' 
#' @param pt Positrack object
#' @param res Time values in electrophysiological samples.
#' @return Dataframe with the x, y and hd at the different time points
#' 
#' @docType methods
#' @rdname getXYHDAtResValues-methods
setGeneric(name="getXYHDAtResValues",
           def=function(pt,res)
           {standardGeneric("getXYHDAtResValues")}
)
#' @rdname getXYHDAtResValues-methods
#' @aliases getXYHDAtResValues,ANY,ANY-method
setMethod(f="getXYHDAtResValues",
          signature="Positrack",
          definition=function(pt,res)
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@x)==0)
              stop("pt@x has length of 0")
            
            if(any(res<0))
              stop("some res<0")
            
            x<-.Call("speed_at_res_values_cwrap", 
                             pt@x,
                             length(pt@x),
                             as.integer(res),
                             length(res),
                             as.integer(pt@resSamplesPerWhlSample))
            y<-.Call("speed_at_res_values_cwrap", 
                     pt@y,
                     length(pt@y),
                     as.integer(res),
                     length(res),
                     as.integer(pt@resSamplesPerWhlSample))
            
            hd<-.Call("spike_head_direction_cwrap",
                      pt@hd,
                      length(pt@hd),
                      as.integer(res),
                      length(res),
                      as.integer(pt@resSamplesPerWhlSample),
                      as.integer(0),
                      as.integer(tail(res,n = 1)+1),
                      as.integer(1))
            
          return(data.frame(x=x,y=y,hd=hd))
      }
  )














#' Shift the position of the animal of a random amount larger than the value of minShiftMs
#' 
#' This is used to do shuffling. The position values (x and y) are shifted in time. 
#' Only valid values are shifted, so if you want to only shift the position only in 
#' one environment or a specific part of the track, first call setInvalidOutsideInterval
#' 
#' @param pt Positrack object
#' @return Positrack object with shifted x and y values
#' 
#' @docType methods
#' @rdname shiftPositionRandom-methods
setGeneric(name="shiftPositionRandom",
           def=function(pt)
           {standardGeneric("shiftPositionRandom")}
)
#' @rdname shiftPositionRandom-methods
#' @aliases shiftPositionRandom,ANY,ANY-method
setMethod(f="shiftPositionRandom",
          signature="Positrack",
          definition=function(pt)
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@x)==0)
              stop("pt@x has length of 0")
            ## only consider valid data points
            xx<-pt@x[!is.na(pt@x)]
            yy<-pt@y[!is.na(pt@x)]
            results<-shiftPositionVectors(x=xx,y=yy,
                                     timePerSampleRes=pt@resSamplesPerWhlSample,
                                     pt@minShiftMs,
                                     pt@samplingRateDat)
            pt@x[!is.na(pt@x)]<-results[[1]]
            pt@y[!is.na(pt@x)]<-results[[2]]
            return(pt)
          }
)



#' Shift the speed of the animal of a random amount larger than the value of minShiftMs
#' 
#' This is used to do shuffling. The speed values are shifted in time. 
#' Only valid values are shifted.
#' 
#' @param pt Positrack object
#' @return Positrack object with shifted speed values
#' 
#' @docType methods
#' @rdname shiftSpeedRandom-methods
setGeneric(name="shiftSpeedRandom",
           def=function(pt)
           {standardGeneric("shiftSpeedRandom")}
)
#' @rdname shiftSpeedRandom-methods
#' @aliases shiftSpeedRandom,ANY,ANY-method
setMethod(f="shiftSpeedRandom",
          signature="Positrack",
          definition=function(pt)
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@speed)==0)
              stop("pt@speed has length of 0")
            ## only consider valid data points
            index<-!is.na(pt@speed)
            pt@speed[index]<-shiftPositionVector(x=pt@speed[index],
                                     timePerSampleRes=pt@resSamplesPerWhlSample,
                                     pt@minShiftMs,
                                     pt@samplingRateDat)
            return(pt)
          }
)


#' Shift the head direction data of the animal of a random amount larger than the value of minShiftMs
#' 
#' This is used to do shuffling. The hd values are shifted in time. 
#' Only valid values are shifted.
#' 
#' @param pt Positrack object
#' @return Positrack object with shifted head direction (hd) values
#' 
#' @docType methods
#' @rdname shiftHdRandom-methods
setGeneric(name="shiftHdRandom",
           def=function(pt)
           {standardGeneric("shiftHdRandom")}
)
#' @rdname shiftHdRandom-methods
#' @aliases shiftHdRandom,ANY,ANY-method
setMethod(f="shiftHdRandom",
          signature="Positrack",
          definition=function(pt)
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@hd)==0)
              stop("pt@hd has length of 0")
            ## only consider valid data points
            index<-!is.na(pt@hd)
            pt@hd[index]<-shiftPositionVector(x=pt@hd[index],
                                               timePerSampleRes=pt@resSamplesPerWhlSample,
                                               pt@minShiftMs,
                                               pt@samplingRateDat)
            return(pt)
          }
)



#' Shift the linear position data of the animal of a random amount larger than the value of minShiftMs
#' 
#' This is used to do shuffling. The lin values are shifted in time. 
#' Only valid values are shifted.
#' 
#' @param pt Positrack object
#' @return Positrack object with shifted linear position (lin) values
#' 
#' @docType methods
#' @rdname shiftLinRandom-methods
setGeneric(name="shiftLinRandom",
           def=function(pt)
           {standardGeneric("shiftLinRandom")}
)
#' @rdname shiftLinRandom-methods
#' @aliases shiftLinRandom,ANY,ANY-method
setMethod(f="shiftLinRandom",
          signature="Positrack",
          definition=function(pt)
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@lin)==0)
              stop("pt@lin has length of 0")
            ## only consider valid data points
            index<-!is.na(pt@lin)
            pt@lin[index]<-shiftPositionVector(x=pt@lin[index],
                                         timePerSampleRes=pt@resSamplesPerWhlSample,
                                         pt@minShiftMs,
                                         pt@samplingRateDat)
            return(pt)
          }
)


#' Linearize the position data
#' 
#' This is usually done if the animal ran on a linear track. 
#' The regression line of the position data is calculated and the 
#' closest point on the line is found for each 2D coordinate.
#' 
#' @param pt Positrack object
#' @return Positrack object with linear position (lin)
#' 
#' @docType methods
#' @rdname linearizeLinearTrack-methods
setGeneric(name="linearzeLinearTrack",
           def=function(pt)
           {standardGeneric("linearzeLinearTrack")}
)
#' @rdname linearizeLinearTrack-methods
#' @aliases linearizeLinearTrack,ANY,ANY-method
setMethod(f="linearzeLinearTrack",
          signature="Positrack",
          definition=function(pt)
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@x)==0)
              stop("pt@x has length of 0")
            
            # only data from linear track should be valid when calling this        
            x<-pt@x
            y<-pt@y
            x[which(x==-1.0)]<-NA
            y[which(y==-1.0)]<-NA
            #plot(x,y)
            
            ###
            ### get the regression lines that best fit the maze
            ###
            ## for some reasons I don't understand today,
            ## I can only get lm to work if I swap x and y
            ## would need to understand this point
            lm1<-lm(x~y)
         #   plot(y,x,xlim=c(0,90),ylim=c(0,90))
         #   abline(lm1,col="red")
    
            ###	find the closest point on the line for each position data
            ### some high school stuff
            ### https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
            ### y = mx + b  
            ### 
            ### b<-lm1$coefficients[1]
            ### m<-lm1$coefficients[2]
            ### rewrite as Ax + By + C = 0; standard form
            ###
            ### -mx + y - b = 0
            A=0-lm1$coefficients[2]
            B=1
            C=0-lm1$coefficients[1]
            ###
            ### given point x0,y0
            ### point on this line which is closest to (x0,y0) has coordinates
            ### x = B(Bx0-Ay0)-AC/(A^2+B^2)
            ### y = A(-Bx0 +Ay0) - BC/ (A^2 + B^2)
            ###
            x.on.line<-function(x,y,A,B,C){
              return( (B*(B*x-A*y)-A*C)/(A^2+B^2))
            }
            y.on.line<-function(x,y,A,B,C){
              return( (A*(-B*x+A*y)-B*C)/(A^2+B^2))
            }
            xl<-x.on.line(y,x,A,B,C)
            yl<-y.on.line(y,x,A,B,C)
            #points(xl,yl,col="blue")
            ### now use the points on the line to find the 
            ### distance between the point with the smallest x and all other points
            ### that is our transformation from 2d to 1d
            X<-xl[which.min(xl)]
            Y<-yl[which.min(xl)]
            #points(X,Y,col="red")
            distance<-function(x1,y1,x2,y2){
              return(sqrt((x1-x2)^2+(y1-y2)^2))
            }
            lin<-distance(xl,yl,X,Y)
            #points(xl,lin,col="green")
            pt@lin<-lin
            pt@lin[is.na(pt@lin)]<- -1.0
            ### smooth the 1d array
            pt@lin<-smoothGaussian(x=pt@lin,sd=5,invalid=-1.0)
            pt@lin[which(pt@lin==-1.0)]<-NA
            
            ### get the direction T or F
            pt@dir<-rep(0,length(pt@x))
            pt@dir[2:length(pt@dir)]<-ifelse(diff(pt@lin)>0,1,0)
            pt@dir[which(pt@lin==-1.0)]<-NA
            #plot(pt@lin,xlim=c(240000,245000),type='l')
            #lines(ifelse(pt@dir>0,50,0),xlim=c(240000,245000),col="red")
            
            return(pt)
          }
)


### show ###
setMethod("show", "Positrack",
          function(object){
            print(paste("session:",object@session))
            print(paste("path:",object@path))
            if(length(object@samplingRateDat)==0){
              print(paste("Oject not yet loaded"))
              print("call loadPositrack")
              return()
            }
            print(paste("samplingRate:",object@samplingRateDat,"Hz"))
            print(paste("resSamplesPerWhlSample:",object@resSamplesPerWhlSample))
            print(paste("pxPerCm:",object@pxPerCm))
            print(paste("number data points:", length(object@xWhl)))
            print(paste("time in data:", (length(object@xWhl)*object@resSamplesPerWhlSample)/object@samplingRateDat,"sec"))
            print(paste("min max x:",min(object@x,na.rm=T),"cm,", max(object@x,na.rm=T),"cm"))
            print(paste("min max y:",min(object@y,na.rm=T),"cm,", max(object@y,na.rm=T),"cm"))
            print(paste("total valid time:", 
                        sum(!is.na(object@x))*object@resSamplesPerWhlSample/object@samplingRateDat,
                        "sec"))
            print(paste("percent of invalid position:", round(object@percentInvalidPosition)))
          })
