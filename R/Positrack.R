#' An S4 class representing the path of an animal during a recording session
#' 
#' This class is used to manipulate and represent the position data of the animal.
#' The first data point is at time resSamplesPerWhlSample and not 0.
#' The data are usually loaded from .whl, .whd, .res_samples_per_whl_sample, .sampling_rate_dat and .px_per_cm
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
#' @minShiftMs Minimum shift of the position vector during a shuffling procedure
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
          minShiftMs="numeric"
          ),  # cell list to limit the analysis to these cells
  prototype = list(session="",path="",defaultXYSmoothing=2,minShiftMs=20000))

### loadPositrack ###
setGeneric(name="loadPositrack",
           def=function(pt)
           {standardGeneric("loadPositrack")}
)
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
            if(!file.exists(paste(pathSession,"whl",sep=".")))
              stop("need",paste(pathSession,"whl",sep="."))
            if(!file.exists(paste(pathSession,"whd",sep=".")))
              stop("need",paste(pathSession,"whd",sep="."))
            if(!file.exists(paste(pathSession,"res_samples_per_whl_sample",sep=".")))
              stop("need",paste(pathSession,"res_samples_per_whl_sample",sep="."))
            if(!file.exists(paste(pathSession,"sampling_rate_dat",sep=".")))
              stop("need",paste(pathSession,"sampling_rate_dat",sep="."))
            if(!file.exists(paste(pathSession,"px_per_cm",sep=".")))
              stop("need",paste(pathSession,"px_per_cm",sep="."))
            ## get sampling rate
            whl<-read.table(paste(pathSession,"whl",sep="."))
            pt@xWhl<-whl$V1
            pt@yWhl<-whl$V2
            whd<-read.table(paste(pathSession,"whd",sep="."))
            pt@hdWhd<-whd$V3
            pt@resSamplesPerWhlSample<-read.table(paste(pathSession,"res_samples_per_whl_sample",sep="."))$V1
            pt@pxPerCm<-read.table(paste(pathSession,"px_per_cm",sep="."))$V1
            pt@samplingRateDat<-read.table(paste(pathSession,"sampling_rate_dat",sep="."))$V1
            
            ## check that things add up
            if(length(pt@xWhl)!=length(pt@hdWhd))
              stop("Problem with length of whl and whd files")
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
                            4, # look ahead
                            4, # look back
                            as.numeric(pt@pxPerCm),
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
            
            # get cm
            pt@x<-pt@x/pt@pxPerCm
            pt@y<-pt@y/pt@pxPerCm
                        
            return(pt)
          }
)


### smoothxy ###
setGeneric(name="smoothxy",
           def=function(pt,sd)
           {standardGeneric("smoothxy")}
)
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

### smoothhd ###
setGeneric(name="smoothhd",
           def=function(pt,sd)
           {standardGeneric("smoothhd")}
)
setMethod(f="smoothhd",
          signature="Positrack",
          definition=function(pt,sd=2)
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@x)==0)
              stop("pt@x has length of 0")
            ## smooth x and y
            pt@hd[which(is.na(pt@hd))]<- -1.0
            
            print("Smoothing of hd data with an function that does not take into account that data are circular")
            pt@hd<-smoothGaussian(pt@hd,sd=sd,invalid=-1.0)
            pt@hd[which(pt@hd==-1.0)]<- NA
            return(pt)
          }
)

### filter for speed ###
setGeneric(name="speedFilter",
           def=function(pt,minSpeed,maxSpeed)
           {standardGeneric("speedFilter")}
)
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

#### setInvalidOutsideInterval
setGeneric(name="setInvalidOutsideInterval",
           def=function(pt,s,e)
           {standardGeneric("setInvalidOutsideInterval")}
)
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
            return(pt)
          }
)

#### getIntervalsAtSpeed
setGeneric(name="getIntervalsAtSpeed",
           def=function(pt,minSpeed,maxSpeed)
           {standardGeneric("getIntervalsAtSpeed")}
)
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
            
            rownames(results)<-c("start","end")
            return(results)
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
          })

#### getSpeedAtResValues
setGeneric(name="getSpeedAtResValues",
           def=function(pt,res)
           {standardGeneric("getSpeedAtResValues")}
)
setMethod(f="getSpeedAtResValues",
          signature="Positrack",
          definition=function(pt,res)
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@x)==0)
              stop("pt@x has length of 0")
            
            if(any(res<0))
              stop("some res<0")
            
            results<-.Call("speed_at_res_values_cwrap", 
                           pt@speed,
                           length(pt@speed),
                           as.integer(res),
                           length(res),
                           as.integer(pt@resSamplesPerWhlSample))
            
            return(results)
          }
)



#### shiftPositionRandom
setGeneric(name="shiftPositionRandom",
           def=function(pt)
           {standardGeneric("shiftPositionRandom")}
)
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


#### shiftSpeedRandom
setGeneric(name="shiftSpeedRandom",
           def=function(pt)
           {standardGeneric("shiftSpeedRandom")}
)
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


#### shiftHdRandom
setGeneric(name="shiftHdRandom",
           def=function(pt)
           {standardGeneric("shiftHdRandom")}
)
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


#### shiftLinRandom
setGeneric(name="shiftLinRandom",
           def=function(pt)
           {standardGeneric("shiftLinRandom")}
)
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


#### linearzeLinearTrack
setGeneric(name="linearzeLinearTrack",
           def=function(pt)
           {standardGeneric("linearzeLinearTrack")}
)
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


