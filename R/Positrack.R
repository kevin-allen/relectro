############################################
#### definition of Positrack Class      ###
############################################
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
          speed="numeric",
          angularSpeed="numeric",
          res="numeric",
          defaultXYSmoothing="numeric"
          ),  # cell list to limit the analysis to these cells
  prototype = list(session="",path="",defaultXYSmoothing=2))

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
            pathSession=paste(rs@path,rs@session,sep="/")
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
            pt@x[which(pt@speed<minSpeed|pt@speed>maxSpeed)]<- NA
            pt@y[which(pt@speed<minSpeed|pt@speed>maxSpeed)]<- NA
            pt@hd[which(pt@speed<minSpeed|pt@speed>maxSpeed)]<- NA
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
          definition=function(pt,s,e)
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
            print(paste("min max y:",min(object@y,na.rm=T),"cm,",max(object@y,na.rm=T),"cm"))
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


