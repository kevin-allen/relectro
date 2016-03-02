############################################
#### definition of Positrack Class      ###
############################################
Positrack <- setClass(
  "Positrack", ## name of the class
  slots=c(session="character",
          px.per.cm="numeric",
          samplingRateDat="numeric",
          res.samples.per.whl.sample="numeric",
          x.whl="numeric",
          y.whl="numeric",
          hd.whd="numeric",
          x="numeric", # in cm
          y="numeric",
          hd="numeric",
          speed="numeric",
          angular.speed="numeric",
          res="numeric",
          default.xy.smoothing="numeric"
          ),  # cell list to limit the analysis to these cells
  prototype = list(session="",default.xy.smoothing=2))

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
            if(!file.exists(paste(pt@session,"whl",sep=".")))
              stop("need",paste(pt@session,"whl",sep="."))
            if(!file.exists(paste(pt@session,"whd",sep=".")))
              stop("need",paste(pt@session,"whd",sep="."))
            if(!file.exists(paste(pt@session,"res_samples_per_whl_sample",sep=".")))
              stop("need",paste(pt@session,"res_samples_per_whl_sample",sep="."))
            if(!file.exists(paste(pt@session,"sampling_rate_dat",sep=".")))
              stop("need",paste(pt@session,"sampling_rate_dat",sep="."))
            if(!file.exists(paste(pt@session,"px_per_cm",sep=".")))
              stop("need",paste(pt@session,"px_per_cm",sep="."))
            
            
            ## get sampling rate
            whl<-read.table(paste(pt@session,"whl",sep="."))
            pt@x.whl<-whl$V1
            pt@y.whl<-whl$V2
            whd<-read.table(paste(pt@session,"whd",sep="."))
            pt@hd.whd<-whd$V3
            pt@res.samples.per.whl.sample<-read.table(paste(pt@session,"res_samples_per_whl_sample",sep="."))$V1
            pt@px.per.cm<-read.table(paste(pt@session,"px_per_cm",sep="."))$V1
            pt@samplingRateDat<-read.table(paste(pt@session,"sampling_rate_dat",sep="."))$V1
            
            ## check that things add up
            if(length(pt@x.whl)!=length(pt@hd.whd))
              stop("Problem with length of whl and whd files")
            if(pt@samplingRateDat<1|pt@samplingRateDat>100000)
              stop(paste("pt@sampingRateDat is out of bound:",pt@samplingRateDat))
            if(pt@px.per.cm<1|pt@px.per.cm>10000)
              stop(paste("pt@px.per.cm is out of bound:", pt@px.per.cm))
            if(max(pt@x.whl)>2000|min(pt@x.whl)< -1.0){
              print(paste("min x:",min(pt@x.whl),"max x:",max(pt@x.whl)))
              stop(paste("values of pt@x.whl are out of bound"))
            }
            if(max(pt@y.whl)>2000|min(pt@y.whl)< -1.0){
              print(paste("min y:",min(pt@y.whl),"max y:",max(pt@y.whl)))
              stop(paste("values of pt@x.whl are out of bound"))
            }
            if(max(pt@hd.whd>360))
              stop(paste("max value of pt@hd.whd > 360:",max(pt@hd.whd)))
            if(min(pt@hd.whd< -1.0))
              stop(paste("min value of pt@hd.whd < -1.0:",min(pt@hd.whd)))
            ## pt@x.whl will never be changed, 
            pt@x<-pt@x.whl
            pt@y<-pt@y.whl
            pt@hd<-pt@hd.whd
            
            ## set the res value associated with eack whl sample
            pt@res<-seq(from=pt@res.samples.per.whl.sample,by=pt@res.samples.per.whl.sample,length.out=length(pt@x))
            
            ## apply a bit a smoothing to the x and y
            pt@x<-smooth.gaussian(x=pt@x,sd=pt@default.xy.smoothing,invalid=-1.0)
            pt@y<-smooth.gaussian(x=pt@y,sd=pt@default.xy.smoothing,invalid=-1.0)
         
            ## apply no smoothing to head direction
               
            ## get the speed from position
            pt@speed<- .Call("speed_from_whl_cwrap",
                            pt@x,
                            pt@y,
                            length(pt@x),
                            4, # look ahead
                            4, # look back
                            as.numeric(pt@px.per.cm),
                            pt@samplingRateDat, 
                            pt@res.samples.per.whl.sample)
         
            ## get the angular speed from position
            pt@angular.speed<- .Call("angular_speed_from_hd_cwrap",
                                     pt@hd,
                                     length(pt@hd),
                                     4,
                                     4,
                                     pt@samplingRateDat,
                                     pt@res.samples.per.whl.sample)
            
            # set -1.0 to NA for calculation in R       
            pt@x[which(pt@x==-1.0)]<-NA
            pt@y[which(pt@y==-1.0)]<-NA
            pt@hd[which(pt@hd==-1.0)]<-NA
            
            # get cm
            pt@x<-pt@x/pt@px.per.cm
            pt@y<-pt@y/pt@px.per.cm
                        
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
            pt@x<-smooth.gaussian(pt@x,sd=sd,invalid=-1.0)
            pt@y<-smooth.gaussian(pt@y,sd=sd,invalid=-1.0)
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
            pt@hd<-smooth.gaussian(pt@hd,sd=sd,invalid=-1.0)
            pt@hd[which(pt@hd==-1.0)]<- NA
            return(pt)
          }
)

### filter for speed ###
setGeneric(name="speed.filter",
           def=function(pt,min.speed,max.speed)
           {standardGeneric("speed.filter")}
)
setMethod(f="speed.filter",
          signature="Positrack",
          definition=function(pt,min.speed=3,max.speed=100)
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@x)==0)
              stop("pt@x has length of 0")
            
            ## set time when outside of speed to NA
            pt@x[which(pt@speed<min.speed|pt@speed>max.speed)]<- NA
            pt@y[which(pt@speed<min.speed|pt@speed>max.speed)]<- NA
            pt@hd[which(pt@speed<min.speed|pt@speed>max.speed)]<- NA
            return(pt)
          }
)

#### set.invalid.outside.interval
setGeneric(name="set.invalid.outside.interval",
           def=function(pt,s,e)
           {standardGeneric("set.invalid.outside.interval")}
)
setMethod(f="set.invalid.outside.interval",
          signature="Positrack",
          definition=function(pt,s,e)
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@x)==0)
              stop("pt@x has length of 0")

            ## if s is a matrix, then e is ignored
            if(class(s)=="matrix"){
              start.interval<-as.numeric(s[,1])
              end.interval<-as.numeric(s[,2])
            }else{
              start.interval<-s
              end.interval<-e
            }
            if(length(start.interval)!=length(end.interval))
              stop("problem with length of start.interval and end.interval")
            if(any(start.interval>end.interval))
              stop("start.interval>end.interval")
            if(any(diff(start.interval)<0))
              stop("problem with chonology within start.interval")
            if(any(diff(end.interval)<0))
              stop("problem with chonology within end.interval")
            if(any(start.interval[-1]-end.interval[-length(end.interval)]<0))
              stop("problem with chronology between intervals, from end to next start")

            ## slow, 1 sec for a list of 2 intervals
            #isin<-sapply(pt@res,function(x)any(x>=start.interval&x<=end.interval)) # takes seconds
            isin<-as.logical(.Call("resWithinIntervals", ## about 1000 times faster
                  length(start.interval),
                  as.integer(start.interval),
                  as.integer(end.interval),
                  length(pt@res),
                  as.integer(pt@res)))
            pt@x[which(!isin)]<-NA
            pt@y[which(!isin)]<-NA
            pt@hd[which(!isin)]<-NA
            return(pt)
          }
)

#### get.intervals.at.speed
setGeneric(name="get.intervals.at.speed",
           def=function(pt,min.speed,max.speed)
           {standardGeneric("get.intervals.at.speed")}
)
setMethod(f="get.intervals.at.speed",
          signature="Positrack",
          definition=function(pt,min.speed,max.speed)
          {
            if(pt@session=="")
              stop("pt@session is empty")
            if(length(pt@x)==0)
              stop("pt@x has length of 0")
            
            if(min.speed>max.speed)
              stop("min.speed>max.speed")
            
            
            results<-.Call("speed_intervals_cwrap", 
                  pt@speed,
                  length(pt@speed),
                  pt@res.samples.per.whl.sample,
                  min.speed,
                  max.speed)
            
            rownames(results)<-c("start","end")
            return(results)
          }
)





### show ###
setMethod("show", "Positrack",
          function(object){
            print(paste("session:",object@session))
            if(length(object@samplingRateDat)==0){
              print(paste("Oject not yet loaded"))
              print("call loadPositrack")
              return()
            }
              
            print(paste("samplingRate:",object@samplingRateDat,"Hz"))
            print(paste("res.samples.per.whl.sample:",object@res.samples.per.whl.sample))
            print(paste("px.per.cm:",object@px.per.cm))
            print(paste("number data points:", length(object@x.whl)))
            print(paste("time in data:", (length(object@x.whl)*object@res.samples.per.whl.sample)/object@samplingRateDat,"sec"))
            print(paste("min max x:",min(object@x,na.rm=T),"cm,", max(object@x,na.rm=T),"cm"))
            print(paste("min max y:",min(object@y,na.rm=T),"cm,",max(object@y,na.rm=T),"cm"))
          })




#### get.speed.at.res.values
setGeneric(name="get.speed.at.res.values",
           def=function(pt,res)
           {standardGeneric("get.speed.at.res.values")}
)
setMethod(f="get.speed.at.res.values",
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
                           as.integer(pt@res.samples.per.whl.sample))
            
            return(results)
          }
)


