###############################################
#### definition of SpikeTrain class         ###
###############################################
SpikeTrain <- setClass(
  "SpikeTrain", ## name of the class
  slots=c(session="character",samplingRate="numeric",res="numeric",time="numeric",clu="numeric",
          nCells="numeric",nSpikes="numeric",nSpikesPerCell="numeric",
          startInterval="numeric",endInterval="numeric", # to limit analysis to these intervals
          startResIndex="numeric",endResIndex="numeric",
          cellList="numeric"),  # cell list to limit the analysis to these cells
  prototype = list(session=""))


### spikeTimeAutocorrelation ###
# function using .Call() to call c wrapper
# It calls the the c function with the cwrap suffix
setGeneric(name="spikeTimeAutocorrelation",
           def=function(st,...)
             {standardGeneric("spikeTimeAutocorrelation")})

setMethod(f="spikeTimeAutocorrelation",
          signature = "SpikeTrain",
          definition=function(st,bin.size.ms=1,window.size.ms=200,probability=FALSE,...)
            {
            n.bins=(window.size.ms*2)/bin.size.ms
            window.size=2*window.size.ms*st@samplingRate/1000 # window size in res value from - to + extrems
            # call cwrapper function
            
            results<- .Call("autocorrelation_cwrap",
                        st@cellList,
                        length(st@cellList),
                        st@clu,
                        st@res,
                        st@nSpikes,
                        n.bins,
                        window.size,
                        st@startResIndex,
                        st@endResIndex,
                        length(st@startResIndex),
                        probability)
            
            # return a data fram with clu time count
            if(probability==F)
              data.frame(clu=rep(st@cellList,each=n.bins),
                         time=seq(-window.size.ms+bin.size.ms,window.size.ms,bin.size.ms)-(bin.size.ms/2),
                         count=results)
            else{
              data.frame(clu=rep(st@cellList,each=n.bins),
                         time=seq(-window.size.ms+bin.size.ms,window.size.ms,bin.size.ms)-(bin.size.ms/2),
                         prob=results)
            }
          }
          )

### loadSpikeTrain ###
setGeneric(name="loadSpikeTrain",
           def=function(st)
           {standardGeneric("loadSpikeTrain")}
)
setMethod(f="loadSpikeTrain",
          signature="SpikeTrain",
          definition=function(st)
          {
            st@res<-read.table(paste(st@session,"res",sep="."))$V1
            st@clu<-read.table(paste(st@session,"clu",sep="."))$V1
            st@samplingRate<-read.table(paste(st@session,"sampling_rate_dat",sep="."))$V1
            if(st@samplingRate<1 | st@samplingRate > 100000)
              stop(paste("samplingRate is out of range:",st@samplingRate))
            st@clu<-st@clu[-1] ## remove first number
            if(length(st@res)!=length(st@clu))
              stop(paste("length of res and clu files not equals:",length(st@res),length(st@clu)))
            # remove noise cluster 1 or 0
            df<-data.frame(res=st@res,clu=st@clu)
            df<-df[which(df$clu>1),]
            st@res<-df$res
            st@clu<-df$clu
            # get time in sec
            st@time<-st@res/st@samplingRate
            st@nSpikes<-length(st@res)
            st@nCells<-length(unique(st@clu))
            st@nSpikesPerCell<-as.numeric(table(st@clu))
            # by default analysis on all cells and all recording period
            st@cellList<-1:st@nCells
            st@startInterval<-0
            st@endInterval<-max(st@res)
            st@startResIndex<-0
            st@endResIndex<-length(st@res)
            return(st)
          }
)



### loadSpikeTrain ###
setGeneric(name="setSpikeTrain",
           def=function(st,...)
           {standardGeneric("setSpikeTrain")}
)


setMethod(f="setSpikeTrain",
          signature="SpikeTrain",
          definition=function(st,res,clu,sampling.rate)
          {
            st@res<-res
            st@clu<-clu
            st@samplingRate<-sampling.rate
          
            if(st@samplingRate<1 | st@samplingRate > 100000)
              stop(paste("samplingRate is out of range:",st@samplingRate))
          
            if(length(st@res)!=length(st@clu))
              stop(paste("length of res and clu files not equals:",length(st@res),length(st@clu)))
          
            # get time in sec
            st@time<-st@res/st@samplingRate
            st@nSpikes<-length(st@res)
            st@nCells<-length(unique(st@clu))
            st@nSpikesPerCell<-as.numeric(table(st@clu))
            # by default analysis on all cells and all recording period
            st@cellList<-1:st@nCells
            st@startInterval<-0
            st@endInterval<-max(st@res)
            st@startResIndex<-0
            st@endResIndex<-length(st@res)
            return(st)
          }
)


### show ###
setMethod("show", "SpikeTrain",
          function(object){
            print(paste("session:",object@session))
            print(paste("samplingRate:",object@samplingRate))
            print(paste("nCells:",object@nCells))
            print(paste("nSpikes:",object@nSpikes))
            print(paste("nSpikesPerCell:"))
            print(object@nSpikesPerCell)
            print(paste("cellList:"))
            print(object@cellList)
            print(paste("nIntervals:",length(object@startInterval))) 
                  })
