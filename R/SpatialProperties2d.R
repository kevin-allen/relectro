#################################################
#### definition of SpatialProperties2d Class  ###
#################################################
SpatialProperties2d<- setClass(
  "SpatialProperties2d", ## name of the class
    slots=c(session="character",
            cm.per.bin="numeric",
            smooth.occupancy.sd="numeric", ## in cm
            smooth.rate.map.sd="numeric", ## in cm
            ncol.map="integer",
            nrow.map="integer",
            ncol.auto="integer",
            nrow.auto="integer",
            x.spikes="numeric",
            y.spikes="numeric",
            maps="array",
            occupancy="matrix",
            auto="array",
            spatial.autocorrelation="numeric",
            cell.list="numeric",
            reduce.size="logical",
            min.valid.bins.auto="numeric",
            peak.rates="numeric",
            info.score="numeric",
            sparsity="numeric"),
      
  prototype = list(session="",cm.per.bin=2,smooth.occupancy.sd=3,smooth.rate.map.sd=3,min.valid.bins.auto=20,reduce.size=T))


### show ###
setMethod("show", "SpatialProperties2d",
          function(object){
            print(paste("session:",object@session))
            print(paste("cm.per.bin:",object@cm.per.bin))
            print(paste("smooth.occupancy.sd:",object@smooth.occupancy.sd))
            print(paste("smooth.rate.map.sd:",object@smooth.rate.map.sd))
            print(paste("reduce.size:",object@reduce.size))
            print(paste("ncol.map:",object@ncol.map))
            print(paste("nrow.map:",object@nrow.map))
            if(length(object@cell.list)!=0){
              print(paste("cell.list:"))
              print(object@cell.list)
            }
              
          })

#### make firing rate maps 
setGeneric(name="firing.rate.map.2d",
           def=function(sp,st,pt)
           {standardGeneric("firing.rate.map.2d")}
)
setMethod(f="firing.rate.map.2d",
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
            
            sp@cell.list<-st@cellList
            
            ## reduce the size of maps and map autocorrelation
            if(sp@reduce.size==T){
              x<-pt@x-min(pt@x,na.rm=T)+sp@cm.per.bin
              y<-pt@y-min(pt@y,na.rm=T)+sp@cm.per.bin
            }else{
              x<-pt@x
              y<-pt@y
            }
            
            #plot(x,y)
            ## use -1 as invalid values in c functions
            x[is.na(x)]<- -1
            y[is.na(y)]<- -1
            
            ## get the dimensions of the map
            sp@ncol.map=as.integer(((max(x)+1)/sp@cm.per.bin)+1) # x 
            sp@nrow.map=as.integer(((max(y)+1)/sp@cm.per.bin)+1) # y
            
            ## get spike position
            results<-.Call("spike_position_cwrap",
                  x,
                  y,
                  length(x),
                  as.integer(st@res),
                  st@nSpikes,
                  pt@res.samples.per.whl.sample,
                  as.integer(st@startInterval),
                  as.integer(st@endInterval),
                  length(st@startInterval))
            sp@x.spikes<-results[1,]
            sp@y.spikes<-results[2,]
            #plot(head(sp@x.spikes[which(sp@x.spikes!=-1.0)],n=20000),head(sp@y.spikes[which(sp@x.spikes!=-1.0)],n=20000))
            
            ## make the occupancy map
            sp@occupancy<-.Call("occupancy_map_cwrap",
                           sp@ncol.map,
                            sp@nrow.map,
                            sp@cm.per.bin,
                            sp@cm.per.bin,
                            x,
                            y,
                            length(x),
                            pt@res.samples.per.whl.sample/pt@samplingRateDat*1000, ## ms per whl samples
                            as.integer(st@startInterval),
                            as.integer(st@endInterval),
                            length(st@startInterval),
                            pt@res.samples.per.whl.sample)
            #image(t(sp@occupancy),zlim=c(0,max(sp@occupancy,na.rm=T)))
            
            ## smooth the occupancy map
            sp@occupancy<- .Call("smooth_double_gaussian_2d_cwrap",
                  as.numeric(sp@occupancy),
                  sp@ncol.map,
                  sp@nrow.map,
                  sp@smooth.occupancy.sd/sp@cm.per.bin,
                  -1.0)
            #image(t(sp@occupancy),zlim=c(0,max(sp@occupancy,na.rm=T)))
            
            ## make the 2d maps
            results<- .Call("firing_rate_map_2d_cwrap",
                            sp@ncol.map,
                            sp@nrow.map,
                            sp@cm.per.bin,
                            sp@cm.per.bin,
                            sp@x.spikes,
                            sp@y.spikes,
                            as.integer(st@clu),
                            st@nSpikes,
                            sp@cell.list,
                            length(sp@cell.list),
                            as.numeric(sp@occupancy),
                            sp@smooth.rate.map.sd/sp@cm.per.bin)
            sp@maps<-array(data=results,dim=(c(sp@nrow.map,sp@ncol.map,length(sp@cell.list))))
            
            
            return(sp)
          }
)



#### get.map.stats
setGeneric(name="get.map.stats",
           def=function(sp)
           {standardGeneric("get.map.stats")}
)
setMethod(f="get.map.stats",
          signature="SpatialProperties2d",
          definition=function(sp)
          {
            
            if(length(sp@maps)==0)
              stop("Need to call firing.rate.map.2d first to run get.map.stats")
            if(length(sp@occupancy)==0)
              stop("sp@occupancy length ==0")
          
            ### get peak rates
            sp@peak.rates<-apply(sp@maps,3,max)

            ### get info scores
            sp@info.score<- .Call("information_score_cwrap",
                as.integer(sp@cell.list),
                length(sp@cell.list),
                as.numeric(sp@maps),
                as.numeric(sp@occupancy),
                as.integer(sp@ncol.map*sp@nrow.map))
            
            ### get sparsity scores
            sp@sparsity<- .Call("sparsity_score_cwrap",
                                as.integer(sp@cell.list),
                                length(sp@cell.list),
                                as.numeric(sp@maps),
                                as.numeric(sp@occupancy),
                                as.integer(sp@ncol.map*sp@nrow.map))
            
          return(sp)
          }
)

setGeneric(name="map.spatial.autocorrelation",
           def=function(sp)
           {standardGeneric("map.spatial.autocorrelation")}
)
setMethod(f="map.spatial.autocorrelation",
          signature="SpatialProperties2d",
          definition=function(sp)
          {
            
            if(length(sp@maps)==0)
              stop("Need to call firing.rate.map.2d first to run get.map.stats")
            if(length(sp@occupancy)==0)
              stop("sp@occupancy length ==0")
            
            sp@ncol.auto = as.integer((sp@ncol.map*2)+1)
            sp@nrow.auto = as.integer((sp@nrow.map*2)+1)
            
            results<- .Call("map_autocorrelation_cwrap",
                  sp@cell.list,
                  length(sp@cell.list),
                  as.numeric(sp@maps),
                  sp@ncol.map,
                  sp@nrow.map,
                  sp@ncol.auto,
                  sp@nrow.auto,
                  sp@min.valid.bins.auto)
                  
            sp@auto<-array(data=results,dim=(c(sp@nrow.auto,sp@ncol.auto,length(sp@cell.list))))
            
            return(sp)
          }
)


dyn.load("~/repo/r_packages/relectro/src/relectro.so")

setwd("~/repo/r_packages/data")
session="jp4298-15022016-0106"
rs<-new("RecSession",session=session) ## info about rec session
rs<-loadRecSession(rs)
pt<-new("Positrack",session=session) ## info about position
pt<-loadPositrack(pt)
st<-new("SpikeTrain",session=session) ## info about spike trains
st<-loadSpikeTrain(st)

sp<-new("SpatialProperties2d",session=session) ## object to get spatial properties
pt<-set.invalid.outside.interval(pt,s=getIntervalsEnvironment(rs,env="sqr70")) ## select position data for one environment
sp<-firing.rate.map.2d(sp,st,pt) ## make firing rate maps
sp<-get.map.stats(sp)
sp<-map.spatial.autocorrelation(sp)

?benchmark

## plot one map
jet.colors = colorRampPalette(c("#00007F", "blue","#007FFF",  "cyan", "#7FFF7F", "yellow", "#FF7F00","red"))
map<-sp@auto[,,7]
image(t(map),zlim=c(-1,max(map,na.rm=T)), col=jet.colors(200),xlab='',ylab='',axes=FALSE)
