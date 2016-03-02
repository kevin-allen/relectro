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
            x.spikes="numeric",
            y.spikes="numeric",
            maps="array",
            occupancy="matrix",
            spatial.autocorrelation="numeric",
            cell.list="numeric",
            reduce.size="logical"),
  prototype = list(session="",cm.per.bin=2,smooth.occupancy.sd=3,smooth.rate.map.sd=3,reduce.size=T))


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
            sp@maps<-array(data=results,dim=(c(sp@nrow.map,sp@ncol.map,length(sp@cell.list))), )
            return(sp)
          }
)
