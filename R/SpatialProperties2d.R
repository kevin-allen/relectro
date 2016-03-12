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
            autos="array",
            spatial.autocorrelation="numeric",
            cell.list="numeric",
            reduce.size="logical",
            min.valid.bins.auto="numeric",
            peak.rates="numeric",
            info.score="numeric",
            border.score="numeric",
            border.cm="numeric",
            border.dm="numeric",
            border.num.field.detected="numeric",
            sparsity="numeric",
            grid.score="numeric",
            grid.orientation="numeric",
            grid.spacing="numeric",
            grid.score.number.fields.to.detect="numeric",
            grid.score.min.num.bins.per.field="numeric",
            grid.score.field.threshold="numeric",
            border.percentage.threshold.field="numeric",
            border.min.bins.in.field="numeric"),
      
  prototype = list(session="",cm.per.bin=2,smooth.occupancy.sd=3,smooth.rate.map.sd=3,min.valid.bins.auto=20,reduce.size=T,
                   grid.score.number.fields.to.detect=40,grid.score.min.num.bins.per.field=10,grid.score.field.threshold=0.1,
                   border.percentage.threshold.field=20,border.min.bins.in.field=10))


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
            sp@autos<-array(data=results,dim=(c(sp@nrow.auto,sp@ncol.auto,length(sp@cell.list))))
            return(sp)
          }
)



setGeneric(name="grid.score",
           def=function(sp)
           {standardGeneric("grid.score")}
)
setMethod(f="grid.score",
          signature="SpatialProperties2d",
          definition=function(sp)
          {
            if(length(sp@autos)==0)
              stop("Need to call map.spatial.autocorrelation first to run grid.score()")
            
            dyn.load("~/repo/r_packages/relectro/src/relectro.so")
            
            sp@grid.score<-.Call("grid_score_cwrap",
                  sp@cell.list,
                  length(sp@cell.list),
                  sp@autos,
                  sp@ncol.auto,
                  sp@nrow.auto,
                  sp@cm.per.bin,
                  sp@grid.score.number.fields.to.detect,
                  sp@grid.score.min.num.bins.per.field,
                  sp@grid.score.field.threshold,
                  -1.0) #invalid
            return(sp)
          }
)
setGeneric(name="grid.orientation",
           def=function(sp)
           {standardGeneric("grid.orientation")}
)
setMethod(f="grid.orientation",
          signature="SpatialProperties2d",
          definition=function(sp)
          {
            if(length(sp@autos)==0)
              stop("Need to call map.spatial.autocorrelation first to run grid.orientation()")
            

            sp@grid.orientation<- .Call("grid_orientation_cwrap",
                                 sp@cell.list,
                                 length(sp@cell.list),
                                 sp@autos,
                                 sp@ncol.auto,
                                 sp@nrow.auto,
                                 sp@cm.per.bin,
                                 sp@grid.score.number.fields.to.detect,
                                 sp@grid.score.min.num.bins.per.field,
                                 sp@grid.score.field.threshold,
                                 -1.0)
            return(sp)
          }
)
setGeneric(name="grid.spacing",
           def=function(sp)
           {standardGeneric("grid.spacing")}
)
setMethod(f="grid.spacing",
          signature="SpatialProperties2d",
          definition=function(sp)
          {
            if(length(sp@autos)==0)
              stop("Need to call map.spatial.autocorrelation first to run grid.spacing()")
            
            sp@grid.spacing<- .Call("grid_spacing_cwrap",
                                        sp@cell.list,
                                        length(sp@cell.list),
                                        sp@autos,
                                        sp@ncol.auto,
                                        sp@nrow.auto,
                                        sp@cm.per.bin,
                                        sp@grid.score.number.fields.to.detect,
                                        sp@grid.score.min.num.bins.per.field,
                                        sp@grid.score.field.threshold,
                                        -1.0)
            return(sp)
          }
)


setGeneric(name="border.score",
           def=function(sp)
           {standardGeneric("border.score")}
)
setMethod(f="border.score",
          signature="SpatialProperties2d",
          definition=function(sp)
          {
            if(length(sp@maps)==0)
              stop("Need to call firing.rate.map.2d first to run border.score()")
            
            results<-.Call("border_score_rectangular_environment_cwrap",
                    sp@cell.list,
                    length(sp@cell.list),
                    sp@ncol.map,
                    sp@nrow.map,
                    sp@occupancy,
                    sp@maps,
                    sp@border.percentage.threshold.field,
                    sp@border.min.bins.in.field)
            
            sp@border.score<-results[1,]
            sp@border.cm<-results[2,]
            sp@border.dm<-results[3,]
            sp@border.num.field.detected<-results[4,]
            return(sp)
          }
)





