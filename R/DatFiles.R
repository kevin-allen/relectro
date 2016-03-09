############################################
#### definition of DatFiles Class      ###
############################################
DatFiles <- setClass(
  "DatFiles", ## name of the class
  slots=c(file.names="character",
          resofs="numeric",
          num.channels="numeric"
          ),
  prototype = list(file.names="",num.channels=0))


### dat.files.set ###
setGeneric(name="dat.files.set",
           def=function(df,file.names,num.channels)
           {standardGeneric("dat.files.set")}
)
setMethod(f="dat.files.set",
          signature="DatFiles",
          definition=function(df,file.names="",num.channels=0)
          {
            if(length(file.names)==0)
              stop("file.names length = 0")
            if(num.channels<=0)
              stop("num.channels should be larger than 0")
            
            df@file.names=file.names
            df@num.channels=num.channels
            return(df)                        
          })

### dat.files.get.ups ###
setGeneric(name="dat.files.get.ups",
           def=function(df)
           {standardGeneric("dat.files.get.ups")}
)
setMethod(f="dat.files.get.ups",
          signature="DatFiles",
          definition=function(df)
          {
            if(length(df@file.names)==0)
              stop("df@file.names length = 0")
            if(num.channels<=0)
              stop("num.channels should be larger than 0")
            
            ## call c code that will read the dat files and return the up times
            
            
            
            
          })



### show ###
setMethod("show", "DatFiles",
          function(object){
            print(paste("file.names:"))
            print(paste(object@file.names))
            print(paste("num.channels:",object@num.channels))
          })


## get session information
rs<-new("RecSession")
setwd("/data/processing/jp4298/jp4298-12022016-0104/")
rs@session="jp4298-12022016-0104"
rs<-loadRecSession(rs)

## get data from files
df<-new("DatFiles")
setwd("/data/bindata/jp4298/jp4298-12022016-0104/")
df<-dat.files.set(df,file.names=paste(rs@trial.names,"dat",sep="."),num.channels=rs@n.channels)

dyn.load("~/repo/r_packages/relectro/src/relectro.so")

