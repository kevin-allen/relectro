############################################
#### definition of CellGroup Class      ###
############################################
CellGroup <- setClass(
  "CellGroup", ## name of the class
  slots=c(session="character",
          path="character", 
          fileBase="character", # path+session
          nTetrodes="numeric",
          nCells="numeric",
          id="character",
          animal="numeric",
          clu="numeric",
          tetrode="numeric",
          tetrodeId="character",
          cluToTetrode="numeric",
          brainRegion="character"),  # cell list to limit the analysis to these cells
  prototype = list(session="",path="",nTetrodes=0))



### show ###
setMethod("show", "CellGroup",
          function(object){
            print(paste("session:",object@session))
            print(paste("path:",object@path))
            print(paste("nTetrodes:",object@nTetrodes))
            print(paste("nCells:",object@nCells))
            print("clu:")
            print(object@clu)
            print("id:")
            print(object@id)
            print("tetrode:")
            print(object@tetrode)
            print("cluToTetrode:")
            print(object@cluToTetrode)
            print("brainRegion")
            print(object@brainRegion)
            print("tetrodeId:")
            print(object@tetrodeId)
          })



### loadCellGroup ###
setGeneric(name="loadCellGroup",
           def=function(cg)
           {standardGeneric("loadCellGroup")}
)
setMethod(f="loadCellGroup",
          signature="CellGroup",
          definition=function(cg)
          {
            if(cg@session=="")
              stop("cg@session is empty")
            if(cg@path=="")
              cg@path=getwd()
            cg@fileBase=paste(cg@path,cg@session,sep="/")
            if(cg@nTetrodes==0)
              stop("cg@nTetrodes==0")
            if(!file.exists(paste(cg@fileBase,"clu",sep="."))) # need the main clu file
              stop("needs ",paste(cg@fileBase,"clu",sep=".")) 
            
            cg@nCells<-as.numeric(readLines(con=paste(cg@fileBase,"clu",sep="."),n=1))-1
            # clu 1 is noise, ignore
            cg@clu=2:(cg@nCells+1)
            cg@id=paste(cg@session,cg@clu,sep="_")
            
            cg@tetrode=vector("numeric")        
            cg@cluToTetrode=vector("numeric")        
            for (t in 1:cg@nTetrodes){
              n=as.numeric(readLines(con=paste(cg@fileBase,"clu",t,sep="."),n=1))-1
              if(n>0){
                cg@tetrode<-c(cg@tetrode,rep(t,n))
                cg@cluToTetrode<-c(cg@cluToTetrode,2:(n+1))
              }
            }
            if(length(cg@cluToTetrode)!=length(cg@clu))
            {
              print("problem with the number of clusters in main clu files and in tetrode clu files")  
              print(paste(cg@session,":",length(cg@clu),"vs",length(cg@cluToTetrode)))
              stop()
            }
            cg@tetrodeId<-paste(cg@session,cg@tetrode,sep="_")
            if(file.exists(paste(cg@fileBase,"desel",sep="."))){
              br<-readLines(con=paste(cg@fileBase,"desel",sep="."))
              cg@brainRegion<-br[cg@tetrode]          
            }
            return(cg)
          }
)
