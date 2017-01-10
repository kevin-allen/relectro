#' An S4 class representing a group of cells
#'
#' This class gets you the tetrode id and brain region associated with each cell of a recording session.
#' It also get the cluster id of the cell on its respective tetrode.
#'
#' @slot session A character vector containing the names of the recording session.
#' @slot path The directory in which the files of the session are located.
#' @slot fileBase Is the path and session
#' @slot nTetrodes Number of tetrodes
#' @slot nCells Number of cells
#' @slot id Cell id
#' @slot animal Name of the animal
#' @slot clu Cluster id of the cells
#' @slot tetrode tetrode number of the cells
#' @slot tetrodeId id of the tetrode for each cell
#' @slot cluToTetrode clu id of the cell on its respective tetrode (used for clustering)
#' @slot brainRegion region in which the cell was recorded
#' @examples
#' df<-new("CellGroup")
#'
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
          brainRegion="character"),
  prototype = list(session="",path="",nTetrodes=0))

#' Load the information regarding a group of cells
#'
#' The object will get the tetrode id and
#' cluster number on the tetrode for each cluster.
#' The object also load the brain region associated with each cell.
#' This infomration is read from the .desel file and the clu file
#' of each tetrode.
#'
#' @param cg A CellGroup object
#' @return A CellGroup object with the information loaded
#'
#' @docType methods
#' @rdname loadCellGroup-methods
setGeneric(name="loadCellGroup",
           def=function(cg)
           {standardGeneric("loadCellGroup")}
)

#' @rdname loadCellGroup-methods
#' @aliases loadCellGroup,ANY,ANY-method
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
              stop(paste("cg@nTetrodes==0", "consider setting it when creating the CellGroup object"))
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





#' Get the brain region from a list of cluNo
#'
#'
#' @param cg A CellGroup object
#' @param cluNo Numeric containing the cluNo of the cluster for which you want the brain region.
#' @return Character vector with be name of the brain region for each cluNo
#'
#' @docType methods
#' @rdname brainRegionFromCluNo-methods
setGeneric(name="brainRegionFromCluNo",
           def=function(cg,cluNo)
           {standardGeneric("brainRegionFromCluNo")}
)
#' @rdname brainRegionFromCluNo-methods
#' @aliases brainRegionFromCluNo,ANY,ANY-method
setMethod(f="brainRegionFromCluNo",
          signature="CellGroup",
          definition=function(cg,cluNo)
          {
            if(cg@session=="")
              stop("cg@session is empty")
            if(length(cluNo)==0)
              return()
            if(any(!cluNo%in%cg@clu)){
              stop("cluNo contains clusters that are not in cg object")
            }
            ## use merge to match the brainRegion based on cluNo
            df1<-data.frame(clu=cg@clu,brainRegion=cg@brainRegion)
            df2<-data.frame(clu=cluNo)
            df3<-merge(df1,df2,by="clu")
            return(as.character(df3$brainRegion))
          }
)



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
