#' Create the directories of a recording session. 
#' 
#' This creates the directories require to add a recording session to the database.
#' It will create a directory for the subject if it does not exists.
#' A session directory is created within the subject directory
#' 
#' It can create a symbolic link to the subject directory
#' 
#' 
#' @param session Character vector giving the session name.
#' @param path Path where the directories will be created
#' @param linkTo A location where to put a symbolic link to the subject directory
databaseDirectories<-function(session,
                              path=NA,
                              link=NA)
  {
  if(rs@session=="")
    stop(paste("whdFromPositrack, rs@session == \"\""))
  
  if(is.na(rs@samplingRate))
    stop(paste("rs@samplingRate is NA"))
  if(file.exists(paste(rs@fileBase,"whd",sep="."))&overwrite==FALSE)
    return()
  
  
}


