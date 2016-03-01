make.pairs<-function(cl1="",cl2=NULL){
  if(is.null(cl2)){
    m<-combn(cl1,m=2)
    data.frame(Var1=m[1,],Var2=m[2,])
  }
  else{
    expand.grid(1:3,4:5)
  }
}

smooth.gaussian<-function(x,sd=2,invalid=-1.0)
{
  if(length(x)==0)
    return
  if(sd==0)
    return
  if(sd<0)
    stop(paste("sd is smaller than 0:",sd))
  results<- .Call("smooth_double_gaussian_cwrap",
                    x, length(x), sd, invalid)
  return(results)
}
smooth.gaussian.degrees<-function(x,sd=2,invalid=-1.0)
{
  if(length(x)==0)
    return
  if(sd==0)
    return
  if(sd<0)
    stop(paste("sd is smaller than 0:",sd))
  if(any(x>360))
    stop(paste("x values larger than 360"))
  if(any(x<0&x!=-1.0))
    stop(paste("x values < 0 & != -1.0"))
  results<- .Call("smooth_double_gaussian_degrees_cwrap",
                  x, length(x), sd, invalid)
  return(results)
}
