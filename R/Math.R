makePairs<-function(cl1="",cl2=NULL){
  if(is.null(cl2)){
    m<-combn(cl1,m=2)
    data.frame(Var1=m[1,],Var2=m[2,])
  }
  else{
    expand.grid(1:3,4:5)
  }
}

smoothGaussian<-function(x,sd=2,invalid=-1.0)
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
smoothGaussianDegrees<-function(x,sd=2,invalid=-1.0)
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
shift <- function (v, places, dir = "right") 
{# shift a vector in place
  vnew <- NULL
  d <- substring(dir, 1, 1)
  if (d == "r" & places < 0) {
    d <- "l"
  }
  else {
    if (d == "l" & places < 0) {
      d <- "r"
    }
  }
  n <- length(v)
  p <- abs(places)
  if (p == 0) {
    vnew <- v
  }
  else {
    if (d == "r") {
      vnew <- c(v[(n - p + 1):n], v[1:(n - p)])
    }
    else {
      vnew <- c(v[(p + 1):n], v[1:p])
    }
  }
  return(vnew)
}

shuffle.vector<-function(x,
                         time.per.sample.res,
                         min.mv.ms,
                         sampling_rate){
  min.mv<- min.mv.ms*(sampling_rate/1000)/time.per.sample.res
  mv<-sample(min.mv:length(x)-min.mv,1)
  return(shift(x,mv))
}

shuffle.vectors<-function(x,y,
                         time.per.sample.res,
                         min.mv.ms,
                         sampling_rate){
  min.mv<- min.mv.ms*(sampling_rate/1000)/time.per.sample.res
  mv<-sample(min.mv:length(x)-min.mv,1)
  x<-shift(x,mv)
  y<-shift(y,mv)
  return(list(x=x,y=y))
}