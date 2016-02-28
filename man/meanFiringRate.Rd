\name{meanFiringRate}
\alias{meanFiringRate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the mean firing rate
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
meanFiringRate(st)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{st}{
    A SpikeTrain object
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector containing the firing rate of neurons
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## generate spikes for 2 neurons  
res1<-cumsum(rpois(n=100,lambda=10))
res2<-cumsum(rpois(n=200,lambda=15))
clu<-c(rep(1,100),rep(2,200))
df<-data.frame(res=c(res1,res2),clu=clu)
df<-df[order(df$res),] # sort according to res values


## create SpikeTrain object from random spikes ###
st<-new("SpikeTrain")
## set the spike trains in the object
st<-setSpikeTrain(st=st,res=df$res,clu=df$clu,sampling.rate=20000)
## get the mean rate
rate<-meanFiringRate(st)
## plot the rate
plot(rate,type='l')
 
 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
