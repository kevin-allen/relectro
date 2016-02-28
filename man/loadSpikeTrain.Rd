\name{loadSpikeTrain}
\alias{loadSpikeTrain}
\title{
Load spike trains from files
}
\description{
Load spike trains from .res, .clu, .sampling_rate_dat files
}
\usage{
loadSpikeTrain(st)
}

\arguments{
  \item{st}{
SpikeTime object
}
}
\details{
}
\value{

\item{st }{SpikeTime object with the spikes loaded}

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

\dontrun{
#### read from file ###
setwd("~/repo/r_packages/data")
session="jp4298-15022016-0106"
st<-new("SpikeTrain") # create SpikeTrain object
st@session=session # set session name
st<-loadSpikeTrain(st) # load res clu and sampling rate
auto<-spikeTimeAutocorrelation(st,bin.size.ms=1,window.size.ms=200,probability=T) # calculate autocorrelation
plot(auto$prob,type='l')
}
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
