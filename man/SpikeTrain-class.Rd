\name{SpikeTrain-class}
\Rdversion{1.1}
\docType{class}
\alias{SpikeTrain-class}
\alias{SpikeTrain}
\alias{loadSpikeTrain,SpikeTrain-method}
\alias{setSpikeTrain,SpikeTrain-method}
\alias{show,SpikeTrain-method}
\alias{spikeTimeAutocorrelation,SpikeTrain-method}

\title{Class \code{"SpikeTrain"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{SpikeTrain(...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{session}:}{Object of class \code{"character"} ~~ }
    \item{\code{samplingRate}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{res}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{time}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{clu}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nCells}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nSpikes}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nSpikesPerCell}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{startInterval}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{endInterval}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{startResIndex}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{endResIndex}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{cellList}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{loadSpikeTrain}{\code{signature(st = "SpikeTrain")}: ... }
    \item{setSpikeTrain}{\code{signature(st = "SpikeTrain")}: ... }
    \item{show}{\code{signature(object = "SpikeTrain")}: ... }
    \item{spikeTimeAutocorrelation}{\code{signature(st = "SpikeTrain")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("SpikeTrain")
}
\keyword{classes}
\usage{
SpikeTrain(...)
}
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
