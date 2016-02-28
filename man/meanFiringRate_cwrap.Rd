\name{meanFiringRate_cwraps}
\alias{meanFiringRate_cwraps}
\title{
c function that calculates the mean firing rate
}
\description{
}
\usage{
 results<- .Call("meanFiringRate_cwrap",
                            st@cellList,
                            length(st@cellList),
                            st@clu,
                            st@res,
                            st@nSpikes,
                            st@startInterval,
                            st@endInterval,
                            st@startResIndex,
                            st@endResIndex,
                            length(st@startResIndex),
                            st@samplingRate)
}
\arguments{
  \item{cell.list}{}
  \item{n.cells}{}
  \item{clu}{}
  \item{res}{}
  \item{nSpikes}{}
  \item{start.interval}{}
  \item{end.interval}{}
  \item{start.res.index}{}
  \item{end.res.index}{}
  \item{lenght.intervals}{}
  \item{sampling.rate}{}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{df }{A vector with the meanFiringRate}

}
\references{
}
\author{
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
