\name{autocorrelation_cwraps}
\alias{autocorrelation_cwraps}
\title{
c function that calculates spike-time autocorrelation
}
\description{
}
\usage{
 results<- .Call("autocorrelation_cwrap",
                        st@cellList,
                        length(st@cellList),
                        st@clu,
                        st@res,
                        st@nSpikes,
                        n.bins,
                        window.size,
                        st@startResIndex,
                        st@endResIndex,
                        length(st@startResIndex),
                        probability)
}
\arguments{
  \item{cell.list}{}
  \item{n.cells}{}
  \item{clu}{}
  \item{res}{}
  \item{sSpikes}{}
  \item{n.bins}{}
  \item{window.size}{}
  \item{start.res.index}{}
  \item{end.res.index}{}
  \item{lenght.intervals}{}
  \item{probability}{}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{df }{A vector with the autocorrelation}

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
