\name{join.intervals}
\alias{join.intervals}
\title{
Join intervals with a AND logic
}
\description{
Join intervals so that the time that is within both list of interval in obtained
}
\usage{
join.intervals(s1, e1, s2, e2)
}
\arguments{
  \item{s1}{
start times of the first set of intervals
}
  \item{e1}{
end times of the first set of intervals
}
  \item{s2}{
start times of the second set of intervals
}
  \item{e2}{
end times of the second set of intervals
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{df }{A data.frame with 2 columns, named s and e for the start and
  end times of the resulting intervals}

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

s1<-c(0,30,50,100)
e1<-c(10,46,55,150)
s2<-c(45,125)
e2<-c(55,155)
join.intervals(s1,e1,s2,e2)
rm(s1,s2,e1,e2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
