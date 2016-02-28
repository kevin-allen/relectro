\name{make.pairs}
\alias{make.pairs}
\title{
Make all combinations of pairs
}
\description{
Takes one or two lists of numbers and makes all possible combination
within or across the lists.
}
\usage{
make.pairs(cl1=1:4)
make.pairs(cl1=1:4,cl2=5:9)
}
\arguments{
  \item{cl1}{
first list of numbers
}
  \item{cl2}{
second list of numbers
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{df }{A data.frame with 2 columns, named Var1 and Var2 containing
  the possible combinations}

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
make.pairs(cl1=1:4)
make.pairs(cl1=1:4,cl2=5:9)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
