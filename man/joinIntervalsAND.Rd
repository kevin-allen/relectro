% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/JoinIntervals.R
\name{joinIntervalsAND}
\alias{joinIntervalsAND}
\title{Join two sets of intervals using an AND logic}
\usage{
joinIntervalsAND(s1, e1, s2, e2)
}
\arguments{
\item{s1}{Numeric vector containing the start times of the first set of intervals}

\item{e1}{Numeric vector containing the end times of the first set of intervals}

\item{s2}{Numeric vector containing the start times of the second set of intervals}

\item{e2}{Numeric vector containing the end times of the second set of intervals}
}
\value{
matrix containing the start and end of the resulting intervals. First col is the start times. Second col is end times.
}
\description{
To be included in the output interval, the time period needs to be included in both set of time intervals.
Adjacent s1-e1 or s2-e2 will be joined together by default
}
\examples{
s1<-c(0,20000)
e1<-c(500,25000)
s2<-c(250,24000)
e2<-c(600,30000)
joinIntervalsAND(s1,e1,s2,e2)
}
