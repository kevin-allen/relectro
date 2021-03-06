% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotsSpatial.R
\name{firingRateMapPlot}
\alias{firingRateMapPlot}
\title{Plot a single firing rate map}
\usage{
firingRateMapPlot(m, name = "", outma = c(2, 2, 2, 2), margin = c(1, 1, 1,
  1), cex.title = 0.5, cex.name = 0.5, xlab = "", ylab = "",
  main.title = "", peak.rate.prefix = "")
}
\arguments{
\item{m}{A data.frame or matrix containing the data of the firing rate map. If a data.frame is given as argument, it should have the column names x, y and rate}

\item{name}{Character vector giving the name of the firing rate map which is displayed before the max firing rate}

\item{outma}{Outer margins of the figure}

\item{margin}{Inner margins of the figure}

\item{cex.title}{Size of the font use for the title}

\item{cex.name}{Size of the font use for the name of the map}

\item{xlab}{Name to display under the x axis}

\item{ylab}{Name to display at the left of the y axis}

\item{main.title}{A title for the figure}

\item{peak.rate.prefix}{Additional information to display before the peak firing rate.}
}
\description{
Plot a 2-dimensional representation of firing rate using the image function included in the graphics package
}
