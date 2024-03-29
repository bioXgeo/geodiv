% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bearingarea.R
\name{height_ba}
\alias{height_ba}
\title{Value of the Bearing Area Curve at a Specified Value}
\usage{
height_ba(x, xval)
}
\arguments{
\item{x}{A raster or matrix.}

\item{xval}{Numeric value along the x-axis.}
}
\value{
A numeric value of the bearing area function
  corresponding to \code{xval}.
}
\description{
Determines the value of the bearing area curve for a
specific value along the x-axis (\code{xval}).
}
\examples{
# import raster image
data(normforest)
normforest <- terra::unwrap(normforest)

# determine the bearing area function value
# corresponding to an x value of 0.4
val <- height_ba(normforest, 0.4)
}
