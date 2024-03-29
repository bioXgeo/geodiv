% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simpsons.R
\name{simpsons}
\alias{simpsons}
\title{Simpson's Rule Empirical Area Under a Curve}
\usage{
simpsons(f, a, b, n = 100)
}
\arguments{
\item{f}{A function.}

\item{a}{Numeric. The left x boundary.}

\item{b}{Numeric. The right x boundary.}

\item{n}{Numeric. The number of subdivisions along the function
line.}
}
\value{
A numeric value representing the area under the curve with
  x bounds \code{a} and \code{b}.
}
\description{
Calculates the area below a curve from
points \code{a} to \code{b}. This function is provided
for general use.
}
\details{
Note that if y-values are negative, this returns the area
above the function line.
}
\examples{
# import raster image
data(normforest)
normforest <- terra::unwrap(normforest)

# basic values
z <- terra::values(normforest)

# calculate cumulative probability density function of surface 'height' (= ndvi)
mod <- ecdf((1 - z))

# calculate integral
int_area <- simpsons(f = mod, b = 1, a = 0.8, n = 500)
}
