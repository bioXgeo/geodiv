% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bearingarea.R
\name{sci}
\alias{sci}
\title{Core Fluid Retention Index}
\usage{
sci(x)
}
\arguments{
\item{x}{A raster or matrix.}
}
\value{
A numeric value representing the core fluid
retention index.
}
\description{
Determines the core fluid retention index (Sci). This
value is the void volume (area under the bearing area
curve) in the 'core' zone. See Figure 2a from Kedron
et al. (2018) for more details.
}
\examples{
# import raster image
data(normforest)
normforest <- terra::unwrap(normforest)

# determine the core fluid retention index
Sci <- sci(normforest)
}
