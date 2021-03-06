% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/focalmetrics.R
\name{.calculate_met_focal}
\alias{.calculate_met_focal}
\title{Calculate Texture Metric for Single Pixel}
\usage{
.calculate_met_focal(landscape, n_row, n_col, points, what, ...)
}
\arguments{
\item{landscape}{A raster or matrix.}

\item{n_row}{Numeric. Number of rows in focal window.}

\item{n_col}{Numeric. Number of columns in focal window.}

\item{points}{Dataframe. Coordinates and values of cells, calculated with the *landscapemetrics*
\code{raster_to_points} function.}

\item{what}{Character. Metric to calculate for each window. Metrics
from the geodiv package are listed below.}

\item{...}{Additional arguments for the metric functions. All applicable arguments
will be applied to the entire list of metrics.}
}
\value{
The metric value over the window.
}
\description{
Calculates the various texture metrics over a window centered
on an individual pixel. This function is modified slightly from the
\code{calculate_lsm_focal} function in the \emph{landscapemetrics} package (Hesselbarth et al. 2019).
}
\details{
Metrics from geodiv package:
\enumerate{
   \item{\code{'sa'}: average surface roughness}
   \item{\code{'sq'}: root mean square roughness}
   \item{\code{'s10z'}: ten-point height}
   \item{\code{'sdq'}: root mean square slope of surface, 2-point method}
   \item{\code{'sdq6'}: root mean square slope of surface, 7-point method}
   \item{\code{'sdr'}: surface area ratio}
   \item{\code{'sbi'}: surface bearing index}
   \item{\code{'sci'}: core fluid retention index}
   \item{\code{'ssk'}: skewness}
   \item{\code{'sku'}: kurtosis}
   \item{\code{'sds'}: summit density}
   \item{\code{'sfd'}: 3d fractal dimension}
   \item{\code{'srw'}: dominant radial wavelength, radial wavelength index, mean half wavelength}
   \item{\code{'std'}: angle of dominating texture, texture direction index}
   \item{\code{'svi'}: valley fluid retention index}
   \item{\code{'stxr'}: texture aspect ratio}
   \item{\code{'ssc'}: mean summit curvature}
   \item{\code{'sv'}: maximum valley depth}
   \item{\code{'sph'}: maximum peak height}
   \item{\code{'sk'}: core roughness depth}
   \item{\code{'smean'}: mean peak height}
   \item{\code{'svk'}: reduced valley depth}
   \item{\code{'spk'}: reduced peak height}
   \item{\code{'scl'}: correlation length}
   \item{\code{'sdc'}: bearing area curve height interval}
}
}
\references{
\enumerate{
\item{Hesselbarth, M.H.K., Sciaini, M., With, K.A., Wiegand, K., Nowosad, J. 2019.
landscapemetrics: an open-source R tool to calculate landscape metrics. - Ecography 42:1648-1657(ver. 0).}
}
}
