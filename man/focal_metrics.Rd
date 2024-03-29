% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/focalmetrics.R
\name{focal_metrics}
\alias{focal_metrics}
\title{Calculate Texture Metrics per Pixel}
\usage{
focal_metrics(x, window, metrics, progress, ...)
}
\arguments{
\item{x}{A raster or matrix. Image over which to apply focal window calculations.}

\item{window}{Matrix. The focal window used to create the image.}

\item{metrics}{List. List of metrics to apply. Function names must be strings.}

\item{progress}{Logical. Display progress through metrics list?}

\item{...}{Additional arguments for the metric functions. All applicable arguments
will be applied to the entire list of metrics.}
}
\value{
A raster of the metric calculated in windows over the raster or matrix.
If the input was a matrix, the function will return a raster with an extent of [0, 1, 0, 1].
}
\description{
Calculates the various texture metrics over windows centered
on individual pixels. This creates a continuous surface of the
texture metric.
This function is a modified version of the \code{window_lsm} function from the
\emph{landscapemetrics} package (Hesselbarth et al. 2019).
}
\details{
Metrics available from geodiv package:
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
\examples{
# import raster image
data(normforest)
normforest <- terra::unwrap(normforest)

# crop raster to smaller area
x <- terra::crop(normforest, terra::ext(normforest[1:100, 1:100, drop = FALSE]))

# get a surface of root mean square roughness
sa_img <- focal_metrics(x = x, window = matrix(1, 5, 5),
                        metrics = list('sa'), progress = TRUE)

# plot the result
terra::plot(sa_img$sa)
}
\references{
\enumerate{
\item{Hesselbarth, M.H.K., Sciaini, M., With, K.A., Wiegand, K., Nowosad, J. 2019.
landscapemetrics: an open-source R tool to calculate landscape metrics. - Ecography 42:1648-1657(ver. 0).}
}
}
