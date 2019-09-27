#' NDVI for a portion of southwestern Oregon, USA.
#'
#' A raster image of Normalized Difference Vegetation Index (NDVI)
#' derived from Landsat data for a small portion of SW Oregon state.
#'
#' NDVI values are derived from Landsat scene path 45,
#' row 30 summarized as the mean NDVI value between June
#' and August 2000 at roughly 30m resolution. Clouds
#' were removed from the Landsat scene before calculating
#' the mean. The image was created using Google Earth
#' Engine in August 2018.
#'
#' @format A raster image with 371 x 371 pixels
#' \describe{
#'   \item{range}{0 -- 1}
#'   \item{bounds}{-123, -122.9, 43.0002, 43.1}
#'   \item{resolution}{30m x 30m}
#'   \item{projection}{WGS84}
#'   \item{scalar}{1}
#'   ...
#' }
"orforest"

#' NDVI errors for a portion of southwestern Oregon, USA.
#'
#' A raster image of Normalized Difference
#' Vegetation Index (NDVI) errors derived from Landsat
#' data for a small portion of SW Oregon state. This
#' raster was calculated by subtracting the best-fit polynomial
#' plane from the \code{orforest} values. The best-fit
#' polynomial plane was calculated using \code{bestfitplane}.
#'
#' NDVI values are derived from Landsat scene path 45,
#' row 30 summarized as the mean NDVI value between June
#' and August 2000 at roughly 30m resolution. Clouds
#' were removed from the Landsat scene before calculating
#' the mean. The image was created using Google Earth
#' Engine in August 2018.
#'
#' @format A raster image with 371 x 371 pixels
#' \describe{
#'   \item{range}{-0.5854638 -- 0.1585918}
#'   \item{bounds}{-123, -122.9, 43.0002, 43.1}
#'   \item{resolution}{30m x 30m}
#'   \item{projection}{WGS84}
#'   \item{scalar}{1}
#'   ...
#' }
"normforest"

#' Raster of Enhanced Vegetation Index (EVI) over Oregon, USA.
#'
#' A raster image of maximum growing season (May-September)
#' EVI over Oregon state, USA derived from MODIS 16-day EVI data.
#' The growing season summary was completed in Google Earth Engine
#' on 30 August, 2019. Data were filtered using the general quality assurance codes
#' included with the MODIS data; only the highest quality (0: VI
#' produced with good quality) data were included.
#'
#' @format A raster image with 1918 x 3670 pixels.
#' \describe{
#'   \item{range}{0 -- 10000}
#'   \item{bounds}{-124.7041, -116.4621, 41.99175, 46.29917}
#'   \item{resolution}{250m x 250m}
#'   \item{projection}{WGS84}
#'   \item{scalar}{0.0001}
#'   ...
#' }
"oregonEVI"

#' Raster of elevation for Oregon, USA.
#'
#' A raster image of Shuttle Radar Topography Mission (SRTM)
#' elevation over Oregon state, USA.
#' The elevation data were prepared (aggregated to 270m resolution)
#' in Google Earth Engine and downloaded on 25 September, 2019.
#'
#' @format A raster image with 1918 x 3670 pixels.
#' \describe{
#'   \item{range}{0 -- 3332}
#'   \item{bounds}{-124.7039, -116.463, 41.99197, 46.29957}
#'   \item{resolution}{270m x 270m}
#'   \item{projection}{WGS84}
#'   \item{scalar}{1}
#'   ...
#' }
"oregonElev"
