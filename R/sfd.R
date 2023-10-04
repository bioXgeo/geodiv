#' Calculate the fractal dimension of a raster.
#'
#' Calculates the 3D fractal dimension of a raster using the
#' triangular prism surface area method.
#'
#' @param x A raster or matrix.
#' @param silent Logical. If \code{FALSE} (default), the function will
#' print warning messages.
#' @return A numeric value representing the fractal dimension of
#' the image.
#' @references Clarke, K.C., 1986. Computation of the fractal dimension
#' of topographic surfaces using the triangular prism surface area method.
#' Computers & Geosciences, 12(5), pp.713-722.
#' @examples
#'
#' # import raster image
#' data(normforest)
#' normforest <- terra::unwrap(normforest)
#'
#' # calculate the fractal dimension
#' Sfd <- sfd(normforest)
#' @importFrom terra rast
#' @import Rcpp
#' @export
sfd <- function(x, silent = FALSE) {
  # check type
  stopifnot('x must be a raster or matrix.' = inherits(x, c('RasterLayer', 'matrix', 'SpatRaster')))

  # if raster, convert to matrix
  if (class(x)[1] %in% c('RasterLayer', 'SpatRaster')) {
    if (silent == FALSE) {
      # tell users that this will always reproject to equal area
      print('Warning: Raster will be converted to matrix format.')
    }
    # matrices are faster for window_metric, so convert to matrix
    mat <- matrix(x[], nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
  } else {
    mat <- x
  }

  # send to C function
  out <- sfd_(mat)

  return(out)
}
