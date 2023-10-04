#' Root Mean Square Slope of Surface
#'
#' Calculates the root mean square slope of a raster or matrix
#' surface using the two-point method.
#'
#' @param x A raster or matrix.
#' @return A numeric value representing the two-point root
#'   mean square slope, Sdq. The units of the returned value
#'   are change in z per one unit (pixel).
#' @references  This function is based
#' on the equations found at
#' https://www.ntmdt-si.ru/data/media/files/manuals/image_analisys_p9_nov12.e.pdf.
#' @examples
#' # import raster image
#' data(normforest)
#' normforest <- terra::unwrap(normforest)
#'
#' # calculate root mean square slope
#' Sdq <- sdq(normforest)
#' @importFrom terra rast
#' @export
sdq <- function(x) {
  stopifnot('x must be a raster or matrix.' = inherits(x, c('RasterLayer', 'matrix', 'SpatRaster')))

  deltax <- 1
  deltay <- 1

  # calculate z with an offset of x + 1, y + 1
  z_xplus <- zshift(x, xdist = 1, ydist = 0, yrm = 1)
  z_yplus <- zshift(x, xdist = 0, ydist = 1, xrm = 1)
  z <- zshift(x, xrm = 1, yrm = 1)

  # calculate two-point slope
  sdq <- sqrt(sum((((z - z_xplus) / deltax) ^ 2) +
                    (((z - z_yplus) / deltay) ^ 2),
                  na.rm = TRUE) / length(z))

  return(sdq)
}

#' Root Area Mean Square Slope of Surface
#'
#' Calculates the area root mean square slope of a raster or matrix
#' surface using the seven-point method.
#'
#' @param x A raster or matrix.
#' @return A numeric value representing the seven-point root
#'   mean square slope, Sdq6. The units of the returned value
#'   are change in z per one unit (pixel).
#' @references  This function is based
#' on the equations found at
#' https://www.ntmdt-si.ru/data/media/files/manuals/image_analisys_p9_nov12.e.pdf.
#' @examples
#' # import raster image
#' data(normforest)
#' normforest <- terra::unwrap(normforest)
#'
#' # calculate area root mean square slope
#' Sdq6 <- sdq6(normforest)
#' @importFrom terra rast
#' @export
sdq6 <- function(x) {
  stopifnot('x must be a raster or matrix.' = inherits(x, c('RasterLayer', 'matrix', 'SpatRaster')))

  deltax <- 1 # per unit, not per degree, etc.
  deltay <- 1

  # get dimensions
  N <- dim(x)[1] # rows
  M <- dim(x)[2] # cols

  # calculate offset vectors (shifts xdist, ydist, then removes rows xrm, yrm)
  z_xps1 <- zshift(x, xdist = 1, ydist = 0, xrm = c(-3, 3), yrm = c(-3, 3))
  z_xps2 <- zshift(x, xdist = 2, ydist = 0, xrm = c(-3, 3), yrm = c(-3, 3))
  z_xps3 <- zshift(x, xdist = 3, ydist = 0, xrm = c(-3, 3), yrm = c(-3, 3))
  z_xmn1 <- zshift(x, xdist = -1, ydist = 0, xrm = c(-3, 3), yrm = c(-3, 3))
  z_xmn2 <- zshift(x, xdist = -2, ydist = 0, xrm = c(-3, 3), yrm = c(-3, 3))
  z_xmn3 <- zshift(x, xdist = -3, ydist = 0, xrm = c(-3, 3), yrm = c(-3, 3))

  # same as above, but for y
  z_yps1 <- zshift(x, xdist = 0, ydist = 1, xrm = c(-3, 3), yrm = c(-3, 3))
  z_yps2 <- zshift(x, xdist = 0, ydist = 2, xrm = c(-3, 3), yrm = c(-3, 3))
  z_yps3 <- zshift(x, xdist = 0, ydist = 3, xrm = c(-3, 3), yrm = c(-3, 3))
  z_ymn1 <- zshift(x, xdist = 0, ydist = -1, xrm = c(-3, 3), yrm = c(-3, 3))
  z_ymn2 <- zshift(x, xdist = 0, ydist = -2, xrm = c(-3, 3), yrm = c(-3, 3))
  z_ymn3 <- zshift(x, xdist = 0, ydist = -3, xrm = c(-3, 3), yrm = c(-3, 3))

  # calculate two-point slope (see https://www.ntmdt-si.ru/data/media/files/manuals/image_analisys_p9_nov12.e.pdf)
  pklx <- (1 / (60 * deltax)) * (z_xps3 - (9 * z_xps2) + (45 * z_xps1) - (45 * z_xmn1) -
                                   (9 * z_xps2) - z_xmn3)
  pkly <- (1 / (60 * deltay)) * (z_yps3 - (9 * z_yps2) + (45 * z_yps1) - (45 * z_ymn1) -
                                   (9 * z_yps2) - z_ymn3)

  # calculate sdq6
  sdq6 <- sqrt((1 / length(x)) * sum(sum(pklx ^ 2, na.rm = TRUE),
                                sum(pkly ^ 2, na.rm = TRUE),
                                na.rm = TRUE))

  return(sdq6)
}
