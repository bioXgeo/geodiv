# functions to find Sdq and Sdq6 using two-point and seven-point slopes

# both functions from the equations here:
# https://www.ntmdt-si.ru/data/media/files/manuals/image_analisys_p9_nov12.e.pdf

#' Root Mean Square Slope of Surface
#'
#' Calculates the root mean square slope of a raster
#' surface using the two-point method. This function is based
#' on the equations found at
#' https://www.ntmdt-si.ru/data/media/files/manuals/image_analisys_p9_nov12.e.pdf.
#'
#' @param x A raster object.
#' @return A numeric value representing the two-point root
#'   mean square slope, Sdq.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # calculate root mean square slope
#' Sdq <- sdq(normforest)
#' @export
sdq <- function(x) {
  if(class(x) != 'RasterLayer') {stop('x must be a raster.')}

  # z values, coordinates, and resolution (change in x, y)
  z <- getValues(x)
  deltax <- res(x)[1]
  deltay <- res(x)[2]

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
#' Calculates the area root mean square slope of a raster
#' surface using the seven-point method. This function is based
#' on the equations found at
#' https://www.ntmdt-si.ru/data/media/files/manuals/image_analisys_p9_nov12.e.pdf.
#'
#' @param x A raster object.
#' @return A numeric value representing the seven-point root
#'   mean square slope, Sdq6.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # calculate area root mean square slope
#' Sdq6 <- sdq6(normforest)
#' @export
sdq6 <- function(x) {
  if(class(x) != 'RasterLayer') {stop('x must be a raster.')}

  # get resolution
  deltax <- res(x)[1]
  deltay <- res(x)[2]

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
