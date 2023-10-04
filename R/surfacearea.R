#' Flattened Surface Area
#'
#' Calculates the surface area of a flat raster or matrix with the
#' same x, y bounds as the study surface.
#'
#' This function scales both x and y to between 0 and 1. This
#' is done because most satellite data have units where the x,
#' y units do not equal the z units and because the flat
#' surface area is usually compared to the actual surface area.
#' Surface area is calculated over the sample area (N-1, M-1).
#'
#' @param x A raster or matrix.
#' @return A numeric value representing the scaled surface
#'   area of a flattened surface with the same x, y bounds.
#' @examples
#' # import raster image
#' data(normforest)
#' normforest <- terra::unwrap(normforest)
#'
#' # calculate flattened surface area
#' flatsa(normforest)
#' @importFrom terra rast res crds setValues crs
#' @export
flatsa <- function(x) {
  stopifnot('x must be a raster or matrix.' = inherits(x, c('RasterLayer', 'matrix', 'SpatRaster')))

  # In case the value area of the raster is an odd shape,
  # calculate the surface area of the flattened raster in the same way
  # as actual surface area
  # get dimensions
  N <- dim(x)[1] # rows
  M <- dim(x)[2] # cols

  # convert matrix to raster if necessary (equal area)
  if (class(x)[1] == 'matrix') {
    x <- rast(x)
    terra::crs(x) <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  }

  # coordinates and resolution (change in x, y)
  deltax <- res(x)[1]
  deltay <- res(x)[2]

  # create flat raster
  z <- x[]
  z[!is.na(z)] <- 0
  fakerast <- x
  fakerast <- setValues(fakerast, z)

  # get shifted z values of flat raster
  z <- zshift(fakerast, xdist = 0, ydist = 0, xrm = 1, yrm = 1, scale = FALSE)
  z_ypl <- zshift(fakerast, xdist = 0, ydist = 1, xrm = 1, scale = FALSE)
  z_xpl <- zshift(fakerast, xdist = 1, ydist = 0, yrm = 1, scale = FALSE)
  z_xplypl <- zshift(fakerast, xdist = 1, ydist = 1, scale = FALSE)

  # normalize deltas
  divide <- mean(deltax, deltay)
  deltax <- deltax / divide
  deltay <- deltay / divide

  # remove any outside NA values
  z_ypl <- z_ypl[!is.na(z)]
  z_xpl <- z_xpl[!is.na(z)]
  z_xplypl <- z_xplypl[!is.na(z)]
  z <- z[!is.na(z)]

  # calculate area
  akl1 <- ((1 / 2) *
             (sqrt((deltay ^ 2) + ((z_ypl - z) ^ 2)) *
                sqrt((deltax ^ 2) + ((z_xpl - z) ^ 2)))) +
    ((1 / 2) *
       (sqrt((deltay ^ 2) + ((z_xplypl - z_xpl) ^ 2)) *
          sqrt((deltax ^ 2) + ((z_xplypl - z_ypl) ^ 2))))

  akl2 <- ((1 / 2) *
             (sqrt((deltay ^ 2) + ((z_ypl - z) ^ 2)) *
                sqrt((deltax ^ 2) + ((z_xplypl - z_ypl) ^ 2)))) +
    ((1 / 2) *
       (sqrt((deltay ^ 2) + ((z_xplypl - z_xpl) ^ 2)) *
          sqrt((deltax ^ 2) + ((z_xpl - z) ^ 2))))

  akl <- (1 / 2) * (akl1 + akl2)

  sa <- sum(akl, na.rm = TRUE)

  return(sa)
}

#' Surface Area
#'
#' Calculates the scaled surface area of a raster or matrix.
#'
#' This function scales both x and y, as well as the surface value (z),
#' to between 0 and 1 to best match their units. This is done because
#' most satellite data have units where the x, y units do not equal the
#' z units. The surface area represents the surface area of the sample
#' area (N-1, M-1).
#'
#' Note that the surface object may have NA values around the edges,
#' but should not have any missing values within the main area.
#'
#' @param x A raster or matrix.
#' @return A numeric value representing the scaled surface area.
#' @examples
#' # import raster image
#' data(normforest)
#' normforest <- terra::unwrap(normforest)
#'
#' # calculate surface area
#' surface_area(normforest)
#' @importFrom terra rast crop res crs
#' @export
surface_area <- function(x) {
  stopifnot('x must be a raster or matrix.' = inherits(x, c('RasterLayer', 'matrix', 'SpatRaster')))

  # get dimensions
  N <- dim(x)[1] # rows
  M <- dim(x)[2] # cols

  # convert matrix to raster if necessary (equal area)
  if (class(x)[1] == 'matrix') {
    x <- rast(x)
    terra::crs(x) <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  }

  # coordinates and resolution (change in x, y)
  deltax <- res(x)[1]
  deltay <- res(x)[2]

  z <- zshift(x, xdist = 0, ydist = 0, xrm = 1, yrm = 1, scale = TRUE)
  z_ypl <- zshift(x, xdist = 0, ydist = 1, xrm = 1, scale = TRUE)
  z_xpl <- zshift(x, xdist = 1, ydist = 0, yrm = 1, scale = TRUE)
  z_xplypl <- zshift(x, xdist = 1, ydist = 1, scale = TRUE)

  # normalize deltas
  divide <- mean(deltax, deltay)
  deltax <- deltax / divide
  deltay <- deltay / divide

  # remove any outside NA values
  z_ypl <- z_ypl[!is.na(z)]
  z_xpl <- z_xpl[!is.na(z)]
  z_xplypl <- z_xplypl[!is.na(z)]
  z <- z[!is.na(z)]

  # calculate area
  akl1 <- ((1 / 2) *
             (sqrt((deltay ^ 2) + ((z_ypl - z) ^ 2)) *
                sqrt((deltax ^ 2) + ((z_xpl - z) ^ 2)))) +
    ((1 / 2) *
       (sqrt((deltay ^ 2) + ((z_xplypl - z_xpl) ^ 2)) *
          sqrt((deltax ^ 2) + ((z_xplypl - z_ypl) ^ 2))))

  akl2 <- ((1 / 2) *
             (sqrt((deltay ^ 2) + ((z_ypl - z) ^ 2)) *
                sqrt((deltax ^ 2) + ((z_xplypl - z_ypl) ^ 2)))) +
    ((1 / 2) *
       (sqrt((deltay ^ 2) + ((z_xplypl - z_xpl) ^ 2)) *
          sqrt((deltax ^ 2) + ((z_xpl - z) ^ 2))))

  akl <- (1 / 2) * (akl1 + akl2)

  sa <- sum(akl, na.rm = TRUE)

  return(sa)
}

#' Surface Area Ratio
#'
#' Calculates the surface area ratio of a raster or matrix. This is the
#' ratio of a flat surface to the actual surface.
#'
#' This function scales both x and y, as well as the surface value (z),
#' to between 0 and 1 to best match their units. This is done because
#' most satellite data have units where the x, y units do not equal the
#' z units. Surface area is calculated over the sample area (N-1, M-1).
#'
#' @param x A raster or matrix.
#' @return A numeric value representing the surface area ratio.
#' @examples
#' # import raster image
#' data(normforest)
#' normforest <- terra::unwrap(normforest)
#'
#' # calculate the surface area ratio
#' Sdr <- sdr(normforest)
#' @importFrom terra rast
#' @export
sdr <- function(x) {
  stopifnot('x must be a raster or matrix.' = inherits(x, c('RasterLayer', 'matrix', 'SpatRaster')))

  # get area of flat plane
  flat_area <- flatsa(x)

  # get surface area of raster
  sa <- surface_area(x)

  # calculate area ratio
  adr <- ((sa - flat_area) / flat_area) * 100

  return(adr)
}
