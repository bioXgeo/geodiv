#' Radian to Degree Conversion
#'
#' Converts radian value(s) to degrees.
#'
#' @param x Numeric. Radian value(s).
#' @return Numeric of degree value(s).
.rad2deg <- function(x) {(x * 180) / (pi)}

#' Degree to Radian Conversion
#'
#' Converts degree value(s) to radians.
#'
#' @param x Numeric. Degree value(s).
#' @return Numeric of degree value(s).
.deg2rad <- function(x) {(x * pi) / (180)}

#' Calculate a Least Squares Polynomial Surface
#'
#' Fits a polynomial surface of order \code{n} to a raster
#' or matrix.
#'
#' @param x A raster or matrix.
#' @param order Numeric. Indicates the polynomial order to be fit.
#' @return A matrix with values predicted from the polynomial fit.
#' @examples
#'
#' # import raster image
#' data(orforest)
#' orforest <- terra::unwrap(orforest)
#'
#' # find the 2nd order least squares polynomial surface
#' polyfit <- fitplane(orforest, order = 2)
#'
#' # create raster of polyfit
#' x <- terra::setValues(orforest, polyfit)
#'
#' # plot the fit
#' terra::plot(x)
#' @importFrom terra rast crop crds
#' @importFrom spatial surf.ls
#' @export
fitplane <- function(x, order) {
  stopifnot('x must be a raster or matrix.' = inherits(x, c('RasterLayer', 'matrix', 'SpatRaster')))
  if(length(order) > 1) {stop('too many values supplied to order.')}
  stopifnot('order must be numeric or integer.' = inherits(order, c('numeric', 'integer')))

  if(order %% 1 > 0) {
    warning('order will be rounded to the nearest integer.')
    order <- as.integer(floor(order))}
  if(order < 0) {stop('order must be >= 0.')}

  order <- as.integer(order)

  if (class(x)[1] %in% c('RasterLayer')) {
    x <- rast(x)
  }

  if (class(x)[1] %in% c('RasterLayer', 'SpatRaster')) {
    # extract coordinates and values
    xcoord <- crds(x, na.rm = FALSE)[, 1]
    ycoord <- crds(x, na.rm = FALSE)[, 2]
    z <- as.numeric(x[])
  } else {
    xcoord <- rep(seq(1, ncol(x)), nrow(x))
    ycoord <- rep(seq(1, nrow(x)), each = ncol(x))
    z <- c(x)
  }

  # fit least squares polynomial with order = order
  surfmod <- spatial::surf.ls(np = order, x = xcoord[!is.na(z) & !is.na(xcoord)], y = ycoord[!is.na(z) & !is.na(xcoord)], z = z[!is.na(z) & !is.na(xcoord)])

  # predict polynomial model over raster
  surfvals <- matrix(predict(surfmod, xcoord, ycoord), nrow = nrow(x), ncol = ncol(x), byrow = TRUE)

  return(surfvals)
}

#' Finds the Best Fit Polynomial Surface
#'
#' Finds the best fit polynomial surface for a raster or matrix. This
#' function tests least squares polynomial fits with orders of
#' 0 - 3 and determines which order minimizes the error when the
#' fit is subtracted from the original image.
#'
#' @param x A raster or matrix.
#' @return A raster or matrix of the same size as the input with values
#'   predicted from the best polynomial fit.
#' @examples
#'
#' # import raster image
#' data(orforest)
#' orforest <- terra::unwrap(orforest)
#'
#' # find the least squares polynomial surface
#' poly <- bestfitplane(orforest)
#'
#' # plot the fit
#' terra::plot(poly)
#' @importFrom terra rast setValues crop
#' @export
bestfitplane <- function(x) {
  stopifnot('x must be a raster or matrix.' = inherits(x, c('RasterLayer', 'matrix', 'SpatRaster')))

  # fit least squares plane for polynomials from orders 0-3
  mods <- lapply(seq(0, 3), FUN = function(i) fitplane(x, order = i))

  # convert raster to matrix
  if (class(x)[1] %in% c('RasterLayer', 'SpatRaster')) {
    xmat <- matrix(x, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
  } else {
    xmat <- x
  }

  # calculate raster errors from best fit planes for each order tested
  errlist <- lapply(seq(1, 4), FUN = function(i) xmat - mods[[(i)]])
  meanerr <- lapply(seq(1, 4), FUN = function(i) abs(mean(errlist[[(i)]], na.rm = TRUE)))

  # determine which order is best
  bestfit <- which(as.numeric(meanerr) == min(as.numeric(meanerr), na.rm = TRUE))

  # if more than one, get first
  if (length(bestfit) > 1) {
    bestfit <- min(bestfit, na.rm = TRUE)
  }

  # fill in raster with best fit values
  if (class(x)[1] == 'RasterLayer' | class(x)[1] == 'SpatRaster'){
    bfx <- terra::setValues(x, mods[[bestfit]])
  } else {
    bfx <- matrix(mods[[bestfit]], nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
  }

  print(paste('Order of polynomial that minimizes errors: ', bestfit - 1, sep = ''))
  return(bfx)
}

#' Removes the Best Fit Polynomial Surface
#'
#' Finds the best fit polynomial surface for a raster or matrix and
#' subtracts it from the actual values. The output image has positive values
#' where the actual values are higher than the surface and negative values
#' where the actual value are lower than the surface.
#'
#' @param x A raster or matrix.
#' @return A raster or matrix of the same size as the input with values
#'   equal to the difference between the original and bestfit
#'   plane.
#' @examples
#' # import raster image
#' data(orforest)
#' orforest <- terra::unwrap(orforest)
#'
#' # remove the least squares polynomial surface
#' new_rast <- remove_plane(orforest)
#'
#' # plot
#' terra::plot(new_rast)
#' @importFrom terra rast
#' @export
remove_plane <- function(x) {
  stopifnot('x must be a raster or matrix.' = inherits(x, c('RasterLayer', 'matrix', 'SpatRaster')))

  bfx <- bestfitplane(x)

  errors <- x - bfx # higher = above plane, lower = below plane

  return(errors)
}

#' Rotates a matrix 180 degrees.
#'
#' Rotates a matrix 180 degrees. Code is from https://stackoverflow.com/questions/16496210/rotate-a-matrix-in-r-by-90-degrees-clockwise.
#'
#' @param x A matrix.
#' @return A matrix rotate 180 degrees.
#' @export
rotate <- function(x) {
  x[] <- rev(x)
  return(x)
}
