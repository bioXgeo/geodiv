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

#' Calculate a Least Squares Polynomial Plane
#'
#' Fits a polynomial plane of order \code{n} to a raster
#' or matrix.
#'
#' @param x A raster or matrix.
#' @param order Numeric. Indicates the polynomial order to be fit.
#' @return A matrix with values predicted from the polynomial fit.
#' @examples
#' library(raster)
#'
#' # import raster image
#' data(orforest)
#'
#' # find the 2nd order least squares polynomial plane
#' polyfit <- fitplane(orforest, order = 2)
#'
#' # create raster of polyfit
#' x <- setValues(orforest, polyfit)
#'
#' # plot the fit
#' plot(x)
#' @export
fitplane <- function(x, order) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}
  if(length(order) > 1) {stop('too many values supplied to order.')}
  if(class(order) != 'integer' & class(order) != 'numeric') {stop('order must be numeric or integer.')}
  if(order %% 1 > 0) {
    warning('order will be rounded to the nearest integer.')
    order <- as.integer(floor(order))}
  if(order < 0) {stop('order must be >= 0.')}

  order <- as.integer(order)

  if (class(x) == 'RasterLayer') {
    # extract coordinates and values
    xcoord <- sp::coordinates(x)[, 1]
    ycoord <- sp::coordinates(x)[, 2]
    z <- getValues(x)
  } else {
    xcoord <- rep(seq(1, ncol(x)), nrow(x))
    ycoord <- rep(seq(1, nrow(x)), each = ncol(x))
    z <- c(x)
  }

  # fit least squares polynomial with order = order
  surfmod <- spatial::surf.ls(np = order, xcoord[!is.na(z)], ycoord[!is.na(z)], z[!is.na(z)])

  # predict polynomial model over raster
  surfvals <- matrix(predict(surfmod, xcoord, ycoord), nrow = nrow(x), ncol = ncol(x), byrow = TRUE)

  return(surfvals)
}

#' Finds the Best Fit Polynomial Plane
#'
#' Finds the best fit polynomial plane for a surface. This
#' function tests least squares polynomial fits with orders of
#' 1 - 3 and determines which order minimizes the error when the
#' fit is subtracted from the original image.
#'
#' @param x A raster or matrix.
#' @return A raster or matrix of the same size as the input with values
#'   predicted from the best polynomial fit.
#' @examples
#' library(raster)
#'
#' # import raster image
#' data(orforest)
#'
#' # find the least squares polynomial plane
#' poly <- bestfitplane(orforest)
#'
#' # plot the fit
#' plot(poly)
#' @export
bestfitplane <- function(x) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}

  # fit least squares plane for polynomials from orders 2 - 3
  mods <- lapply(seq(0, 3), FUN = function(i) fitplane(x, order = i))

  # convert raster to matrix
  if (class(x) == 'RasterLayer') {
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
  if (class(x) == 'RasterLayer'){
    bfx <- setValues(x, mods[[bestfit]])
  } else {
    bfx <- matrix(mods[[bestfit]], nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
  }

  print(paste('Order of polynomial that minimizes errors: ', bestfit - 1, sep = ''))
  return(bfx)
}

#' Removes the Best Fit Polynomial Plane from a Surface
#'
#' Finds the best fit polynomial plane for a surface and
#' subtracts it from the actual values. The remaining
#' surface has positive values where the actual values are higher
#' than the plane and negative values where the actual value
#' are lower than the plane.
#'
#' @param x A raster or matrix.
#' @return A raster or matrix of the same size as the input with values
#'   equal to the difference between the original and bestfit
#'   plane.
#' @examples
#' library(raster)
#'
#' # import raster image
#' data(orforest)
#'
#' # remove the least squares polynomial plane
#' new_rast <- remove_plane(orforest)
#'
#' # plot
#' plot(new_rast)
#' @export
remove_plane <- function(x) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}

  bfx <- bestfitplane(x)

  errors <- x - bfx # higher = above plane, lower = below plane

  return(errors)
}
