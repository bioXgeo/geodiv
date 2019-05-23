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
#' image.
#'
#' @param x A raster.
#' @param order Numeric. Indicates the polynomial order to be fit.
#' @return A matrix of the same size as the raster with values
#'   predicted from the polynomial fit.
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
  if(class(x) != 'RasterLayer') {stop('x must be a raster.')}
  if(class(order) == 'numeric') {
    warning('numeric order, converting to integer.')
    order <- as.integer(floor(order))}
  if(class(order) != 'integer') {stop('order must be integer.')}
  if(length(order) > 1) {stop('too many values supplied to order.')}
  if(order <= 0) {stop('order must be >= 1.')}

  # extract coordinates and values
  xcoord <- sp::coordinates(x)[, 1]
  ycoord <- sp::coordinates(x)[, 2]
  z <- getValues(x)

  # fit least squares polynomial with order = order
  surfmod <- spatial::surf.ls(np = order, xcoord[!is.na(z)], ycoord[!is.na(z)], z[!is.na(z)])

  # predict polynomial model over raster
  surfvals <- matrix(predict(surfmod, xcoord, ycoord), nrow = nrow(x), ncol = ncol(x), byrow = TRUE)

  return(surfvals)
}

#' Finds the Best Fit Polynomial Plane
#'
#' Finds the best fit polynomial plane for a raster image. This
#' function tests least squares polynomial fits with orders of
#' 2 - 3 and determines which order minimizes the error when the
#' fit is subtracted from the original image.
#'
#' @param x A raster.
#' @return A raster of the same size as the input with values
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
  if(class(x) != 'RasterLayer') {stop('x must be a raster.')}

  # fit least squares plane for polynomials from orders 2 - 3
  mods <- lapply(seq(2, 3), FUN = function(i) fitplane(x, order = i))

  # convert raster to matrix
  xmat <- matrix(x, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)

  # calculate raster errors from best fit planes for each order tested
  errlist <- lapply(seq(2, 3), FUN = function(i) xmat - mods[[(i - 1)]])
  meanerr <- lapply(seq(2, 3), FUN = function(i) mean(errlist[[(i - 1)]], na.rm = TRUE))

  # determine which order is best
  bestfit <- which(as.numeric(meanerr) == min(as.numeric(meanerr), na.rm = TRUE))

  # if more than one, get first
  if (length(bestfit) > 1) {
    bestfit <- min(bestfit, na.rm = TRUE)
  }

  # fill in raster with best fit values
  bfx <- setValues(x, mods[[bestfit]])

  print(paste('Order of polynomial that minimizes errors: ', bestfit, sep = ''))
  return(bfx)
}

#' Removes the Best Fit Polynomial Plane from a Raster
#'
#' Finds the best fit polynomial plane for a raster image and
#' subtracts it from the actual raster values. The remaining
#' raster has positive values where the actual values are higher
#' than the plane and negative values where the actual value
#' are lower than the plane.
#'
#' @param x A raster.
#' @return A raster of the same size as the input with values
#'   equal to the difference between the original and bestfit
#'   plane rasters.
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
  if(class(x) != 'RasterLayer') {stop('x must be a raster.')}

  bfx <- bestfitplane(x)

  errors <- x - bfx # higher = above plane, lower = below plane

  return(errors)
}
