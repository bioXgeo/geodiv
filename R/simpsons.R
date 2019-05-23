#' Area Above the Bearing Area Curve
#'
#' Calculates the area above the bearing area curve from
#' points \code{a} to \code{b}. If a box is drawn around
#' a function with the upper-left at \code{a} and the
#' bottom-right at \code{b}, this function extracts the area
#' above the function within the box.
#'
#' The area under the curve used to calculate area above the
#' curve is calculated as the numerical integral
#' of the Bearing Area function from \code{a} to \code{b}
#' using the trapezoid rule with n subdivisions. Assume
#' \code{a < b} and \code{n} is a positive integer.
#'
#' @param f The function for the Bearing Area curve produced by
#'   \code{stats::ecdf()}.
#' @param a Numeric. The left x boundary.
#' @param b Numeric. The right x boundary.
#' @param n Numeric. The number of subdivisions along the function
#'   line.
#' @return A numeric value representing the area above the curve with
#'   x bounds \code{a} and \code{b}.
#' @examples
#' library(raster)
#'
#' # import raster image
#' data(normforest)
#'
#' # basic values
#' z <- getValues(normforest)
#'
#' # calculate cumulative probability density function of surface 'height' (= ndvi)
#' mod <- ecdf((1 - z))
#'
#' # valley fluid retention index = void volume in 'valley' zone
#' Svi <- area_above(f = mod, b = 1, a = 0.8, n = 500)
#' @export
area_above <- function(f, a, b, n = 100) {
  if(('function' %in% class(f)) != TRUE) {stop('f must be a function.')}
  if(class(a) != 'numeric') {stop('a must be numeric.')}
  if(class(b) != 'numeric') {stop('b must be numeric.')}
  if(class(n) != 'numeric') {stop('n must be numeric.')}
  if(length(a) > 1) {stop('too many values supplied to a.')}
  if(length(b) > 1) {stop('too many values supplied to b.')}
  if(length(n) > 1) {stop('too many values supplied to n.')}
  if(n <= 0) {stop('n must be greater than 0.')}
  if(a >= b) {stop('b must be greater than a.')}

  h <- (b - a) / n # sub-interval width
  x <- seq(a, b, by = h)
  y <- (1 - stats::quantile(f, probs = x)) # get y-values of inverse cdf function

  # if y is negative, shift the whole thing up so that min(y) = 0
  if (sum(y < 0) >= 1) {
    y <- y + abs(min(y, na.rm = TRUE))
  }

  # area under curve from simpson's rule
  s <- (h / 3) * (y[[1]] + sum(4 * y[seq(2, n - 1, by = 2)]) +
                                sum(2 * y[seq(3, n - 1, by = 2)]) +
                                y[[n]])

  # get inverse of s for actual area above curve
  area_above <- ((max(y) - min(y)) * (max(x) - min(x))) - s

  return(area_above)
}

#' Simpson's Rule Empirical Area Under a Curve
#'
#' Calculates the area below a curve from
#' points \code{a} to \code{b}. This function is provided
#' for general use.
#'
#' Note that if y-values are negative, this returns the area
#' above the function line.
#'
#' @param f A function.
#' @param a Numeric. The left x boundary.
#' @param b Numeric. The right x boundary.
#' @param n Numeric. The number of subdivisions along the function
#'   line.
#' @return A numeric value representing the area under the curve with
#'   x bounds \code{a} and \code{b}.
#' @examples
#' library(raster)
#'
#' # import raster image
#' data(normforest)
#'
#' # basic values
#' z <- getValues(normforest)
#'
#' # calculate cumulative probability density function of surface 'height' (= ndvi)
#' mod <- ecdf((1 - z))
#'
#' # calculate integral
#' int_area <- simpsons(f = mod, b = 1, a = 0.8, n = 500)
#' @export
simpsons <- function(f, a, b, n = 100) {
  if(('function' %in% class(f)) != TRUE) {stop('f must be a function.')}
  if(class(a) != 'numeric') {stop('a must be numeric.')}
  if(class(b) != 'numeric') {stop('b must be numeric.')}
  if(class(n) != 'numeric') {stop('n must be numeric.')}
  if(length(a) > 1) {stop('too many values supplied to a.')}
  if(length(b) > 1) {stop('too many values supplied to b.')}
  if(length(n) > 1) {stop('too many values supplied to n.')}
  if(n <= 0) {stop('n must be greater than 0.')}
  if(a >= b) {stop('b must be greater than a.')}

  h <- (b - a) / n # sub-interval width
  x <- seq(a, b, by = h)

  # get y-values of inverse cdf function
  y <- (1 - stats::quantile(f, probs = x))

  # if y is negative, shift the whole thing up so that min(y) = 0
  if (sum(y < 0) >= 1) {
    y <- y + abs(min(y, na.rm = TRUE))
  }

  s <- (h / 3) * (y[[1]] + sum(4 * y[seq(2, n - 1, by = 2)]) +
                                sum(2 * y[seq(3, n - 1, by = 2)]) +
                                y[[n]])
  return(s)
}
