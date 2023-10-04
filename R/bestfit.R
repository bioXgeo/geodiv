# functions to find the best fit line for forty percent of slope with lowest slope

#' Determines the Slopes Along the Bearing Area Curve
#'
#' Calculates the slopes along the bearing area curve
#' of a raster or matrix. Slopes are determined at points x,
#' from point x - h to x + h.
#'
#' @param x A vector of x values.
#' @param h Spacing before and after each point.
#' 2h is the distance over which slopes are calculated.
#' @param f Bearing area function as calculated with
#' bearing_area.
#' @return A dataframe with the slope for each segment
#' with centerpoint x.
#' @examples
#' # import raster image
#' data(normforest)
#' normforest <- terra::unwrap(normforest)
#'
#' # find the slopes along the bearing area curve
#' ba <- bearing_area(normforest)
#' x <- seq(0, 1, length.out = 100000)
#' slopes <- slopecalc(x = x, h = 0.01, f = ba)
#' @importFrom terra rast
#' @importFrom stats quantile
#' @export
slopecalc <- function(x, h, f) {
  stopifnot('x must be numeric.' = inherits(x, 'numeric'))
  stopifnot('h must be numeric.' = inherits(h, 'numeric'))

  if(sum(class(f) %in% c('ecdf', 'stepfun', 'function')) != 3) {stop('f was not produced with bearing_area function.')}
  if(length(h) > 1) {stop('too many values for h.')}
  if(h >= 1 | h <= 0) {stop('h must be less than 1 and greater than 0.')}

  xplus <- x + h
  xminus <- x - h
  x <- x

  # figure out ends (need distance on both sides of x)
  space <- 1 / length(x)
  end_length <- h / space
  if (end_length < 1) { end_length = 1}
  begin <- ceiling(1 + end_length)
  end <- floor(length(x) - end_length)
  begin1 <- begin - 1
  end1 <- end + 1

  fxh_pos <- seq(1, length(x))
  fxh_neg <- seq(1, length(x))
  fxh_pos[1:begin1] <- (1 - stats::quantile(f, probs = xplus[1:begin1]))
  fxh_neg[1:begin1] <- (1 - stats::quantile(f, probs = x[1:begin1]))
  # variation on newton's difference quotient at far end
  fxh_pos[end1:length(x)] <- (1 - stats::quantile(f, probs = x[end1:length(x)]))
  fxh_neg[end1:length(x)] <- (1 - stats::quantile(f, probs = xminus[end1:length(x)]))
  # symmetric difference quotient everywhere else (99800 points)
  fxh_pos[begin:end] <- (1 - stats::quantile(f, probs = xplus[begin:end]))
  fxh_neg[begin:end] <- (1 - stats::quantile(f, probs = xminus[begin:end]))

  diff_quo <- (fxh_pos - fxh_neg) / (2 * h)
  slopes <- data.frame(slope = diff_quo, x = x)

  return(slopes)
}

#' Determines the Average Slope Along Larger Segments of
#' the Bearing Area Curve
#'
#' Calculates the average slope over every segment
#' of a specified percentage length of the total bearing
#' area curve.
#'
#' @param slopes A dataframe containing all slopes along
#' the bearing area curve, calculated using the slopecalc
#' function.
#' @param l Percentage of the curve over which to calculate
#' mean slope.
#' @return A dataframe with the average slope over segments
#' beginning at specified x locations along the bearing area
#' curve. 'slope' represents the mean slope over the segment,
#' 'xstart' is the beginning x location of the segment, and
#' 'xend' is the concluding x location of the segment.
#' @examples
#' # import raster image
#' data(normforest)
#' normforest <- terra::unwrap(normforest)
#'
#' # find the average slope of segments of the bearing area
#' # curve.
#' ba <- bearing_area(normforest)
#' x <- seq(0, 1, length.out = 10000)
#' slopes <- slopecalc(x = x, h = 0.01, f = ba)
#' slopes_forty <- slopemeans(slopes = slopes, l = 0.4)
#' @importFrom terra rast
#' @export
slopemeans <- function(slopes, l = 0.4) {
  stopifnot('slopes must be a dataframe.' = inherits(slopes, 'data.frame'))
  stopifnot('l must be numeric.' = inherits(l, 'numeric'))

  if(length(l) > 1) {stop('too many values for l.')}
  if(l >= 1 | l <= 0) {stop('l must be less than 1 and greater than 0.')}
  if(sum(names(slopes) %in% c('slope', 'x')) != 2) {stop('incorrect column names for slopes dataframe -- need slope and x.')}

  x <- slopes$x
  slope <- slopes$slope
  length <- length(x) * l
  end_ind <- length(x) - length
  begin_ind <- 1 + length
  xstart <- x[1:end_ind]
  istart <- seq(1, end_ind)
  xend <- x[begin_ind:length(x)]
  iend <- seq(begin_ind, length(x))

  # create matrix of slopes
  xmat <- matrix(c(xstart, xend), ncol = 2)
  mean <- sapply(seq(1, end_ind), function(i) {mean(abs(slope[istart[i]:iend[i]]))})

  data <- data.frame(slope = mean, xstart = xstart, xend = xend)
  return(data)
}
