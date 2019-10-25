#' Calculates the Rotated Bearing Area Curve
#'
#' Finds a rotated version of the Bearing Area (Abbott-Firestone)
#' curve from a raster or matrix. The resulting function should be
#' rotated 90 degrees clockwise to get the actual Bearing
#' Area curve.
#'
#' @param x A raster or matrix.
#' @return A function describing the rotated Bearing Area curve.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # find the rotated Bearing Area curve.
#' ba_func <- bearing_area(normforest)
#'
#' # rotate the values and re-plot
#' xval <- environment(ba_func)$y
#' yval <- (1 - environment(ba_func)$x)
#' plot(yval ~ xval)
#' @export
bearing_area <- function(x) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}

  if (class(x) == 'RasterLayer') {
    z <- getValues(x)
  } else {
    z <- x
  }

  # basic values
  N <- length(z)
  s <- stats::sd(z)
  zbar <- mean(z, na.rm = TRUE)

  f <- stats::ecdf(1 - z)

  return(f)
}

#' Plots the Bearing Area Curve
#'
#' Calculates and plots the Bearing Area curve for a raster
#' or matrix using the \code{bearing_area()} function (with correctly
#' rotated results).
#'
#' If \code{divisions = TRUE}, the lines representing the
#' best fit line to the flattest 40 percent of the curve will be
#' shown, as well as both the x and y interception points
#' of that line.
#'
#' @param x A raster or matrix.
#' @param divisions Logical, defaults to \code{FALSE}. If
#'   \code{TRUE}, divisions of the curve will be plotted.
#'   See details section for more information.
#' @return Plots the Bearing Area curve.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # plot the bearing area curve
#' plot_ba_curve(normforest, divisions = TRUE)
#' @export
plot_ba_curve <- function(x, divisions = FALSE) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}
  if(class(divisions) != 'logical') {stop('divisions argument must be TRUE/FALSE.')}

  f <- bearing_area(x)

  xval <- environment(f)$y
  yval <- (1 - environment(f)$x)

  plot(yval ~ xval)

  if (divisions == TRUE) {
    line_fit <- find_flat(x, perc = 0.4)

    graphics::abline(line_fit[[1]]$coefficients[[1]], line_fit[[1]]$coefficients[[2]], col = 'blue')
    graphics::abline(h = line_fit[[3]], col = 'red')
    graphics::abline(h = line_fit[[4]], col = 'red')
    graphics::abline(v = line_fit[[5]], col = 'green')
    graphics::abline(v = line_fit[[6]], col = 'green')
  }
}

#' Finds the Flattest Part of the Bearing Area Curve
#'
#' Locates the flattest x percentage of the Bearing Area
#' curve. Meant to locate the flattest 40 percent of the
#' Bearing Area curve as used in several roughness parameter
#' calculations.
#'
#' @param x A raster or matrix.
#' @param perc Numeric between 0 and 1. The percentage of
#'   the curve over which to fit the line.
#' @return A list containing the equation for the best fit
#'   line, the predicted values from that line, the high
#'   and low y-intercept values for the intersection points
#'   of the line with the Bearing Area curve, and the high
#'   and low x-intercept values for the intersection points
#'   of the line with the Bearing Area curve.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # locate the flattest 40% of the bearing area curve
#' line_data <- find_flat(normforest, perc = 0.4)
#'
#' # extract the equation of the line
#' bf_line <- line_data[[1]]
#' @export
find_flat <- function(x, perc = 0.4) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}
  if(class(perc) != 'numeric') {stop('perc must be numeric.')}
  if(length(perc) > 1) {stop('too many values supplied to perc.')}
  if(perc > 1 | perc < 0) {stop('perc must be between 0 and 1.')}

  f <- bearing_area(x)

  xval <- environment(f)$y
  yval <- (1 - environment(f)$x)

  # find 40% of curve with least decline
  # use symmetric difference quotient to estimate the derivative at evenly spaced points
  # then find 40% consecutive section with lowest mean slope
  even_x <- seq(0, 1, length.out = 10000)
  even_y <- (1 - stats::quantile(f, probs = even_x))
  forty_length <- perc * length(even_x)
  h <- 0.001
  slopes <- slopecalc(even_x, h, f = f) # calculate slope at every point
  means <- slopemeans(slopes) # calculate averages for each 40% segment

  # x value of start of 40% section with smallest decline
  slope_min <- means[means$slope == min(means$slope),]

  # calculate least-squares line for 40% of curve with smallest decline (lowest slope)
  lm_data <- data.frame(x = xval[xval >= slope_min$xstart & xval <= slope_min$xend],
                        y = yval[xval >= slope_min$xstart & xval <= slope_min$xend])
  ls_line <- stats::lm(y ~ x, data = lm_data)

  # get value of ls line between 0 and 1
  pred_data <- tibble::remove_rownames(data.frame(x = even_x, y = even_y))
  pred_data$y <- stats::predict(ls_line, newdata = pred_data)

  # what is the ls line y-value at x = 0, x = 1?
  ls_int_high <- pred_data$y[pred_data$x == 0]
  ls_int_low <- pred_data$y[pred_data$x == 1]

  # Smr1/Smr2 = x values that correspond to cdf y values at ls_int_high/low
  Smr1 <- f(1 - ls_int_high)
  Smr2 <- f(1 - ls_int_low)

  return(list(ls_line, pred_data, ls_int_high, ls_int_low, Smr1, Smr2))
}

#' Value of the Bearing Area Curve at a Specified Value
#'
#' Determines the value of the bearing area curve for a
#' specific value along the x-axis (\code{xval}).
#'
#' @param x A raster or matrix.
#' @param xval Numeric value along the x-axis.
#' @return A numeric value of the bearing area function
#'   corresponding to \code{xval}.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # determine the bearing area function value
#' # corresponding to an x value of 0.4
#' val <- height_ba(normforest, 0.4)
#' @export
height_ba <- function(x, xval) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}
  if(class(xval) != 'numeric') {stop('xval must be numeric.')}
  if(length(xval) > 1) {stop('too many values supplied to xval.')}
  if(xval > 1 | xval < 0) {stop('xval must be between 0 and 1.')}

  f <- bearing_area(x)

  val <- (1 - stats::quantile(f, probs = c(xval))[[1]])

  return(val)
}

#' Height Intervals of the Bearing Area Curve
#'
#' Determines the height interval (height distance) for
#' points along the bearing area curve as defined by
#' their x values.
#'
#' @param x A raster or matrix.
#' @param low Numeric value along the x-axis corresponding
#'   to the lowest value of interest along the x-axis.
#' @param high Numeric value along the y-axis corresponding
#'   to the highest value of interest along the x-axis.
#' @return A numeric value of the difference in height of
#'   the y values along the bearing area curve corresponding
#'   to the specified x values.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # determine the 10-40% height interval of the
#' # bearing area curve
#' val <- sdc(normforest, 0.1, 0.4)
#' @export
sdc <- function(x, low, high) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}
  if(class(low) != 'numeric') {stop('low value must be numeric.')}
  if(length(low) > 1) {stop('too many values supplied to low.')}
  if(low > 1 | low < 0) {stop('low value must be between 0 and 1.')}
  if(class(high) != 'numeric') {stop('high value must be numeric.')}
  if(length(high) > 1) {stop('too many values supplied to high.')}
  if(high > 1 | high < 0) {stop('high value must be between 0 and 1.')}
  if(high <= low) {stop('high value must be greater than low value.')}


  val_low <- height_ba(x, low)
  val_high <- height_ba(x, high)

  val <- val_low - val_high

  return(val)
}

#' Surface Bearing Index
#'
#' Determines the surface bearing index (Sbi), calculated as the ratio
#' of root mean square roughness (Sq) to height at 5\%
#' of bearing area (z05).
#'
#' @param x A raster or matrix.
#' @return A numeric value representing the surface bearing index.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # determine the surface bearing index
#' Sbi <- sbi(normforest)
#' @export
sbi <- function(x) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}

  Sq <- sq(x)
  z05 <- height_ba(x, 0.05)

  val <- Sq / z05

  return(val)
}

#' Valley Fluid Retention Index
#'
#' Determines the valley fluid retention index (Svi). This
#' value is the void volume (area under the bearing area
#' curve) in the 'valley' zone. See Figure 2a from Kedron
#' et al. (2018) for more details.
#'
#' @param x A raster or matrix.
#' @return A numeric value representing the valley fluid
#' retention index.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # determine the valley fluid retention index
#' Svi <- svi(normforest)
#' @export
svi <- function(x) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}

  f <- bearing_area(x)

  val <- area_above(f = f, b = 1, a = 0.8, n = 500)

  return(val)
}

#' Core Fluid Retention Index
#'
#' Determines the core fluid retention index (Sci). This
#' value is the void volume (area under the bearing area
#' curve) in the 'core' zone. See Figure 2a from Kedron
#' et al. (2018) for more details.
#'
#' @param x A raster or matrix.
#' @return A numeric value representing the core fluid
#' retention index.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # determine the core fluid retention index
#' Sci <- sci(normforest)
#' @export
sci <- function(x) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}

  f <- bearing_area(x)

  core_above <- area_above(f = f, b = 1, a = 0.05, n = 1000)

  # remove the valley zone to get the core zone
  Svi <- svi(x)
  val <- core_above - Svi

  return(val)
}

#' Core Roughness Depth
#'
#' Determines the core roughness depth (Sk), the
#' height difference between y values of the
#' intersection points of the least mean square line
#' fit to the flattest 40\% of the bearing area curve.
#' See Figure 2a from Kedron et al. (2018) for more details.
#'
#' @param x A raster.
#' @return A numeric value representing the core roughness
#'   depth of the image.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # determine the core roughness depth
#' Sk <- sk(normforest)
#' @export
sk <- function(x) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}

  line_info <- find_flat(x, perc = 0.4)

  ls_int_high <- line_info[[3]]
  ls_int_low <- line_info[[4]]

  val <- abs(ls_int_high - ls_int_low)

  return(val)
}

#' Reduced Valley Depth
#'
#' Determines the reduced valley depth (Svk), the
#' height difference between y value of the lowest
#' intersection point of the least mean square line
#' fit to the flattest 40\% of the bearing area curve and
#' the minimum y value of the bearing area curve.
#' See Figure 2a from Kedron et al. (2018) for more details.
#'
#' @param x A raster or matrix.
#' @return A numeric value representing the reduced valley depth.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # determine the reduced valley depth
#' Svk <- svk(normforest)
#' @export
svk <- function(x) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}

  # find the bearing area curve
  f <- bearing_area(x)

  # find the flattest 40% of the bearing area curve
  line_info <- find_flat(x, perc = 0.4)

  smr2 <- line_info[[6]]

  val <- abs((1 - stats::quantile(f, probs = 1)) - (1 - stats::quantile(f, probs = smr2)))[[1]]

  return(val)
}

#' Reduced Peak Height
#'
#' Determines the reduced peak height (Spk), the
#' height difference between the maximum y value of the
#' bearing area curve and the y value of the highest
#' intersection point of the least mean square line
#' fit to the flattest 40\% of the bearing area curve.
#' See Figure 2a from Kedron et al. (2018) for more details.
#'
#' @param x A raster or matrix.
#' @return A numeric value representing the reduced peak height.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # determine the reduced peak height
#' Spk <- spk(normforest)
#' @export
spk <- function(x) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}

  # find the bearing area curve
  f <- bearing_area(x)

  # find the flattest 40% of bearing area curve
  line_info <- find_flat(x, perc = 0.4)

  smr1 <- line_info[[5]]

  val <- abs((1 - stats::quantile(f, probs = 0)) - (1 - stats::quantile(f, probs = smr1)))[[1]]

  return(val)
}
