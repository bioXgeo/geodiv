#' Estimate the Areal Autocorrelation Function
#'
#' Calculates the areal autocorrelation function (AACF) as the
#' inverse of the Fourier power spectrum. \code{aacf(x)} returns
#' the AACF in both matrix and image format.
#'
#' @param x An n x n raster object.
#' @return A list containing matrix and image representations
#'   of the AACF. Both matrix and image values are normalized
#'   so that the maximum is equal to 1.
#' @examples
#' library(raster)
#'
#' # import raster image
#' data(normforest)
#'
#' # calculate aacf img and matrix
#' aacf_list <- aacf(normforest)
#'
#' # plot resulting aacf image
#' plot(aacf_list[[2]])
#' @export
aacf <- function(x) {
  if(class(x) != 'RasterLayer') {stop('x must be a raster.')}

  # get raster dimensions
  M <- ncol(x)
  N <- nrow(x)

  # get matrix of values
  zmat <- matrix(getValues(x), ncol = M, nrow = N, byrow = TRUE)

  # if irregular non-na area, cut to biggest square possible
  if (sum(is.na(zmat)) != 0) {
    origin <- c(mean(sp::coordinates(x)[, 1]), mean(sp::coordinates(x)[, 2]))
    potentials <- data.frame(xmin = rep(origin[1], floor(N / 2)),
                             xmax = origin[1],
                             ymin = origin[2],
                             ymax = origin[2])
    potentials$xmin <- origin[1] - (res(x)[1] * seq(1, round(N / 2)))
    potentials$xmax <- origin[1] + (res(x)[1] * seq(1, round(N / 2)))
    potentials$ymin <- origin[2] - (res(x)[2] * seq(1, round(N / 2)))
    potentials$ymax <- origin[2] + (res(x)[2] * seq(1, round(N / 2)))

    potentials$na <- sapply(seq(1, nrow(potentials)), FUN = function(i) {
      xmin <-
        newrast <- crop(x, extent(potentials$xmin[i], potentials$xmax[i], potentials$ymin[i], potentials$ymax[i]))
      return(sum(is.na(getValues(newrast))))
    })

    max_dim <- potentials[max(which(potentials$na <= 0)),]

    x <- crop(x, extent(max_dim$xmin, max_dim$xmax, max_dim$ymin, max_dim$ymax))

    # get raster dimensions
    M <- ncol(x)
    N <- nrow(x)

    # get matrix of values
    zmat <- matrix(getValues(x), ncol = M, nrow = N, byrow = TRUE)
  }

  # create windows to prevent leakage
  wc <- e1071::hanning.window(M)
  wr <- e1071::hanning.window(N)

  # create matrix of weights
  w <- pracma::meshgrid(wc, wr)
  w <- w[[1]] * w[[2]]

  # apply window weights
  zmatw <- zmat * w

  # perform Fourier transform
  ft <- stats::fft(zmatw)

  # calculate power spectrum
  ps <- (abs(ft) ^ 2) / (M * N)

  # autocorrelation function
  af <- Re(stats::fft(ps, inverse = TRUE) / (M * N))
  af_shift <- fftshift(af) # this should be symmetric!

  # normalize to max 1
  af_norm <- af_shift / max(as.numeric(af_shift), na.rm = TRUE)

  # set values of new raster
  af_img <- setValues(x, af_norm)

  return (list(af_norm, af_img))
}


#' Calculate Correlation Length
#'
#' Calculates the smallest and largest distances to specified autocorrelation
#' values (e.g., 0.2) of the areal autocorrelation function (AACF). All 180
#' degrees from the origin of the AACF image are considered for the calculation.
#'
#' @param x A raster.
#' @param threshold A numeric vector containing values between 0 and 1. Indicates
#'   the autocorrelation values to which the rates of decline are measured.
#' @param plot Logical. Defaults to \code{FALSE}. If \code{TRUE}, the AACF and
#'   lines showing the considered directions of autocorrelation from the origin
#'   will be plotted.
#' @return A list containing the minimum and maximum distances from an
#'   autocorrelation value of 1 to the specified autocorrelation values < 1.
#'   Distances are in the units of the x, y coordinates of the raster image. If more
#'   than one threshold value is specified, the order of this list will be
#'   \code{[minval(t1), minval(t2), maxval(t1), maxval(t2)]}.
#' @examples
#' library(raster)
#'
#' # import raster image
#' data(normforest)
#'
#' # crop raster to much smaller area
#' x <- crop(normforest, extent(-123, -122.99, 43, 43.01))
#'
#' # calculate aacf img and matrix
#' aacf_list <- aacf(x)
#'
#' # estimate the fastest/slowest declines to 0.20 and 0.37 (1/e) autocorrelation
#' sclvals <- scl(aacf_list[[2]])
#'
#' # calculate Scl20, the minimum distance to an autocorrelation value of 0.2 in the AACF
#' Scl20 <- sclvals[[1]]
#' @export
scl <- function(x, threshold = c(0.20, 1 / exp(1)), plot = FALSE) {
  if(class(x) != 'RasterLayer') {stop('x must be a raster.')}
  if(class(plot) != 'logical') {stop('plot argument must be TRUE/FALSE.')}
  if(class(threshold) != 'numeric') {stop('threshold must be numeric.')}
  if(sum(threshold < 0) >= 1) {stop('threshold values cannot be less than 0.')}

  # get aacf img
  aacfimg <- aacf(x)[[2]]

  # take amplitude image, cut in half (y direction)
  half_dist <- (ymax(aacfimg) - ymin(aacfimg)) / 2
  ymin <- ymax(aacfimg) - half_dist
  aacfimg <- crop(aacfimg, c(xmin(aacfimg), xmax(aacfimg), ymin, ymax(aacfimg)))

  # get origin of image (actually bottom center)
  origin <- c(mean(sp::coordinates(aacfimg)[, 1]), ymin(aacfimg))

  ### line calculations are taken from the plotrix function draw.radial.line
  # calculate rays extending from origin
  M <- 180
  j <- seq(0, (M - 1))
  alpha <- (pi * j) / M # angles
  px <- c(0, half_dist) # line length
  linex <- unlist(lapply(seq(1, length(alpha)), function(x) origin[1] + px * cos(alpha[x])))
  liney <- unlist(lapply(seq(1, length(alpha)), function(x) origin[2] + px * sin(alpha[x])))
  linelist <- lapply(seq(1, length(linex), 2),
                     FUN = function(i) sp::Lines(sp::Line(cbind(linex[i:(i + 1)], liney[i:(i + 1)])),
                                             ID = paste('l', i, sep = '')))
  lines <- sp::SpatialLines(linelist, proj4string = sp::CRS(sp::proj4string(aacfimg)))

  # plot and calculate amplitude sums along rays
  if(plot == TRUE) {
    plot(aacfimg)
    lines(lines)
  }

  # get values for all points along line
  Aalpha <- list()
  for (i in 1:length(lines)) {
    Aalpha[[i]] <- extract(aacfimg, lines[i], along = TRUE, cellnumbers = TRUE)
  }

  # each line has length = half_dist, with each point approx. 1 pixel apart
  fast_dists <- list()
  for (i in 1:length(threshold)) {
    fast_dists[[i]] <- suppressWarnings(.mindist(threshold[i], Aalpha, aacfimg))
  }

  slow_dists <- list()
  for (i in 1:length(threshold)) {
    slow_dists[[i]] <- suppressWarnings(.maxdist(threshold[i], Aalpha, aacfimg))
  }

  return(c(fast_dists, slow_dists))
}

#' Estimate Minimum Correlation Length
#'
#' Internal function to calculates the minimum distances to specified
#' autocorrelation values (e.g., 0.2) of the areal autocorrelation
#' function (AACF). All 180 degrees from the origin of the AACF image
#' are considered for the calculation.
#'
#' @param threshold A number with a value between 0 and 1. Indicates
#'   the autocorrelation value to which the rates of decline are measured.
#' @param Aalpha An list of dataframes produced by \code{scl()} that contain
#'   the AACF values along lines extending in multiple directions from the
#'   AACF origin (autocorrelation = 1).
#' @param aacfimg A raster of the areal autocorrelation function. This
#'   is the AACF raster split in two in terms of height.
#' @return A list containing the minimum distances from an
#'   autocorrelation value of 1 to the specified autocorrelation value < 1.
#'   Distances are in the units of the x, y coordinates of the raster image.
.mindist <- function(threshold, Aalpha, aacfimg) {
  # get index of minimum <= threshold value
  decay_ind <- list()
  for (j in 1:length(Aalpha)) {
    decay_ind[[j]] <- lapply(Aalpha[[j]], FUN = function(x)
      min(which(x[, 2] <= threshold), na.rm = TRUE))
  }
  decay_ind <- unlist(decay_ind)

  # get distance to minimum index below threshold value
  x <- sp::coordinates(aacfimg)[, 1]
  y <- sp::coordinates(aacfimg)[, 2]
  origin <- c(mean(sp::coordinates(aacfimg)[, 1]), ymin(aacfimg))
  decay_celln <- list()
  for (j in 1:length(Aalpha)) {
    decay_celln[[j]] <- unlist(lapply(Aalpha[[j]], FUN = function(x) x[decay_ind[j], 1]))
  }
  decay_celln <- unlist(decay_celln)
  decay_coords <- data.frame(x = x[decay_celln], y = y[decay_celln])
  decay_coords <- decay_coords[decay_ind != Inf,]
  decay_dist <- min(spatstat::crossdist.default(X = decay_coords$x[!is.na(decay_coords$x)], Y = decay_coords$y[!is.na(decay_coords$x)],
                                      x2 = origin[1], y2 = origin[2]))

  return(decay_dist)
}

#' Estimate Maximum Correlation Length
#'
#' Internal function to calculates the maximum distances to specified
#' autocorrelation values (e.g., 0.2) of the areal autocorrelation
#' function (AACF). All 180 degrees from the origin of the AACF image
#' are considered for the calculation.
#'
#' @param threshold A number with a value between 0 and 1. Indicates
#'   the autocorrelation value to which the rates of decline are measured.
#' @param Aalpha An list of dataframes produced by \code{scl()} that contain
#'   the AACF values along lines extending in multiple directions from the
#'   AACF origin (autocorrelation = 1).
#' @param aacfimg A raster of the areal autocorrelation function. This
#'   is the AACF raster split in two in terms of height.
#' @return A list containing the maximum distances from an
#'   autocorrelation value of 1 to the specified autocorrelation value < 1.
#'   Distances are in the units of the x, y coordinates of the raster image.
.maxdist <- function(threshold, Aalpha, aacfimg) {
  # get index of minimum <= threshold value
  decay_ind <- list()
  for (j in 1:length(Aalpha)) {
    decay_ind[[j]] <- lapply(Aalpha[[j]], FUN = function(x)
      min(which(x[, 2] <= threshold), na.rm = TRUE))
  }
  decay_ind <- unlist(decay_ind)

  # get distance to minimum index below threshold value
  x <- sp::coordinates(aacfimg)[, 1]
  y <- sp::coordinates(aacfimg)[, 2]
  origin <- c(mean(sp::coordinates(aacfimg)[, 1]), ymin(aacfimg))
  decay_celln <- list()
  for (j in 1:length(Aalpha)) {
    decay_celln[[j]] <- unlist(lapply(Aalpha[[j]], FUN = function(x) x[decay_ind[j], 1]))
  }
  decay_celln <- unlist(decay_celln)
  decay_coords <- data.frame(x = x[decay_celln], y = y[decay_celln])
  decay_coords <- decay_coords[decay_ind != Inf,]
  decay_dist <- max(spatstat::crossdist.default(X = decay_coords$x[!is.na(decay_coords$x)], Y = decay_coords$y[!is.na(decay_coords$x)],
                                      x2 = origin[1], y2 = origin[2]))

  return(decay_dist)
}

#' Estimate Texture Aspect Ratio
#'
#' Calculates the texture aspect ratio (Str) at defined autocorrelation
#' values. The texture aspect ratio is the ratio of the fastest to
#' the slowest decay lengths of the autocorrelation function to the
#' defined autocorrelation values.
#'
#' @param x A raster.
#' @param threshold A vector of autocorrelation values with values
#'   between 0 and 1. Indicates the autocorrelation value(s) to
#'   which the rates of decline are measured.
#' @return A vector with length equal to that of \code{threshold}
#'   containing the texture aspect ratio(s) for the input autocorrelation
#'   value(s).
#' @examples
#' library(raster)
#'
#' # import raster image
#' data(normforest)
#'
#' # crop raster to much smaller area
#' x <- crop(normforest, extent(-123, -122.99, 43, 43.01))
#'
#' # estimate the texture aspect ratio for autocorrelation
#' # thresholds of 0.20 and 0.37 (1/e)
#' strvals <- str(x, threshold = c(0.20, 1 / exp(1)))
#'
#' # calculate Str20, the texture aspect ratio for
#' # autocorrelation value of 0.2 in the AACF
#' Str20 <- strvals[[1]]
#' @export
str <- function(x, threshold = c(0.20, 1 / exp(1))) {
  if(class(x) != 'RasterLayer') {stop('x must be a raster.')}
  if(class(threshold) != 'numeric') {stop('threshold must be numeric.')}
  if(sum(threshold < 0) >= 1) {stop('threshold values cannot be less than 0.')}

  sclvals <- scl(x, threshold = threshold, plot = FALSE)

  vals <- list()
  # because the list contains both min/max vals, need double the length
  for (i in 1:length(threshold)) {
    minval <- sclvals[[i]]
    j <- length(threshold) + i
    maxval <- sclvals[[j]]
    vals[[i]] <- minval / maxval
  }

  return(vals)
}
