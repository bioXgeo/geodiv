#' Estimate the Areal Autocorrelation Function
#'
#' Calculates the areal autocorrelation function (AACF) as the
#' inverse of the Fourier power spectrum. \code{aacf(x)} returns
#' the AACF in both matrix and raster format.
#'
#' @param x An n x n raster or matrix.
#' @return A raster or matrix representation
#'   of the AACF. Both raster and matrix values are normalized
#'   so that the maximum is equal to 1.
#' @examples
#' # import raster image
#' data(normforest)
#' normforest <- terra::unwrap(normforest)
#'
#' # calculate aacf img and matrix
#' aacf_out <- aacf(normforest)
#'
#' # plot resulting aacf image
#' terra::plot(aacf_out)
#' @importFrom terra rast crds crop setValues crs
#' @importFrom e1071 hanning.window
#' @importFrom pracma meshgrid
#' @importFrom stats fft
#' @export
aacf <- function(x) {
  stopifnot('x must be a raster or matrix.' = inherits(x, c('RasterLayer', 'matrix', 'SpatRaster')))

  # get raster dimensions
  M <- ncol(x)
  N <- nrow(x)

  data_type <- if(class(x)[1] == 'matrix') {'matrix'} else if (class(x)[1] == 'RasterLayer') {'RasterLayer'} else {'SpatRaster'}

  # convert matrix to raster if necessary (equal area)
  if (data_type == 'matrix') {
    x <- rast(x)
    terra::crs(x) <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  }

  # convert RasterLayer to SpatRaster
  if (data_type == 'RasterLayer') {
    x <- rast(x)
  }

  # get matrix of values
  zmat <- matrix(x[], ncol = M, nrow = N, byrow = TRUE)

  # if irregular non-na area, cut to biggest square possible
  if (sum(is.na(zmat)) != 0) {
    origin <- c(mean(crds(x, na.rm = FALSE)[, 1]), mean(crds(x, na.rm = FALSE)[, 2]))
    potentials <- data.frame(xmin = rep(origin[1], floor(N / 2)),
                             xmax = origin[1],
                             ymin = origin[2],
                             ymax = origin[2])
    potentials$xmin <- origin[1] - (res(x)[1] * seq(1, floor(N / 2)))
    potentials$xmax <- origin[1] + (res(x)[1] * seq(1, floor(N / 2)))
    potentials$ymin <- origin[2] - (res(x)[2] * seq(1, floor(N / 2)))
    potentials$ymax <- origin[2] + (res(x)[2] * seq(1, floor(N / 2)))

    potentials$na <- sapply(seq(1, nrow(potentials)), FUN = function(i) {
      xmin <-
      newrast <- terra::crop(x, ext(potentials$xmin[i], potentials$xmax[i], potentials$ymin[i], potentials$ymax[i]))
      return(sum(is.na(newrast[])))
    })

    max_dim <- potentials[max(which(potentials$na <= 0)),]

    if (sum(is.na(max_dim$xmin)) != 0) {
      zmat <- zmat
    } else if (sum(is.na(max_dim$xmin)) == 0) {
      x <- terra::crop(x, ext(max_dim$xmin, max_dim$xmax, max_dim$ymin, max_dim$ymax))

      # get raster dimensions
      M <- ncol(x)
      N <- nrow(x)

      # get matrix of values
      zmat <- matrix(x[], ncol = M, nrow = N, byrow = TRUE)
    }
  }

  if (sum(is.na(zmat)) == length(zmat)) {
    return(NA)
  } else if (sum(is.na(zmat)) != length(zmat)) {

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
    af <- base::as.matrix(base::Re(stats::fft(ps, inverse = TRUE) / (M * N)))
    af_shift <- geodiv::fftshift(af) # this should be symmetric!

    # normalize to max 1
    af_norm <- af_shift / max(as.numeric(af_shift), na.rm = TRUE)

    if (data_type == 'RasterLayer' | data_type == 'SpatRaster') {
      # set values of new raster
      af_img <- terra::setValues(x, af_norm)
      return(af_img)
    } else {
      return(af_norm)
    }
  }
}

#' Calculate Correlation Length
#'
#' Calculates the smallest and largest distances to specified autocorrelation
#' values (e.g., 0.2) of the areal autocorrelation function (AACF). All 180
#' degrees from the origin of the AACF image are considered for the calculation.
#'
#' @param x A raster or matrix.
#' @param threshold A numeric vector containing values between 0 and 1. Indicates
#'   the autocorrelation values to which the rates of decline are measured.
#' @param create_plot Logical. Defaults to \code{FALSE}. If \code{TRUE}, the AACF and
#'   lines showing the considered directions of autocorrelation from the origin
#'   will be plotted.
#' @return A list containing the minimum and maximum distances from an
#'   autocorrelation value of 1 to the specified autocorrelation values < 1.
#'   Distances are in the units of the x, y coordinates of the raster image. If more
#'   than one threshold value is specified, the order of this list will be
#'   \code{[minval(t1), minval(t2), maxval(t1), maxval(t2)]}.
#' @examples
#' # import raster image
#' data(normforest)
#' normforest <- terra::unwrap(normforest)
#'
#' # calculate Scl20, the minimum distance to an autocorrelation value of 0.2 in the AACF
#' Scl20 <- scl(normforest)[1]
#' @importFrom dplyr %>% summarize group_by
#' @importFrom terra crop rast xmin xmax ymin ymax unwrap crs values
#' @importFrom sf st_geometry_type
#' @importFrom rlang .data
#' @export
scl <- function(x, threshold = c(0.20, 1 / exp(1)), create_plot = FALSE) {
  stopifnot('x must be a raster or matrix.' = inherits(x, c('RasterLayer', 'matrix', 'SpatRaster')))
  stopifnot('create_plot argument must be TRUE/FALSE.' = inherits(create_plot, 'logical'))
  stopifnot('threshold must be numeric.' = inherits(threshold, 'numeric'))

  if(sum(threshold < 0) >= 1) {stop('threshold values cannot be less than 0.')}

  # get aacf img
  aacfimg <- aacf(x)

  if (!(class(aacfimg)[1] %in% c('matrix', 'RasterLayer', 'SpatRaster')) | sum(is.na(aacfimg)[]) == length(aacfimg)) {
    return(c(NA, NA, NA, NA))
  } else if (class(aacfimg)[1] %in% c('matrix', 'RasterLayer', 'SpatRaster')) {

    data_type <- if(class(x)[1] == 'matrix') {'matrix'} else if (class(x)[1] == 'RasterLayer') {'RasterLayer'} else {'SpatRaster'}

    if (data_type %in% c('RasterLayer', 'SpatRaster')) {
      # take amplitude image, cut in half (y direction)
      half_dist <- (ymax(aacfimg) - ymin(aacfimg)) / 2
      ymin <- ymax(aacfimg) - half_dist
      aacfimg <- terra::crop(aacfimg, c(xmin(aacfimg), xmax(aacfimg), ymin, ymax(aacfimg)))

      # get origin of image (actually bottom center)
      origin <- c(mean(crds(aacfimg, na.rm = FALSE)[, 1]), ymin(aacfimg))

    } else {
      # take amplitude image, cut in half (y direction)
      half_dist <- nrow(aacfimg) / 2
      ymin <- round(half_dist)
      aacfimg <- aacfimg[0:ymin, 0:ncol(aacfimg)]

      # get origin of image (actually bottom center)
      origin <- c(nrow(aacfimg), round(ncol(aacfimg) / 2))
    }

    # convert matrix to raster if necessary
    if (data_type == 'matrix') {
      aacf_rast <- rast(aacfimg)
      terra::crs(aacf_rast) <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
      aacfimg <- aacf_rast
    }

    ### line calculations are taken from the plotrix function draw.radial.line
    # calculate rays extending from origin
    if (create_plot == TRUE) {
      M <- 180
      j <- seq(0, (M - 1))
      alpha <- (pi * j) / M # angles
      px <- c(0, half_dist) # line length
      linex <- unlist(lapply(seq(1, length(alpha)), function(x) origin[1] + px * cos(alpha[x])))
      liney <- unlist(lapply(seq(1, length(alpha)), function(x) origin[2] + px * sin(alpha[x])))
      linelist <- lapply(seq(1, length(linex), 2), FUN = function(i) {
        data.frame(x = linex[i:(i + 1)], y = liney[i:(i + 1)],
                                                           id = paste('l', i, sep = ''))})
      linelist <- bind_rows(linelist)
      lines <- linelist %>%
        st_as_sf(coords = c("x", "y"), na.fail = FALSE, crs = terra::crs(aacfimg)) %>%
        group_by(.data$id) %>%
        summarize()
      multi_inds <- which(st_geometry_type(lines$geometry) == 'MULTIPOINT')
      lines <- st_cast(lines[multi_inds, ], "LINESTRING")

      # plot and calculate amplitude sums along rays
      plot(aacfimg)
      lines(lines)
    }

    # calculate distances from center to all other points
    dist_rast <- aacfimg
    terra::values(dist_rast) <- NA
    center <- ceiling(dim(x) / 2)
    nce <- ifelse(ncol(aacfimg) / 2 == round(ncol(aacfimg) / 2), 1, 0)
    dist_rast[nrow(aacfimg), center[2] + nce] <- 1
    dist_rast <- distance(dist_rast)

    # each line has length = half_dist, with each point approx. 1 pixel apart
    fast_dists <- list()
    for (i in 1:length(threshold)) {
      fast_dists[[i]] <- suppressWarnings(.mindist(threshold[i], aacfimg, dist_rast))
    }

    slow_dists <- list()
    for (i in 1:length(threshold)) {
      slow_dists[[i]] <- suppressWarnings(.maxdist(threshold[i], aacfimg, dist_rast))
    }

    return(c(unlist(fast_dists), unlist(slow_dists)))
  }
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
#' @param aacfimg A raster of the areal autocorrelation function. This
#'   is the AACF raster split in two in terms of height.
#' @param distimg A raster of distances to all pixels from the center of the
#'   original image. Distances are in meters if original raster was
#'   unprojected, and are in map units (usually meters) if raster was projected
#'   (see raster::distance documentation for more details).
#' @return A list containing the minimum distances from an
#'   autocorrelation value of 1 to the specified autocorrelation value < 1.
#'   Distances are meters if original raster was unprojected, and are in
#'   map units (usually meters) if raster was projected (see
#'   raster::distance documentation for more details).
.mindist <- function(threshold, aacfimg, distimg) {
  # get indices where value <= threshold value
  decay_ind <- which(aacfimg[] <= threshold)

  # find distances associated with aacf values less than threshold
  decay_dist <- distimg[decay_ind]

  # find minimum distance where aacf is less than threshold
  min_decay <- min(decay_dist, na.rm = TRUE)

  return(min_decay)
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
#' @param aacfimg A raster of the areal autocorrelation function. This
#'   is the AACF raster split in two in terms of height.
#' @param distimg A raster of distances to all pixels from the center of the
#'   original image. Distances are in meters if original raster was
#'   unprojected, and are in map units (usually meters) if raster was projected
#'   (see raster::distance documentation for more details).
#' @return A list containing the maximum distances from an
#'   autocorrelation value of 1 to the specified autocorrelation value < 1.
#'   Distances are meters if original raster was unprojected, and are in
#'   map units (usually meters) if raster was projected (see
#'   raster::distance documentation for more details).
.maxdist <- function(threshold, aacfimg, distimg) {
  # get indices where value <= threshold value
  decay_ind <- which(aacfimg[] <= threshold)

  # find distances associated with aacf values less than threshold
  decay_dist <- distimg[decay_ind]

  # find maximum distance where aacf is less than threshold
  max_decay <- max(decay_dist, na.rm = TRUE)

  return(max_decay)
}

#' Estimate Texture Aspect Ratio
#'
#' Calculates the texture aspect ratio (Str) at defined autocorrelation
#' values. The texture aspect ratio is the ratio of the fastest to
#' the slowest decay lengths of the autocorrelation function to the
#' defined autocorrelation values.
#'
#' @param x A raster or matrix.
#' @param threshold A vector of autocorrelation values with values
#'   between 0 and 1. Indicates the autocorrelation value(s) to
#'   which the rates of decline are measured.
#' @return A vector with length equal to that of \code{threshold}
#'   containing the texture aspect ratio(s) for the input autocorrelation
#'   value(s).
#' @examples
#' # import raster image
#' data(normforest)
#' normforest <- terra::unwrap(normforest)
#'
#' # estimate the texture aspect ratio for autocorrelation
#' # thresholds of 0.20 and 0.37 (1/e)
#' strvals <- stxr(normforest, threshold = c(0.20, 1 / exp(1)))
#'
#' # calculate Str20, the texture aspect ratio for
#' # autocorrelation value of 0.2 in the AACF
#' Str20 <- strvals[1]
#' @importFrom terra rast
#' @export
stxr <- function(x, threshold = c(0.20, 1 / exp(1))) {
  stopifnot('x must be a raster or matrix.' = inherits(x, c('RasterLayer', 'matrix', 'SpatRaster')))
  stopifnot('threshold must be numeric.' = inherits(threshold, 'numeric'))

  if(sum(threshold < 0) >= 1) {stop('threshold values cannot be less than 0.')}

  sclvals <- scl(x, threshold = threshold, create_plot = FALSE)

  vals <- list()
  # because the list contains both min/max vals, need double the length
  for (i in 1:length(threshold)) {
    minval <- sclvals[[i]]
    j <- length(threshold) + i
    maxval <- sclvals[[j]]
    vals[[i]] <- minval / maxval
  }

  return(unlist(vals))
}
