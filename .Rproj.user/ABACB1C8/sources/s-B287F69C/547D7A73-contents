#' Texture Direction Metrics
#'
#' Calculates the angle of dominating texture and the texture
#' direction index of the Fourier spectrum image calculated
#' from a raster image (see Kedron et al. 2018).
#'
#' @param x A raster or matrix.
#' @param plot Logical. If \code{TRUE}, returns a plot of the
#'   amplitude spectrum with lines showing directions in which
#'   amplitude is summed for the Std and Stdi calculations.
#'   Plotting is not possible when input is a matrix.
#' @return A vector containing numeric values for the angle of
#'   dominating texture and the texture direction index.
#' @examples
#' library(raster)
#'
#' # import raster image
#' data(normforest)
#'
#' # calculate Std and Stdi
#' stdvals <- std(normforest)
#'
#' # extract each value
#' Std <- stdvals[1]
#' Stdi <- stdvals[2]
#' @export
std <- function(x, plot = FALSE) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}
  if(class(plot) != 'logical') {stop('plot must be logical.')}

  # get raster dimensions
  M <- ncol(x)
  N <- nrow(x)

  data_type <- if(class(x) == 'matrix') {'matrix'} else {'RasterLayer'}

  # convert matrix to raster if necessary (equal area)
  if (data_type == 'matrix') {
    x <- raster(x)
    extent(x) <- c(0, ncol(x), 0, nrow(x))
    crs(x) <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  }

  # get matrix of values
  zmat <- matrix(getValues(x), ncol = M, nrow = N, byrow = TRUE)

  # if irregular non-na area, cut to biggest square possible
  if (sum(is.na(zmat)) != 0) {
    origin <- c(mean(sp::coordinates(x)[, 1]), mean(sp::coordinates(x)[, 2]))
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

  # get fourier transform
  # complex spectrum from fast fourier transform
  ft <- fft(zmat)
  ft_shift <- geodiv::fftshift(ft)

  # amplitude spectrum
  amplitude <- sqrt((Re(ft_shift) ^ 2) + (Im(ft_shift) ^ 2))

  # create amplitude image
  amp_img <- setValues(x, amplitude)

  # take amplitude image, cut in half (y direction)
  half_dist <- (ymax(amp_img) - ymin(amp_img)) / 2
  ymin <- ymax(amp_img) - half_dist
  amp_img <- crop(amp_img, c(xmin(amp_img), xmax(amp_img), ymin, ymax(amp_img)))

  # get origin of image (actually bottom center)
  origin <- c(mean(sp::coordinates(amp_img)[, 1]), ymin(amp_img))

  if(plot == TRUE) {
    if (data_type == 'matrix') {
      print("cannot draw lines for object of class 'matrix.'")
    }

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
    lines <- sp::SpatialLines(linelist, proj4string = sp::CRS(sp::proj4string(amp_img)))

    # plot and calculate amplitude sums along rays
    plot(amp_img)
    lines(lines)
  }

  # replace: find distances and angle to each pixel, sum all values at given angle

  # calculate distances from center to all other points
  center <- ceiling(dim(x) / 2)
  nce <- ifelse(ncol(amp_img) / 2 == round(ncol(amp_img) / 2), 1, 0)
  dist_rast <- amp_img
  values(dist_rast) <- NA
  dist_rast[nrow(amp_img), center[2] + nce] <- 1
  dist_rast <- distance(dist_rast)
  min_dist <- min(max(dist_rast[nrow(dist_rast),]), max(dist_rast[center[2] + nce]))

  # calculate angles from center to all points
  angle_rast <- amp_img
  values(angle_rast) <- NA
  nce <- ifelse(ncol(amp_img) / 2 == round(ncol(amp_img) / 2), 1, 0)
  center_ind <- cellFromRowCol(angle_rast, nrow(amp_img), center[2])
  coords <- coordinates(amp_img)
  angles <- atan2(coords[, 2] - coords[center_ind, 2], coords[, 1] - coords[center_ind, 1])
  angles <- pracma::rad2deg(angles)
  angle_rast <- setValues(angle_rast, angles)

  # angles at which you want to sum values
  alpha <- seq(0, 180, 0.5)

  # sum values along lines
  Aalpha <- list()
  for (i in 1:length(alpha)) {
    # get indices where angle is within range of desired angle
    angle_ind <- which(dplyr::near(getValues(angle_rast), alpha[i], 0.25))
    # cut off at shortest distance to edge (so all rays are the same length)
    good_cells <- which(getValues(dist_rast) <= min_dist)
    angle_ind <- angle_ind[angle_ind %in% good_cells]
    angle_sum <- sum(amp_img[angle_ind], na.rm = TRUE)
    Aalpha[i] <- angle_sum
  }

  # find direction where sum of values is highest
  std <- min(alpha[which(unlist(Aalpha) == max(unlist(Aalpha), na.rm = TRUE))], na.rm = TRUE)
  # how dominant is the largest sum of values
  stdi <- mean(unlist(Aalpha), na.rm = TRUE) / max(unlist(Aalpha), na.rm = TRUE)

  return(c(std, stdi))
}

#' Radial Wavelength Metrics
#'
#' Calculates the dominant radial wavelength, radial wavelength
#' index, and mean half wavelength of the radial Fourier spectrum.
#' See Kedron et al. (2018) for more detailed description.
#'
#' @param x A raster or matrix.
#' @param plot Logical. If \code{TRUE}, returns a plot of the
#'   amplitude spectrum with lines showing the radii at which
#'   Srw, Srwi, and Shw are calculated.
#' @return A vector containing numeric values for the dominant
#'   radial wavelength, radial wavelength index, and mean half
#'   wavelength.
#' @examples
#' library(raster)
#'
#' # import raster image
#' data(normforest)
#'
#' # calculate metrics
#' srwvals <- srw(normforest)
#'
#' # extract each value
#' Srw <- srwvals[1]
#' Srwi <- srwvals[2]
#' Shw <- srwvals[3]
#' @export
srw <- function(x, plot = FALSE) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}
  if(class(plot) != 'logical') {stop('plot must be logical.')}

  # get raster dimensions
  M <- ncol(x)
  N <- nrow(x)

  data_type <- if(class(x) == 'matrix') {'matrix'} else {'RasterLayer'}

  # convert matrix to raster if necessary (equal area)
  if (data_type == 'matrix') {
    x <- raster(x)
    extent(x) <- c(0, ncol(x), 0, nrow(x))
    crs(x) <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  }

  # get matrix of values
  zmat <- matrix(getValues(x), ncol = M, nrow = N, byrow = TRUE)

  # if irregular non-na area, cut to biggest square possible
  if (sum(is.na(zmat)) != 0) {
    origin <- c(mean(sp::coordinates(x)[, 1]), mean(sp::coordinates(x)[, 2]))
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

  # complex spectrum from fast fourier transform
  ft <- fft(zmat)
  ft_shift <- geodiv::fftshift(ft)

  # amplitude spectrum
  amplitude <- sqrt((Re(ft_shift) ^ 2) + (Im(ft_shift) ^ 2))

  # take amplitude image, cut in half (y direction)
  amp_img <- setValues(x, amplitude)
  half_dist <- (ymax(amp_img) - ymin(amp_img)) / 2
  ymin <- ymax(amp_img) - half_dist
  amp_img <- crop(amp_img, c(xmin(amp_img), xmax(amp_img), ymin, ymax(amp_img)))

  # get origin of image (actually bottom center)
  origin <- c(mean(sp::coordinates(amp_img)[,1]), ymin(amp_img))

  # figure out number of circles
  if ((0.5 * ncol(amp_img)) <= 100) {
    ncircles <- floor(0.5 * ncol(amp_img))
  } else {
    ncircles <- 100
  }

  if (plot == TRUE) {
    if (data_type == 'matrix') {
      print("cannot draw lines for object of class 'matrix.'")
    }
    nv <- 100
    angle.inc <- 2 * pi / nv
    angles <- seq(0, 2 * pi - angle.inc, by = angle.inc)
    radius <- seq(0, half_dist, length.out = ncircles)
    linex <- unlist(lapply(seq(1, length(radius)), function(x) origin[1] + radius[x] * cos(angles)))
    liney <- unlist(lapply(seq(1, length(radius)), function(x) origin[2] + radius[x] * sin(angles)))
    # number points per line
    npts <- length(linex) / ncircles
    linelist <- lapply(round(seq(1, length(linex) - npts, length.out = ncircles)),
                       FUN = function(i) sp::Lines(list(sp::Line(cbind(linex[i:(i + (npts - 1))], liney[i:(i + (npts - 1))]))), ID = paste('p', i, sep = '')))
    lines <- sp::SpatialLines(linelist, proj4string = sp::CRS(sp::proj4string(amp_img)))

    # plot and get amplitude sums within each radius
    plot(amp_img)
    plot(lines, add = TRUE)
  }

  # calculate distances from center to all other points
  dist_rast <- amp_img
  values(dist_rast) <- NA
  center <- ceiling(dim(x) / 2)
  nce <- ifelse(ncol(amp_img) / 2 == round(ncol(amp_img) / 2), 1, 0)
  dist_rast[nrow(amp_img), center[2] + nce] <- 1
  dist_rast <- distance(dist_rast)

  # get maximum distance from center in x direction and radii to analyze
  half_dist <- max(dist_rast[nrow(dist_rast), ])
  radius <- seq(0, half_dist, length.out = ncircles)

  # get all points some distance away from center, then
  # sum all values at that distance and add to Br list
  Br <- list()
  Br[1] <- 0
  tolerance <- (radius[3] - radius[2]) / 2
  for (i in 2:length(radius)) {
    radius_inds <- which(dplyr::near(getValues(dist_rast), radius[i], tol = tolerance))
    Br[i] <- sum(amp_img[radius_inds], na.rm = TRUE)
  }

  if (plot == TRUE) {
    plot(unlist(Br) ~ (1 / radius), type = 'l')
  }

  Srw <- radius[which(unlist(Br) == max(unlist(Br), na.rm = TRUE))]
  Srwi <- mean(unlist(Br), na.rm = TRUE) / max(unlist(Br), na.rm = TRUE)

  half_val <- 0.5 * sum(unlist(Br), na.rm = TRUE)
  vals <- as.numeric()
  for (i in 1:length(unlist(Br))) {
    vals[i] <- sum(unlist(Br)[1:i], na.rm = TRUE)
  }
  shw_ind <- vals - half_val
  Shw <- radius[which(abs(shw_ind) == min(abs(shw_ind), na.rm = TRUE))]

  return(c(Srw, Srwi, Shw))
}
