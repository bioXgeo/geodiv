#' Texture Direction Metrics
#'
#' Calculates the angle of dominating texture and the texture
#' direction index of the Fourier spectrum image calculated
#' from a raster image (see Kedron et al. 2018).
#'
#' @param x A raster or matrix.
#' @param create_plot Logical. If \code{TRUE}, returns a plot of the
#'   amplitude spectrum with lines showing directions in which
#'   amplitude is summed for the Std and Stdi calculations.
#'   Plotting is not possible when input is a matrix.
#' @param option Numeric. Code for which output metric(s) to return. 1 = Std, 2 = Stdi.
#' @return A vector containing numeric values for the angle of
#'   dominating texture and the texture direction index.
#' @examples
#' # import raster image
#' data(normforest)
#' normforest <- terra::unwrap(normforest)
#'
#' # calculate Std and Stdi
#' stdvals <- std(normforest)
#'
#' # extract each value
#' Std <- stdvals[1]
#' Stdi <- stdvals[2]
#' @importFrom terra crds rast crop ext xmin xmax ymin ymax res setValues crs distance values cellFromRowCol
#' @importFrom dplyr %>% filter group_by summarize near bind_rows
#' @importFrom sf st_as_sf st_cast st_geometry_type
#' @importFrom rlang .data
#' @export
std <- function(x, create_plot = FALSE, option = c(1, 2)) {
  stopifnot('x must be a raster or matrix.' = inherits(x, c('RasterLayer', 'matrix', 'SpatRaster')))
  stopifnot('create_plot must be logical.' = inherits(create_plot, 'logical'))

  # get raster dimensions
  M <- ncol(x)
  N <- nrow(x)

  data_type <- if(class(x)[1] == 'matrix') {'matrix'} else if (class(x)[1] == 'RasterLayer') {'RasterLayer'} else {'SpatRaster'}

  # convert matrix to raster if necessary (equal area)
  if (data_type == 'matrix') {
    x <- rast(x)
    terra::crs(x) <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  }

  if (data_type == 'RasterLayer') {
    x <- rast(x)
  }

  # get matrix of values
  zmat <- matrix(x[], ncol = M, nrow = N, byrow = TRUE)

  # if irregular non-na area, cut to biggest square possible
  if (sum(is.na(zmat)) != 0) {
    coords <- crds(x, na.rm = FALSE)
    origin <- c(mean(coords[, 1], na.rm = TRUE), mean(coords[, 2], na.rm = TRUE))

    potentials <- data.frame(xmin = rep(origin[1], floor(N / 2)),
                             xmax = origin[1],
                             ymin = origin[2],
                             ymax = origin[2])

    change1 <- (res(x)[1] * seq(1, floor(N / 2)))
    change2 <- (res(x)[2] * seq(1, floor(N / 2)))
    potentials$xmin <- origin[1] - change1
    potentials$xmax <- origin[1] + change1
    potentials$ymin <- origin[2] - change2
    potentials$ymax <- origin[2] + change2

    potentials$na <- sapply(seq(1, nrow(potentials)), FUN = function(i) {

      xmin <-
      newrast <- terra::crop(x, ext(potentials$xmin[i], potentials$xmax[i], potentials$ymin[i], potentials$ymax[i]))
      return(sum(is.na(newrast[])))

    })

    max_dim <- potentials[max(which(potentials$na <= 0)),]

    if (sum(is.na(max_dim$xmin)) != 0) {

      zmat <- zmat

    } else {#if (sum(is.na(max_dim$xmin)) == 0) {

      x <- terra::crop(x, ext(max_dim$xmin, max_dim$xmax, max_dim$ymin, max_dim$ymax))

      # get raster dimensions
      M <- ncol(x)
      N <- nrow(x)

      # get matrix of values
      zmat <- matrix(x[], ncol = M, nrow = N, byrow = TRUE)

    }
  }

  if (sum(is.na(zmat)) == length(zmat)) {

    return(c(NA, NA))

  } else {#if (sum(is.na(zmat)) != length(zmat)) {

    # get fourier transform
    # complex spectrum from fast fourier transform
    ft <- fft(zmat)
    ft_shift <- geodiv::fftshift(ft)

    # amplitude spectrum
    amplitude <- sqrt((Re(ft_shift) ^ 2) + (Im(ft_shift) ^ 2))

    # create amplitude image
    amp_img <- terra::setValues(x, amplitude)

    # take amplitude image, cut in half (y direction)
    ymin_amp <- ymin(amp_img)
    ymax_amp <- ymax(amp_img)
    xmin_amp <- xmin(amp_img)
    xmax_amp <- xmax(amp_img)
    half_dist <- (ymax_amp - ymin_amp) / 2
    ymin <- ymax_amp - half_dist
    amp_img <- terra::crop(amp_img, c(xmin_amp, xmax_amp, ymin, ymax_amp))

    # get origin of image (actually bottom center)
    coords <- crds(amp_img, na.rm = FALSE)

    origin <- c(mean(coords[, 1], na.rm = TRUE), ymin)

    if(create_plot == TRUE) {

      if (data_type == 'matrix') {

        print("cannot draw lines for object of class 'matrix.'")

      }

      ### line calculations are taken from the plotrix function draw.radial.line
      # calculate rays extending from origin
      M <- 180
      j <- seq(0, (M - 1))
      alpha <- (pi * j) / M # angles
      px <- c(0, half_dist) # line length
      a_len <- length(alpha)
      linex <- unlist(lapply(seq(1, a_len), function(x) origin[1] + px * cos(alpha[x])))
      liney <- unlist(lapply(seq(1, a_len), function(x) origin[2] + px * sin(alpha[x])))
      linelist <- lapply(seq(1, length(linex), 2), FUN = function(i) {
        data.frame(x = linex[i:(i + 1)], y = liney[i:(i + 1)],
                   id = paste('l', i, sep = ''))})
      linelist <- bind_rows(linelist)
      lines <- linelist %>%
        st_as_sf(coords = c("x", "y"), na.fail = FALSE, crs = terra::crs(amp_img)) %>%
        group_by(.data$id) %>%
        summarize()
      multi_inds <- which(st_geometry_type(lines$geometry) == 'MULTIPOINT')
      lines <- st_cast(lines[multi_inds, ], "LINESTRING")

      # plot and calculate amplitude sums along rays
      terra::plot(amp_img)
      terra::lines(lines)

    }

    # replace: find distances and angle to each pixel, sum all values at given angle

    # calculate distances from center to all other points
    center <- ceiling(dim(x) / 2)
    ncol_amp2 <- ncol(amp_img) / 2
    nce <- ifelse(ncol_amp2 == round(ncol_amp2), 1, 0)
    dist_rast <- amp_img
    terra::values(dist_rast) <- NA
    dist_rast[nrow(amp_img), center[2] + nce] <- 1
    dist_rast <- distance(dist_rast)
    min_dist <- min(max(dist_rast[nrow(dist_rast),]), max(dist_rast[center[2] + nce]))

    # calculate angles from center to all points
    angle_rast <- amp_img
    terra::values(angle_rast) <- NA
    nce <- ifelse(ncol_amp2 == round(ncol_amp2), 1, 0)
    center_ind <- cellFromRowCol(angle_rast, nrow(amp_img), center[2])
    angles <- atan2(coords[, 2] - coords[center_ind, 2], coords[, 1] - coords[center_ind, 1])
    angles <- pracma::rad2deg(angles)
    angle_rast <- terra::setValues(angle_rast, angles)

    # angles at which you want to sum values
    alpha <- seq(0, 180, 0.5)

    # sum values along lines
    vals <- angle_rast[]
    vals_dist <- dist_rast[]
    # cut off at shortest distance to edge (so all rays are the same length)
    good_cells <- which(vals_dist <= min_dist)
    Aalpha <- sapply(alpha, FUN = function(i) {
      angle_ind <- which(dplyr::near(vals[good_cells], i, 0.25))
      angle_sum <- sum(amp_img[good_cells[angle_ind]], na.rm = TRUE)
      return(angle_sum)
    })

    # find direction where sum of values is highest
    std <- min(alpha[which(Aalpha == max(Aalpha, na.rm = TRUE))], na.rm = TRUE)
    # how dominant is the largest sum of values
    stdi <- mean(Aalpha, na.rm = TRUE) / max(Aalpha, na.rm = TRUE)

    out <- c(std, stdi)
    return(out[option])
  }
}

#' Radial Wavelength Metrics
#'
#' Calculates the dominant radial wavelength, radial wavelength
#' index, and mean half wavelength of the radial Fourier spectrum.
#' See Kedron et al. (2018) for more detailed description.
#'
#' @param x A raster or matrix.
#' @param create_plot Logical. If \code{TRUE}, returns a plot of the
#'   amplitude spectrum with lines showing the radii at which
#'   Srw, Srwi, and Shw are calculated.
#' @param option Numeric. Code for which output metric(s) to return. 1 = Srw, 2 = Srwi, 3 = Shw.
#' @return A vector containing numeric values for the dominant
#'   radial wavelength, radial wavelength index, and mean half
#'   wavelength.
#' @examples
#' # import raster image
#' data(normforest)
#' normforest <- terra::unwrap(normforest)
#'
#' # calculate metrics
#' srwvals <- srw(normforest)
#'
#' # extract each value
#' Srw <- srwvals[1]
#' Srwi <- srwvals[2]
#' Shw <- srwvals[3]
#' @importFrom dplyr %>% filter group_by summarize near bind_rows
#' @importFrom sf st_as_sf st_cast st_geometry_type
#' @importFrom terra crds rast crop ext xmin xmax ymin ymax res setValues values crs distance
#' @importFrom rlang .data
#' @export
srw <- function(x, create_plot = FALSE, option = c(1, 2, 3)) {
  stopifnot('x must be a raster or matrix.' = inherits(x, c('RasterLayer', 'matrix', 'SpatRaster')))
  stopifnot('create_plot must be logical.' = inherits(create_plot, 'logical'))

  # get raster dimensions
  M <- ncol(x)
  N <- nrow(x)

  data_type <- if(class(x)[1] == 'matrix') {'matrix'} else if (class(x)[1] == 'RasterLayer') {'RasterLayer'} else {'SpatRaster'}

  # convert matrix to raster if necessary (equal area)
  if (data_type == 'matrix') {
    x <- rast(x)
    terra::crs(x) <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  }

  # get matrix of values
  zmat <- matrix(x[], ncol = M, nrow = N, byrow = TRUE)

  # if irregular non-na area, cut to biggest square possible
  if (sum(is.na(zmat)) != 0) {

    coords <- crds(x, na.rm = FALSE)
    origin <- c(mean(coords[, 1], na.rm = TRUE), mean(coords[, 2], na.rm = TRUE))
    potentials <- data.frame(xmin = rep(origin[1], floor(N / 2)),
                             xmax = origin[1],
                             ymin = origin[2],
                             ymax = origin[2])

    change1 <- (res(x)[1] * seq(1, floor(N / 2)))
    change2 <- (res(x)[2] * seq(1, floor(N / 2)))
    potentials$xmin <- origin[1] - change1
    potentials$xmax <- origin[1] + change1
    potentials$ymin <- origin[2] - change2
    potentials$ymax <- origin[2] + change2

    potentials$na <- sapply(seq(1, nrow(potentials)), FUN = function(i) {
      xmin <-
        newrast <- terra::crop(x, ext(potentials$xmin[i], potentials$xmax[i], potentials$ymin[i], potentials$ymax[i]))
      return(sum(is.na(newrast[])))
    })

    max_dim <- potentials[max(which(potentials$na <= 0)),]

    if (sum(is.na(max_dim$xmin)) != 0) {

      zmat <- zmat

    } else {#if (sum(is.na(max_dim$xmin)) == 0) {

      x <- terra::crop(x, ext(max_dim$xmin, max_dim$xmax, max_dim$ymin, max_dim$ymax))

      # get raster dimensions
      M <- ncol(x)
      N <- nrow(x)

      # get matrix of values
      zmat <- matrix(x[], ncol = M, nrow = N, byrow = TRUE)
    }
  }

  if (sum(is.na(zmat)) == length(zmat)) {

    return(c(NA, NA, NA))

  } else {#if (sum(is.na(zmat)) != length(zmat)) {

    # complex spectrum from fast fourier transform
    ft <- fft(zmat)
    ft_shift <- geodiv::fftshift(ft)

    # amplitude spectrum
    amplitude <- sqrt((Re(ft_shift) ^ 2) + (Im(ft_shift) ^ 2))

    # take amplitude image, cut in half (y direction)
    amp_img <- terra::setValues(x, amplitude)
    ymax_amp <- ymax(amp_img)
    ymin_amp <- ymin(amp_img)
    half_dist <- (ymax_amp - ymin_amp) / 2
    ymin <- ymax(amp_img) - half_dist
    amp_img <- terra::crop(amp_img, c(xmin(amp_img), xmax(amp_img), ymin, ymax_amp))

    # get origin of image (actually bottom center)
    coords_amp <- terra::crds(amp_img, na.rm = FALSE)
    origin <- c(mean(coords_amp[, 1], na.rm = TRUE), ymin_amp)

    # figure out number of circles
    if ((0.5 * ncol(amp_img)) <= 100) {

      ncircles <- floor(0.5 * ncol(amp_img))

    } else {

      ncircles <- 100

    }

    if (create_plot == TRUE) {

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
      linelist <- lapply(round(seq(1, length(linex) - npts, length.out = ncircles)), FUN = function(i) {
        data.frame(x = linex[i:(i + (npts - 1))], y = liney[i:(i + (npts - 1))],
                   id = paste('p', i, sep = ''))})
      linelist <- bind_rows(linelist)
      lines <- linelist %>%
        st_as_sf(coords = c("x", "y"), na.fail = FALSE, crs = terra::crs(amp_img)) %>%
        group_by(.data$id) %>%
        summarize()
      multi_inds <- which(st_geometry_type(lines$geometry) == 'MULTIPOINT')
      lines <- st_cast(lines[multi_inds, ], "LINESTRING")

      # plot and get amplitude sums within each radius
      terra::plot(amp_img)
      terra::plot(lines, add = TRUE)

    }

    # calculate distances from center to all other points
    dist_rast <- amp_img
    terra::values(dist_rast) <- NA
    center <- ceiling(dim(x) / 2)
    nce <- ifelse(ncol(amp_img) / 2 == round(ncol(amp_img) / 2), 1, 0)
    dist_rast[nrow(amp_img), center[2] + nce] <- 1
    dist_rast <- distance(dist_rast)

    # get maximum distance from center in x direction and radii to analyze
    half_dist <- max(dist_rast[nrow(dist_rast), ])
    radius <- seq(0, half_dist, length.out = ncircles)

    # get all points some distance away from center, then
    # sum all values at that distance and add to Br list
    tolerance <- (radius[3] - radius[2]) / 2
    dist_vals <- dist_rast[]
    Br <- sapply(radius, FUN = function(i) {
      radius_inds <- which(dplyr::near(dist_vals, i, tol = tolerance))
      return(sum(amp_img[radius_inds], na.rm = TRUE))
    })
    Br[1] <- 0

    if (create_plot == TRUE) {

      plot(Br ~ (1 / radius), type = 'l')

    }

    Srw <- radius[which(Br == max(Br, na.rm = TRUE))]
    Srwi <- mean(Br, na.rm = TRUE) / max(Br, na.rm = TRUE)

    half_val <- 0.5 * sum(Br, na.rm = TRUE)
    vals <- as.numeric()
    for (i in 1:length(Br)) {
      vals[i] <- sum(Br[1:i], na.rm = TRUE)
    }
    shw_ind <- vals - half_val
    Shw <- radius[which(abs(shw_ind) == min(abs(shw_ind), na.rm = TRUE))]

    # break for tiny matrices that this can't be calculated for
    if (length(Srw) > 1) {
      Srw <- NA
    }
    if (length(Srwi) > 1) {
      Srwi <- NA
    }
    if (length(Shw) > 1) {
      Shw <- NA
    }

    out <- c(Srw, Srwi, Shw)
    return(out[option])
  }
}
