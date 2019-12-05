# figure out peaks and valleys of surface

# from Image Metrology:
# local minimum = points where all 8 surrounding points are higher & below zero
# local maximum = points where all 8 surrounding points are lower & above zero

#' Find Local Peaks
#'
#' Locates local peaks on a raster or matrix. A peak is defined as any pixel where
#' all 8 surrounding pixels have lower values, and the center pixel
#' has a positive value.
#'
#' @param x A raster or matrix.
#' @return A dataframe of local peak locations (\code{x, y}) and
#'   values (\code{val}). The raster or matrix location index (\code{ind}),
#'   row (\code{row}), and column (\code{col}) are also listed.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # locate peaks
#' peaks <- findpeaks(normforest)
#'
#' # calculate the summit density (# peaks/area)
#' N <- ncol(normforest)
#' M <- nrow(normforest)
#' Sds <- nrow(peaks) / ((N - 1) * (M - 1))
#' @export
findpeaks <- function(x) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}

  N <- dim(x)[1] # rows
  M <- dim(x)[2] # cols

  peaks <- data.frame(x = NA, y = NA, val = NA, ind = NA, row = NA, col = NA)

  # convert matrix to raster if necessary (equal area)
  if (class(x) == 'matrix') {
    x <- raster(x)
    extent(x) <- c(0, ncol(x), 0, nrow(x))
    crs(x) <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  }

  # center values, indices, and coordinates
  centers <- getValues(x)
  xcoords <- sp::coordinates(x)[, 1]
  ycoords <- sp::coordinates(x)[, 2]

  # create matrix of centers to get surrounding from
  zmat <- matrix(centers, nrow = N, ncol = M, byrow = TRUE)

  # row/col of each center
  rows <- rep(1:N, each = M)
  cols <- rep(rep(1:M), N)

  # get rid of edge points
  rm_inds <- which(rows < 2 | rows == max(rows) | cols < 2 | cols == max(cols))
  centers <- centers[-rm_inds]
  xcoords <- xcoords[-rm_inds]
  ycoords <- ycoords[-rm_inds]
  rows <- rows[-rm_inds]
  cols <- cols[-rm_inds]

  # gather surrounding points
  xmin <- rows - 1
  xmax <- rows + 1
  ymin <- cols - 1
  ymax <- cols + 1
  ind <- seq(1, length(centers))
  surrounding <- lapply(ind, function(i) {zmat[xmin[i]:xmax[i], ymin[i]:ymax[i]][-5]})

  # check for peak requirements
  check <- lapply(ind, function(i) {((sum(centers[i] > surrounding[[i]]) ==
                                        length(surrounding[[i]])) & centers[i] > 0)})

  # create dataframe, limit to actual peaks
  peaks <- data.frame(x = xcoords, y = ycoords, val = centers, ind = ind,
                      row = rows, col = cols, check = unlist(check))
  peaks <- peaks[peaks$check ==  TRUE,]
  peaks <- peaks[, 1:6]
  peaks <- peaks[!is.na(peaks$val),]

  return(peaks)
}

#' Find Local Valleys
#'
#' Locates local valleys on a raster or matrix. A valley is defined as any pixel where
#' all 8 surrounding pixels have higher values, and the center pixel
#' has a negative value.
#'
#' @param x A raster or matrix.
#' @return A dataframe of local valley locations (\code{x, y}) and
#'   values (\code{val}). The raster or matrix location index (\code{ind}),
#'   row (\code{row}), and column (\code{col}) are also listed.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # locate peaks and valleys
#' peaks <- findpeaks(normforest)
#' valleys <- findvalleys(normforest)
#'
#' # find top 5 peaks, valleys
#' top_peaks <- peaks[order(peaks$val, decreasing = TRUE)[1:5],]
#' bottom_valleys <- valleys[order(valleys$val)[1:5],]
#'
#' # calculate the ten-point height
#' S10z <- (sum(top_peaks$val) + sum(abs(bottom_valleys$val))) / 5
#' @export
findvalleys <- function(x) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}

  N <- dim(x)[1] # rows
  M <- dim(x)[2] # cols

  peaks <- data.frame(x = NA, y = NA, val = NA, ind = NA, row = NA, col = NA)

  # convert matrix to raster if necessary (equal area)
  if (class(x) == 'matrix') {
    x <- raster(x)
    extent(x) <- c(0, ncol(x), 0, nrow(x))
    crs(x) <- "+proj=aea +la
    t_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  }

  # center values, indices, and coordinates
  centers <- getValues(x)
  ind <- seq(1, length(centers))
  xcoords <- sp::coordinates(x)[, 1]
  ycoords <- sp::coordinates(x)[, 2]

  # create matrix of centers to get surrounding from
  zmat <- matrix(centers, nrow = N, ncol = M)

  # row/col of each center
  rows <- rep(1:N, each = M)
  cols <- rep(rep(1:M), N)

  # get rid of edge points
  rm_inds <- which(rows < 2 | rows == max(rows) | cols < 2 | cols == max(cols))
  centers <- centers[-rm_inds]
  xcoords <- xcoords[-rm_inds]
  ycoords <- ycoords[-rm_inds]
  rows <- rows[-rm_inds]
  cols <- cols[-rm_inds]

  # gather surrounding points
  xmin <- rows - 1
  xmax <- rows + 1
  ymin <- cols - 1
  ymax <- cols + 1
  ind <- seq(1, length(centers))
  surrounding <- lapply(ind, function(i) {zmat[xmin[i]:xmax[i], ymin[i]:ymax[i]][-5]})

  # check for valley requirements
  check <- lapply(ind, function(i) {((sum(centers[i] < surrounding[[i]]) ==
                                        length(surrounding[[i]])) & centers[i] < 0)})

  # create dataframe, limit to actual valleys
  valleys <- data.frame(x = xcoords, y = ycoords, val = centers, ind = ind,
                      row = rows, col = cols, check = unlist(check))
  valleys <- valleys[valleys$check ==  TRUE,]
  valleys <- valleys[, 1:6]
  valleys <- valleys[!is.na(valleys$val),]

  return(valleys)
}

#' Mean Summit Curvature
#'
#' Calculates the mean summit curvature of a raster or matrix. Mean summit
#' curvature is the average principle curvature of local maximas
#' on the surface.
#'
#' @param x A raster or matrix.
#' @return A numeric value representing the average curvature of
#'   surface peaks.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # calculate mean summit curvature
#' Ssc <- ssc(normforest)
#' @export
ssc <- function(x) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}

  # convert matrix to raster if necessary (equal area)
  if (class(x) == 'matrix') {
    x <- raster(x)
    extent(x) <- c(0, ncol(x), 0, nrow(x))
    crs(x) <- "+proj=aea +la
    t_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  }

  # z values, coordinates, and resolution (change in x, y)
  z <- getValues(x)
  xcoords <- sp::coordinates(x)[, 1]
  ycoords <- sp::coordinates(x)[, 2]
  deltax <- res(x)[1] / mean(res(x)[1], res(x)[2])
  deltay <- res(x)[2] / mean(res(x)[1], res(x)[2])

  # zmat
  zmat <- zshift(x, xdist = 0, ydist = 0, scale = TRUE)
  zmat <- matrix(zmat, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)

  # number of peaks
  peaks <- findpeaks(x)
  n <- nrow(peaks)

  # zshift of 1
  z_xpl <- zshift(x, xdist = 1, ydist = 0, scale = TRUE)
  z_xpl <- matrix(z_xpl, nrow = nrow(x), ncol = ncol(x) - 1, byrow = TRUE)

  # yshift of 1
  z_ypl <- zshift(x, xdist = 0, ydist = 1, scale = TRUE)
  z_ypl <- matrix(z_ypl, nrow = nrow(x) - 1, ncol = ncol(x), byrow = TRUE)

  # zshift of 1
  z_xmn <- zshift(x, xdist = -1, ydist = 0, scale = TRUE)
  z_xmn <- matrix(z_xmn, nrow = nrow(x), ncol = ncol(x))

  # zshift of 1
  z_ymn <- zshift(x, xdist = 0, ydist = -1, scale = TRUE)
  z_ymn <- matrix(z_ymn, nrow = nrow(x), ncol = ncol(x))

  if (n != 0) {
    # get z_xpl, z_ypl at peaks (add to df)
    peaks$val_xpl <- unlist(lapply(seq(1, nrow(peaks)),
                                   FUN = function(i) z_xpl[peaks$row[i], peaks$col[i]]))
    peaks$val_ypl <- unlist(lapply(seq(1, nrow(peaks)),
                                   FUN = function(i) z_ypl[peaks$row[i], peaks$col[i]]))
    peaks$val_xmn <- unlist(lapply(seq(1, nrow(peaks)),
                                   FUN = function(i) z_xmn[peaks$row[i], peaks$col[i]]))
    peaks$val_ymn <- unlist(lapply(seq(1, nrow(peaks)),
                                   FUN = function(i) z_ymn[peaks$row[i], peaks$col[i]]))

    # new center val
    peaks$val <- unlist(lapply(seq(1, nrow(peaks)),
                               FUN = function(i) zmat[peaks$row[i], peaks$col[i]]))

    # calculate curvature with normalized values
    ssc <- -(1 / (2 * n)) * sum(((((peaks$val - peaks$val_xpl) + (peaks$val - peaks$val_xmn)) / 2) / deltax),
                                ((((peaks$val - peaks$val_ypl) + (peaks$val - peaks$val_ymn)) / 2) / deltay))
  } else {
    print('No peaks to analyze for curvature.')
    ssc <- NA
  }

  return(ssc)
}

#' Summit Density
#'
#' Calculates the summit density of a raster or matrix. Summit density is the number of local
#' peaks per unit area.
#'
#' @param x A raster or matrix.
#' @return A numeric value representing the summit density.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # calculate summit density.
#' Sds <- sds(normforest)
#' @export
sds <- function(x) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}

  M <- nrow(x)
  N <- ncol(x)

  peaks <- findpeaks(x)

  val <- nrow(peaks) / ((N - 1) * (M - 1))

  return(val)
}

#' Ten-Point Height
#'
#' Calculates the average height above the mean surface for the five highest local maxima
#' plus the average height below the mean surface for the five lowest local minima.
#'
#' @param x A raster or matrix.
#' @return A numeric value representing the ten-point height.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # calculate ten-point height.
#' S10z <- s10z(normforest)
#' @export
s10z <- function(x) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}

  peaks <- findpeaks(x)
  valleys <- findvalleys(x)

  # find top 5 peaks, valleys
  top_peaks <- peaks[order(peaks$val, decreasing = TRUE)[1:5],]
  bottom_valleys <- valleys[order(valleys$val)[1:5],]

  val <- (sum(top_peaks$val, na.rm = TRUE) + sum(abs(bottom_valleys$val), na.rm = TRUE)) / 5

  return(val)
}
