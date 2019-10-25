#### !!!! srw and the like should not be done on very small windows (e.g., size = 3)

#' Calculate Texture Metrics per Pixel
#'
#' Calculates the various texture metrics over windows centered
#' on individual pixels. This creates a continuous surface of the
#' texture metric.
#'
#' @param x A raster or matrix.
#' @param window_type Character. Type of window, either circular or square.
#' @param size Numeric. Size of window, in number of pixels extra on each side,
#' or distance from center (in meters) for circular windows.
#' @param epsg_proj Numeric. Appropriate equal area EPSG code used to
#' crop raster to each circular window.
#' @param metric Character. Metric to calculate for each window. Metrics
#' from the geodiv package are listed below.
#' @param args List. Arguments from function to be applied over each window
#' (e.g., list(threshold = 0.2)).
#' @param parallel Logical. Option to run the calculations in
#' parallel on available cores.
#' @param ncores Numeric. If parallel is TRUE, number of cores on which to
#' run the calculations. Defaults to all available, minus 1.
#' @param nclumps Numeric. Number of clumps to split the raster or matrix into.
#' @return A raster or list of rasters (if function results in multiple outputs)
#' with pixel values representative of the metric value for the window
#' surrounding that pixel.
#' @note The total window size for square windows will be (size * 2) + 1.
#' @details Metrics available from geodiv package:
#' \enumerate{
#'    \item{\code{'sa'}: average surface roughness}
#'    \item{\code{'sq'}: root mean square roughness}
#'    \item{\code{'s10z'}: ten-point height}
#'    \item{\code{'sdq'}: root mean square slope of surface, 2-point method}
#'    \item{\code{'sdq6'}: root mean square slope of surface, 7-point method}
#'    \item{\code{'sdr'}: surface area ratio}
#'    \item{\code{'sbi'}: surface bearing index}
#'    \item{\code{'sci'}: core fluid retention index}
#'    \item{\code{'ssk_adj'}: adjusted skewness}
#'    \item{\code{'ssk'}: skewness}
#'    \item{\code{'sku_exc'}: excess kurtosis}
#'    \item{\code{'sku'}: kurtosis}
#'    \item{\code{'sds'}: summit density}
#'    \item{\code{'sfd'}: 3d fractal dimension}
#'    \item{\code{'srw'}: dominant radial wavelength}
#'    \item{\code{'srwi'}: radial wavelength index}
#'    \item{\code{'shw'}: mean half wavelength}
#'    \item{\code{'std'}: angle of dominating texture}
#'    \item{\code{'stdi'}: texture direction index}
#'    \item{\code{'svi'}: valley fluid retention index}
#'    \item{\code{'str'}: texture aspect ratio}
#'    \item{\code{'ssc'}: mean summit curvature}
#'    \item{\code{'sv'}: maximum valley depth}
#'    \item{\code{'sph'}: maximum peak height}
#'    \item{\code{'sk'}: core roughness depth}
#'    \item{\code{'smean'}: mean peak height}
#'    \item{\code{'svk'}: reduced valley depth}
#'    \item{\code{'spk'}: reduced peak height}
#'    \item{\code{'scl'}: correlation length}
#'    \item{\code{'sdc'}: bearing area curve height interval}
#' }
#' @examples
#' library(raster)
#'
#' # import raster image
#' data(normforest)
#'
#' # crop raster to much smaller area
#' x <- crop(normforest, extent(-123, -122.99, 43, 43.01))
#'
#' # get a surface of root mean square roughness
#' sq_img <- texture_image(x = x, window = 'square',
#' size = 11, metric = 'sq',
#' parallel = TRUE, ncores = 2)
#'
#' # plot the result
#' plot(sq_img)
#' @export
texture_image <- function(x, window_type = 'square', size = 11, epsg_proj = 5070,
                          metric, args = NULL, parallel = TRUE, ncores = NULL, nclumps = 100){

  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}
  if(class(window_type) != 'character') {stop('window_type must be a string.')}
  if(class(size) != 'numeric') {stop('size must be numeric.')}
  if(class(epsg_proj) != 'numeric') {stop('epsg_proj must be numeric.')}
  if(class(metric) != 'character') {stop('metric must be a character.')}

  if(length(metric) > 1) {stop('too many values provided for metric.')}
  if(length(epsg_proj) > 1) {stop('too many values provided to epsg_proj.')}
  if(length(size) > 1) {stop('too many values provided to size.')}
  if(length(window_type) > 1) {stop('too many values provided to window_type.')}

  if(missing(ncores)) {ncores <- parallel::detectCores() - 1}

  if(Sys.info()['sysname'][[1]] == 'Windows') {
    warning('mclapply is not supported on Windows, using parLapply instead (much slower).')
    os_type = 'windows'
  } else {
    os_type = 'other'
  }

  if (class(x) == 'matrix') {
    # convert to equal area raster
    x <- raster(x)
    extent(x) <- c(0, ncol(x), 0, nrow(x))
    crs(x) <- st_crs(epsg_proj)$proj4string
  }

  # get list of # total pixels, break up into smaller lists (by number of cores available)
  pixlist <- seq(1, length(x), 1)
  seg_length <- ceiling(length(x) / nclumps)
  segment <- ceiling(length(x) / seg_length)
  new_pixlist <- vector('list', segment)
  i <- 1
  while(i <= segment) {
    newi <- ceiling(seq(1, length(x), seg_length))[i]
    if (i < segment) {
      new_pixlist[[i]] <- seq(newi, (newi - 1) + seg_length, 1)
    } else {
      new_pixlist[[i]] <- seq(newi, length(x), 1)
    }
    i <- i + 1
  }
  cat('There are ', segment, ' clumps to analyze.', '\n', sep = '')

  # output raster
  out <- x

  # reset size so it's a common variable (whether square or circular)
  if (window_type == 'square') {
    # change size to distance out from center
    size <- size
  } else {
    # convert to new proj
    projx <- projectRaster(x, crs = sp::CRS(sf::st_crs(epsg_proj)$proj4string))

    # get equivalent # pixels of size
    pixeq_size <- ceiling(size / res(projx))[1]
  }

  # add padding to raster
  if (window_type == 'square') {
    ext_x <- pad_edges(x, window_type = 'square', size = size)

    # data frame of x, y locations
    coords <- data.frame(xyFromCell(x, 1:ncell(x)))
    rownum <- rowFromCell(x, pixlist) + size
    colnum <- colFromCell(x, pixlist) + size
  } else {
    ext_x <- pad_edges(x, window_type = 'circular', size = pixeq_size)

    # data frame of x, y locations
    coords <- data.frame(xyFromCell(x, 1:ncell(x)))
    rownum <- rowFromCell(x, pixlist) + pixeq_size
    colnum <- colFromCell(x, pixlist) + pixeq_size
  }

  # collect arguments
  if (!is.null(args)) {
    input_args <- args
  } else {
    input_args <- NULL
  }

  if (parallel == FALSE) {
    result <- c()
    print('Beginning calculation of metrics over windows...')
    start <- Sys.time()
    for (i in pixlist) {
      outval <- window_metric(ext_x, i, window_type = window_type, size = size, epsg_proj = epsg_proj,
                                coords = coords, rownum = rownum, colnum = colnum,
                                metric = metric, args = input_args)

      result <- c(result, outval)
    }
    end <- Sys.time()
    cat('Total time to calculate metrics: ', end - start, '\n', sep = '')
  } else if (os_type != 'windows' & parallel == TRUE) {
    # get smaller rasters
    print('Beginning calculation of metrics over windows...')
    start <- Sys.time()
    result <- parallel::mclapply(new_pixlist, FUN = function(l) {
      lapply(l, FUN = function(i) {window_metric(x = ext_x, i = i, window_type = window_type,
                                                 size = size, epsg_proj = epsg_proj,
                                                 coords = coords, rownum = rownum,
                                                 colnum = colnum, metric = metric,
                                                 args = input_args)})
    }, mc.cores = ncores, mc.cleanup = TRUE)
    end <- Sys.time()
    cat('Total time to calculate metrics: ', end - start, '\n', sep = '')

    result <- do.call(unlist, args = list(result, recursive = TRUE))
  } else {
    print('Beginning calculation of metrics over windows...')
    start <- Sys.time()
    # make and start cluster
    try(stopCluster(cl), silent = TRUE)
    cl <- makeCluster(ncores, type = 'SOCK')
    doSNOW::registerDoSNOW(cl)
    snow::clusterExport(cl = cl, list = list('ext_x', 'coords', 'size',
                                             'window_type', 'epsg_proj',
                                             'rownum', 'colnum',
                                             'new_pixlist', 'metric', 'input_args'),
                        envir = environment())
    # for each list in new_pixlist, run lapply
    result <- snow::parLapply(cl, new_pixlist, fun = function(l) {
      lapply(l, FUN = function(i) {geodiv::window_metric(x = ext_x, i = i, window_type = window_type,
                                                 size = size, epsg_proj = epsg_proj,
                                                 coords = coords, rownum = rownum,
                                                 colnum = colnum, metric = metric,
                                                 args = input_args)})
    })
    stopCluster(cl)
    end <- Sys.time()
    cat('Total time to calculate metrics: ', end - start, '\n', sep = '')

    result <- do.call(unlist, args = list(result, recursive = TRUE))
  }

  # deal with functions that have multiple outputs
  outfinal <- list()
  nresult <- length(result) / length(x)
  for (i in 1:nresult) {
    temp <- result[seq(i, length(result), nresult)]
    outfinal[[i]] <- setValues(out, temp)
  }

  return(outfinal)
}

#' Calculate Texture Metric for Single Pixel
#'
#' Calculates the various texture metrics over a window centered
#' on an individual pixel.
#'
#' @param x A raster or matrix.
#' @param i Index of cell at which to calculate the metric.
#' @param window_type Character. Type of window, either circular or square.
#' @param size Numeric. Size of window, in number of pixels on each
#' side for square windows (must be an odd value), or distance from
#' center (in meters) for circular windows.
#' @param epsg_proj Numeric. Appropriate equal area EPSG code used to
#' crop raster to each circular window. Only used for circular windows.
#' @param coords Matrix of coordinates for the input raster or matrix.
#' x-coordinates should be in the first column, and y-coordinates should
#' be in the second column.
#' @param rownum Vector of row numbers at which to calculate the metric.
#' @param colnum Vector of column numbers at which to calculate the metric.
#' @param metric Character. Metric to calculate for each window. Metrics
#' from the geodiv package are listed below.
#' @param args List. Arguments from function to be applied over each window
#' (e.g., list(threshold = 0.2)).
#' @return A raster with pixel values representative of the metric
#' value for the window surrounding that pixel.
#' @note Note that if calculating the metric at the edge of a raster or matrix,
#' the input raster/matrix must be padded. This can be done using the \code{pad_edges}
#' function.
#' @details Metrics from geodiv package:
#' \enumerate{
#'    \item{\code{'sa'}: average surface roughness}
#'    \item{\code{'sq'}: root mean square roughness}
#'    \item{\code{'s10z'}: ten-point height}
#'    \item{\code{'sdq'}: root mean square slope of surface, 2-point method}
#'    \item{\code{'sdq6'}: root mean square slope of surface, 7-point method}
#'    \item{\code{'sdr'}: surface area ratio}
#'    \item{\code{'sbi'}: surface bearing index}
#'    \item{\code{'sci'}: core fluid retention index}
#'    \item{\code{'ssk_adj'}: adjusted skewness}
#'    \item{\code{'ssk'}: skewness}
#'    \item{\code{'sku_exc'}: excess kurtosis}
#'    \item{\code{'sku'}: kurtosis}
#'    \item{\code{'sds'}: summit density}
#'    \item{\code{'sfd'}: 3d fractal dimension}
#'    \item{\code{'srw'}: dominant radial wavelength}
#'    \item{\code{'srwi'}: radial wavelength index}
#'    \item{\code{'shw'}: mean half wavelength}
#'    \item{\code{'std'}: angle of dominating texture}
#'    \item{\code{'stdi'}: texture direction index}
#'    \item{\code{'svi'}: valley fluid retention index}
#'    \item{\code{'str'}: texture aspect ratio}
#'    \item{\code{'ssc'}: mean summit curvature}
#'    \item{\code{'sv'}: maximum valley depth}
#'    \item{\code{'sph'}: maximum peak height}
#'    \item{\code{'sk'}: core roughness depth}
#'    \item{\code{'smean'}: mean peak height}
#'    \item{\code{'svk'}: reduced valley depth}
#'    \item{\code{'spk'}: reduced peak height}
#'    \item{\code{'scl'}: correlation length}
#'    \item{\code{'sdc'}: bearing area curve height interval}
#' }
#' @examples
#' library(raster)
#'
#' # import raster image
#' data(normforest)
#'
#' # crop raster to much smaller area if on a smaller computer
#' x <- crop(normforest, extent(-123, -122.99, 43, 43.01))
#'
#' # get coordinates, rownums, cellnums
#' ext_x <- pad_edges(x, window_type = 'square', size = 4)
#' coords <- data.frame(xyFromCell(x, 1:ncell(x)))
#' rownum <- rowFromCell(x, pixlist) + 4
#' colnum <- colFromCell(x, pixlist) + 4
#'
#' # get a surface of root mean square roughness
#' sq_img <- window_metric(x = x, i = 40, window = 'square',
#' size = 4, epsg_proj = 5070, coords = coords,
#' rownum = rownum, colnum = colnum, metric = 'sq')
#' @export
window_metric <- function(x, i, window_type = 'square', size = 11, epsg_proj = 5070,
                          coords, rownum, colnum, metric, args = NULL) {

  # row and column number
  rownum <- rownum[i]
  colnum <- colnum[i]
  coords <- coords[i, ]

  if (class(x) == 'matrix') {
    # convert to equal area raster
    x <- raster(x)
    extent(x) <- c(0, ncol(x), 0, nrow(x))
    crs(x) <- st_crs(epsg_proj)$proj4string
  }

  if (window_type == 'square') {
    # crop to square
    y1 <- rownum - size
    y2 <- rownum + size
    x1 <- colnum - size
    x2 <- colnum + size
    cropped_x <- raster::crop(x, extent(x, y1, y2, x1, x2))
  } else {
    # crop to circle
    pt_sf <- st_as_sf(coords, coords = c("x", "y"), crs = st_crs(x))
    pt_sf <- st_transform(pt_sf, epsg_proj)
    poly_circ <- st_buffer(pt_sf, size)
    poly_circ <- st_transform(poly_circ, st_crs(x))
    poly_circ <- as_Spatial(poly_circ)
    cropped_x <- raster::crop(x, poly_circ)
    cropped_x <- mask(cropped_x, poly_circ)
  }
  # append cropped_x to args
  args$x <- cropped_x

  # calculate metric
  out_val <- do.call(metric, args)

  return(out_val)
}

#' Extend edges of a raster or matrix.
#'
#' Extends edge values of a raster or matrix by a specified number of pixels.
#'
#' @param x A raster or matrix.
#' @param window_type Character. Type of window, either circular or square.
#' @param size Numeric. Size of window, in number of pixels on each
#' side. For circular windows, this will be the number of pixels of the radius.
#' @param epsg_proj Numeric. Appropriate equal area EPSG code used to
#' crop raster to each circular window. Only used for circular windows.
#' @return A raster with edges padded \code{size} number of pixels on each edge.
#' @examples
#' library(raster)
#'
#' # import raster image
#' data(normforest)
#'
#' # crop raster to much smaller area
#' x <- pad_edges(normforest, 'square', 11)
#' @export
pad_edges <- function(x, window_type = 'square', size = 11, epsg_proj = 5070) {
  if (class(x) == 'matrix') {
    # convert to equal area raster
    x <- raster(x)
    extent(x) <- c(0, ncol(x), 0, nrow(x))
    crs(x) <- st_crs(epsg_proj)$proj4string
  }

  # add padding to raster or matrix
  if (window_type == 'square') {
    # continue values to edges to account for edge effect (# pixels radius/edge)

    # first, get edge values that will be extended
    firstrow_vals <- x[1, ]
    firstcol_vals <- x[, 1]
    lastrow_vals <- x[nrow(x), ]
    lastcol_vals <- x[, ncol(x)]

    # add pixels on all sides, increasing the extent of the raster as well
    ext_x <- x
    dim(ext_x) <- c(nrow(x) + (2 * size), ncol(x) + (2 * size))
    extra_x <- size * res(x)[2]
    extra_y <- size * res(x)[1]
    extent(ext_x) <- extent(c(xmin(ext_x) - extra_x, xmax(ext_x) + extra_x,
                              ymin(ext_x) - extra_y, ymax(ext_x) + extra_y))

    # fill in top rows
    ext_x[1:size, (size + 1):(ncol(ext_x) - size)] <- rep(firstrow_vals, size)
    # fill in bottom rows
    ext_x[(nrow(ext_x) - (size - 1)):nrow(ext_x), (size + 1):(ncol(ext_x) - size)] <- rep(lastrow_vals, size)
    # fill in left columns
    ext_x[(size + 1):(nrow(ext_x) - size), 1:size] <- t(firstcol_vals * matrix(1, nrow = nrow(x), ncol = size))
    # fill in right columns
    ext_x[(size + 1):(nrow(ext_x) - size), (ncol(ext_x) - (size - 1)):ncol(ext_x)] <- t(lastcol_vals * matrix(1, nrow = nrow(x), ncol = size))
    # fill in middle
    ext_x[(size + 1):(nrow(ext_x) - size), (size + 1):(ncol(ext_x) - size)] <- getValues(x)
    # fill in corners with nearest point value (always the same)
    ext_x_mat <- zoo::na.approx(matrix(ext_x, nrow = nrow(ext_x), ncol = ncol(ext_x)), rule = 2)
    ext_x <- setValues(ext_x, t(ext_x_mat))
  } else {
    # convert to new proj
    projx <- projectRaster(x, crs = sp::CRS(sf::st_crs(epsg_proj)$proj4string))

    # extend...
    # continue values to edges to account for edge effect (# pixels radius/edge)
    firstrow_vals <- x[1, ]
    firstcol_vals <- x[, 1]
    lastrow_vals <- x[nrow(x), ]
    lastcol_vals <- x[, ncol(x)]

    # add pixels on all sides, increasing the extent of the raster as well
    ext_x <- x
    dim(ext_x) <- c(nrow(x) + (2 * size), ncol(x) + (2 * size))
    extra_x <- size * res(x)[2]
    extra_y <- size * res(x)[1]
    extent(ext_x) <- extent(c(xmin(ext_x) - extra_x, xmax(ext_x) + extra_x,
                              ymin(ext_x) - extra_y, ymax(ext_x) + extra_y))

    # fill in top rows
    ext_x[1:size, (size + 1):(ncol(ext_x) - size)] <- rep(firstrow_vals, size)
    # fill in bottom rows
    ext_x[(nrow(ext_x) - (size - 1)):nrow(ext_x), (size + 1):(ncol(ext_x) - size)] <- rep(lastrow_vals, size)
    # fill in left columns
    ext_x[(size + 1):(nrow(ext_x) - size), 1:size] <- t(firstcol_vals * matrix(1, nrow = nrow(x), ncol = size))
    # fill in right columns
    ext_x[(size + 1):(nrow(ext_x) - size), (ncol(ext_x) - (size - 1)):ncol(ext_x)] <- t(lastcol_vals * matrix(1, nrow = nrow(x), ncol = size))
    # fill in middle
    ext_x[(size + 1):(nrow(ext_x) - size), (size + 1):(ncol(ext_x) - size)] <- getValues(x)
    # fill in corners with nearest point value (always the same)
    ext_x_mat <- zoo::na.approx(matrix(ext_x, nrow = nrow(ext_x), ncol = ncol(ext_x)), rule = 2)
    ext_x <- setValues(ext_x, t(ext_x_mat))
  }
  if (class(x) == 'matrix') {
    return(as.matrix(ext_x))
  } else {
    return(ext_x)
  }
}
