#' Calculate Texture Metrics per Pixel
#'
#' Calculates the various texture metrics over windows centered
#' on individual pixels. This creates a continuous surface of the
#' texture metric.
#'
#' @param x A raster or matrix. If a raster is given, it will be projected to
#' an equal area projection (given by \code{epsg_proj} argument).
#' @param window_type Character. Type of window, either circular or square.
#' @param size Numeric. Size of window, in number of pixels extra on each side,
#' or radius (in meters).
#' @param in_meters Logical. Is the size given in meters?
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
#' surrounding that pixel. Note that the raster will always be projected to an
#' equal area projection because calculations are done on matrices with a
#' radius of number of pixels.
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
#'    \item{\code{'stxr'}: texture aspect ratio}
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
#' sa_img <- texture_image(x = x, window = 'square',
#' size = 11, metric = 'sa',
#' parallel = FALSE)
#'
#' # plot the result
#' plot(sa_img)
#' @export
texture_image <- function(x, window_type = 'square', size = 5, in_meters = FALSE, epsg_proj = 5070,
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
    print('mclapply is not supported on Windows, using parLapply instead.')
    os_type = 'windows'
  } else {
    os_type = 'other'
  }


  if ('matrix' %in% base::class(x)) {
    # convert to equal area raster
    x <- raster(x)
    extent(x) <- c(0, ncol(x), 0, nrow(x))
    #crs(x) <- paste0('+init=EPSG:', epsg_proj)
    crs(x) <- st_crs(5070)$proj4string
  } else {
    # tell users that this will always reproject to equal area
    print('Warning: This function will reproject rasters to equal area! Make sure that the epsg_code argument contains an appropriate projection for your data.')
  }

  # reset size so it's a common variable (always number pixels)
  if (in_meters == TRUE) {
    # convert to new proj
    tempproj <- raster::projectRaster(x, crs = sp::CRS(sf::st_crs(epsg_proj)$proj4string))


    # get equivalent # pixels of size
    size <- ceiling(size / res(tempproj))[1]
  }

  # add padding to raster
  ext_x <- pad_edges(x, size = (size * 2), val = NULL)
  dummy_ext_x <- pad_edges(x, size = (size * 2), val = -88888)

  # make sure that internal na values are substituted until later
  dummy_ext_x[is.na(dummy_ext_x)] <- Inf

  # convert to new proj
  if (st_crs(x)$proj4string != st_crs(epsg_proj)$proj4string) {
    projx <- raster::projectRaster(ext_x, crs = sp::CRS(sf::st_crs(epsg_proj)$proj4string))
    noext_projx <- raster::projectRaster(from = dummy_ext_x, to = projx)
  } else {
    projx <- ext_x
    noext_projx <- dummy_ext_x
  }
  # get projected unpadded raster coordinates
  noext_coords <- data.frame(xyFromCell(noext_projx, 1:ncell(noext_projx)))

  # remove NA value locations from noext_coords
  na_ind <- which((is.na(getValues(noext_projx)) | getValues(noext_projx) < -9999))
  noext_coords <- noext_coords[-na_ind,]

  # data frame of x, y locations
  coords <- data.frame(xyFromCell(projx, 1:ncell(projx)))
  coords$ind <- seq(1, nrow(coords))
  good_coords <- dplyr::inner_join(coords, noext_coords, by = c('x', 'y'))$ind

  # get list of # total pixels, break up into smaller lists (by number of cores available)
  pixlist <- good_coords
  seg_length <- ceiling(length(pixlist) / nclumps)
  segment <- ceiling(length(pixlist) / seg_length)
  new_pixlist <- vector('list', segment)
  i <- 1
  while(i <= segment) {
    newi <- ceiling(seq(pixlist[1], max(pixlist), seg_length))[i]
    if (i < segment) {
      new_pixlist[[i]] <- seq(newi, (newi - 1) + seg_length, 1)
      new_pixlist[[i]] <- new_pixlist[[i]][new_pixlist[[i]] %in% pixlist]
    } else {
      new_pixlist[[i]] <- seq(newi, max(pixlist), 1)
      new_pixlist[[i]] <- new_pixlist[[i]][new_pixlist[[i]] %in% pixlist]
    }
    i <- i + 1
  }

  if  (parallel == TRUE) {
    cat('There are ', segment, ' clumps to analyze.', '\n', sep = '')
  }

  rownum <- rowFromCell(projx, seq(1, length(projx)))
  colnum <- colFromCell(projx, seq(1, length(projx)))

  # collect arguments
  if (!is.null(args)) {
    input_args <- args
  } else {
    input_args <- NULL
  }

  # matrices are faster for window_metric, so convert to matrix
  ext_mat <- as.matrix(projx, nrow = nrow(projx), ncol = ncol(projx), byrow = TRUE)

  if (parallel == FALSE) {
    print('Beginning calculation of metrics over windows...')
    start <- Sys.time()
    result <- lapply(pixlist, FUN = function(i) {window_metric(x = ext_mat, i = i,
                                                               window_type = window_type,
                                                               size = size,
                                                               rownum = rownum,
                                                               colnum = colnum,
                                                               metric = metric,
                                                               args = input_args)})
    end <- Sys.time()
    cat('Total time to calculate metrics: ', end - start, '\n', sep = '')
    result <- unlist(result)
  } else if (os_type != 'windows' & parallel == TRUE) {
    print('Beginning calculation of metrics over windows...')
    start <- Sys.time()
    result <- parallel::mclapply(new_pixlist, FUN = function(l) {
      lapply(l, FUN = function(i) {window_metric(x = ext_mat, i = i, window_type = window_type,
                                                 size = size,
                                                 rownum = rownum,
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
    # doSNOW::registerDoSNOW(cl)
    parallel::clusterExport(cl = cl, list = list('ext_mat', 'coords', 'size',
                                             'window_type',
                                             'rownum', 'colnum',
                                             'new_pixlist', 'metric', 'input_args'),
                        envir = environment())
    # for each list in new_pixlist, run lapply
    result <- parallel::parLapply(cl, new_pixlist, fun = function(l) {
      lapply(l, FUN = function(i) {geodiv::window_metric(x = ext_mat, i = i, window_type = window_type,
                                                 size = size,
                                                 rownum = rownum,
                                                 colnum = colnum, metric = metric,
                                                 args = input_args)})
    })
    stopCluster(cl)
    end <- Sys.time()
    cat('Total time to calculate metrics: ', end - start, '\n', sep = '')

    result <- do.call(unlist, args = list(result, recursive = TRUE))
  }

  # create clean out raster
  out <- noext_projx
  out[na_ind] <- NA

  # deal with functions that have multiple outputs
  outfinal <- list()
  nresult <- length(result) / length(good_coords)
  for (i in 1:nresult) {
    temp <- data.frame(newvals = result[seq(i, length(result), nresult)], ind = pixlist)
    temprast <- out
    temprast[temp$ind] <- temp$newvals
    outfinal[[i]] <- temprast
    outfinal[[i]] <- trim(outfinal[[i]])
  }

  if (length(outfinal) == 1) {
    outfinal <- outfinal[[1]]
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
#' @param size Numeric. Radius of window in number of pixels.
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
#'    \item{\code{'stxr'}: texture aspect ratio}
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
#' pixlist <- seq(1, length(x), 1)
#' ext_x <- pad_edges(x, size = 4)
#' rownum <- rowFromCell(x, pixlist) + 4
#' colnum <- colFromCell(x, pixlist) + 4
#'
#' # get a surface of root mean square roughness
#' sq_val <- window_metric(x = x, i = 40, window = 'square',
#' size = 4, rownum = rownum, colnum = colnum, metric = 'sq')
#' @export
window_metric <- function(x, i, window_type = 'square', size = 11,
                          rownum, colnum, metric, args = NULL) {

  # row and column number
  rownum <- rownum[i]
  colnum <- colnum[i]

  if ('RasterLayer' %in% base::class(x)) {
    # tell users that this will always reproject to equal area
    print('Warning: This function assumes an equal area raster!')

    # convert to matrix
    x <- raster::as.matrix(x, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
  }

  # crop to square
  x1 <- rownum - size
  x2 <- rownum + size
  y1 <- colnum - size
  y2 <- colnum + size
  cropped_x <- matrix(x[x1:x2, y1:y2], nrow = (size * 2) + 1, ncol = (size * 2) + 1)

  # crop to circle if necessary
  if (window_type == 'circle') {
    center <- c(round(ncol(cropped_x) / 2), round(nrow(cropped_x) / 2))
    row_ind <- matrix(rep(seq(1, nrow(cropped_x)), each = ncol(cropped_x)), nrow = nrow(cropped_x), ncol = ncol(cropped_x),
                      byrow = TRUE)
    col_ind <- matrix(rep(seq(1, ncol(cropped_x)), nrow(cropped_x)), nrow = nrow(cropped_x), ncol = ncol(cropped_x),
                      byrow = TRUE)
    mat_dists_x <- col_ind - col_ind[center[2], center[1]]
    mat_dists_y <- row_ind - row_ind[center[2], center[1]]
    mat_dists <- sqrt((mat_dists_x ^ 2) + (mat_dists_y ^ 2))
    cropped_x[mat_dists > size] <- NA
  }

  args[[length(args) + 1]] <- cropped_x

  # make sure that the matrix is assigned the name of the first function argument
  names(args)[[length(args)]] <- names(as.list(args(metric)))[1]

  # calculate metric
  out_val <- do.call(metric, args)

  return(out_val)
}

#' Extend edges of a raster or matrix.
#'
#' Extends edge values of a raster or matrix by a specified number of pixels.
#'
#' @param x A raster or matrix.
#' @param size Numeric. Number of pixels to add to each side.
#' @param val Numeric. If NULL (default), this extends the edge values
#' out. If not null, this value will be used for the extended cells.
#' @return A raster with edges padded \code{size} number of pixels on each edge.
#' @examples
#' library(raster)
#'
#' # import raster image
#' data(normforest)
#'
#' # crop raster to much smaller area
#' x <- pad_edges(normforest, 11)
#' @export
pad_edges <- function(x, size = 11, val = NULL) {
  if (class(x) == 'matrix') {
    # convert to equal area raster
    x <- raster(x)
    extent(x) <- c(0, ncol(x), 0, nrow(x))
    crs(x) <- st_crs(5070)$proj4string
  }

  # add padding to raster or matrix
    # continue values to edges to account for edge effect (# pixels radius/edge)
  x[is.na(x)] <- -9999999 # change NA values for now

  if (is.null(val)) {
    # first, get edge values that will be extended
    firstrow_vals <- x[1, ]
    firstcol_vals <- x[, 1]
    lastrow_vals <- x[nrow(x), ]
    lastcol_vals <- x[, ncol(x)]
  } else {
    # first, get edge values that will be extended
    firstrow_vals <- val
    firstcol_vals <- val
    lastrow_vals <- val
    lastcol_vals <- val
  }

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
  vals <- c(ext_x_mat)
  vals[vals == -9999999] <- NA
  ext_x_mat <- matrix(data = vals, nrow = dim(ext_x)[1], ncol = dim(ext_x)[2], byrow = TRUE)
  ext_x <- setValues(ext_x, ext_x_mat)

  if (class(x) == 'matrix') {
    return(as.matrix(ext_x))
  } else {
    return(ext_x)
  }
}
