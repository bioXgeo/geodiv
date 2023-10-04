#' Calculate Texture Metrics per Pixel
#'
#' Calculates the various texture metrics over windows centered
#' on individual pixels. This creates a continuous surface of the
#' texture metric.
#'
#' Note that this function is meant to work on rasters with an equal area
#' projection.
#'
#' @param x A raster or matrix.
#' @param window_type Character. Type of window, either circular or square.
#' @param size Numeric. Size of window (edge length) or diameter (in meters).
#' @param in_meters Logical. Is the size given in meters?
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
#'    \item{\code{'ssk'}: skewness}
#'    \item{\code{'sku'}: kurtosis}
#'    \item{\code{'sds'}: summit density}
#'    \item{\code{'sfd'}: 3d fractal dimension}
#'    \item{\code{'srw'}: dominant radial wavelength, radial wavelength index, mean half wavelength}
#'    \item{\code{'std'}: angle of dominating texture, texture direction index}
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
#' # import raster image
#' data(normforest)
#' normforest <- terra::unwrap(normforest)
#'
#' # crop raster to smaller area
#' x <- terra::crop(normforest, terra::ext(normforest[1:100, 1:100, drop = FALSE]))
#'
#' # get a surface of root mean square roughness
#' sa_img <- texture_image(x = x, window = 'square',
#' size = 5, metric = 'sa',
#' parallel = TRUE, ncores = 1, nclumps = 20)
#'
#' # plot the result
#' terra::plot(sa_img)
#' @importFrom terra rast crop crds
#' @importFrom parallel mclapply clusterExport clusterEvalQ parLapply makeCluster stopCluster
#' @export
texture_image <- function(x, window_type = 'square', size = 5, in_meters = FALSE,
                          metric, args = NULL, parallel = TRUE, ncores = NULL, nclumps = 100){

  stopifnot('x must be a raster or matrix.' = inherits(x, c('RasterLayer', 'matrix', 'SpatRaster')))
  stopifnot('window_type must be a string.' = inherits(window_type, 'character'))
  stopifnot('size must be numeric.' = inherits(size, 'numeric'))
  stopifnot('metric must be a character.' = inherits(metric, 'character'))

  if(length(metric) > 1) {stop('too many values provided for metric.')}
  if(length(size) > 1) {stop('too many values provided to size.')}
  if(length(window_type) > 1) {stop('too many values provided to window_type.')}

  if(missing(ncores)) {ncores <- parallel::detectCores() - 1}

  if(Sys.info()['sysname'][[1]] == 'Windows') {
    print('mclapply is not supported on Windows, using parLapply instead.')
    os_type = 'windows'
  } else {
    os_type = 'other'
  }

  # Is size argument is in meters, convert it to a number of pixels
  if (in_meters == TRUE) {

    # get equivalent # pixels of size (this assumes an equal area projection)
    size <- ceiling(size / res(x))[1]

  }

  # change to matrix for faster processing
  new_x <- as.matrix(x)

  # add padding to matrix and create matrix of same size (all NA) for finding
  # which values are extra later
  shift <- floor(size / 2)
  ext_x <- dummy_ext_x <- pad_edges(new_x, size = shift, val = NA)
  dummy_ext_x[] <- NA

  # find indices of center, left, right, top, and bottom of each window
  start_ind <- 1 + shift
  center <- left <- right <- top <- bottom <- dummy_ext_x
  # first, set those values to 1 for each side
  center[start_ind:(shift + nrow(new_x)), start_ind:(shift + ncol(new_x))] <- 1
  left[start_ind:(shift + nrow(new_x)), 1:(ncol(new_x))] <- 1
  right[start_ind:(shift + nrow(new_x)), (start_ind + shift):((shift * 2) + ncol(new_x))] <- 1
  top[1:nrow(new_x), (shift + 1):(shift + ncol(new_x))] <- 1
  bottom[(start_ind + shift):((shift * 2) + nrow(new_x)), start_ind:(shift + ncol(new_x))] <- 1
  # get those indices
  center_inds <- which(center == 1)
  left_inds <- which(left == 1, arr.ind = TRUE)
  right_inds <- which(right == 1, arr.ind = TRUE)
  top_inds <- which(top == 1, arr.ind = TRUE)
  bottom_inds <- which(bottom == 1, arr.ind = TRUE)

  # combine window centers and edges into data table
  coord_list <- data.frame(center = center_inds, left = left_inds[, 2], right = right_inds[, 2],
                           top = top_inds[, 1], bottom = bottom_inds[, 1])

  # collect arguments
  if (!is.null(args)) {

    input_args <- args

  } else {

    input_args <- NULL

  }

  # sequence of coord rows to run
  pixlist <- seq(1, nrow(coord_list))

  # get list of # total pixels, break up into smaller lists (by number of cores available)
  seg_length <- ceiling(length(pixlist) / nclumps)
  segment <- ceiling(length(pixlist) / seg_length)
  new_pixlist <- vector('list', segment)
  i <- 1
  while(i <= segment) {
    newi <- ceiling(seq(pixlist[1], max(pixlist), seg_length))[i]
    if (i < segment) {
      new_pixlist[[i]] <- seq(newi, (newi - 1) + seg_length, 1)
    } else {
      new_pixlist[[i]] <- seq(newi, max(pixlist), 1)
    }
    i <- i + 1
  }

  if (parallel == FALSE) {

    print('Beginning calculation of metrics over windows...')
    start <- Sys.time()
    result <- lapply(pixlist, FUN = function(i) {window_metric(x = ext_x, coords = coord_list[i, ],
                                                                   window_type = window_type,
                                                                   size = size,
                                                                   metric = metric,
                                                                   args = input_args)})
    end <- Sys.time()
    cat('Total time to calculate metrics: ', end - start, '\n', sep = '')
    result <- unlist(result)

  } else if (os_type != 'windows' & parallel == TRUE) {

    print('Beginning calculation of metrics over windows...')
    start <- Sys.time()
    result <- parallel::mclapply(pixlist, FUN = function(l) {
      lapply(l, FUN = function(i) {window_metric(x = ext_x, coords = coord_list[i, ], window_type = window_type,
                                                     size = size, metric = metric,
                                                     args = input_args)})
    }, mc.cores = ncores, mc.cleanup = TRUE)
    end <- Sys.time()
    cat('Total time to calculate metrics: ', end - start, '\n', sep = '')

    result <- do.call(unlist, args = list(result, recursive = TRUE))

  } else {

    print('Beginning calculation of metrics over windows...')
    start <- Sys.time()
    # make and start cluster
    try(parallel::stopCluster(cl), silent = TRUE)
    cl <- parallel::makeCluster(ncores, type = 'PSOCK')
    parallel::clusterExport(cl = cl, list('ext_x', 'coord_list', 'size',
                                          'window_type',
                                          'new_pixlist', 'metric', 'input_args',
                                          'window_metric'),
                            envir = environment())
    parallel::clusterEvalQ(cl, library('geodiv'))
    parallel::clusterEvalQ(cl, library('terra'))
    # for each list in new_pixlist, run lapply
    result <- parallel::parLapply(cl, new_pixlist, fun = function(l) {
      lapply(l, FUN = function(i) {window_metric(x = ext_x, coords = coord_list[i, ],
                                                     window_type = window_type,
                                                     size = size, metric = metric,
                                                     args = input_args)})
    })
    parallel::stopCluster(cl)
    end <- Sys.time()
    cat('Total time to calculate metrics: ', end - start, '\n', sep = '')

    result <- do.call(unlist, args = list(result, recursive = TRUE))

  }

  # deal with functions that have multiple outputs
  outfinal <- list()
  nresult <- length(result) / nrow(coord_list)
  for (i in 1:nresult) {
    temp <- data.frame(newvals = result[seq(i, length(result), nresult)], ind = pixlist)
    outfinal[[i]] <- x
    if (class(x)[1] %in% c('RasterLayer', 'SpatRaster')) {
      outfinal[[i]] <- setValues(x, matrix(temp$newvals, nrow = nrow(x), ncol = ncol(x)))
    } else {
      outfinal[[i]] <- matrix(temp$newvals, nrow = nrow(x), ncol = ncol(x))
    }
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
#' @param coords Dataframe. Coordinates of window edges.
#' @param window_type Character. Type of window, either circular or square.
#' @param size Numeric. Edge length or diameter of window in number of pixels.
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
#'    \item{\code{'ssk'}: skewness}
#'    \item{\code{'sku'}: kurtosis}
#'    \item{\code{'sds'}: summit density}
#'    \item{\code{'sfd'}: 3d fractal dimension}
#'    \item{\code{'srw'}: dominant radial wavelength, radial wavelength index, mean half wavelength}
#'    \item{\code{'std'}: angle of dominating texture, texture direction index}
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
#' @importFrom terra crds crop rast
#' @export
window_metric <- function(x, coords, window_type = 'square', size = 11, metric, args = NULL) {

  # crop to square
  cropped_x <- x[coords$top:coords$bottom, coords$left:coords$right]

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

#' Extend edges of a matrix.
#'
#' Extends edge values of a raster or matrix by a specified number of pixels.
#'
#' @param x A matrix.
#' @param size Numeric. Number of pixels to add to each side.
#' @param val Numeric. If NULL (default), this extends the edge values
#' out. If not null, this value will be used for the extended cells.
#' @return A raster with edges padded \code{size} number of pixels on each edge.
#' @examples
#' # import raster image
#' data(normforest)
#' normforest <- terra::unwrap(normforest)
#'
#' # crop raster to much smaller area
#' x <- pad_edges(as.matrix(normforest), 3, val = NA)
#' @importFrom zoo na.approx
#' @export
pad_edges <- function(x, size, val = NULL) {

  # add padding to matrix
  # continue values to edges to account for edge effect (# pixels radius/edge)
  x[is.na(x)] <- -999999999 # change NA values for now

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

  # Create new matrix to be filled. This will have an increased extent (by
  # 2*size on each side).
  ext_x <- matrix(NA, nrow = nrow(x) + (2 * size), ncol = ncol(x) + (2 * size))

  # fill in top rows
  ext_x[1:size, (size + 1):(ncol(ext_x) - size)] <- matrix(rep(firstrow_vals, size), nrow = size, byrow = TRUE)
  # fill in bottom rows
  ext_x[(nrow(ext_x) - (size - 1)):nrow(ext_x), (size + 1):(ncol(ext_x) - size)] <- matrix(rep(lastrow_vals, size), nrow = size, byrow = TRUE)
  # fill in left columns
  ext_x[(size + 1):(nrow(ext_x) - size), 1:size] <- matrix(rep(firstcol_vals, size), ncol = size, byrow = FALSE)
  # fill in right columns
  ext_x[(size + 1):(nrow(ext_x) - size), (ncol(ext_x) - (size - 1)):ncol(ext_x)] <- matrix(rep(lastcol_vals, size), ncol = size, byrow = FALSE)
  # fill in middle
  ext_x[(size + 1):(nrow(ext_x) - size), (size + 1):(ncol(ext_x) - size)] <- x

  # fill in corners with nearest point value (always the same)
  if (is.null(val)) {
    ext_x_mat <- zoo::na.approx(matrix(ext_x, nrow = nrow(ext_x), ncol = ncol(ext_x)), rule = 2)
  } else {
    ext_x_mat <- ext_x
    ext_x_mat[is.na(ext_x_mat)] <- val
  }

  # If the values were NA before (-9999999 now), convert them back to NA.
  ext_x_mat[ext_x_mat == -999999999] <- NA

  return(ext_x_mat)
}
