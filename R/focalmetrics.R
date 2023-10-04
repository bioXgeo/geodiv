#' Calculate Texture Metrics per Pixel
#'
#' Calculates the various texture metrics over windows centered
#' on individual pixels. This creates a continuous surface of the
#' texture metric.
#' This function is a modified version of the \code{window_lsm} function from the
#' \emph{landscapemetrics} package (Hesselbarth et al. 2019).
#'
#' @param x A raster or matrix. Image over which to apply focal window calculations.
#' @param window Matrix. The focal window used to create the image.
#' @param metrics List. List of metrics to apply. Function names must be strings.
#' @param progress Logical. Display progress through metrics list?
#' @param ... Additional arguments for the metric functions. All applicable arguments
#' will be applied to the entire list of metrics.
#' @return A raster of the metric calculated in windows over the raster or matrix.
#' If the input was a matrix, the function will return a raster with an extent of [0, 1, 0, 1].
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
#' @references
#' \enumerate{
#' \item{Hesselbarth, M.H.K., Sciaini, M., With, K.A., Wiegand, K., Nowosad, J. 2019.
#' landscapemetrics: an open-source R tool to calculate landscape metrics. - Ecography 42:1648-1657(ver. 0).}
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
#' sa_img <- focal_metrics(x = x, window = matrix(1, 5, 5),
#'                         metrics = list('sa'), progress = TRUE)
#'
#' # plot the result
#' terra::plot(sa_img$sa)
#' @importFrom terra rast crds ext focal
#' @export
focal_metrics <- function(x,
                          window,
                          metrics,
                          progress,
                          ...) {

  # check if window has uneven sides
  if (any(dim(window) %% 2 == 0)) {

    stop("The window must have uneven sides.", call. = FALSE)
  }

  if (inherits(x, "matrix") == TRUE) {

    x <- rast(x)

  }

  if (class(x)[1] == 'RasterLayer') {

    x <- rast(x)

  }

  number_metrics <- length(metrics)

  # get coordinates of cells
  points <- cbind(crds(x, na.rm = FALSE), x[])
  colnames(points)[3] <- 'z'

  # get dimensions of window
  n_row <- nrow(window)
  n_col <- ncol(window)

  # create object for warning messages
  warning_messages <- character(0)

  result <- withCallingHandlers(expr = {lapply(seq_along(metrics), function(current_metric) {

    # print progess using the non-internal name
    if (progress) {
      cat("\r> Progress metrics: ", current_metric, "/", number_metrics)
    }

    terra::focal(x = x, w = window, fun = function(x) {

      .calculate_met_focal(landscape = x,
                           n_row = n_row,
                           n_col = n_col,
                           points = points,
                           what = metrics[[current_metric]],
                           ...)})
  })},
  warning = function(cond) {

    warning_messages <<- c(warning_messages, conditionMessage(cond))

    invokeRestart("muffleWarning")})

  names(result) <- metrics

  if (progress) {cat("\n")}

  # warnings present
  if (length(warning_messages) > 0) {

    # only unique warnings
    warning_messages <- unique(warning_messages)

    # print warnings
    lapply(warning_messages, function(x){ warning(x, call. = FALSE)})
  }

  return(result)
}


#' Calculate Texture Metric for Single Pixel
#'
#' Calculates the various texture metrics over a window centered
#' on an individual pixel. This function is modified slightly from the
#' \code{calculate_lsm_focal} function in the \emph{landscapemetrics} package (Hesselbarth et al. 2019).
#'
#' @param landscape A raster or matrix.
#' @param n_row Numeric. Number of rows in focal window.
#' @param n_col Numeric. Number of columns in focal window.
#' @param points Dataframe. Coordinates and values of cells, calculated with the *landscapemetrics*
#' \code{raster_to_points} function.
#' @param what Character. Metric to calculate for each window. Metrics
#' from the geodiv package are listed below.
#' @param ... Additional arguments for the metric functions. All applicable arguments
#' will be applied to the entire list of metrics.
#' @return The metric value over the window.
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
#' @references
#' \enumerate{
#' \item{Hesselbarth, M.H.K., Sciaini, M., With, K.A., Wiegand, K., Nowosad, J. 2019.
#' landscapemetrics: an open-source R tool to calculate landscape metrics. - Ecography 42:1648-1657(ver. 0).}
#' }
#' @importFrom terra rast crds ext focal
#' @export
.calculate_met_focal <- function(landscape,
                                n_row,
                                n_col,
                                points,
                                what,
                                ...) {

  # convert focal window to matrix
  raster_window <- matrix(landscape, n_row, n_col)

  # match function name
  foo <- get(what, mode = "function")

  # get argument
  arguments <- names(formals(foo))[-1]

  # get provided arguments
  arguments_provided <- substitute(...())

  # landscape argument
  arguments_values <- list(raster_window)

  # sort alphabetically to match later with defaults
  if (!is.null(arguments_provided)) {

    arguments_provided <- arguments_provided[order(names(arguments_provided))]

    # exchange arguments
    arguments[arguments %in% names(arguments_provided)] <- arguments_provided

    # replace general argument option with specific if na.rm is present
    if ("..." %in% arguments & "na.rm" %in% names(arguments_provided)) {
      arguments <- c(arguments, arguments_provided[which(names(arguments_provided) == "na.rm")])
    }
    # remove general dots argument
    if ("..." %in% arguments) {
      dots_ind <- which(arguments == "...")
      arguments <- arguments[-dots_ind]
    }

    # combine input raster with argument values
    arguments_values <- c(arguments_values, arguments)
  }

  # run function
  result <- do.call(what = foo,
                    args = arguments_values)

  return(result)
}
