#' Offset Raster or Matrix Values
#'
#' Calculates a matrix of values with a negative
#' or positive, x or y, offset.
#'
#' @param r A raster or matrix.
#' @param xdist Numeric indicating the number and direction (+, -)
#'   of columns for the offset.
#' @param ydist Numeric indicating the number and direction (+, -)
#'   of rows for the offset.
#' @param xrm Numeric value or vector indicating the number of
#'   columns to be removed from the final matrix. If not set,
#'   this value defaults to \code{xdist}. Positive values remove
#'   columns from the right, while negative values remove columns
#'   from the left. The absolute value of \code{xrm} must be
#'   \code{>= abs(xdist)}.
#' @param yrm Numeric value or vector indicating the number
#'   of rows to be removed from the final matrix. If not set,
#'   this value defaults to \code{ydist}. Positive values remove
#'   rows from the bottom, while negative values remove rows from
#'   the top. The absolute value must be \code{>= abs(ydist)}.
#' @param scale Logical. Indicates whether or not to scale the values of
#'   the raster.
#' @return A numeric vector of values created from a matrix of the values
#'   with the specified offset. The vector is created from a matrix with
#'    \code{xrm} fewer columns and \code{yrm} fewer rows than the original
#'   raster value matrix.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # remove right and bottom borders 2 deep
#' noborder <- zshift(normforest, xdist = 2, ydist = 2)
#' @export
zshift <- function(r, xdist = 0, ydist = 0, xrm, yrm, scale = FALSE) {
  # xdist is distance away in x direction
  # ydist is distance away in y direction
  # xrm, yrm can be specified if you want matrix clipped more than xdist/ydist
  # xrm, yrm cannot be less than corresponding dist, or of a different sign
  # need to provide at least one of xdist or ydist
  try(if(missing(xrm)) (xrm = xdist))
  try(if(missing(yrm)) (yrm = ydist))

  if(class(r) != 'RasterLayer' & class(r) != 'matrix') {stop('r must be a raster or matrix.')}
  if(class(xdist) != 'numeric') {stop('xdist must be numeric.')}
  if(class(ydist) != 'numeric') {stop('ydist must be numeric.')}
  if(class(xrm) != 'numeric') {stop('xrm must be numeric.')}
  if(class(yrm) != 'numeric') {stop('yrm must be numeric.')}
  if(class(scale) != 'logical') {stop('scale must be logical.')}
  if(length(xdist) > 1) {stop('too many values supplied to xdist.')}
  if(length(ydist) > 1) {stop('too many values supplied to ydist.')}

  # get dimensions
  N <- dim(r)[1] # rows
  M <- dim(r)[2] # cols

  # calculate zmat and coordinates
  if (class(r) == 'RasterLayer') {
    z <- getValues(r)
  } else if (class(r) == 'matrix') {
    z <- as.numeric(r)
  }

  if (scale == TRUE) {
    zmat <- matrix(((z - min(z, na.rm = TRUE)) / (max(z, na.rm = TRUE) - min(z, na.rm = TRUE))),
                   nrow = N, ncol = M, byrow = TRUE)
  } else {
    zmat <- matrix(z, nrow = N, ncol = M, byrow = TRUE)
  }

  # row/col of each center
  rows <- rep(1:N, each = M)
  cols <- rep(rep(1:M), N)

  # get rid of edge points x distance away
  rm_inds <- rm_indsx <- rm_indsy <- numeric(0)
  for (i in 1:length(xrm)){
    if (xrm[i] > 0) {
      posx_rm <- which(cols > (max(cols) - xrm[i]))
      rm_indsx <- c(rm_indsx, posx_rm)
    } else if (xrm[i] < 0) {
      negx_rm <- which(cols < (abs(xrm[i]) + 1))
      rm_indsx <- c(rm_indsx, negx_rm)
    } else {
      posx_rm <- NULL
      negx_rm <- NULL
      rm_indsx <- c(rm_indsx, posx_rm, negx_rm)
    }
    rm_inds <- c(rm_inds, rm_indsx)
  }
  for (i in 1:length(yrm)){
    if (yrm[i] > 0) {
      posy_rm <- which(rows > (max(rows) - yrm[i]))
      rm_indsy <- c(rm_indsy, posy_rm)
    } else if (yrm[i] < 0) {
      negy_rm <- which(rows < (abs(yrm[i]) + 1))
      rm_indsy <- c(rm_indsy, negy_rm)
    } else {
      posy_rm <- NULL
      negy_rm <- NULL
      rm_indsy <- c(rm_indsy, posy_rm, negy_rm)
    }
    rm_inds <- c(rm_inds, rm_indsy)
  }

  if (length(rm_inds) < 1) {
    z <- z
    rows <- rows
    cols <- cols
  } else {
    z <- z[-rm_inds]
    rows <- rows[-rm_inds]
    cols <- cols[-rm_inds]
  }

  yshift <- rows + ydist
  xshift <- cols + xdist
  ind <- seq(1, length(z))
  if (xdist != 0 & ydist != 0) {
    z_shift <- unlist(lapply(ind, function(i) {zmat[yshift[i], xshift[i]]}))
  } else if (xdist != 0 & ydist == 0){
    z_shift <- unlist(lapply(ind, function(i) {zmat[rows[i], xshift[i]]}))
  } else if (xdist == 0 & ydist != 0) {
    z_shift <- unlist(lapply(ind, function(i) {zmat[yshift[i], cols[i]]}))
  } else {
    z_shift = unlist(lapply(ind, function(i) {zmat[rows[i], cols[i]]}))
  }

  if (xdist < 0 & length(xrm) <= 1) {
    z_shift <- matrix(z_shift, nrow = length(unique(yshift)),
                      ncol = length(unique(xshift)), byrow = TRUE)
    z_shift <- cbind(rep(rep(NA, abs(xdist)), nrow(z_shift)), z_shift)
    z_shift <- as.numeric(z_shift)
  }
  if (ydist < 0 & length(yrm) <= 1) {
    z_shift <- matrix(z_shift, nrow = length(unique(yshift)),
                      ncol = length(unique(xshift)), byrow = TRUE)
    z_shift <- rbind(rep(rep(NA, abs(ydist)), ncol(z_shift)), z_shift)
    z_shift <- as.numeric(z_shift)
  }

  return(z_shift)
}
