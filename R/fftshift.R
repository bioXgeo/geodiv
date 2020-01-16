#' Fourier Transform Shift
#'
#' This function serves to shift the zero-frequency component of the
#' Fourier transform to the center of the matrix.
#'
#' @param x An n x n Fourier transform matrix.
#' @param dim Which dimension to shift the matrix. -1 swaps up/down and
#' left/right. 1 swaps up/down. 2 swaps left/right.
#' @return An n x n matrix with the zero-frequency component of
#'   the Fourier transform in the center.
#' @references #' This function was created from code posted by rayryeng at:
#' https://stackoverflow.com/questions/38230794/how-to-write-fftshift-and-ifftshift-in-r.
#' @examples
#' library(raster)
#'
#' # import raster image
#' data(normforest)
#'
#' # convert to matrix form
#' M <- ncol(normforest)
#' N <- nrow(normforest)
#' zmat <- matrix(raster::getValues(normforest), ncol = M, nrow = N, byrow = TRUE)
#'
#' # calculate fourier transform and shift
#' ftmat <- fft(zmat)
#' ftshift <- fftshift(ftmat)
#'
#' # plot real component
#' r <- setValues(normforest, Re(ftshift))
#' plot(r)
#' @export
fftshift <- function(x, dim = -1) {
  if(length(base::class(x)) > 1) {
    if(!('matrix' %in% class(x))) {stop('x must be a matrix.')}
  } else {
    if(base::class(x) != 'matrix') {stop('x must be a matrix.')}
  }
  if(length(dim) > 1) {stop('too many values provided for dim.')}
  if(class(dim) != 'numeric') {stop('dim must be numeric.')}
  if(dim != -1 & dim != 1 & dim != 2) {stop('invalid value for dim -- must be -1, 1, or 2.')}

  rows <- dim(x)[1]
  cols <- dim(x)[2]

  swap_up_down <- function(x) {
    rows_half <- ceiling(rows / 2)
    return(rbind(x[((rows_half + 1):rows), (1:cols)], x[(1:rows_half), (1:cols)]))
  }

  swap_left_right <- function(x) {
    cols_half <- ceiling(cols / 2)
    return(cbind(x[1:rows, ((cols_half + 1):cols)], x[1:rows, 1:cols_half]))
  }

  if (dim == -1) {
    x <- swap_up_down(x)
    return(swap_left_right(x))
  }
  else if (dim == 1) {
    return(swap_up_down(x))
  }
  else if (dim == 2) {
    return(swap_left_right(x))
  }
  else {
    stop("Invalid dimension parameter")
  }
}
