#' Calculates the Average Roughness of a Surface
#'
#' Finds the average roughness of a surface (Sa) as the absolute
#' deviation of surface heights from the mean surface height.
#' Height is measured as the value of a raster and may not
#' necessarily represent actual height.
#'
#' @param x A raster or matrix.
#' @return A value of average roughness in the units of the
#'   original raster or matrix.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # find the surface roughness
#' roughness <- sa(normforest)
#' @export
sa <- function(x) {
  if(!('RasterLayer' %in% class(x)) & !('matrix' %in% class(x))) {stop('x must be a raster or matrix.')}

  if ('RasterLayer' %in% class(x)) {
    z <- getValues(x)
  } else {
    z <- x
  }

  zbar <- mean(z, na.rm = TRUE)
  N <- length(z)

  val <- sum(abs(z - zbar), na.rm = TRUE) / N

  return(val)
}

#' Calculates the Root Mean Square Roughness of a Surface
#'
#' Finds the root mean square roughness of a surface
#' (Sq) as the standard deviation of surface heights
#' from the mean surface height. Height is measured as
#' the value of a raster and may not necessarily
#' represent actual height.
#'
#' @param x A raster or matrix.
#' @return A value of root mean square roughness in
#'   the units of the original raster or matrix.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # find the surface roughness
#' roughness <- sq(normforest)
#' @export
sq <- function(x) {
  if(!('RasterLayer' %in% class(x)) & !('matrix' %in% class(x))) {stop('x must be a raster or matrix.')}

  if ('RasterLayer' %in% class(x)) {
    z <- getValues(x)
  } else {
    z <- x
  }

  val <- sd(z, na.rm = TRUE)

  return(val)
}

#' Calculates the Skewness of Raster Values
#'
#' Finds the Fisher-Pearson coefficient of skewness
#' for raster or matrix values (Ssk). Skewness represents the
#' asymmetry of the surface height distribution.
#' Height is measured as the value of a raster/matrix and
#' may not necessarily represent actual height.
#'
#' @param x A raster or matrix.
#' @param adj Logical, defaults to \code{TRUE}. If \code{TRUE},
#'   the adjusted Fisher-Pearson coefficient of skewness
#'   is calculated. Otherwise, the standard coefficient is
#'   calculated.
#' @return A numeric value representing skewness.
#' # import raster image
#' data(normforest)
#'
#' # find the adjusted coefficient of skewness
#' Ssk <- ssk(normforest, adj = TRUE)
#' @export
ssk <- function(x, adj = TRUE) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}
  if(class(adj) != 'logical') {stop('adj argument must be TRUE/FALSE.')}

  if (class(x) == 'RasterLayer') {
    z <- getValues(x)
  } else {
    z <- x
  }

  zbar <- mean(z, na.rm = TRUE)
  s <- sd(z, na.rm = TRUE)
  N <- length(z)

  val_unadj <- (sum((z - zbar) ^ 3, na.rm = TRUE) / N) / (s ^ 3)

  if (adj == TRUE) {
    val <- (sqrt((N * (N - 1))) / (N - 2)) * val_unadj # adjusted Fisher-Pearson coefficient of skewness
  } else {
    val <- val_unadj # Fisher-Pearson coefficient of skewness
  }

  return(val)
}

#' Calculates the Kurtosis of Raster Values
#'
#' Finds the kurtosis for a distribution of raster or matrix
#' values (Sku). Kurtosis represents the peakedness
#' of the raster surface height distribution. Height
#' is measured as the value of a raster/matrix and may not
#' necessarily represent actual height.
#'
#' @param x A raster or matrix.
#' @param excess Logical, defaults to \code{TRUE}. If
#'   \code{TRUE}, excess kurtosis is calculated. If \code{FALSE},
#'   kurtosis is calculated as the difference from the
#'   normal distribution.
#' @return A numeric value representing kurtosis.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # find the excess kurtosis of the raster distribution
#' Sku <- sku(normforest, excess = TRUE)
#' @export
sku <- function(x, excess = TRUE) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}
  if(class(excess) != 'logical') {stop('excess argument must be TRUE/FALSE.')}

  if (class(x) == 'RasterLayer') {
    z <- getValues(x)
  } else {
    z <- x
  }

  zbar <- mean(z, na.rm = TRUE)
  s <- sd(z, na.rm = TRUE)
  N <- length(z)

  val_unadj <- (sum((z - zbar) ^ 4, na.rm = TRUE) / N) / (s ^ 4)

  if (excess == TRUE) {
    val <- val_unadj - 3 # excess kurtosis (i.e., diff from normal distribution kurtosis)
  } else {
    val <- val_unadj
  }

  return(val)
}

#' Calculates the Maximum Valley Depth of a Surface Image
#'
#' Finds the absolute value of the lowest value in the
#' landscape (maximum valley depth; Sv) for a raster or matrix
#' representing a surface.
#'
#' @param x A raster or matrix.
#' @return A numeric value of maximum valley depth.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # find the maximum valley depth
#' Sv <- sv(normforest)
#' @export
sv <- function(x) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}

  if (class(x) == 'RasterLayer') {
    z <- getValues(x)
  } else {
    z <- x
  }

  val <- abs(min(z, na.rm = TRUE))

  return(val)
}

#' Calculates the Maximum Peak Height of a Surface Image
#'
#' Finds the absolute value of the highest value in the
#' landscape (maximum peak height; Sph) for a raster or matrix
#' representing a surface.
#'
#' @param x A raster or matrix.
#' @return A numeric value of maximum peak height.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # find the maximum peak height
#' Sph <- sph(normforest)
#' @export
sph <- function(x) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}

  if (class(x) == 'RasterLayer') {
    z <- getValues(x)
  } else {
    z <- x
  }

  val <- abs(max(z, na.rm = TRUE))

  return(val)
}

#' Calculates the Mean Peak Height of a Surface Image
#'
#' Finds the mean height of positive values in the
#' landscape (mean peak height; Smean) for a raster or matrix
#' representing a surface.
#'
#' @param x A raster or matrix.
#' @return A numeric value of mean peak height.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # find the maximum peak height
#' Smean <- smean(normforest)
#' @export
smean <- function(x) {
  if(class(x) != 'RasterLayer' & class(x) != 'matrix') {stop('x must be a raster or matrix.')}

  if (class(x) == 'RasterLayer') {
    z <- getValues(x)
  } else {
    z <- x
  }

  val <- mean(z[z > 0], na.rm = TRUE)

  return(val)
}
