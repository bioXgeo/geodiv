#' Texture Direction Metrics
#'
#' Calculates the angle of dominating texture and the texture
#' direction index of the Fourier spectrum image calculated
#' from a raster image (see Kedron et al. 2018).
#'
#' @param x A raster.
#' @param plot Logical. If \code{TRUE}, returns a plot of the
#'   amplitude spectrum with lines showing directions in which
#'   amplitude is summed for the Std and Stdi calculations.
#' @return A list containing numeric values for the angle of
#'   dominating texture and the texture direction index.
#' @examples
#' library(raster)
#'
#' # import raster image
#' data(normforest)
#'
#' # crop raster to much smaller area
#' x <- crop(normforest, extent(-123, -122.99, 43, 43.01))
#'
#' # calculate Std and Stdi
#' stdvals <- std(x)
#'
#' # extract each value
#' Std <- stdvals[[1]]
#' Stdi <- stdvals[[2]]
#' @export
std <- function(x, plot = FALSE) {
  if(class(x) != 'RasterLayer') {stop('x must be a raster.')}
  if(class(plot) != 'logical') {stop('plot must be logical.')}

  # get raster dimensions
  M <- ncol(x)
  N <- nrow(x)

  # get matrix of values
  zmat <- matrix(getValues(x), ncol = M, nrow = N, byrow = TRUE)

  # if irregular non-na area, cut to biggest square possible
  if (sum(is.na(zmat)) != 0) {
    origin <- c(mean(sp::coordinates(x)[, 1]), mean(sp::coordinates(x)[, 2]))
    potentials <- data.frame(xmin = rep(origin[1], floor(N / 2)),
                             xmax = origin[1],
                             ymin = origin[2],
                             ymax = origin[2])
    potentials$xmin <- origin[1] - (res(x)[1] * seq(1, round(N / 2)))
    potentials$xmax <- origin[1] + (res(x)[1] * seq(1, round(N / 2)))
    potentials$ymin <- origin[2] - (res(x)[2] * seq(1, round(N / 2)))
    potentials$ymax <- origin[2] + (res(x)[2] * seq(1, round(N / 2)))

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
  ft_shift <- fftshift(ft)

  # amplitude spectrum
  amplitude <- sqrt((Re(ft_shift) ^ 2) + (Im(ft_shift) ^ 2))

  # take amplitude image, cut in half (y direction)
  amp_img <- setValues(x, amplitude)
  half_dist <- (ymax(amp_img) - ymin(amp_img)) / 2
  ymin <- ymax(amp_img) - half_dist
  amp_img <- crop(amp_img, c(xmin(amp_img), xmax(amp_img), ymin, ymax(amp_img)))

  # get origin of image (actually bottom center)
  origin <- c(mean(sp::coordinates(amp_img)[,1]), ymin(amp_img))

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
  if(plot == TRUE) {
    plot(amp_img)
    lines(lines)
  }

  Aalpha <- list()
  for (i in 1:length(lines)) {
    Aalpha[i] <- extract(amp_img, lines[i], fun = sum)
  }

  std <- .rad2deg(alpha[which(unlist(Aalpha) == max(unlist(Aalpha), na.rm = TRUE))])
  stdi <- mean(unlist(Aalpha), na.rm = TRUE) / max(unlist(Aalpha), na.rm = TRUE)

  return(list(std, stdi))
}

#' Radial Wavelength Metrics
#'
#' Calculates the dominant radial wavelength, radial wavelength
#' index, and mean half wavelength of the radial Fourier spectrum.
#' See Kedron et al. (2018) for more detailed description.
#'
#' @param x A raster.
#' @param plot Logical. If \code{TRUE}, returns a plot of the
#'   amplitude spectrum with lines showing the radii at which
#'   Srw, Srwi, and Shw are calculated.
#' @return A list containing numeric values for the dominant
#'   radial wavelength, radial wavelength index, and mean half
#'   wavelength.
#' @examples
#' library(raster)
#'
#' # import raster image
#' data(normforest)
#'
#' # crop raster to much smaller area
#' x <- crop(normforest, extent(-123, -122.99, 43, 43.01))
#'
#' # calculate metrics
#' srwvals <- srw(x)
#'
#' # extract each value
#' Srw <- srwvals[[1]]
#' Srwi <- srwvals[[2]]
#' Shw <- srwvals[[3]]
#' @export
srw <- function(x, plot = FALSE) {
  if(class(x) != 'RasterLayer') {stop('x must be a raster.')}
  if(class(plot) != 'logical') {stop('plot must be logical.')}

  # get raster dimensions
  M <- ncol(x)
  N <- nrow(x)

  # get matrix of values
  zmat <- matrix(getValues(x), ncol = M, nrow = N, byrow = TRUE)

  # if irregular non-na area, cut to biggest square possible
  if (sum(is.na(zmat)) != 0) {
    origin <- c(mean(sp::coordinates(x)[, 1]), mean(sp::coordinates(x)[, 2]))
    potentials <- data.frame(xmin = rep(origin[1], floor(N / 2)),
                             xmax = origin[1],
                             ymin = origin[2],
                             ymax = origin[2])
    potentials$xmin <- origin[1] - (res(x)[1] * seq(1, round(N / 2)))
    potentials$xmax <- origin[1] + (res(x)[1] * seq(1, round(N / 2)))
    potentials$ymin <- origin[2] - (res(x)[2] * seq(1, round(N / 2)))
    potentials$ymax <- origin[2] + (res(x)[2] * seq(1, round(N / 2)))

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
  ft_shift <- fftshift(ft)

  # amplitude spectrum
  amplitude <- sqrt((Re(ft_shift) ^ 2) + (Im(ft_shift) ^ 2))

  # take amplitude image, cut in half (y direction)
  amp_img <- setValues(x, amplitude)
  half_dist <- (ymax(amp_img) - ymin(amp_img)) / 2
  ymin <- ymax(amp_img) - half_dist
  amp_img <- crop(amp_img, c(xmin(amp_img), xmax(amp_img), ymin, ymax(amp_img)))

  # get origin of image (actually bottom center)
  origin <- c(mean(sp::coordinates(amp_img)[,1]), ymin(amp_img))

  # calculate half circles extending from origin
  if ((0.5 * ncol(amp_img)) <= 100) {
    ncircles <- floor(0.5 * ncol(amp_img))
  } else {
    ncircles <- 100
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
  if (plot == TRUE){
    plot(amp_img)
    plot(lines, add = TRUE)
  }
  Br <- list()
  Br[1] <- 0
  for (i in 2:length(lines)) {
    templine <- crop(lines[i], extent(xmin(amp_img), xmax(amp_img), ymin, ymax(amp_img)))
    Br[i] <- extract(amp_img, templine, fun = sum)
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

  return(list(Srw, Srwi, Shw))
}

#' Calculates the Fractal Dimension
#'
#' Calculates the 2D fractal dimension of a raster using the
#' Fourier transform. This function is a modified version, adapted
#' to R, of the Matlab function 'fdsurfft' version 1.0.0.0
#' developed by Jianbo Zhang.
#'
#' @param x A raster.
#' @return A numeric value representing the 2D fractal dimension of
#'   the image.
#' @examples
#' library(raster)
#'
#' # import raster image
#' data(normforest)
#'
#' # crop raster to much smaller area
#' x <- crop(normforest, extent(-123, -122.99, 43, 43.01))
#'
#' # calculate the fractal dimension
#' Sfd <- sfd_old(x)
#' @export
sfd_old <- function(x) {
  if(class(x) != 'RasterLayer') {stop('x must be a raster.')}

  # this has been checked against the matlab version and produces the same results

  # get raster dimensions
  M <- ncol(x)
  N <- nrow(x)

  # get matrix of values
  zmat <- matrix(getValues(x), ncol = M, nrow = N, byrow = TRUE)

  # if irregular non-na area, cut to biggest square possible
  if (sum(is.na(zmat)) != 0) {
    origin <- c(mean(sp::coordinates(x)[, 1]), mean(sp::coordinates(x)[, 2]))
    potentials <- data.frame(xmin = rep(origin[1], floor(N / 2)),
                             xmax = origin[1],
                             ymin = origin[2],
                             ymax = origin[2])
    potentials$xmin <- origin[1] - (res(x)[1] * seq(1, round(N / 2)))
    potentials$xmax <- origin[1] + (res(x)[1] * seq(1, round(N / 2)))
    potentials$ymin <- origin[2] - (res(x)[2] * seq(1, round(N / 2)))
    potentials$ymax <- origin[2] + (res(x)[2] * seq(1, round(N / 2)))

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
  ft_shift <- fftshift(ft)

  # amplitude spectrum
  amplitude <- sqrt((Re(ft_shift) ^ 2) + (Im(ft_shift) ^ 2))

  # take amplitude image, cut in half (y direction)
  amp_img <- setValues(x, amplitude)
  half_dist <- (ymax(amp_img) - ymin(amp_img)) / 2
  ymin <- ymax(amp_img) - half_dist
  amp_img <- crop(amp_img, c(xmin(amp_img), xmax(amp_img), ymin, ymax(amp_img)))

  # get origin of image (actually bottom center)
  origin <- c(mean(sp::coordinates(amp_img)[,1]), ymin(amp_img))

  # from matlab fdsurfft function
  num_dir <- 24 # number of directions that the frequency space is uniformally divided
  num_rad <- 30 # number of points that the radius is uniformally divided

  M <- nrow(ft_shift)
  N <- ncol(ft_shift)
  xcoords <- sp::coordinates(x)[, 1]
  ycoords <- sp::coordinates(x)[, 2]
  xdist <- matrix(xcoords, nrow = M, ncol = N, byrow = TRUE)[1, ] - origin[1]
  ydist <- matrix(ycoords, nrow = M, ncol = N, byrow = TRUE)[, 1] - origin[2]
  xctr <- which(xdist == min(abs(xdist))) # should actually be index of x that = origin x
  yctr <- which(ydist == min(abs(ydist))) # same as xctr
  fim <- ft_shift

  # calculate power spectrum
  mag <- Re(log(fim * fim + 10 ^ (-6)))
  sumBrite <- pracma::zeros(num_dir, num_rad) # accumulation magnitude for each direction and radius
  nCount <- pracma::zeros(num_dir, num_rad) # number of magnitude
  radius <- pracma::zeros(2 * num_rad, 1) # accumulation magnitude for all directions
  radCount <- pracma::zeros(2 * num_rad, 1) # number of magnitude for all directions

  # Compute phase image and phase histogram
  phaseim <- pracma::zeros(M, N)
  phase <- phonTools::zeros(180)
  for (j  in 1:M) {
    for (i in 1:N) {
      realv <- Re(fim[j, i])
      imagv <- Im(fim[j, i])
      if (realv == 0) {
        value = pi / 2
      } else {
        value <- atan((imagv / realv))
        phaseim[j, i] <- value
        ang <- floor(180 * (pi / 2 + value) / pi)
      }
      if (ang < 0) {
        ang <- 0
      }
      if (ang > 179) {
        ang <- 179
      }
      phase[ang + 1] <- phase[ang + 1] + 1
    }
  }

  maxphase <- max(phase)

  # accumulation of magnitude for each direction and radius
  rmax <- log(min(M, N) / 2) # maximum radius (log scale)
  for (j in 1:M) {
    if (j != yctr) {
      yval <- yctr - j
      y2 <- yval * yval
    for (i in 1:N) {
      if (i != xctr) {
        xval <- i - xctr
        rho <- log(sqrt(y2 + xval * xval))
        if (rho > 0 & rho <= rmax) {
          mval <- mag[j, i]
          temp <- yval / xval
          theta <- atan(temp)
          if (xval < 0) {
            theta <- theta + pi
          }
          if (theta < 0) {
            theta <- theta + 2 * pi
          }
          ang <- floor(num_dir * theta / (2 * pi))
          if (ang > num_dir - 1 | ang < 0) {
            ang <- num_dir - 1
          }
          k <- floor(2 * num_rad * rho / rmax)
          h <- floor(k / 2)
          if (k > 2 * num_rad - 1) {
            h <- num_rad - 1
            k <- 2 * num_rad - 1
          }
          if (h >= 5) {
            sumBrite[ang + 1, h + 1] <- sumBrite[ang + 1, h + 1] + mval
            nCount[ang + 1, h + 1] <- nCount[ang + 1, h + 1] + 1
          }
          if (k >= 5) {
            radius[k + 1] <- radius[k + 1] + mval
            radCount[k + 1] <- radCount[k + 1] + 1
          }
        }
      }
    }
    }
  }

  # linear regression
  slope <- as.numeric()
  intercept <- as.numeric()
  for (ang in 1:num_dir) {
    sumx <- 0
    sumy <- 0
    sumx2 <- 0
    sumxy <- 0
    sumn <- 0
    for (range in 6:num_rad) {
      if (nCount[ang, range] > 0) {
        yval <- sumBrite[ang, range] / nCount[ang, range]
        xval <- (range -1) * rmax / num_rad
        sumx <- sumx + xval
        sumy <- sumy + yval
        sumx2 <- sumx2 + xval * xval
        sumxy <- sumxy + xval * yval
        sumn <- sumn + 1
      }
    }
    slope[ang] <- (sumn * sumxy - sumx * sumy) / (sumn * sumx2 - sumx * sumx)
    intercept[ang] <- (sumy - slope[ang] * sumx) / sumn
  }

  # compute average slope over all directions and scales
  sumn <- 0
  yval <- as.numeric()
  tempr <- as.numeric()
  for (k in 6:(2 * num_rad)) {
    if (radCount[k] > 0) {
      sumn <- sumn + 1
      yval[sumn] <- radius[k] / radCount[k]
      tempr[sumn] <- (k - 1) * rmax / (2 * num_rad)
    }
  }
  model_data <- data.frame(yval = yval, tempr = tempr)
  p <- lm(yval ~ tempr, data = model_data)
  averslope <- p$coefficients[2]
  averIC <- p$coefficients[1]

  fitln <- predict(p, newdata = model_data)
  slope[num_dir + 1] <- slope[1]
  intercept[num_dir + 1] <- intercept[1]

  # draw rose plot of slope and intercept
  ang <- seq(1, (num_dir + 1))

  return(abs(averslope)[[1]])
}
