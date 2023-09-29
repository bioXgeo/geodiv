//sfd.cpp

#include <stdio.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

/* =================================================== */
  /* FRACT3D : Computation of Fractal Dimension for */
  /* digital surfaces using the triangular prism */
  /* surface area method. Keith Clarke 8-13-85 */
  /* Modified Annie C. Smith 5-22-2019 */
/* =================================================== */
#include <stdlib.h>
#include <math.h>
#define MAXCOL 2000000

//' Calculate the fractal dimension of a raster (C function).
//'
//' Calculates the 3D fractal dimension of a raster using the
//' triangular prism surface area method.
//'
//' @param mat A matrix.
//' @return A numeric value representing the fractal dimension of
//' the image.
//' @references Clarke, K.C., 1986. Computation of the fractal dimension of topographic
//' surfaces using the triangular prism surface area method. Computers & Geosciences,
//' 12(5), pp.713-722.
//' @examples
//'
//' # import raster image
//' data(normforest)
//' normforest <- terra::unwrap(normforest)
//'
//' # convert to matrix
//' mat <- matrix(normforest[], ncol = ncol(normforest), nrow = nrow(normforest))
//'
//' # calculate the fractal dimension
//' Sfd <- sfd_(mat)
//' @export
// [[Rcpp::export]]
double sfd_(NumericMatrix mat){
  /* find center area (if rectangular) */
  int begin_row = 1;
  int begin_col = 1;
  int rowmax =  mat.nrow();
  int colmax = mat.ncol();
  int msize;
  if (rowmax > colmax) {
    // rowmax is larger than colmax
    msize = rowmax;
  } else if (colmax > rowmax) {
    // colmax is larger than rowmax
    msize = colmax;
  } else {
    // rowmax and colmax are equal
    msize = rowmax;
  }
  double* area = new double[msize];
  double* resolution = new double[msize];
  double crossr;
  /* normalize matrix values */
  NumericMatrix newmat(rowmax, colmax);
  double minimum = 99999;
  for (int x=0; x<rowmax; x++) {
    for (int y=0; y<colmax; y++) {
       minimum = std::min(mat(x, y), minimum);
    }
  }
  double maximum = -99999;
  for (int x=0; x<rowmax; x++) {
    for (int y=0; y<colmax; y++) {
      maximum = std::max(mat(x, y), maximum);
    }
  }
  for (int x=0; x<rowmax; x++) {
    for (int y=0; y<colmax; y++) {
      newmat(x,y) = (mat(x,y) - minimum) / (maximum - minimum);
    }
  }
  int n = 1, size, slop;
  /* Select short side of array */
  /* if rowmax > colmax, size = colmax. otherwise, size = rowmax. */
  size = (rowmax > colmax) ? colmax : rowmax;
  int steps = 1;
  /* Find power of two which is less than short side */
  do {
    steps++;
    n *= 2;
  } while (n < size);
  n /= 2;
  steps--;
  /* Calculate begin and end rows & cols for processing */
  slop = std::floor((rowmax - n) / 2.);
  begin_row = slop + 1;
  int end_row = n + slop + 1;
  slop = std::floor((colmax - n) / 2.);
  begin_col = slop + 1;
  int end_col = n + slop + 1;
  /* compute the fractal dimension for each step out from center */
  int row, col, step = 1;
  double side, diag;
  double a, b, c, d, e, w, x, y, z, o, p, q, r, sa, sb, sc, sd, aa, ab, ac, ad, surface_area;
  /* Repeat for area sequence 1,4,16,64,255 etc. */
  int time;
  for (time=1; time<=steps; time++) {
    surface_area = 0.0;
    /* Set length of sides of triangles */
    side = (double) step;
    diag = (double) step * sqrt (2.0) / 2.0;
    /* Process whole array at this size */
    for (row = begin_row; row < end_row; row += step) {
      for (col = begin_col; col < end_col; col += step) {
        a = newmat(row,col);
        b = newmat(row,col + step);
        c = newmat(row + step,col + step);
        d = newmat(row + step,col);
        /* e is the center point of four pixel values */
        e = 0.25 * (a + b + c + d);
        /* w,x,y,z are external sides of the square */
        w = sqrt ((a - b) * (a - b) + side * side);
        x = sqrt ((b - e) * (b - c) + side * side);
        y = sqrt ((c - d) * (c - d) + side * side);
        z = sqrt ((a - d) * (a - d) + side * side);
        /* o,p,q,r are internal sides of triangles */
        o = sqrt ((a - e) * (a - e) + diag * diag);
        p = sqrt ((b - e) * (b - e) + diag * diag);
        q = sqrt ((c - e) * (c - e) + diag * diag);
        r = sqrt ((d - e) * (d - e) + diag * diag);
        /* Compute values for Heron's formula */
        sa = 0.5 * (w + p + o);
        sb = 0.5 * (x + p + q);
        sc = 0.5 * (y + q + r);
        sd = 0.5 * (z + o + r);
        /* Solve area~ from Heron's formula */
        aa = sqrt (fabs(sa * (sa - w) * (sa - p) * (sa - o)));
        ab = sqrt (fabs(sb * (sb - x) * (sb - p) * (sb - q)));
        ac = sqrt (fabs(sc * (sc - y) * (sc - q) * (sc - r)));
        ad = sqrt (fabs(sd * (sd - z) * (sd - o) * (sd - r)));
        /* Add to total surface area */
        surface_area += aa + ab + ac + ad;
      }
    }
    /* Save area and resolution, increment step size */
    area[time] = surface_area;
    resolution[time] = step * step;
    step *= 2;
  }
  /* find final fractal dimension */
  double resavg = 0.0, areaavg = 0.0, cross = 0.0, sumres = 0.0,
    sumarea = 0.0, fd, beta;
  /* Do log transform and compute means */
  for (n = 1; n <= steps; n++) {
    resolution[n] = log(resolution[n]);
    area[n] = log(area[n]);
    resavg += resolution[n];
    areaavg += area[n];
  }
  if (steps < 3) {
    errno = 5;
    perror("Error: Matrix is too small");
  }
  resavg /= (double) (steps);
  areaavg /= (double) (steps);
  /* Compute sums of squares */
  for (n = 1; n <= steps; n++) {
    cross += ((resolution[n] - resavg) * (area[n] - areaavg));
    sumres += ((resolution[n] - resavg) * (resolution[n] - resavg));
    sumarea += ((area[n] - areaavg) * (area[n] - areaavg));
  }
  if (sumres == 0.0) {
    sumres = 1;
  }
  if (sumarea == 0.0) {
    sumarea = 1;
  }
  /* Compute correlation coefficient and fractal dimension */
  crossr = cross / sqrt(sumres * sumarea);
  beta = crossr * sqrt(sumarea) / sqrt(sumres);
  fd = 2.0 - beta;
  delete[] area;
  delete[] resolution;
  return fd;
}
