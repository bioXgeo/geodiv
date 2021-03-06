
<!-- README.md is generated from README.Rmd. Please edit that file -->
geodiv
======

<!-- badges: start -->
<!-- badges: end -->
*geodiv* calculates gradient surface metrics in R. These metrics are applied to continuous spatial data (i.e., rasters or matrices) and represent spatial heterogeneity.

Publications
------------

*geodiv* is featured in the following publications:

Smith, AC, Dahlin KM, Record S, Costanza JK, Wilson AM, Zarnetske PL. 2021. The geodiv R package: tools for creating gradient surface metrics. Accepted in *Methods in Ecology and Evolution*.

Installation
------------

You can install the released version of *geodiv* from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("geodiv")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("bioXgeo/geodiv")
```

On Mac OS X, you may need to install the development tools here to get the package to install:

<https://cran.r-project.org/bin/macosx/tools/>

Example
-------

This is a basic example which shows you how to calculate several metrics over an entire image. *geodiv* may also be applied with moving windows over an entire image using the 'texture\_image' function.

``` r
library(geodiv)

# import example raster
data(normforest)

# apply metrics
sa(normforest) # average surface roughness
#> [1] 0.0442945
svk(normforest) # reduced valley depth
#> [1] 0.5449867
ssc(normforest) # mean summit curvature
#> [1] -0.02192238
```

Vignette Data
-------------

The intermediate data for the vignette is available on figshare: https://doi.org/10.6084/m9.figshare.12834896.v5
