---
title: "NEWS.md"
author: "Annie C. Smith"
date: "September 1, 2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
## geodiv 0.1.0

Initial release of geodiv to CRAN.

## geodiv 0.1.1

Updated sfd.cpp to fix an installation issue. Changed to use 'std::floor()' instead of 'double floor(double).'
Decreased version requirement for parallel.

## geodiv 0.2.0

Added a vignette that shows more complex examples and demonstrates the correlation among metrics.
Fixed a new issue with NA autocorrelation images in the scl function.

## geodiv 1.0.0

Added an alternative function for applying focal metrics to create a texture image (focal_metrics), added the vignette with example uses and figures demonstrating relationships among metrics to the R package, fixed several minor bugs resulting in erroneous warnings in most metric functions (issues with checking input classes), and increased the speed of the 'std' and 'srw' functions.

Additionally, the 'texture_image' function no longer requires a projection as input, and will not reproject any rasters. Users should now reproject rasters on their own and then apply the function.

## geodiv 1.0.1

Fixed a minor bug in the 'texture_image' and 'pad_edges' functions. Added text to clarify the vignette objectives. Fixed a small error in the DESCRIPTION.

## geodiv 1.0.2

Fixed errors in function descriptions and examples. Simplified the vignette.

## geodiv 1.0.3

Added a longer timeout option before the download.file arguments in the vignette.

## geodiv 1.0.4

Added a tryCatch so that download.file in the vignette fails gracefully.
