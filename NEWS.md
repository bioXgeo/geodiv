---
title: "NEWS.md"
author: "Annie C. Smith"
date: "September 3, 2020"
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
