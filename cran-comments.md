---
title: "cran-comments.md"
author: "Annie C. Smith"
date: "March 29, 2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
## Test environments
* local Kubuntu 20.04 install, R 4.0.4
* local Windows 10, R 4.0.2
* local Windows 10, R 3.6.1

## R CMD check results
There were no ERRORS or WARNINGS.

There were 2 NOTEs: 

"Suggests orphaned package: 'ggmap'." This package is only used for a figure in the vignette and does not influence anything in the geodiv package functions.

"Uses the superseded package: snow." The tests failed without including 'snow' as an import, although all 'snow' functions are actually accessed through the package 'parallel.'

This update fixes minor bugs resulting in warnings, adds a new focal windowing function, updates the texture_image function, speeds up the 'std' and 'srw' functions, and adds a vignette.

## Downstream dependencies
There are currently no downstream dependencies for this package.
