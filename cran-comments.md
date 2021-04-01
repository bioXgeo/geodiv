---
title: "cran-comments.md"
author: "Annie C. Smith"
date: "April 1, 2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
## Test environments
* local Kubuntu 20.04 install, R 4.0.4
* local Windows 10, R 4.0.2
* local Windows 10, R 3.6.1
* win-builder, R-devel

## R CMD check results
There were no ERRORS or WARNINGS.

There were 2 NOTEs: 

"Suggests orphaned package: 'ggmap'." This package is only used for a figure in the vignette and does not influence anything in the geodiv package functions.

"Uses the superseded package: snow." The tests failed without including 'snow' as an import, although all 'snow' functions are actually accessed through the package 'parallel.'

This update fixes minor a minor bug in the texture_image and pad_edges functions, adds text to the vignette, and fixes an error in the DESCRIPTION.

## Downstream dependencies
There are currently no downstream dependencies for this package.
