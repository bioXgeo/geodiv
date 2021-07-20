---
title: "cran-comments.md"
author: "Annie C. Smith"
date: "July 20, 2021"
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

"Uses the superseded package: snow." We use 'parallel' for all parallel operations, but 'parallel' depends on 'snow for some functions. We now include 'snow' because a user reported that it showed up as a warning during installation.

This update fixes an issue with a failed download in the vignette by adding a longer timeout argument.

## Downstream dependencies
There are currently no downstream dependencies for this package.
