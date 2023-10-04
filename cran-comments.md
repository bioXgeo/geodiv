---
title: "cran-comments.md"
author: "Annie C. Smith"
date: "October 2, 2023"
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

There was 1 NOTE: 

"Uses the superseded package: snow." We use 'parallel' for all parallel operations, but 'parallel' depends on 'snow for some functions. We now include 'snow' because a user reported that it showed up as a warning during installation.

There is 1 WARNING.

Two functions of "parallel" are superceded by "snow." The functions are specified (as snow::makeCluster) in the geodiv functions in which they are called. This warning does not impact the functionality 
of geodiv.

## Downstream dependencies
There are currently no downstream dependencies for this package.
