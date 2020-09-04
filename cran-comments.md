---
title: "cran-comments.md"
author: "Annie C. Smith"
date: "September 3, 2020"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
## Test environments
* local Kubuntu 20.04 install, R 4.0.2
* local Windows 10, R 4.0.2
* R CMD CHECK Windows 10, R development (4.0.0)

## R CMD check results
There were no ERRORs or WARNINGS. 

There was 1 NOTE, "unable to verify current time," which was only produced using devtools::check() on Linux and appears to be due to a reference to an outdated website, not the R package.

This update fixes a bug caused by a bad is.na() call, and adds a vignette.

## Downstream dependencies
There are currently no downstream dependencies for this package.
