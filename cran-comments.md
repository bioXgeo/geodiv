---
title: "cran-comments.md"
author: "Annie C. Smith"
date: "March 19, 2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
## Test environments
* local Kubuntu 20.04 install, R 4.0.2
* local Windows 10, R 4.0.2
* local Windows 10, R 3.6.1

## R CMD check results
There were no ERRORs, 0 WARNINGS, or NOTES. 

This update fixes warnings caused by bad class checks in many of the functions, adds a new focal windowing function, removes reprojections from all functions, speeds up the 'std' and 'srw' functions, and adds a vignette to demonstrate correlations among included spatial heterogeneity metrics.

## Downstream dependencies
There are currently no downstream dependencies for this package.
