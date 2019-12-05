---
title: "cran-comments.md"
author: "Annie C. Smith"
date: "November 22, 2019"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
## Test environments
* local Ubuntu 18.04.3 install, R 3.6.1
* local OS X 10.14.6 install, R 3.6.1
* local Windows 10, R 3.6.1

## R CMD check results
There were no ERRORs, or WARNINGS.

With the check_rhub() function, there was one NOTE (below).

Found the following files/directories:
    'examples_i386' 'examples_x64' 'geodiv-Ex_i386.Rout'
    'geodiv-Ex_x64.Rout' 'tests_i386' 'tests_x64'
    
These files are not included in the R package, and this note
does not appear with the devtools::check() function in any environment.

This is the initial submission of this package.

## Downstream dependencies
There are currently no downstream dependencies for this package.
