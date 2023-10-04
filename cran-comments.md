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
* local Ubuntu 20.04 install, R 4.3.1
* local Windows 10, R 4.0.2
* local Windows 10, R 3.6.1
* win-builder, R-devel

## R CMD check results
There were no ERRORS or WARNINGS.

There was 1 NOTE: 

"Uses the superseded package: snow." We use 'parallel' for all parallel operations, but 'parallel' depends on 'snow for some functions. We now include 'snow' because a user reported that it showed up as a warning during installation.


## Downstream dependencies
There are currently no downstream dependencies for this package.
