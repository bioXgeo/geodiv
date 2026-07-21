---
title: "cran-comments.md"
author: "Kyla M. Dahlin"
date: "July 21, 2026"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
## Test environments
* local Windows 11 install, R 4.6.0
* win-builder (devel and release)
* R-hub v2 (Linux, Windows, macOS)
* GitHub Actions (ubuntu-latest, macOS-latest, windows-latest): R release, devel, 
oldrel-1

## R CMD check results
There were no ERRORS or WARNINGS, and one note re: ghostscript size reduction.

## Downstream dependencies
This package has one reverse dependency, `prior3D`. I have 
checked that it passes R CMD check with this updated version installed 
and confirmed no new problems.

## Comments
This is an update of contact information for authors and a change in the 
maintainer in response to a flag from CRAN about the maintainer's email address 
no longer working.
