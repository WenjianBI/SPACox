# SPACox
A semiparametric empirical SPA approach based on a Cox regression model fitting (SPACox), that is scalable for a genome-wide single-variant survival analysis in large cohorts

### How to install and load this package

```{r}      
library(devtools)  # author version: 2.1.0
install_github("WenjianBi/SPACox")
library(SPACox)
?SPACox  # manual of SPACox package
```
Current version is 0.1.1. For older version and version update information, plesase refer to OldVersions/

Please do not hesitate to contact me (wenjianb@umich.edu) if you meet any problem. Suggestions or comments are also welcome.

### Reference

Wenjian Bi, Lars G. Fritsche, Bhramar Mukherjee, Sehee Kim, Seunggeun Lee, A Fast and Accurate Method for Genome-Wide Time-to-Event Data Analysis and Its Application to UK Biobank. American Journal of Human Genetics (2020), https://doi.org/10.1016/j.ajhg.2020.06.003

## coxphf
This directory contains some discussions about R package coxphf, an R package to conduct Firth's correction for Cox regression.
