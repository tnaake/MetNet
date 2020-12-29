# MetNet

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
![R-CMD-check-bioc](https://github.com/tnaake/MetNet/workflows/R-CMD-check-bioc/badge.svg)
[![codecov.io](http://codecov.io/github/tnaake/MetNet/coverage.svg?branch=master)](http://codecov.io/github/tnaake/MetNet?branch=master)
[![license](http://img.shields.io/badge/license-GPL%20%28%3E=%203%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
[![BioC checks](http://bioconductor.org/checkResults/devel/bioc-LATEST/MetNet/malbec1-checksrc.html)](http://bioconductor.org/checkResults/devel/bioc-LATEST/MetNet/malbec1-checksrc.html)

Inferring metabolic networks from untargeted high-resolution mass spectrometry data

## Description
Please visit the Bioconductor page of 
[MetNet](https://bioconductor.org/packages/MetNet) for further information. 

## Contact 

You are welcome to 

 * write a mail to <thomasnaake@googlemail.com> 
 * submit suggestions and issues: <https://github.com/tnaake/MetNet/issues>
 * send a pull request: <https://github.com/tnaake/MetNet/issues> 

## Install
To install MetNet, please use the stable version available via Bioconductor. 
To install, enter 

```r 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MetNet")
``` 

to your console. The installation via BiocManager requires R version 3.6. 


If you would like to install the development version of MetNet, you will first
have to install [devtools](http://cran.r-project.org/web/packages/devtools/index.html) package: 

```r
install.packages("devtools")
library("devtools")
install_github("tnaake/MetNet")
```


