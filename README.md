![alt text](https://github.com/mtmcgowan/GGcall/blob/master/GGcall_oval.png "GGcall Logo")

# GGcall: Gaussian Genotype Calling

An R package for calling polyploid SNP chip data using Gaussian mixture models.

Author: Matthew McGowan

## Overview
GGcall is an R pipeline for clustering complex signals generated by the Illumina Infinium assays when used in polyploid species. The only input requirement is an R/Theta export generated from a Genome Studio Project. After clustering, users have the option to interactively filter, plot, and call genotypes. Final calls are output in a simple numerical format for further imputation and use in genome-wide-associaton-studies or genomic-selection analysis.

## Requirements
R/Rstudio (version 3.4 or higher)

A spreadsheet editor (Excel, OpenOffice-Calc, Numbers, etc.)

## Installation

The easiest way to get GGcall is to install the package using devtools.
If you don't have devtools already installed:
```{r, eval = FALSE}
install.packages("devtools")
```
Load devtools and then install GGcall from Github:
```{r, eval = FALSE}
library(devtools)
devtools::install_github("mtmcgowan/GGcall")
```

## Usage
The easiest way to learn to use GGcall is to go through the tutorial provided here:
[Guided Tutorial](https://github.com/mtmcgowan/GGcall/wiki/Beginner-Tutorial)

Depending on your population size and marker number, it may be useful to use a high-performance cluster to speed up the clustering step. 

## Getting help

The best way to get help is to submit an issue on Github. Before submitting, be sure to clearly describe the problem and provide a reproducible example.

## More Information

An earlier version of GGcall was presented at the 2018 Plant and Animal Genome Conference. Here is a poster that provides more information for how it works:

[](https://github.com/mtmcgowan/GGcall/PAG2018_McGowanPoster.png)
