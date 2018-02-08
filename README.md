# GGcall
An R package for calling polyploid SNP chip data using Gaussian mixture models.

Author: Matthew McGowan

## Overview
GGcall is an R pipeline for clustering complex signals generated by the Illumina Infinium assay when used in polyploid species. The only input requirement is an R/Theta export generated from a Genome Studio Project. After clustering, users have the option to interactively filter markers and call genotypes for biallelic markers and output them in a simple numerical format for further imputation and use in genome-wide-associaton-study or genomic-selection analysis.

## Requirements
R/Rstudio (version 3.4 or higher)

A spreadsheet editor (Excel, OpenOffice-Calc, Numbers, etc.)

## Installation

```{r, eval = FALSE}
# The easiest way to get GGcall is to install the package using devtools:

# If you don't have devtools already installed:
install.packages("devtools")

# Then load devtools and then install GGcall from Github:
library(devtools)
devtools::install_github("mtmcgowan/GGcall")
```
## Usage
The easiest way to learn to use GGcall is to go through the tutorial provided here:
[Guided Tutorial](https://github.com/mtmcgowan/GGcall/wiki/Guided-Tutorial)

Depending on your population size and marker number, it may be useful to use a high-performance cluster to speed up the clustering step. Here is an example of how this was done for a cluster using the SLURM job scheduler:

LINK HERE (Not implemented yet)

## Getting help

The best way to get help is to submit an issue on Github. Before submitting, be sure to clearly describe the problem and provide a reproducible example.
