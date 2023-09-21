# Epistasis Analysis

This repository represents the analysis hub for the epistatic analysis script. It contains both the bulk analysis script of the raw data, as well as the statistical analysis of the output files.

## Instructions

### Library requirements

All scripts were run using R version 4.1.2

All scripts were run using the following packages:

tidyverse 2.0.0
✔ dplyr     1.1.3     ✔ readr     2.1.4
✔ forcats   1.0.0     ✔ stringr   1.5.0
✔ ggplot2   3.4.3     ✔ tibble    3.2.1
✔ lubridate 1.9.2     ✔ tidyr     1.3.0
✔ purrr     1.0.2 

✔ splines 4.1.2     
✔ mgcv  1.8.42
✔ knitr 1.43
✔ minpack.lm 1.2.3
✔ ggrepel 0.9.3
✔ ggpubr 0.6.0
✔ readxl 1.4.3

All scripts must be run by opening the R project file (epistasis_analysis_v2.Rproj).

### Script order

First it is necessary to run Pre non-linear transformation scripts as the non-linear transformation depends on these outputs. Next, non-lin-trans.Rmd must be run to export the four-parameter transformsed data. Finally, non-lin_epistasis_model_bulk.R must be run to generate the all outputs for the analyses, followed by non-lin_statistical_epistasis_analysis_rev.Rmd and Supplementary File 2.Rmd.

---

## File overview

### Project file

The file epistasis_analysis_v2.Rproj is the project file for this repository.

### Epistasis model bulk

The non-lin_epistasis_model_bulk.R code requires an Input directory with *.csv* files which specify the WT string, the mutated positions within the protein, then a list of genotype-phenotype pairs with WT-normalized and log10 transformed data in the Phenotype column. The code will generate subdirectories in the Output directory, all of which contain files representative of various analyses for each combinatorial landscape. The approximate run time for all input files is ~ 30 min.

### Manuscript output

The statistical analysis is provided as an R markdown (Rmd) document called non-linwith a knitted html file that summarizes all the statistical code on the data in the Output directory. All of the data in the manuscript can be found in this document.

### Statistical analysis

The statistical analysis is provided as an R markdown (Rmd) document with a knitted html file that summarizes all the statistical code on the data in the Output directory. This document has all exploratory statistical analysis performed in this project - not all of this data is present in the publication.

### Non-linear transformation (Supplementary file 2)

After the first round of revisions for the manuscript, all manuscript data were analysed using a non-linear transformation, where both a monotonic spline fit and four parameter fit were evaluated. These data are found in the Pre Non-linear Transformation directory. The analysis of the data can be found in the Supplementary File 2, while the export of the data for the analysis can be found in non-lin_trans.Rmd.

### Pre non-linear transformation

Prior to the first round of revisions for the manuscript, all manuscript data were analysed without using a non-linear transformation. These data are found in the Pre Non-linear Transformation directory.

### Supplementary

This directory contains the 3 supplementary data files. Pre non-linear transformation data can be found within the "old" directory. The pos-neg_combos.Rmd explores the proportion of positive and negative SMEs and EEs for all positions and combinations marked as positive-negative (see Supplementary Fig. 1 and 2 for details).

---

For more information regarding the script and analyses contact Nobuhiko Tokuriki at tokuriki@msl.ubc.ca
