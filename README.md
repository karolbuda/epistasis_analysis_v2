# Epistasis Analysis

This repository represents the analysis hub for the epistatic analysis script. It contains both the bulk analysis script of the raw data, as well as the statistical analysis of the output files.

## Epistais Model Bulk

The epistasis_model_bulk.R code requires an Input folder with *.csv* files which specify the WT string, the mutated positions within the protein, then a list of genotype-phenotype pairs with WT-normalized and log10 transformed data in the Phenotype column. The code will generate subdirectories in the Output folder, all of which contain files representative of various analyses for each combinatorial landscape.

## Statistical Analysis

The statistical analysis is provided as an R markdown (Rmd) document which can be Knit into an html summarizing the statistical exploration on all data in the Output folder.

## Nonlinear Transformation Exploration

The nonlin_trans.Rmd file explores the effect of a non-linear transformation and how it can affect interpretability of epistatic data. It was used as a supplementary document in the manuscript.

For more information regarding the script and analyses contact Nobuhiko Tokuriki at tokuriki@msl.ubc.ca
