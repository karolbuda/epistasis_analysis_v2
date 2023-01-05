# Epistasis Analysis

This repository represents the analysis hub for the epistatic analysis script. It contains both the bulk analysis script of the raw data, as well as the statistical analysis of the output files.

## Epistasis Model Bulk

The epistasis_model_bulk.R code requires an Input folder with *.csv* files which specify the WT string, the mutated positions within the protein, then a list of genotype-phenotype pairs with WT-normalized and log10 transformed data in the Phenotype column. The code will generate subdirectories in the Output folder, all of which contain files representative of various analyses for each combinatorial landscape.

## Manuscript output

The statistical analysis is provided as an R markdown (Rmd) document with a knitted html file that summarizes all the statistical code on the data in the Output folder. All of the data in the manuscript can be found in this document.

## Statistical Analysis

The statistical analysis is provided as an R markdown (Rmd) document with a knitted html file that summarizes all the statistical code on the data in the Output folder. This document has all exploratory statistical analysis performed in this project - not all of this data is present in the publication.

## Nonlinear Transformation Exploration

The nonlin_trans.Rmd file explores the effect of a non-linear transformation and how it can affect interpretability of epistatic data. This data was not explicitly included in the manuscript.

## Structural Exploration

The distance_analysis.R script measures the correlation between distances of alpha-carbon atoms in the PDB files, and the associate pairwise epistasis between those positions. It was used to generate the extended data in the manuscript. This data was not explicitly included in the manuscript.

---

For more information regarding the script and analyses contact Nobuhiko Tokuriki at tokuriki@msl.ubc.ca
