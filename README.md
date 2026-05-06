# RVF risk mapping in Cameroon – WLC-MCDA workflow

This repository contains the R scripts used to reproduce the weighted linear combination (WLC-MCDA) modelling, sensitivity analysis and validation presented in the manuscript.

## Scope of the repository
The repository starts from preprocessed and standardized risk-factor layers.  
Scripts related to the upstream geoprocessing and preparation of the environmental and livestock predictor layers are not included.

## Included workflow
- weighted linear combination for dry season
- weighted linear combination for wet season
- sensitivity analysis
- model validation
- figure and table generation

## Inputs
The scripts require preprocessed raster layers prepared on a common grid and projection, as described in the manuscript and supplementary materials.

## Requirements
R version X.X.X

Main packages:
- terra
- sf
- tmap
- ggplot2
- pROC
- dplyr
- writexl

## Reproducibility
Run scripts in the following order:
1. `01_risk_model_dry.R`
2. `02_risk_model_wet.R`
3. `03_sensitivity_analysis.R`
4. `04_validation_auc.R`

## Data availability
Due to size and/or ownership restrictions, the full raster datasets are not deposited in this repository. Processed inputs or metadata needed to reproduce the workflow are described in the manuscript and supplementary files.
