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
R version 2024.12.1-563

Main packages:
- terra
- sf
- dplyr
- tidyr
- ggplot2
- tmap
- pROC
- writexl
- openxlsx
- forcats
- viridis

## Reproducibility
Run scripts in the following order:
1. `01_wlc_experts.R`
2. `02_wlc_literature.R`
3. `03_combined_models.R`
4. `04_sensitivity_analysis.R`
5. `05_validation_and_figures.R`
6. `06_exposure_and_risk_drivers.R`

06_exposure_and_risk_drivers.R provides two complementary analyses based on the combined seasonal models:
- estimation of the number of domestic ruminants located in high-risk divisions and in high + medium-risk divisions
- identification of the three dominant risk factors per division and mapping of the dominant factor
- 
## Data availability
Due to size and/or ownership restrictions, the full raster datasets are not deposited in this repository. Processed inputs or metadata needed to reproduce the workflow are described in the manuscript and supplementary files.
