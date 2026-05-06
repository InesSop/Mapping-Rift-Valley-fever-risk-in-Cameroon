# Script: 03_combined_models.R
# Purpose: build dry- and wet-season combined RVF risk maps by averaging
#          expert-based and literature-based WLC-MCDA outputs
# Inputs: expert- and literature-based seasonal risk rasters
# Outputs: combined dry- and wet-season risk rasters in outputs/

library(terra)
library(here)

# -------------------------------------------------------------------
# 1. Define projection and paths
# -------------------------------------------------------------------

utm33n <- "EPSG:32633"

output_dir <- here("outputs")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -------------------------------------------------------------------
# 2. DRY SEASON - combined model
# -------------------------------------------------------------------

message("Building dry-season combined risk map...")

# Load input rasters
risk_lit_dry <- rast(here(output_dir, "risk_index_RVF_lit_dry.tif"))
risk_expert_dry <- rast(here(output_dir, "risk_index_RVF_expert_dry.tif"))

# Check spatial alignment
if (!compareGeom(risk_lit_dry, risk_expert_dry)) {
  stop("Dry-season literature and expert rasters are not aligned.")
}

message("Dry-season literature and expert rasters are correctly aligned.")

# Build combined model (equal-weight average: 50/50)
combined_risk_dry <- (0.5 * risk_lit_dry) + (0.5 * risk_expert_dry)

# Export final raster
writeRaster(
  combined_risk_dry,
  here(output_dir, "Combined_risk_index_RVF_dry.tif"),
  overwrite = TRUE
)

message("Dry-season combined risk raster exported successfully.")

# Optional quick visualization
plot(combined_risk_dry, main = "Combined RVF risk index - dry season")

# -------------------------------------------------------------------
# 3. WET SEASON - combined model
# -------------------------------------------------------------------

message("Building wet-season combined risk map...")

# Load input rasters
risk_lit_wet <- rast(here(output_dir, "risk_index_RVF_lit_wet.tif"))
risk_expert_wet <- rast(here(output_dir, "risk_index_RVF_expert_wet.tif"))

# Check spatial alignment
if (!compareGeom(risk_lit_wet, risk_expert_wet)) {
  stop("Wet-season literature and expert rasters are not aligned.")
}

message("Wet-season literature and expert rasters are correctly aligned.")

# Build combined model (equal-weight average: 50/50)
combined_risk_wet <- (0.5 * risk_lit_wet) + (0.5 * risk_expert_wet)

# Export final raster
writeRaster(
  combined_risk_wet,
  here(output_dir, "Combined_risk_index_RVF_wet.tif"),
  overwrite = TRUE
)

message("Wet-season combined risk raster exported successfully.")

# Optional quick visualization
plot(combined_risk_wet, main = "Combined RVF risk index - wet season")
