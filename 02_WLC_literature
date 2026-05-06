# Script: 02_wlc_literature.R
# Purpose: build dry- and wet-season RVF risk maps using WLC-MCDA with literature-derived risk factors and weights
# Inputs: standardized raster layers produced during data processing
# Outputs: literature-driven dry- and wet-season risk rasters in outputs/

library(terra)
library(here)

# -------------------------------------------------------------------
# 1. Define projection and paths
# -------------------------------------------------------------------

utm33n <- "EPSG:32633"

input_dir_lit_dry <- here("inputs", "standardized_rasters", "literature", "dry")
input_dir_lit_wet <- here("inputs", "standardized_rasters", "literature", "wet")
input_dir_risk_factors <- here("inputs", "risk_factors")
output_dir <- here("outputs")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -------------------------------------------------------------------
# 2. Define helper functions
# -------------------------------------------------------------------

align_to_ref <- function(r, ref) {
  project(r, ref, method = "bilinear")
}

check_alignment <- function(ref, raster_list) {
  checks <- sapply(raster_list, function(x) compareGeom(ref, x))
  all(checks)
}

# -------------------------------------------------------------------
# 3. DRY SEASON - literature-based WLC model
# -------------------------------------------------------------------

message("Building dry-season literature-based risk map...")

# Step 1: Load reference raster
ref_dry <- rast(here(input_dir_lit_dry, "rainfall_dry_std.tif"))

# Assign CRS only if missing or incorrectly defined in the source file
crs(ref_dry) <- utm33n

# Step 2: Load and align all rasters to the same reference grid
rainfall <- ref_dry

protected_areas <- align_to_ref(
  rast(here(input_dir_lit_dry, "protected_areas_dist_std.tif")),
  ref_dry
)

cattle_density <- align_to_ref(
  rast(here(input_dir_lit_dry, "cattle_density_std.tif")),
  ref_dry
)

sheep_density <- align_to_ref(
  rast(here(input_dir_lit_dry, "sheep_density_std.tif")),
  ref_dry
)

water_bodies <- align_to_ref(
  rast(here(input_dir_lit_dry, "water_bodies_perm_std.tif")),
  ref_dry
)

irrigation <- align_to_ref(
  rast(here(input_dir_lit_dry, "irrigation_std.tif")),
  ref_dry
)

NDVI_dry <- align_to_ref(
  rast(here(input_dir_risk_factors, "ndvi", "Standardized_NDVI_dry_UTM33N_1km_filled.tif")),
  ref_dry
)

rivers_streams <- align_to_ref(
  rast(here(input_dir_lit_dry, "river_streams_dry_std.tif")),
  ref_dry
)

elevation <- align_to_ref(
  rast(here(input_dir_lit_dry, "elevation_std.tif")),
  ref_dry
)

Aedes_density <- align_to_ref(
  rast(here(input_dir_lit_dry, "mosquito_score_UTM33N_1km_cameroon_F32.tif")),
  ref_dry
)

# Step 3: Check alignment
rasters_dry <- list(
  rainfall = rainfall,
  protected_areas = protected_areas,
  cattle_density = cattle_density,
  sheep_density = sheep_density,
  water_bodies = water_bodies,
  irrigation = irrigation,
  NDVI_dry = NDVI_dry,
  rivers_streams = rivers_streams,
  elevation = elevation,
  Aedes_density = Aedes_density
)

if (!check_alignment(ref_dry, rasters_dry)) {
  stop("Dry-season literature rasters are not correctly aligned.")
}
message("Dry-season literature rasters are correctly aligned.")

# Step 4: Define literature-derived weights
weights_dry <- c(
  rainfall = 0.152173914,
  protected_areas = 0.065217392,
  cattle_density = 0.065217392,
  sheep_density = 0.065217392,
  water_bodies = 0.108695653,
  irrigation = 0.065217392,
  NDVI_dry = 0.152173914,
  rivers_streams = 0.108695653,
  elevation = 0.152173914,
  Aedes_density = 0.065217392
)

if (abs(sum(weights_dry) - 1) > 1e-6) {
  stop("Dry-season literature weights do not sum to 1.")
}

# Step 5: Compute weighted risk index
risk_index_lit_dry <- rainfall        * weights_dry["rainfall"] +
                      protected_areas * weights_dry["protected_areas"] +
                      cattle_density  * weights_dry["cattle_density"] +
                      sheep_density   * weights_dry["sheep_density"] +
                      water_bodies    * weights_dry["water_bodies"] +
                      irrigation      * weights_dry["irrigation"] +
                      NDVI_dry        * weights_dry["NDVI_dry"] +
                      rivers_streams  * weights_dry["rivers_streams"] +
                      elevation       * weights_dry["elevation"] +
                      Aedes_density   * weights_dry["Aedes_density"]

# Step 6: Export final raster
writeRaster(
  risk_index_lit_dry,
  here(output_dir, "risk_index_RVF_lit_dry.tif"),
  overwrite = TRUE
)

message("Dry-season literature-based risk raster exported successfully.")

# Optional quick visualization
plot(risk_index_lit_dry, main = "RVF risk index - dry season (literature weights)")

# -------------------------------------------------------------------
# 4. WET SEASON - literature-based WLC model
# -------------------------------------------------------------------

message("Building wet-season literature-based risk map...")

# Step 1: Load reference raster
ref_wet <- rast(here(input_dir_lit_wet, "rainfall_wet_std.tif"))

# Assign CRS only if missing or incorrectly defined in the source file
crs(ref_wet) <- utm33n

# Step 2: Load and align all rasters to the same reference grid
rainfall <- ref_wet

protected_areas <- align_to_ref(
  rast(here(input_dir_lit_wet, "protected_areas_dist_std.tif")),
  ref_wet
)

cattle_density <- align_to_ref(
  rast(here(input_dir_lit_wet, "cattle_density_std.tif")),
  ref_wet
)

sheep_density <- align_to_ref(
  rast(here(input_dir_lit_wet, "sheep_density_std.tif")),
  ref_wet
)

water_bodies <- align_to_ref(
  rast(here(input_dir_lit_wet, "water_bodies_perm_temp_std.tif")),
  ref_wet
)

irrigation <- align_to_ref(
  rast(here(input_dir_lit_wet, "irrigation_std.tif")),
  ref_wet
)

NDVI_wet <- align_to_ref(
  rast(here(input_dir_risk_factors, "ndvi", "Standardized_NDVI_wet_UTM33N_1km_filled.tif")),
  ref_wet
)

rivers_streams <- align_to_ref(
  rast(here(input_dir_lit_wet, "river_streams_wet_std.tif")),
  ref_wet
)

elevation <- align_to_ref(
  rast(here(input_dir_lit_wet, "elevation_std.tif")),
  ref_wet
)

Aedes_density <- align_to_ref(
  rast(here(input_dir_lit_dry, "mosquito_score_UTM33N_1km_cameroon_F32.tif")),
  ref_wet
)

# Step 3: Check alignment
rasters_wet <- list(
  rainfall = rainfall,
  protected_areas = protected_areas,
  cattle_density = cattle_density,
  sheep_density = sheep_density,
  water_bodies = water_bodies,
  irrigation = irrigation,
  NDVI_wet = NDVI_wet,
  rivers_streams = rivers_streams,
  elevation = elevation,
  Aedes_density = Aedes_density
)

if (!check_alignment(ref_wet, rasters_wet)) {
  stop("Wet-season literature rasters are not correctly aligned.")
}
message("Wet-season literature rasters are correctly aligned.")

# Step 4: Define literature-derived weights
weights_wet <- c(
  rainfall = 0.152173914,
  protected_areas = 0.065217392,
  cattle_density = 0.065217392,
  sheep_density = 0.065217392,
  water_bodies = 0.108695653,
  irrigation = 0.065217392,
  NDVI_wet = 0.152173914,
  rivers_streams = 0.108695653,
  elevation = 0.152173914,
  Aedes_density = 0.065217392
)

if (abs(sum(weights_wet) - 1) > 1e-6) {
  stop("Wet-season literature weights do not sum to 1.")
}

# Step 5: Compute weighted risk index
risk_index_lit_wet <- rainfall        * weights_wet["rainfall"] +
                      protected_areas * weights_wet["protected_areas"] +
                      cattle_density  * weights_wet["cattle_density"] +
                      sheep_density   * weights_wet["sheep_density"] +
                      water_bodies    * weights_wet["water_bodies"] +
                      irrigation      * weights_wet["irrigation"] +
                      NDVI_wet        * weights_wet["NDVI_wet"] +
                      rivers_streams  * weights_wet["rivers_streams"] +
                      elevation       * weights_wet["elevation"] +
                      Aedes_density   * weights_wet["Aedes_density"]

# Step 6: Export final raster
writeRaster(
  risk_index_lit_wet,
  here(output_dir, "risk_index_RVF_lit_wet.tif"),
  overwrite = TRUE
)

message("Wet-season literature-based risk raster exported successfully.")

# Optional quick visualization
plot(risk_index_lit_wet, main = "RVF risk index - wet season (literature weights)")
