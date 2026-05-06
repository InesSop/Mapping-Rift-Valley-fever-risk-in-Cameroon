# Script: 01_wlc_experts.R
# Purpose: build dry- and wet-season RVF risk maps using WLC-MCDA with risk factors and weights derived from experts elicitation
# Inputs: standardized raster layers produced during data processing
# Outputs: expert-driven dry- and wet-season risk rasters in outputs/

library(terra)
library(here)

# -------------------------------------------------------------------
# 1. Define projection and paths
# -------------------------------------------------------------------

utm33n <- "EPSG:32633"

input_dir_experts_dry <- here("inputs", "standardized_rasters", "experts", "dry")
input_dir_experts_wet <- here("inputs", "standardized_rasters", "experts", "wet")
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
# 3. DRY SEASON - expert-based WLC model
# -------------------------------------------------------------------

message("Building dry-season expert-based risk map...")

# Step 1: Load reference raster
ref_dry <- rast(here(input_dir_experts_dry, "rainfall_dry_std.tif"))

# Assign CRS only if missing or incorrectly defined in the source file
crs(ref_dry) <- utm33n

# Step 2: Load and align all rasters to the same reference grid
rainfall <- ref_dry

protected_areas <- align_to_ref(
  rast(here(input_dir_risk_factors, "wild_ruminant_density", "standardised_distance_to_wildlife_R.tif")),
  ref_dry
)

cattle_density <- align_to_ref(
  rast(here(input_dir_risk_factors, "livestock_density", "cattle_suitability_0_1_UTM33N_1km.tif")),
  ref_dry
)

small_ruminant_density <- align_to_ref(
  rast(here(input_dir_risk_factors, "livestock_density", "small_ruminant_density_std_UTM33N_1km.tif")),
  ref_dry
)

water_bodies <- align_to_ref(
  rast(here(input_dir_experts_dry, "water_bodies_perm_std.tif")),
  ref_dry
)

irrigation <- align_to_ref(
  rast(here(input_dir_experts_dry, "irrigation_std.tif")),
  ref_dry
)

Aedes_density <- align_to_ref(
  rast(here(input_dir_experts_dry, "mosquito_score_UTM33N_1km_cameroon_F32.tif")),
  ref_dry
)

livestock_trade <- align_to_ref(
  rast(here(input_dir_experts_dry, "distance_livestock_markets_standardized.tif")),
  ref_dry
)

# Step 3: Check alignment
rasters_dry <- list(
  rainfall = rainfall,
  protected_areas = protected_areas,
  cattle_density = cattle_density,
  small_ruminant_density = small_ruminant_density,
  water_bodies = water_bodies,
  irrigation = irrigation,
  Aedes_density = Aedes_density,
  livestock_trade = livestock_trade
)

if (!check_alignment(ref_dry, rasters_dry)) {
  stop("Dry-season rasters are not correctly aligned.")
}
message("Dry-season rasters are correctly aligned.")

# Step 4: Define expert-derived weights
weights_dry <- c(
  rainfall = 0.066666667,
  protected_areas = 0.170416667,
  cattle_density = 0.133333333,
  small_ruminant_density = 0.175833300,
  water_bodies = 0.113750000,
  irrigation = 0.062083333,
  Aedes_density = 0.148750000,
  livestock_trade = 0.129166667
)

if (abs(sum(weights_dry) - 1) > 1e-6) {
  stop("Dry-season weights do not sum to 1.")
}

# Step 5: Compute weighted risk index
risk_index_expert_dry <- rainfall               * weights_dry["rainfall"] +
                         protected_areas        * weights_dry["protected_areas"] +
                         cattle_density         * weights_dry["cattle_density"] +
                         small_ruminant_density * weights_dry["small_ruminant_density"] +
                         water_bodies           * weights_dry["water_bodies"] +
                         irrigation             * weights_dry["irrigation"] +
                         Aedes_density          * weights_dry["Aedes_density"] +
                         livestock_trade        * weights_dry["livestock_trade"]

# Step 6: Export final raster
writeRaster(
  risk_index_expert_dry,
  here(output_dir, "risk_index_RVF_expert_dry.tif"),
  overwrite = TRUE
)

message("Dry-season expert-based risk raster exported successfully.")

# Optional quick visualization
plot(risk_index_expert_dry, main = "RVF risk index - dry season (expert weights)")

# -------------------------------------------------------------------
# 4. WET SEASON - expert-based WLC model
# -------------------------------------------------------------------

message("Building wet-season expert-based risk map...")

# Step 1: Load reference raster
ref_wet <- rast(here(input_dir_experts_wet, "rainfall_wet_std.tif"))

# Assign CRS only if missing or incorrectly defined in the source file
crs(ref_wet) <- utm33n

# Step 2: Load and align all rasters to the same reference grid
rainfall <- ref_wet

protected_areas <- align_to_ref(
  rast(here(input_dir_experts_wet, "protected_areas_dist_std.tif")),
  ref_wet
)

cattle_density <- align_to_ref(
  rast(here(input_dir_experts_wet, "cattle_density_std.tif")),
  ref_wet
)

small_ruminant_density <- align_to_ref(
  rast(here(input_dir_risk_factors, "livestock_density", "small_ruminant_density_std_UTM33N_1km.tif")),
  ref_wet
)

water_bodies <- align_to_ref(
  rast(here(input_dir_experts_wet, "water_bodies_perm_temp_std.tif")),
  ref_wet
)

irrigation <- align_to_ref(
  rast(here(input_dir_experts_wet, "irrigation_std.tif")),
  ref_wet
)

Aedes_density <- align_to_ref(
  rast(here(input_dir_experts_dry, "mosquito_score_UTM33N_1km_cameroon_F32.tif")),
  ref_wet
)

livestock_trade <- align_to_ref(
  rast(here(input_dir_experts_dry, "distance_livestock_markets_standardized.tif")),
  ref_wet
)

# Step 3: Check alignment
rasters_wet <- list(
  rainfall = rainfall,
  protected_areas = protected_areas,
  cattle_density = cattle_density,
  small_ruminant_density = small_ruminant_density,
  water_bodies = water_bodies,
  irrigation = irrigation,
  Aedes_density = Aedes_density,
  livestock_trade = livestock_trade
)

if (!check_alignment(ref_wet, rasters_wet)) {
  stop("Wet-season rasters are not correctly aligned.")
}
message("Wet-season rasters are correctly aligned.")

# Step 4: Define expert-derived weights
weights_wet <- c(
  rainfall = 0.15771,
  protected_areas = 0.13042,
  cattle_density = 0.09750,
  small_ruminant_density = 0.14406,
  water_bodies = 0.09406,
  irrigation = 0.05219,
  Aedes_density = 0.19771,
  livestock_trade = 0.12635
)

if (abs(sum(weights_wet) - 1) > 1e-6) {
  stop("Wet-season weights do not sum to 1.")
}

# Step 5: Compute weighted risk index
risk_index_expert_wet <- rainfall               * weights_wet["rainfall"] +
                         protected_areas        * weights_wet["protected_areas"] +
                         cattle_density         * weights_wet["cattle_density"] +
                         small_ruminant_density * weights_wet["small_ruminant_density"] +
                         water_bodies           * weights_wet["water_bodies"] +
                         irrigation             * weights_wet["irrigation"] +
                         Aedes_density          * weights_wet["Aedes_density"] +
                         livestock_trade        * weights_wet["livestock_trade"]

# Step 6: Export final raster
writeRaster(
  risk_index_expert_wet,
  here(output_dir, "risk_index_RVF_expert_wet.tif"),
  overwrite = TRUE
)

message("Wet-season expert-based risk raster exported successfully.")

# Optional quick visualization
plot(risk_index_expert_wet, main = "RVF risk index - wet season (expert weights)")
