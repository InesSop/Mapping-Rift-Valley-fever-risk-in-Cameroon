# Script: 06_exposure_and_risk_drivers.R
# Purpose:
#   1. Estimate the number of domestic ruminants located in high-risk divisions
#      and in high + medium-risk divisions for dry and wet seasons
#   2. Identify the three dominant risk factors per division in the combined model
#      and map the dominant factor for each season
#
# Inputs:
#   - standardized expert- and literature-based raster layers
#   - combined seasonal risk rasters
#   - division boundaries
#   - livestock density rasters (cattle, sheep, goat)
#
# Outputs:
#   - Excel tables of livestock exposure by risk class
#   - Excel tables of top 3 risk factors by division
#   - dominant factor maps for dry and wet seasons
#   - optional shapefiles with dominant factor attributes

library(here)
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(tmap)
library(writexl)

# =====================================================================
# 1. BACKGROUND
# =====================================================================

# This script performs two complementary analyses based on the combined RVF model:
#
# 1) Livestock exposure analysis
#    The script estimates the number of domestic ruminants located in:
#    - high-risk divisions only
#    - high + medium-risk divisions
#
#    Division-level livestock totals are approximated from mean livestock density
#    multiplied by division area.
#
# 2) Division-level risk driver analysis
#    The combined model was obtained as:
#      combined_risk = 0.5 * expert_model + 0.5 * literature_model
#
#    To identify which variables most strongly contribute to risk in each division,
#    the contribution of each harmonized factor was reconstructed as the weighted
#    sum of its expert-side and literature-side components.
#
#    Variable harmonization:
#    - The expert-based model uses 'small_ruminant_density'
#    - The literature-based model uses 'sheep_density'
#    - For the purpose of combined factor attribution, both were harmonized under
#      the common variable name 'small_ruminant_density'
#
#    Therefore, in the combined contribution maps:
#      contribution_small_ruminant_density =
#        0.5 * (expert weight for small_ruminant_density * expert raster) +
#        0.5 * (literature weight for sheep_density * literature raster)
#
#    Variables present in only one sub-model are retained with contribution from
#    that sub-model only:
#    - NDVI, rivers_streams, elevation: literature-only
#    - livestock_trade: expert-only
#
# Terminology:
# - The administrative unit NAME_2 is referred to here as "division"
#   (not "department").

# =====================================================================
# 2. PATHS AND SETTINGS
# =====================================================================

utm33n <- "EPSG:32633"

input_boundaries <- here("inputs", "boundaries", "gadm41_CMR_2.shp")

input_dir_lit_dry <- here("inputs", "standardized_rasters", "literature", "dry")
input_dir_lit_wet <- here("inputs", "standardized_rasters", "literature", "wet")
input_dir_exp_dry <- here("inputs", "standardized_rasters", "experts", "dry")
input_dir_exp_wet <- here("inputs", "standardized_rasters", "experts", "wet")
input_dir_risk_factors <- here("inputs", "risk_factors")

output_models_dir <- here("outputs")
output_dir <- here("outputs", "exposure_and_risk_drivers")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Final classification thresholds derived in 05_validation_and_figures.R
# These thresholds are used here to classify divisions into Low / Medium / High risk.
thresholds <- list(
  dry = list(
    low = 0.4040356,
    high = 0.4732773
  ),
  wet = list(
    low = 0.4047360,
    high = 0.4878685
  )
)

# Color palette for dominant factor maps
factor_colors <- c(
  rainfall = "#1f77b4",
  protected_areas = "#98df8a",
  cattle_density = "#C98A2E",
  small_ruminant_density = "#bcbd22",
  water_bodies = "#17becf",
  irrigation = "#e377c2",
  Aedes_density = "#F08080",
  livestock_trade = "#7A0177",
  NDVI = "#2ca02c",
  rivers_streams = "#9edae5",
  elevation = "#FDD0A2"
)

# =====================================================================
# 3. HELPER FUNCTIONS
# =====================================================================

align_to_ref <- function(r, ref) {
  project(r, ref, method = "bilinear")
}

check_alignment <- function(ref, raster_list) {
  checks <- sapply(raster_list, function(x) compareGeom(ref, x))
  all(checks)
}

extract_mean_to_divisions <- function(divisions_sf, rast_obj, var_name) {
  div_vect <- terra::vect(divisions_sf)
  values <- terra::extract(rast_obj, div_vect, fun = mean, na.rm = TRUE)[, 2]
  divisions_sf[[var_name]] <- values
  divisions_sf
}

classify_risk <- function(x, low_thresh, high_thresh) {
  cut(
    x,
    breaks = c(-Inf, low_thresh, high_thresh, Inf),
    labels = c("Low", "Medium", "High"),
    right = TRUE
  )
}

get_top3_table <- function(mean_contrib_df) {
  mean_contrib_long <- mean_contrib_df %>%
    pivot_longer(-division, names_to = "factor", values_to = "mean_value")
  
  top3_factors <- mean_contrib_long %>%
    group_by(division) %>%
    arrange(desc(mean_value), .by_group = TRUE) %>%
    slice_head(n = 3) %>%
    mutate(rank = row_number()) %>%
    ungroup()
  
  top3_factors %>%
    select(division, rank, factor, mean_value) %>%
    pivot_wider(
      names_from = rank,
      values_from = c(factor, mean_value),
      names_glue = "top{rank}_{.value}"
    )
}

prepare_livestock_summary <- function(divisions_sf) {
  high_risk <- divisions_sf %>%
    filter(risk_cat == "High")
  
  high_medium_risk <- divisions_sf %>%
    filter(risk_cat %in% c("High", "Medium"))
  
  summary_high <- high_risk %>%
    summarise(
      n_divisions = n(),
      cattle_total = sum(cattle_total, na.rm = TRUE),
      sheep_total = sum(sheep_total, na.rm = TRUE),
      goat_total = sum(goat_total, na.rm = TRUE),
      livestock_total = sum(livestock_total, na.rm = TRUE)
    ) %>%
    mutate(risk_group = "High risk only")
  
  summary_high_medium <- high_medium_risk %>%
    summarise(
      n_divisions = n(),
      cattle_total = sum(cattle_total, na.rm = TRUE),
      sheep_total = sum(sheep_total, na.rm = TRUE),
      goat_total = sum(goat_total, na.rm = TRUE),
      livestock_total = sum(livestock_total, na.rm = TRUE)
    ) %>%
    mutate(risk_group = "High + Medium risk")
  
  national_total <- sum(divisions_sf$livestock_total, na.rm = TRUE)
  
  bind_rows(summary_high, summary_high_medium) %>%
    select(
      risk_group,
      n_divisions,
      cattle_total,
      sheep_total,
      goat_total,
      livestock_total
    ) %>%
    mutate(
      percent_national_livestock = 100 * livestock_total / national_total
    ) %>%
    mutate(
      across(
        c(cattle_total, sheep_total, goat_total, livestock_total, percent_national_livestock),
        ~ round(.x, 2)
      )
    )
}

# =====================================================================
# 4. LOAD COMMON INPUTS
# =====================================================================

divisions <- st_read(input_boundaries, quiet = TRUE)
divisions <- st_transform(divisions, crs = utm33n)

# Livestock density rasters used to estimate exposed domestic ruminants
cattle_density_real <- rast(
  here(input_dir_risk_factors, "livestock_density", "cattle_density_adjusted_2021_UTM33N_1km.tif")
)

sheep_density_real <- rast(
  here(input_dir_risk_factors, "livestock_density", "sheep_density_adjusted_2021_UTM33N_1km.tif")
)

goat_density_real <- rast(
  here(input_dir_risk_factors, "livestock_density", "goat_density_adjusted_2021_UTM33N_1km.tif")
)

# =====================================================================
# 5. SEASONAL PROCESSING FUNCTION
# =====================================================================

process_season <- function(season = c("dry", "wet")) {
  season <- match.arg(season)
  message("Processing ", season, " season...")
  
  season_dir <- here(output_dir, season)
  if (!dir.exists(season_dir)) {
    dir.create(season_dir, recursive = TRUE)
  }
  
  # -------------------------------------------------------------------
  # 5.1 Load season-specific combined risk raster
  # -------------------------------------------------------------------
  
  if (season == "dry") {
    combined_risk <- rast(here(output_models_dir, "Combined_risk_index_RVF_dry.tif"))
    ref_raster <- rast(here(input_dir_lit_dry, "rainfall_dry_std.tif"))
    
    # Expert weights (dry)
    w_expert <- c(
      rainfall = 0.066666667,
      protected_areas = 0.170416667,
      cattle_density = 0.133333333,
      small_ruminant_density = 0.175833300,
      water_bodies = 0.113750000,
      irrigation = 0.062083333,
      Aedes_density = 0.148750000,
      livestock_trade = 0.129166667
    )
    
    # Literature weights (dry)
    w_lit <- c(
      rainfall = 0.152173914,
      protected_areas = 0.065217392,
      cattle_density = 0.065217392,
      sheep_density = 0.065217392,
      water_bodies = 0.108695653,
      irrigation = 0.065217392,
      NDVI = 0.152173914,
      rivers_streams = 0.108695653,
      elevation = 0.152173914,
      Aedes_density = 0.065217392
    )
    
    # Expert rasters (dry)
    r_expert <- list(
      rainfall = rast(here(input_dir_exp_dry, "rainfall_dry_std.tif")),
      protected_areas = rast(here(input_dir_exp_dry, "protected_areas_dist_std.tif")),
      cattle_density = rast(here(input_dir_exp_dry, "cattle_density_std.tif")),
      small_ruminant_density = rast(here(input_dir_risk_factors, "livestock_density", "small_ruminant_density_std_UTM33N_1km.tif")),
      water_bodies = rast(here(input_dir_exp_dry, "water_bodies_perm_std.tif")),
      irrigation = rast(here(input_dir_exp_dry, "irrigation_std.tif")),
      Aedes_density = rast(here(input_dir_exp_dry, "mosquito_score_UTM33N_1km_cameroon_F32.tif")),
      livestock_trade = rast(here(input_dir_exp_dry, "distance_livestock_markets_standardized.tif"))
    )
    
    # Literature rasters (dry)
    r_lit <- list(
      rainfall = rast(here(input_dir_lit_dry, "rainfall_dry_std.tif")),
      protected_areas = rast(here(input_dir_lit_dry, "protected_areas_dist_std.tif")),
      cattle_density = rast(here(input_dir_lit_dry, "cattle_density_std.tif")),
      sheep_density = rast(here(input_dir_lit_dry, "sheep_density_std.tif")),
      water_bodies = rast(here(input_dir_lit_dry, "water_bodies_perm_std.tif")),
      irrigation = rast(here(input_dir_lit_dry, "irrigation_std.tif")),
      NDVI = rast(here(input_dir_risk_factors, "ndvi", "Standardized_NDVI_dry_UTM33N_1km_filled.tif")),
      rivers_streams = rast(here(input_dir_lit_dry, "river_streams_dry_std.tif")),
      elevation = rast(here(input_dir_lit_dry, "elevation_std.tif")),
      Aedes_density = rast(here(input_dir_lit_dry, "mosquito_score_UTM33N_1km_cameroon_F32.tif"))
    )
    
    low_thresh <- thresholds$dry$low
    high_thresh <- thresholds$dry$high
    season_label <- "Dry season"
    
  } else {
    combined_risk <- rast(here(output_models_dir, "Combined_risk_index_RVF_wet.tif"))
    ref_raster <- rast(here(input_dir_lit_wet, "rainfall_wet_std.tif"))
    
    # Expert weights (wet)
    w_expert <- c(
      rainfall = 0.157710000,
      protected_areas = 0.130420000,
      cattle_density = 0.097500000,
      small_ruminant_density = 0.144060000,
      water_bodies = 0.094060000,
      irrigation = 0.052190000,
      Aedes_density = 0.197710000,
      livestock_trade = 0.126350000
    )
    
    # Literature weights (wet)
    w_lit <- c(
      rainfall = 0.152173914,
      protected_areas = 0.065217392,
      cattle_density = 0.065217392,
      sheep_density = 0.065217392,
      water_bodies = 0.108695653,
      irrigation = 0.065217392,
      NDVI = 0.152173914,
      rivers_streams = 0.108695653,
      elevation = 0.152173914,
      Aedes_density = 0.065217392
    )
    
    # Expert rasters (wet)
    r_expert <- list(
      rainfall = rast(here(input_dir_exp_wet, "rainfall_wet_std.tif")),
      protected_areas = rast(here(input_dir_exp_wet, "protected_areas_dist_std.tif")),
      cattle_density = rast(here(input_dir_exp_wet, "cattle_density_std.tif")),
      small_ruminant_density = rast(here(input_dir_risk_factors, "livestock_density", "small_ruminant_density_std_UTM33N_1km.tif")),
      water_bodies = rast(here(input_dir_exp_wet, "water_bodies_perm_temp_std.tif")),
      irrigation = rast(here(input_dir_exp_wet, "irrigation_std.tif")),
      Aedes_density = rast(here(input_dir_exp_dry, "mosquito_score_UTM33N_1km_cameroon_F32.tif")),
      livestock_trade = rast(here(input_dir_exp_dry, "distance_livestock_markets_standardized.tif"))
    )
    
    # Literature rasters (wet)
    r_lit <- list(
      rainfall = rast(here(input_dir_lit_wet, "rainfall_wet_std.tif")),
      protected_areas = rast(here(input_dir_lit_wet, "protected_areas_dist_std.tif")),
      cattle_density = rast(here(input_dir_lit_wet, "cattle_density_std.tif")),
      sheep_density = rast(here(input_dir_lit_wet, "sheep_density_std.tif")),
      water_bodies = rast(here(input_dir_lit_wet, "water_bodies_perm_temp_std.tif")),
      irrigation = rast(here(input_dir_lit_wet, "irrigation_std.tif")),
      NDVI = rast(here(input_dir_risk_factors, "ndvi", "Standardized_NDVI_wet_UTM33N_1km_filled.tif")),
      rivers_streams = rast(here(input_dir_lit_wet, "river_streams_wet_std.tif")),
      elevation = rast(here(input_dir_lit_wet, "elevation_std.tif")),
      Aedes_density = rast(here(input_dir_lit_dry, "mosquito_score_UTM33N_1km_cameroon_F32.tif"))
    )
    
    low_thresh <- thresholds$wet$low
    high_thresh <- thresholds$wet$high
    season_label <- "Wet season"
  }
  
  crs(ref_raster) <- utm33n
  
  # -------------------------------------------------------------------
  # 5.2 Align factor rasters to a common reference
  # -------------------------------------------------------------------
  
  r_expert_aligned <- lapply(r_expert, function(r) align_to_ref(r, ref_raster))
  r_lit_aligned <- lapply(r_lit, function(r) align_to_ref(r, ref_raster))
  
  # -------------------------------------------------------------------
  # 5.3 Reconstruct harmonized factor contributions to the combined model
  # -------------------------------------------------------------------
  
  # The combined model is:
  # combined_risk = 0.5 * expert_model + 0.5 * literature_model
  #
  # For each harmonized factor, contribution maps are reconstructed as:
  # contribution_factor = 0.5 * expert_weight * expert_raster
  #                     + 0.5 * literature_weight * literature_raster
  #
  # For factors absent from one sub-model, contribution comes only from
  # the available side.
  
  contribution_maps <- list(
    rainfall =
      0.5 * w_expert["rainfall"] * r_expert_aligned$rainfall +
      0.5 * w_lit["rainfall"] * r_lit_aligned$rainfall,
    
    protected_areas =
      0.5 * w_expert["protected_areas"] * r_expert_aligned$protected_areas +
      0.5 * w_lit["protected_areas"] * r_lit_aligned$protected_areas,
    
    cattle_density =
      0.5 * w_expert["cattle_density"] * r_expert_aligned$cattle_density +
      0.5 * w_lit["cattle_density"] * r_lit_aligned$cattle_density,
    
    # Harmonized factor:
    # expert side = small_ruminant_density
    # literature side = sheep_density
    small_ruminant_density =
      0.5 * w_expert["small_ruminant_density"] * r_expert_aligned$small_ruminant_density +
      0.5 * w_lit["sheep_density"] * r_lit_aligned$sheep_density,
    
    water_bodies =
      0.5 * w_expert["water_bodies"] * r_expert_aligned$water_bodies +
      0.5 * w_lit["water_bodies"] * r_lit_aligned$water_bodies,
    
    irrigation =
      0.5 * w_expert["irrigation"] * r_expert_aligned$irrigation +
      0.5 * w_lit["irrigation"] * r_lit_aligned$irrigation,
    
    Aedes_density =
      0.5 * w_expert["Aedes_density"] * r_expert_aligned$Aedes_density +
      0.5 * w_lit["Aedes_density"] * r_lit_aligned$Aedes_density,
    
    livestock_trade =
      0.5 * w_expert["livestock_trade"] * r_expert_aligned$livestock_trade,
    
    NDVI =
      0.5 * w_lit["NDVI"] * r_lit_aligned$NDVI,
    
    rivers_streams =
      0.5 * w_lit["rivers_streams"] * r_lit_aligned$rivers_streams,
    
    elevation =
      0.5 * w_lit["elevation"] * r_lit_aligned$elevation
  )
  
  # Check alignment
  if (!check_alignment(ref_raster, contribution_maps)) {
    stop("Contribution rasters are not correctly aligned for the ", season, " season.")
  }
  
  contribution_stack <- rast(contribution_maps)
  
  # -------------------------------------------------------------------
  # 5.4 Extract mean factor contribution by division
  # -------------------------------------------------------------------
  
  mean_contrib <- terra::extract(
    contribution_stack,
    terra::vect(divisions),
    fun = mean,
    na.rm = TRUE
  )
  
  mean_contrib <- cbind(
    division = divisions$NAME_2,
    mean_contrib[, -1]
  )
  
  mean_contrib <- as.data.frame(mean_contrib)
  
  # Top 3 table
  top3_table <- get_top3_table(mean_contrib)
  
  write_xlsx(
    top3_table,
    path = here(season_dir, paste0("top3_factors_by_division_", season, ".xlsx"))
  )
  
  # Dominant factor
  dominant_factor <- apply(mean_contrib[, -1], 1, function(x) names(x)[which.max(x)])
  
  divisions_season <- divisions
  divisions_season$DominantFactor <- dominant_factor
  
  # -------------------------------------------------------------------
  # 5.5 Classify divisions using validated combined-model thresholds
  # -------------------------------------------------------------------
  
  divisions_season <- extract_mean_to_divisions(
    divisions_season,
    combined_risk,
    "combined_risk_mean"
  )
  
  divisions_season$risk_cat <- classify_risk(
    divisions_season$combined_risk_mean,
    low_thresh = low_thresh,
    high_thresh = high_thresh
  )
  
  # -------------------------------------------------------------------
  # 5.6 Estimate livestock totals per division
  # -------------------------------------------------------------------
  
  cattle_aligned <- align_to_ref(cattle_density_real, combined_risk)
  sheep_aligned <- align_to_ref(sheep_density_real, combined_risk)
  goat_aligned <- align_to_ref(goat_density_real, combined_risk)
  
  if (!all(
    compareGeom(combined_risk, cattle_aligned),
    compareGeom(combined_risk, sheep_aligned),
    compareGeom(combined_risk, goat_aligned)
  )) {
    stop("Livestock rasters are not correctly aligned for the ", season, " season.")
  }
  
  div_vect <- terra::vect(divisions_season)
  
  divisions_season$cattle_mean <- terra::extract(cattle_aligned, div_vect, fun = mean, na.rm = TRUE)[, 2]
  divisions_season$sheep_mean <- terra::extract(sheep_aligned, div_vect, fun = mean, na.rm = TRUE)[, 2]
  divisions_season$goat_mean <- terra::extract(goat_aligned, div_vect, fun = mean, na.rm = TRUE)[, 2]
  
  divisions_season$area_km2 <- terra::expanse(div_vect, unit = "km")
  
  divisions_season$cattle_total <- divisions_season$cattle_mean * divisions_season$area_km2
  divisions_season$sheep_total <- divisions_season$sheep_mean * divisions_season$area_km2
  divisions_season$goat_total <- divisions_season$goat_mean * divisions_season$area_km2
  
  divisions_season$livestock_total <- divisions_season$cattle_total +
    divisions_season$sheep_total +
    divisions_season$goat_total
  
  # Detailed table for High and Medium divisions
  divisions_at_risk <- divisions_season %>%
    st_drop_geometry() %>%
    filter(risk_cat %in% c("High", "Medium")) %>%
    transmute(
      Division = NAME_2,
      Region = NAME_1,
      Risk_category = as.character(risk_cat),
      Area_km2 = round(area_km2, 2),
      Cattle_mean_density = round(cattle_mean, 2),
      Sheep_mean_density = round(sheep_mean, 2),
      Goat_mean_density = round(goat_mean, 2),
      Cattle_total = round(cattle_total, 2),
      Sheep_total = round(sheep_total, 2),
      Goat_total = round(goat_total, 2),
      Livestock_total = round(livestock_total, 2)
    ) %>%
    arrange(factor(Risk_category, levels = c("High", "Medium")), desc(Livestock_total))
  
  livestock_summary <- prepare_livestock_summary(divisions_season)
  
  write_xlsx(
    list(
      Summary = livestock_summary,
      Divisions_High_Medium = divisions_at_risk
    ),
    path = here(season_dir, paste0("livestock_exposure_", season, ".xlsx"))
  )
  
  # -------------------------------------------------------------------
  # 5.7 Dominant factor map
  # -------------------------------------------------------------------
  
  divisions_labels <- st_point_on_surface(divisions_season)
  
  tmap_mode("plot")
  
  tm_dom_factor <- tm_shape(divisions_season) +
    tm_polygons(
      "DominantFactor",
      palette = factor_colors,
      title = "Dominant risk factor",
      border.col = "gray40",
      lwd = 0.3
    ) +
    tm_shape(divisions_labels) +
    tm_text(
      "NAME_2",
      size = 0.7,
      col = "black",
      shadow = FALSE,
      remove.overlap = TRUE
    ) +
    tm_layout(
      title = season_label,
      legend.outside = TRUE,
      legend.outside.position = "right",
      frame = FALSE
    ) +
    tm_compass(
      type = "arrow",
      position = c("left", "top"),
      size = 2
    ) +
    tm_scale_bar(
      position = c("left", "bottom"),
      text.size = 0.7,
      breaks = c(0, 50, 100),
      width = 0.25
    )
  
  tmap_save(
    tm = tm_dom_factor,
    filename = here(season_dir, paste0("dominant_factor_by_division_", season, ".tiff")),
    width = 2000,
    height = 1600,
    units = "px",
    dpi = 300
  )
  
  # -------------------------------------------------------------------
  # 5.8 Optional shapefile export
  # -------------------------------------------------------------------
  
  dominant_factor_shp <- divisions_season %>%
    select(
      GID_0,
      COUNTRY,
      GID_1,
      NAME_1,
      NAME_2,
      DominantFactor,
      combined_risk_mean,
      risk_cat
    )
  
  st_write(
    dominant_factor_shp,
    here(season_dir, paste0("dominant_factor_by_division_", season, ".shp")),
    delete_layer = TRUE,
    quiet = TRUE
  )
  
  # -------------------------------------------------------------------
  # 5.9 Return objects
  # -------------------------------------------------------------------
  
  list(
    divisions = divisions_season,
    livestock_summary = livestock_summary,
    top3_table = top3_table
  )
}

# =====================================================================
# 6. RUN DRY AND WET ANALYSES
# =====================================================================

dry_results <- process_season("dry")
wet_results <- process_season("wet")

message("Exposure and dominant risk driver analyses completed successfully.")
