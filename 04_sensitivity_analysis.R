# Script: 04_sensitivity_analysis.R
# Purpose: perform sensitivity analysis of the combined RVF risk model
# Inputs: standardized raster layers produced during data processing
# Outputs: surrogate combined risk rasters, uncertainty rasters, average relative change (ARC) rasters,
#          stability tables, and sensitivity plots in outputs/sensitivity_analysis/

library(terra)
library(here)
library(dplyr)
library(ggplot2)
library(forcats)
library(viridis)
library(writexl)
library(tmap)

# =====================================================================
# 1. BACKGROUND AND RATIONALE
# =====================================================================

# The final combined RVF model used in the manuscript was produced by averaging
# the expert-based and literature-based seasonal risk maps with equal weights
# (50:50). Because this final combination was performed at the map level,
# direct perturbation of the original variable-level weights was not possible.
#
# To enable sensitivity analysis of the combined framework, a surrogate
# combined model was therefore constructed. For each harmonized variable i,
# a surrogate combined weight was calculated as:
#
#   w_i_combined = 0.5 * w_i_expert + 0.5 * w_i_literature
#
# This surrogate weighting scheme was designed specifically for sensitivity
# analysis and should be interpreted as an approximation of the integrated
# model structure, consistent with the equal contribution of the expert-based
# and literature-based approaches in the final combined risk map.
#
# Variable harmonization:
# The expert-based and literature-based models did not use exactly the same
# variable names. To make variable-level sensitivity analysis possible,
# variable names were harmonized across the two modelling schemes.
#
# In particular:
# - the harmonized variable name 'small_ruminant_density' was used in the
#   surrogate combined model;
# - on the expert side, this variable corresponds to the original
#   'small_ruminant_density';
# - on the literature side, the corresponding available variable was
#   'sheep_density';
# - therefore, the literature-based weight assigned to 'sheep_density' was
#   mapped to the harmonized variable 'small_ruminant_density'.
#
# This naming harmonization was necessary to ensure comparability between the
# two model structures and to allow calculation of surrogate combined weights.
#
# Likewise, some variables were present in one model but absent from the other:
# - NDVI, rivers_streams, and elevation were included in the literature-based
#   model but not in the expert-based model;
# - livestock_trade was included in the expert-based model but not in the
#   literature-based model.
#
# For the purpose of computing surrogate combined weights, variables absent from
# one of the two models were assigned a weight of 0 in that model.
#
# This script performs:
# 1. construction of surrogate combined dry- and wet-season risk maps;
# 2. sensitivity analysis by perturbing each variable weight from -25% to +25%;
# 3. uncertainty mapping using the pixel-wise standard deviation across
#    perturbed maps;
# 4. leave-one-out sensitivity analysis;
# 5. variable stability metrics;
# 6. average relative change (ARC) mapping.

# =====================================================================
# 2. DEFINE PROJECTION AND PATHS
# =====================================================================

utm33n <- "EPSG:32633"

input_dir_lit_dry <- here("inputs", "standardized_rasters", "literature", "dry")
input_dir_lit_wet <- here("inputs", "standardized_rasters", "literature", "wet")
input_dir_experts_dry <- here("inputs", "standardized_rasters", "experts", "dry")
input_dir_risk_factors <- here("inputs", "risk_factors")

output_dir <- here("outputs", "sensitivity_analysis")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

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

create_risk_map <- function(weights, rasters) {
  risk <- rasters[[1]] * weights[1]
  for (i in 2:length(rasters)) {
    risk <- risk + rasters[[i]] * weights[i]
  }
  return(risk)
}

run_weight_perturbation <- function(base_weights, rasters, deltas) {
  risk_maps_list <- list()
  factors <- names(base_weights)
  
  for (factor in factors) {
    for (delta in deltas) {
      new_weight <- base_weights[factor] * (1 + delta)
      
      other_factors <- setdiff(factors, factor)
      remaining <- 1 - new_weight
      prop_others <- base_weights[other_factors] / sum(base_weights[other_factors])
      
      adjusted_weights <- setNames(numeric(length(base_weights)), factors)
      adjusted_weights[factor] <- new_weight
      adjusted_weights[other_factors] <- prop_others * remaining
      
      risk_temp <- create_risk_map(adjusted_weights, rasters)
      risk_maps_list[[paste(factor, delta, sep = "_")]] <- risk_temp
    }
  }
  
  return(risk_maps_list)
}

run_leave_one_out <- function(base_weights, rasters) {
  risk_maps_exclusion <- list()
  factors <- names(base_weights)
  
  for (var in factors) {
    remaining_factors <- setdiff(factors, var)
    adjusted_weights <- base_weights[remaining_factors] / sum(base_weights[remaining_factors])
    
    new_weights <- base_weights
    new_weights[remaining_factors] <- adjusted_weights
    new_weights[var] <- 0
    
    risk_new <- create_risk_map(new_weights, rasters)
    risk_maps_exclusion[[paste0("removed_", var)]] <- risk_new
  }
  
  return(risk_maps_exclusion)
}

compute_stability_table <- function(risk_maps_list, base_weights, risk_ref) {
  stability_table <- data.frame(
    Variable = character(),
    Mean_change = numeric(),
    SD_change = numeric(),
    CV_change = numeric(),
    SSD_change = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (var in names(base_weights)) {
    var_maps <- risk_maps_list[grep(paste0("^", var, "_"), names(risk_maps_list))]
    
    diff_values <- c()
    ssd_total <- 0
    
    for (m in var_maps) {
      diff_raster <- m - risk_ref
      vals <- values(diff_raster, na.rm = TRUE)
      diff_values <- c(diff_values, vals)
      ssd_total <- ssd_total + sum(vals^2, na.rm = TRUE)
    }
    
    mean_change <- mean(abs(diff_values), na.rm = TRUE)
    sd_change <- sd(abs(diff_values), na.rm = TRUE)
    cv_change <- sd_change / mean_change
    
    stability_table <- rbind(
      stability_table,
      data.frame(
        Variable = var,
        Mean_change = mean_change,
        SD_change = sd_change,
        CV_change = cv_change,
        SSD_change = ssd_total
      )
    )
  }
  
  stability_table <- stability_table[order(stability_table$CV_change), ]
  return(stability_table)
}

compute_arc_map <- function(risk_ref, risk_maps_list) {
  arc_fun <- function(x) {
    r0 <- x[1]
    ri <- x[-1]
    mean(abs((ri - r0) / r0), na.rm = TRUE)
  }
  
  risk_stack_arc <- c(risk_ref, rast(risk_maps_list))
  arc_map <- app(risk_stack_arc, fun = arc_fun)
  return(arc_map)
}

save_cv_plot <- function(stability_table, season_label, filename) {
  plot_data <- stability_table %>%
    arrange(desc(CV_change)) %>%
    mutate(Variable = fct_reorder(Variable, CV_change))
  
  p <- ggplot(plot_data, aes(x = Variable, y = CV_change, fill = CV_change)) +
    geom_col(width = 0.7) +
    scale_fill_viridis(option = "plasma", direction = -1) +
    coord_flip() +
    theme_minimal(base_size = 14) +
    labs(
      title = season_label,
      x = "Variables",
      y = "Coefficient of variation – relative instability",
      fill = "Coefficient of variation"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      legend.position = "right"
    )
  
  ggsave(
    filename = filename,
    plot = p,
    width = 10,
    height = 8,
    units = "in",
    dpi = 600
  )
}

# =====================================================================
# 4. DRY SEASON
# =====================================================================

message("Running sensitivity analysis for the dry season...")

# ---------------------------------------------------------------------
# 4.1 Load reference raster
# ---------------------------------------------------------------------

ref_dry <- rast(here(input_dir_lit_dry, "rainfall_dry_std.tif"))

# Assign CRS only if missing or incorrectly defined in the source file
crs(ref_dry) <- utm33n

# ---------------------------------------------------------------------
# 4.2 Load and align all rasters
# ---------------------------------------------------------------------

rainfall <- ref_dry

protected_areas <- align_to_ref(
  rast(here(input_dir_lit_dry, "protected_areas_dist_std.tif")),
  ref_dry
)

cattle_density <- align_to_ref(
  rast(here(input_dir_lit_dry, "cattle_density_std.tif")),
  ref_dry
)

small_ruminant_density <- align_to_ref(
  rast(here(input_dir_risk_factors, "livestock_density", "small_ruminant_density_std_UTM33N_1km.tif")),
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

NDVI <- align_to_ref(
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

livestock_trade <- align_to_ref(
  rast(here(input_dir_experts_dry, "distance_livestock_markets_standardized.tif")),
  ref_dry
)

rasters_dry <- list(
  rainfall = rainfall,
  protected_areas = protected_areas,
  cattle_density = cattle_density,
  small_ruminant_density = small_ruminant_density,
  water_bodies = water_bodies,
  irrigation = irrigation,
  NDVI = NDVI,
  rivers_streams = rivers_streams,
  elevation = elevation,
  Aedes_density = Aedes_density,
  livestock_trade = livestock_trade
)

if (!check_alignment(ref_dry, rasters_dry)) {
  stop("Dry-season rasters are not correctly aligned.")
}

# ---------------------------------------------------------------------
# 4.3 Define expert-based and literature-based weights
# ---------------------------------------------------------------------

# Expert-based dry-season weights
expert_weights_dry <- c(
  rainfall = 0.066666667,
  protected_areas = 0.170416667,
  cattle_density = 0.133333333,
  small_ruminant_density = 0.175833300,
  water_bodies = 0.113750000,
  irrigation = 0.062083333,
  NDVI = 0,
  rivers_streams = 0,
  elevation = 0,
  Aedes_density = 0.148750000,
  livestock_trade = 0.129166667
)

# Literature-based dry-season weights
# The literature-based model originally used 'sheep_density'.
# For harmonization, that weight is assigned here to the common variable
# name 'small_ruminant_density'.
literature_weights_dry <- c(
  rainfall = 0.152173914,
  protected_areas = 0.065217392,
  cattle_density = 0.065217392,
  small_ruminant_density = 0.065217392,  # harmonized from sheep_density
  water_bodies = 0.108695653,
  irrigation = 0.065217392,
  NDVI = 0.152173914,
  rivers_streams = 0.108695653,
  elevation = 0.152173914,
  Aedes_density = 0.065217392,
  livestock_trade = 0
)

# ---------------------------------------------------------------------
# 4.4 Calculate surrogate combined weights
# ---------------------------------------------------------------------

# For each harmonized variable:
# surrogate_combined_weight = 0.5 * expert_weight + 0.5 * literature_weight
weights_dry <- 0.5 * expert_weights_dry + 0.5 * literature_weights_dry

if (abs(sum(weights_dry) - 1) > 1e-6) {
  stop("Dry-season surrogate combined weights do not sum to 1.")
}

# ---------------------------------------------------------------------
# 4.5 Build baseline surrogate combined risk map
# ---------------------------------------------------------------------

risk_index_dry <- rainfall               * weights_dry["rainfall"] +
                  protected_areas        * weights_dry["protected_areas"] +
                  cattle_density         * weights_dry["cattle_density"] +
                  small_ruminant_density * weights_dry["small_ruminant_density"] +
                  water_bodies           * weights_dry["water_bodies"] +
                  irrigation             * weights_dry["irrigation"] +
                  NDVI                   * weights_dry["NDVI"] +
                  rivers_streams         * weights_dry["rivers_streams"] +
                  elevation              * weights_dry["elevation"] +
                  Aedes_density          * weights_dry["Aedes_density"] +
                  livestock_trade        * weights_dry["livestock_trade"]

writeRaster(
  risk_index_dry,
  here(output_dir, "combined_weights_model_dry.tif"),
  overwrite = TRUE
)

# ---------------------------------------------------------------------
# 4.6 Sensitivity analysis by weight perturbation
# ---------------------------------------------------------------------

deltas <- seq(-0.25, 0.25, by = 0.05)

risk_maps_list_dry <- run_weight_perturbation(
  base_weights = weights_dry,
  rasters = rasters_dry,
  deltas = deltas
)

risk_stack_dry <- rast(risk_maps_list_dry)
uncertainty_map_dry <- app(risk_stack_dry, fun = sd, na.rm = TRUE)

writeRaster(
  uncertainty_map_dry,
  here(output_dir, "uncertainty_surface_combined_dry.tif"),
  overwrite = TRUE
)

# ---------------------------------------------------------------------
# 4.7 Leave-one-out sensitivity analysis
# ---------------------------------------------------------------------

risk_maps_exclusion_dry <- run_leave_one_out(
  base_weights = weights_dry,
  rasters = rasters_dry
)

risk_stack_exclusion_dry <- rast(risk_maps_exclusion_dry)
uncertainty_leave_one_out_dry <- app(risk_stack_exclusion_dry, fun = sd, na.rm = TRUE)

writeRaster(
  uncertainty_leave_one_out_dry,
  here(output_dir, "uncertainty_surface_leave_one_out_dry.tif"),
  overwrite = TRUE
)

# ---------------------------------------------------------------------
# 4.8 Stability indicators
# ---------------------------------------------------------------------

stability_table_dry <- compute_stability_table(
  risk_maps_list = risk_maps_list_dry,
  base_weights = weights_dry,
  risk_ref = risk_index_dry
)

write_xlsx(
  stability_table_dry,
  path = here(output_dir, "stability_table_dry.xlsx")
)

# ---------------------------------------------------------------------
# 4.9 Coefficient of variation (CV) barplot
# ---------------------------------------------------------------------

save_cv_plot(
  stability_table = stability_table_dry,
  season_label = "Dry season",
  filename = here(output_dir, "variable_sensitivity_dry.tiff")
)

# ---------------------------------------------------------------------
# 4.10 Average Relative Change (ARC)
# ---------------------------------------------------------------------

arc_map_dry <- compute_arc_map(
  risk_ref = risk_index_dry,
  risk_maps_list = risk_maps_list_dry
)

writeRaster(
  arc_map_dry,
  here(output_dir, "average_relative_change_dry.tif"),
  overwrite = TRUE
)

# =====================================================================
# 5. WET SEASON
# =====================================================================

message("Running sensitivity analysis for the wet season...")

# ---------------------------------------------------------------------
# 5.1 Load reference raster
# ---------------------------------------------------------------------

ref_wet <- rast(here(input_dir_lit_wet, "rainfall_wet_std.tif"))
crs(ref_wet) <- utm33n

# ---------------------------------------------------------------------
# 5.2 Load and align all rasters
# ---------------------------------------------------------------------

rainfall <- ref_wet

protected_areas <- align_to_ref(
  rast(here(input_dir_lit_wet, "protected_areas_dist_std.tif")),
  ref_wet
)

cattle_density <- align_to_ref(
  rast(here(input_dir_lit_wet, "cattle_density_std.tif")),
  ref_wet
)

small_ruminant_density <- align_to_ref(
  rast(here(input_dir_risk_factors, "livestock_density", "small_ruminant_density_std_UTM33N_1km.tif")),
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

NDVI <- align_to_ref(
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

livestock_trade <- align_to_ref(
  rast(here(input_dir_experts_dry, "distance_livestock_markets_standardized.tif")),
  ref_wet
)

rasters_wet <- list(
  rainfall = rainfall,
  protected_areas = protected_areas,
  cattle_density = cattle_density,
  small_ruminant_density = small_ruminant_density,
  water_bodies = water_bodies,
  irrigation = irrigation,
  NDVI = NDVI,
  rivers_streams = rivers_streams,
  elevation = elevation,
  Aedes_density = Aedes_density,
  livestock_trade = livestock_trade
)

if (!check_alignment(ref_wet, rasters_wet)) {
  stop("Wet-season rasters are not correctly aligned.")
}

# ---------------------------------------------------------------------
# 5.3 Define expert-based and literature-based weights
# ---------------------------------------------------------------------

expert_weights_wet <- c(
  rainfall = 0.157710000,
  protected_areas = 0.130420000,
  cattle_density = 0.097500000,
  small_ruminant_density = 0.144060000,
  water_bodies = 0.094060000,
  irrigation = 0.052190000,
  NDVI = 0,
  rivers_streams = 0,
  elevation = 0,
  Aedes_density = 0.197710000,
  livestock_trade = 0.126350000
)

# The literature-based model originally used 'sheep_density'.
# For harmonization, that weight is assigned here to the common variable
# name 'small_ruminant_density'.
literature_weights_wet <- c(
  rainfall = 0.152173914,
  protected_areas = 0.065217392,
  cattle_density = 0.065217392,
  small_ruminant_density = 0.065217392,  # harmonized from sheep_density
  water_bodies = 0.108695653,
  irrigation = 0.065217392,
  NDVI = 0.152173914,
  rivers_streams = 0.108695653,
  elevation = 0.152173914,
  Aedes_density = 0.065217392,
  livestock_trade = 0
)

# ---------------------------------------------------------------------
# 5.4 Calculate surrogate combined weights
# ---------------------------------------------------------------------

weights_wet <- 0.5 * expert_weights_wet + 0.5 * literature_weights_wet

if (abs(sum(weights_wet) - 1) > 1e-6) {
  stop("Wet-season surrogate combined weights do not sum to 1.")
}

# ---------------------------------------------------------------------
# 5.5 Build baseline surrogate combined risk map
# ---------------------------------------------------------------------

risk_index_wet <- rainfall               * weights_wet["rainfall"] +
                  protected_areas        * weights_wet["protected_areas"] +
                  cattle_density         * weights_wet["cattle_density"] +
                  small_ruminant_density * weights_wet["small_ruminant_density"] +
                  water_bodies           * weights_wet["water_bodies"] +
                  irrigation             * weights_wet["irrigation"] +
                  NDVI                   * weights_wet["NDVI"] +
                  rivers_streams         * weights_wet["rivers_streams"] +
                  elevation              * weights_wet["elevation"] +
                  Aedes_density          * weights_wet["Aedes_density"] +
                  livestock_trade        * weights_wet["livestock_trade"]

writeRaster(
  risk_index_wet,
  here(output_dir, "combined_weights_model_wet.tif"),
  overwrite = TRUE
)

# ---------------------------------------------------------------------
# 5.6 Sensitivity analysis by weight perturbation
# ---------------------------------------------------------------------

risk_maps_list_wet <- run_weight_perturbation(
  base_weights = weights_wet,
  rasters = rasters_wet,
  deltas = deltas
)

risk_stack_wet <- rast(risk_maps_list_wet)
uncertainty_map_wet <- app(risk_stack_wet, fun = sd, na.rm = TRUE)

writeRaster(
  uncertainty_map_wet,
  here(output_dir, "uncertainty_surface_combined_wet.tif"),
  overwrite = TRUE
)

# ---------------------------------------------------------------------
# 5.7 Leave-one-out sensitivity analysis
# ---------------------------------------------------------------------

risk_maps_exclusion_wet <- run_leave_one_out(
  base_weights = weights_wet,
  rasters = rasters_wet
)

risk_stack_exclusion_wet <- rast(risk_maps_exclusion_wet)
uncertainty_leave_one_out_wet <- app(risk_stack_exclusion_wet, fun = sd, na.rm = TRUE)

writeRaster(
  uncertainty_leave_one_out_wet,
  here(output_dir, "uncertainty_surface_leave_one_out_wet.tif"),
  overwrite = TRUE
)

# ---------------------------------------------------------------------
# 5.8 Stability indicators
# ---------------------------------------------------------------------

stability_table_wet <- compute_stability_table(
  risk_maps_list = risk_maps_list_wet,
  base_weights = weights_wet,
  risk_ref = risk_index_wet
)

write_xlsx(
  stability_table_wet,
  path = here(output_dir, "stability_table_wet.xlsx")
)

# ---------------------------------------------------------------------
# 5.9 CV barplot
# ---------------------------------------------------------------------

save_cv_plot(
  stability_table = stability_table_wet,
  season_label = "Wet season",
  filename = here(output_dir, "variable_sensitivity_wet.tiff")
)

# ---------------------------------------------------------------------
# 5.10 Average Relative Change
# ---------------------------------------------------------------------

arc_map_wet <- compute_arc_map(
  risk_ref = risk_index_wet,
  risk_maps_list = risk_maps_list_wet
)

writeRaster(
  arc_map_wet,
  here(output_dir, "average_relative_change_wet.tif"),
  overwrite = TRUE
)

# =====================================================================
# 6. OPTIONAL PUBLICATION-READY COMBINED FIGURES
# =====================================================================

tmap_mode("plot")

# ---------------------------------------------------------------------
# 6.1 Uncertainty maps
# ---------------------------------------------------------------------

tm_uncertainty_dry <- tm_shape(uncertainty_map_dry) +
  tm_raster(
    palette = "-RdBu",
    style = "cont",
    legend.show = FALSE
  ) +
  tm_compass(type = "arrow", position = c("left", "top"), size = 2) +
  tm_scale_bar(
    position = c("left", "bottom"),
    text.size = 0.9,
    breaks = c(0, 50, 100),
    width = 0.2
  ) +
  tm_layout(
    main.title = "Dry season",
    main.title.position = c("center", "top"),
    frame = FALSE,
    inner.margins = c(0.02, 0.02, 0.015, 0.02)
  )

tm_uncertainty_wet <- tm_shape(uncertainty_map_wet) +
  tm_raster(
    palette = "-RdBu",
    style = "cont",
    title = "Standard deviation\n(combined weights)",
    legend.show = TRUE
  ) +
  tm_layout(
    main.title = "Wet season",
    main.title.position = c("center", "top"),
    frame = FALSE,
    legend.outside = TRUE,
    legend.outside.position = "right"
  )

uncertainty_maps_combined <- tmap_arrange(
  tm_uncertainty_dry,
  tm_uncertainty_wet,
  ncol = 2,
  sync = TRUE
)

tmap_save(
  uncertainty_maps_combined,
  filename = here(output_dir, "RVF_uncertainty_maps_combined_dry_wet.tiff"),
  width = 4500,
  height = 2200,
  dpi = 300
)

# ---------------------------------------------------------------------
# 6.2 ARC maps
# ---------------------------------------------------------------------

tm_arc_dry <- tm_shape(arc_map_dry) +
  tm_raster(
    palette = "YlGnBu",
    style = "cont",
    legend.show = FALSE
  ) +
  tm_compass(type = "arrow", position = c("left", "top"), size = 2) +
  tm_scale_bar(
    position = c("left", "bottom"),
    text.size = 0.9,
    breaks = c(0, 50, 100),
    width = 0.2
  ) +
  tm_layout(
    main.title = "Dry season",
    main.title.position = c("center", "top"),
    frame = FALSE,
    inner.margins = c(0.02, 0.02, 0.015, 0.02)
  )

tm_arc_wet <- tm_shape(arc_map_wet) +
  tm_raster(
    palette = "YlGnBu",
    style = "cont",
    title = "Average relative change\n(combined weights)",
    legend.show = TRUE
  ) +
  tm_layout(
    main.title = "Wet season",
    main.title.position = c("center", "top"),
    frame = FALSE,
    legend.outside = TRUE,
    legend.outside.position = "right"
  )

arc_maps_combined <- tmap_arrange(
  tm_arc_dry,
  tm_arc_wet,
  ncol = 2,
  sync = TRUE
)

tmap_save(
  arc_maps_combined,
  filename = here(output_dir, "RVF_average_relative_change_combined_dry_wet.tiff"),
  width = 4500,
  height = 2200,
  dpi = 300
)

message("Sensitivity analysis completed successfully.")
