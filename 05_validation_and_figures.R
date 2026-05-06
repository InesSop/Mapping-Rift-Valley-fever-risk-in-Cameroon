# Script: 05_validation_and_figures.R
# Purpose: validate the combined RVF risk model at division level and generate
#          the main validation tables and figures for dry and wet seasons
# Inputs:
#   - seasonal expert-based risk rasters
#   - seasonal literature-based risk rasters
#   - seasonal combined risk rasters
#   - administrative boundaries (division level)
#   - division-level seroprevalence table
# Outputs:
#   - ROC summary tables
#   - final AUC comparison tables
#   - division-level validation tables
#   - average combined risk maps by division
#   - risk category maps
#   - combined-model ROC figures
#   - optional shapefiles for mapped outputs

library(here)
library(terra)
library(sf)
library(dplyr)
library(stringi)
library(pROC)
library(ggplot2)
library(tmap)
library(writexl)

# =====================================================================
# 1. DEFINE PATHS AND GLOBAL OPTIONS
# =====================================================================

utm33n <- "EPSG:32633"

input_boundaries <- here("inputs", "boundaries", "gadm41_CMR_2.shp")
input_serology <- here("inputs", "validation", "global_seroprev.csv")

output_models_dir <- here("outputs")
output_dir <- here("outputs", "validation_and_figures")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Thresholds explored in ROC analyses
biological_cutoffs <- c(10, 15)
extra_cutoffs_fixed <- c(5, 7.5, 12.5, 17.5, 20, 25, 30)
quantile_probs <- c(0.50, 0.60, 0.70, 0.75, 0.80, 0.85)

# Minimum class size required to retain a cutoff for interpretation
min_class_size <- 5

# Specificity target used to derive the upper threshold for final 3-class maps
target_specificity_for_high <- 0.75

# =====================================================================
# 2. HELPER FUNCTIONS
# =====================================================================

normalize_division_names <- function(x) {
  toupper(stri_trans_general(x, "Latin-ASCII"))
}

roc_summary <- function(data, risk_var, cutoff, boot_n = 2000) {
  sero_bin <- ifelse(data$prevalence_global >= cutoff, 1, 0)
  
  n_pos <- sum(sero_bin == 1, na.rm = TRUE)
  n_neg <- sum(sero_bin == 0, na.rm = TRUE)
  
  if (length(unique(sero_bin[!is.na(sero_bin)])) < 2) {
    return(data.frame(
      cutoff = cutoff,
      model = risk_var,
      auc = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      n_pos = n_pos,
      n_neg = n_neg
    ))
  }
  
  roc_obj <- roc(sero_bin, data[[risk_var]], quiet = TRUE)
  ci_obj <- ci.auc(roc_obj, method = "bootstrap", boot.n = boot_n)
  
  data.frame(
    cutoff = cutoff,
    model = risk_var,
    auc = as.numeric(auc(roc_obj)),
    ci_low = as.numeric(ci_obj[1]),
    ci_high = as.numeric(ci_obj[3]),
    n_pos = n_pos,
    n_neg = n_neg
  )
}

compute_all_roc_results <- function(data, risk_vars, all_cutoffs) {
  results <- data.frame()
  
  for (c in all_cutoffs) {
    for (r in risk_vars) {
      res <- roc_summary(data, r, c)
      results <- rbind(results, res)
    }
  }
  
  results
}

derive_high_threshold <- function(roc_coords, low_thresh, target_specificity = 0.75) {
  candidates <- roc_coords %>%
    filter(
      is.finite(threshold),
      threshold > low_thresh,
      specificity >= target_specificity
    ) %>%
    arrange(threshold)
  
  if (nrow(candidates) > 0) {
    return(candidates$threshold[1])
  }
  
  fallback <- roc_coords %>%
    filter(is.finite(threshold), threshold > low_thresh) %>%
    arrange(desc(specificity), threshold)
  
  if (nrow(fallback) > 0) {
    return(fallback$threshold[1])
  }
  
  return(max(roc_coords$threshold[is.finite(roc_coords$threshold)], na.rm = TRUE))
}

save_roc_step_plot <- function(roc_obj, title_text, output_file) {
  roc_coords <- coords(
    roc_obj,
    x = "all",
    input = "threshold",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  roc_coords <- as.data.frame(roc_coords)
  roc_coords <- roc_coords[order(roc_coords$threshold), ]
  roc_jumps <- roc_coords[c(TRUE, diff(roc_coords$sensitivity) != 0), ]
  
  tiff(
    filename = output_file,
    width = 6,
    height = 6,
    units = "in",
    res = 300,
    compression = "lzw"
  )
  
  par(mar = c(5, 5, 6, 2))
  
  plot(
    1 - roc_coords$specificity,
    roc_coords$sensitivity,
    type = "s",
    lwd = 2,
    lty = 2,
    xlab = "1 - Specificity",
    ylab = "Sensitivity",
    xlim = c(0, 1.05),
    ylim = c(0, 1.05),
    main = title_text,
    col = "black"
  )
  
  abline(0, 1, lty = 2, col = "grey")
  
  text(
    x = 1 - roc_jumps$specificity + 0.03,
    y = roc_jumps$sensitivity,
    labels = round(roc_jumps$threshold, 3),
    pos = 3,
    cex = 0.75
  )
  
  auc_val <- auc(roc_obj)
  text(
    x = 0.72,
    y = 0.45,
    labels = paste0("AUC = ", round(auc_val, 3)),
    cex = 1.2,
    font = 2
  )
  
  dev.off()
}

extract_division_means <- function(zones_sf, rast_obj, output_name) {
  zones_vect <- terra::vect(zones_sf)
  values <- terra::extract(rast_obj, zones_vect, fun = mean, na.rm = TRUE)[, 2]
  zones_sf[[output_name]] <- values
  zones_sf
}

# =====================================================================
# 3. LOAD COMMON INPUTS
# =====================================================================

zones <- st_read(input_boundaries, quiet = TRUE)
zones <- st_transform(zones, crs = utm33n)

seroprev <- read.csv(input_serology, stringsAsFactors = FALSE)

required_cols <- c("Departement", "prevalence_global")
missing_cols <- setdiff(required_cols, names(seroprev))
if (length(missing_cols) > 0) {
  stop("Missing required columns in serology table: ", paste(missing_cols, collapse = ", "))
}

seroprev <- seroprev %>%
  mutate(dept_name = normalize_division_names(Departement))

# =====================================================================
# 4. SEASON PROCESSING FUNCTION
# =====================================================================

process_season <- function(season = c("dry", "wet")) {
  season <- match.arg(season)
  season_dir <- here(output_dir, season)
  if (!dir.exists(season_dir)) {
    dir.create(season_dir, recursive = TRUE)
  }
  
  message("Processing ", season, " season...")
  
  # -------------------------------------------------------------------
  # 4.1 Define season-specific files
  # -------------------------------------------------------------------
  
  if (season == "dry") {
    risk_lit_file <- here(output_models_dir, "risk_index_RVF_lit_dry.tif")
    risk_expert_file <- here(output_models_dir, "risk_index_RVF_expert_dry.tif")
    risk_combined_file <- here(output_models_dir, "Combined_risk_index_RVF_dry.tif")
    map_title <- "Dry season"
  } else {
    risk_lit_file <- here(output_models_dir, "risk_index_RVF_lit_wet.tif")
    risk_expert_file <- here(output_models_dir, "risk_index_RVF_expert_wet.tif")
    risk_combined_file <- here(output_models_dir, "Combined_risk_index_RVF_wet.tif")
    map_title <- "Wet season"
  }
  
  # -------------------------------------------------------------------
  # 4.2 Load rasters
  # -------------------------------------------------------------------
  
  risk_lit <- rast(risk_lit_file)
  risk_expert <- rast(risk_expert_file)
  risk_combined <- rast(risk_combined_file)
  
  if (!compareGeom(risk_lit, risk_expert) || !compareGeom(risk_lit, risk_combined)) {
    stop("Raster geometries are not aligned for the ", season, " season.")
  }
  
  # -------------------------------------------------------------------
  # 4.3 Extract mean risk per division
  # -------------------------------------------------------------------
  
  zones_season <- zones
  
  zones_season <- extract_division_means(zones_season, risk_lit, "risk_lit")
  zones_season <- extract_division_means(zones_season, risk_expert, "risk_expert")
  zones_season <- extract_division_means(zones_season, risk_combined, "risk_combined_mean")
  
  zones_season <- zones_season %>%
    mutate(dept_name = normalize_division_names(NAME_2)) %>%
    left_join(
      seroprev[, c("dept_name", "prevalence_global")],
      by = "dept_name"
    )
  
  # -------------------------------------------------------------------
  # 4.4 Export division-level validation table
  # -------------------------------------------------------------------
  
  validation_table <- zones_season %>%
    st_drop_geometry() %>%
    select(
      dept_name,
      NAME_1,
      NAME_2,
      risk_lit,
      risk_expert,
      risk_combined_mean,
      prevalence_global
    ) %>%
    rename(
      Region = NAME_1,
      Division = NAME_2,
      Combined_risk = risk_combined_mean,
      Literature_risk = risk_lit,
      Expert_risk = risk_expert,
      Seroprevalence = prevalence_global
    )
  
  write_xlsx(
    validation_table,
    path = here(season_dir, paste0("validation_table_", season, ".xlsx"))
  )
  
  # -------------------------------------------------------------------
  # 4.5 Explore ROC performance across seroprevalence cutoffs
  # -------------------------------------------------------------------
  
  quantile_cutoffs <- as.numeric(
    quantile(zones_season$prevalence_global, probs = quantile_probs, na.rm = TRUE)
  )
  
  all_cutoffs <- sort(unique(c(biological_cutoffs, extra_cutoffs_fixed, quantile_cutoffs)))
  risk_vars <- c("risk_lit", "risk_expert", "risk_combined_mean")
  
  results_all_cutoffs <- compute_all_roc_results(
    data = zones_season,
    risk_vars = risk_vars,
    all_cutoffs = all_cutoffs
  )
  
  results_valid <- results_all_cutoffs %>%
    filter(n_pos >= min_class_size, n_neg >= min_class_size)
  
  results_all_export <- results_all_cutoffs %>%
    mutate(
      analysis_type = case_when(
        cutoff %in% biological_cutoffs ~ "biological_candidate",
        TRUE ~ "exploratory"
      ),
      auc = round(auc, 3),
      ci_low = round(ci_low, 3),
      ci_high = round(ci_high, 3)
    )
  
  results_valid_export <- results_valid %>%
    mutate(
      analysis_type = case_when(
        cutoff %in% biological_cutoffs ~ "biological_candidate",
        TRUE ~ "exploratory"
      ),
      auc = round(auc, 3),
      ci_low = round(ci_low, 3),
      ci_high = round(ci_high, 3)
    )
  
  write_xlsx(
    list(
      all_cutoffs = results_all_export,
      valid_cutoffs = results_valid_export
    ),
    path = here(season_dir, paste0("roc_results_all_cutoffs_", season, ".xlsx"))
  )
  
  # -------------------------------------------------------------------
  # 4.6 Compare 10% and 15% for the combined model
  # -------------------------------------------------------------------
  
  results_biological_combined <- results_valid %>%
    filter(
      model == "risk_combined_mean",
      cutoff %in% biological_cutoffs
    ) %>%
    arrange(desc(auc))
  
  # -------------------------------------------------------------------
  # 4.7 Select final cutoff for the combined model
  # -------------------------------------------------------------------
  
  # The final cutoff is selected as the best valid cutoff for the combined model.
  # In the manuscript, this was used to support the final surveillance-oriented
  # validation of the combined seasonal model.
  results_best_combined <- results_valid %>%
    filter(model == "risk_combined_mean") %>%
    arrange(desc(auc))
  
  final_cutoff <- results_best_combined$cutoff[1]
  
  # -------------------------------------------------------------------
  # 4.8 Final ROC models using the retained cutoff
  # -------------------------------------------------------------------
  
  zones_season$sero_binary <- ifelse(
    zones_season$prevalence_global >= final_cutoff,
    1, 0
  )
  
  roc_lit <- roc(zones_season$sero_binary, zones_season$risk_lit, quiet = TRUE)
  roc_expert <- roc(zones_season$sero_binary, zones_season$risk_expert, quiet = TRUE)
  roc_combined <- roc(zones_season$sero_binary, zones_season$risk_combined_mean, quiet = TRUE)
  
  final_auc_table <- data.frame(
    model = c("literature", "expert", "combined"),
    auc = c(
      as.numeric(auc(roc_lit)),
      as.numeric(auc(roc_expert)),
      as.numeric(auc(roc_combined))
    ),
    ci_low = c(
      as.numeric(ci.auc(roc_lit)[1]),
      as.numeric(ci.auc(roc_expert)[1]),
      as.numeric(ci.auc(roc_combined)[1])
    ),
    ci_high = c(
      as.numeric(ci.auc(roc_lit)[3]),
      as.numeric(ci.auc(roc_expert)[3]),
      as.numeric(ci.auc(roc_combined)[3])
    ),
    seroprevalence_cutoff = final_cutoff
  ) %>%
    mutate(
      auc = round(auc, 3),
      ci_low = round(ci_low, 3),
      ci_high = round(ci_high, 3)
    )
  
  write_xlsx(
    list(
      final_auc = final_auc_table,
      biological_comparison = results_biological_combined,
      best_combined_cutoffs = results_best_combined
    ),
    path = here(season_dir, paste0("final_auc_summary_", season, ".xlsx"))
  )
  
  # -------------------------------------------------------------------
  # 4.9 Derive low and high thresholds for combined risk classification
  # -------------------------------------------------------------------
  
  roc_coords <- coords(
    roc_combined,
    x = "all",
    input = "threshold",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  roc_coords <- as.data.frame(roc_coords)
  roc_coords <- roc_coords[is.finite(roc_coords$threshold), ]
  roc_coords$youden <- roc_coords$sensitivity + roc_coords$specificity - 1
  
  low_thresh <- roc_coords$threshold[which.max(roc_coords$youden)]
  high_thresh <- derive_high_threshold(
    roc_coords = roc_coords,
    low_thresh = low_thresh,
    target_specificity = target_specificity_for_high
  )
  
  zones_season$risk_cat <- cut(
    zones_season$risk_combined_mean,
    breaks = c(-Inf, low_thresh, high_thresh, Inf),
    labels = c("Low", "Medium", "High"),
    right = TRUE
  )
  
  thresholds_table <- data.frame(
    season = season,
    final_seroprevalence_cutoff = final_cutoff,
    low_threshold = low_thresh,
    high_threshold = high_thresh,
    target_specificity_for_high = target_specificity_for_high
  )
  
  write_xlsx(
    thresholds_table,
    path = here(season_dir, paste0("classification_thresholds_", season, ".xlsx"))
  )
  
  # -------------------------------------------------------------------
  # 4.10 Figure: average combined risk by division
  # -------------------------------------------------------------------
  
  p_avg_risk <- ggplot(zones_season) +
    geom_sf(aes(fill = risk_combined_mean), color = "gray30", linewidth = 0.2) +
    scale_fill_gradient(
      low = "white",
      high = "red",
      name = "Combined risk index"
    ) +
    labs(
      title = paste("Average combined risk by division (", map_title, ")", sep = "")
    ) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    filename = here(season_dir, paste0("average_combined_risk_by_division_", season, ".tiff")),
    plot = p_avg_risk,
    width = 10,
    height = 8,
    units = "in",
    dpi = 300,
    device = "tiff"
  )
  
  # -------------------------------------------------------------------
  # 4.11 Figure: risk categories map
  # -------------------------------------------------------------------
  
  zones_labels <- st_point_on_surface(zones_season)
  
  tmap_mode("plot")
  
  tm_risk <- tm_shape(zones_season) +
    tm_polygons(
      "risk_cat",
      palette = c("#fef0d9", "#fc8d59", "#b30045"),
      title = "Risk level",
      style = "cat"
    ) +
    tm_shape(zones_labels) +
    tm_text(
      "NAME_2",
      size = 0.65,
      col = "black",
      fontface = "bold",
      shadow = TRUE,
      remove.overlap = TRUE
    ) +
    tm_layout(
      title = map_title,
      title.size = 1.4,
      legend.outside = TRUE,
      legend.title.size = 1.2,
      legend.text.size = 1.0,
      frame = FALSE,
      inner.margins = c(0.02, 0.02, 0.08, 0.02)
    ) +
    tm_credits(
      paste0(
        "Seroprev. cut-off = ", round(final_cutoff, 2),
        "%; Low risk <= ", round(low_thresh, 3),
        "; High risk > ", round(high_thresh, 3)
      ),
      position = c("right", "BOTTOM"),
      size = 1.0
    ) +
    tm_compass(
      type = "arrow",
      position = c("left", "top"),
      size = 2
    ) +
    tm_scale_bar(
      position = c("left", "bottom"),
      text.size = 1.0,
      breaks = c(0, 50, 100),
      width = 0.2
    )
  
  tmap_save(
    tm = tm_risk,
    filename = here(season_dir, paste0("risk_categories_by_division_combined_", season, ".tiff")),
    width = 2000,
    height = 1600,
    units = "px",
    dpi = 300
  )
  
  # -------------------------------------------------------------------
  # 4.12 Figure: ROC curve for final combined model
  # -------------------------------------------------------------------
  
  save_roc_step_plot(
    roc_obj = roc_combined,
    title_text = paste0(
      "ROC curve (seroprevalence >= ",
      round(final_cutoff, 2),
      "%) - ",
      map_title
    ),
    output_file = here(season_dir, paste0("ROC_combined_", season, ".tiff"))
  )
  
  # -------------------------------------------------------------------
  # 4.13 Optional shapefile exports
  # -------------------------------------------------------------------
  
  risk_categories_shp <- zones_season %>%
    select(
      GID_0,
      COUNTRY,
      GID_1,
      NAME_1,
      NAME_2,
      dept_name,
      prevalence_global,
      risk_lit,
      risk_expert,
      risk_combined_mean,
      risk_cat
    )
  
  st_write(
    risk_categories_shp,
    here(season_dir, paste0("risk_categories_by_division_combined_", season, ".shp")),
    delete_layer = TRUE,
    quiet = TRUE
  )
  
  average_risk_shp <- zones_season %>%
    select(
      GID_0,
      COUNTRY,
      GID_1,
      NAME_1,
      NAME_2,
      dept_name,
      prevalence_global,
      risk_lit,
      risk_expert,
      risk_combined_mean
    )
  
  st_write(
    average_risk_shp,
    here(season_dir, paste0("average_combined_risk_by_division_", season, ".shp")),
    delete_layer = TRUE,
    quiet = TRUE
  )
  
  # -------------------------------------------------------------------
  # 4.14 Return objects
  # -------------------------------------------------------------------
  
  list(
    zones = zones_season,
    final_cutoff = final_cutoff,
    low_thresh = low_thresh,
    high_thresh = high_thresh,
    final_auc_table = final_auc_table,
    results_all_cutoffs = results_all_cutoffs,
    results_valid = results_valid
  )
}

# =====================================================================
# 5. RUN DRY AND WET VALIDATION
# =====================================================================

dry_results <- process_season("dry")
wet_results <- process_season("wet")

# =====================================================================
# 6. EXPORT COMBINED SEASON SUMMARY TABLE
# =====================================================================

combined_auc_summary <- bind_rows(
  dry_results$final_auc_table %>% mutate(season = "dry"),
  wet_results$final_auc_table %>% mutate(season = "wet")
) %>%
  select(season, everything())

write_xlsx(
  combined_auc_summary,
  path = here(output_dir, "combined_auc_summary_dry_wet.xlsx")
)

message("Validation and figure generation completed successfully.")
