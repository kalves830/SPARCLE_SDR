
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(purrr)

# =========================================================
# 1. INPUTS
# =========================================================

suit_paths <- list(
  Baseline = "C:/CI_Alves/suitability/rasters-20260706T203255Z-3-001/rasters/suitability_agroforestry.tif",
  SSP126 = "C:/CI_Alves/suitability/rasters-20260706T203255Z-3-001/rasters/future_feasibility_agroforestry_SSP126.tif",  # 1-5 scale
  SSP370 = "C:/CI_Alves/suitability/rasters-20260706T203255Z-3-001/rasters/future_feasibility_agroforestry_SSP370.tif",
  SSP585 = "C:/CI_Alves/suitability/rasters-20260706T203255Z-3-001/rasters/future_feasibility_agroforestry_SSP585.tif"
)

lulc_path <- "C:/CI_Alves/analysis/ssa_invest_runs/Lulc/lulc_ssa.tif"

out_dir <- "outputs/baseline_vs_ssp_suitability"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

top_pct_thresholds <- c(10, 15, 20, 25)

# pixel area only valid if raster is in a projected CRS with constant cell size
# change this if your resolution is different
pixel_area_km2 <- 0.09

lulc_labels <- tibble(
  lulc = c(-128, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95),
  lulc_name = c(
    "No data", "Tree cover", "Shrubland", "Grassland", "Cropland",
    "Sparse vegetation", "Water", "Bare", "Wetland / flooded",
    "Mangrove", "Urban"
  ),
  af_appropriate = c(
    FALSE, FALSE, TRUE, TRUE, TRUE,
    TRUE, FALSE, FALSE, FALSE, FALSE, FALSE
  )
)

# =========================================================
# 2. HELPERS
# =========================================================

extract_suit_lulc <- function(suit_r, lulc_r, scenario_name) {
  
  if (!compareGeom(suit_r, lulc_r, stopOnError = FALSE)) {
    lulc_r <- resample(lulc_r, suit_r, method = "near")
  }
  
  stk <- c(suit_r, lulc_r)
  names(stk) <- c("suitability", "lulc")
  
  as.data.frame(stk, na.rm = TRUE) |>
    left_join(lulc_labels, by = "lulc") |>
    mutate(scenario = scenario_name) |>
    filter(!is.na(lulc_name))
}

analyse_thresholds <- function(df, top_pcts = top_pct_thresholds) {
  
  suit_quantiles <- quantile(
    df$suitability,
    probs = 1 - top_pcts / 100,
    na.rm = TRUE
  )
  
  out <- vector("list", length(top_pcts))
  
  for (i in seq_along(top_pcts)) {
    
    pct <- top_pcts[i]
    suit_cut <- unname(suit_quantiles[i])
    
    subset_df <- df |>
      filter(suitability >= suit_cut)
    
    total_px <- nrow(subset_df)
    total_area_km2 <- total_px * pixel_area_km2
    
    lulc_summary <- subset_df |>
      count(lulc_name, af_appropriate, name = "n_pixels") |>
      mutate(
        area_km2 = n_pixels * pixel_area_km2,
        pct_of_total = 100 * n_pixels / total_px,
        top_pct_threshold = pct,
        suit_cutoff = suit_cut,
        total_area_km2 = total_area_km2
      )
    
    af_summary <- subset_df |>
      filter(af_appropriate) |>
      summarise(
        n_af_px = n(),
        af_area_km2 = n() * pixel_area_km2,
        pct_af_of_total = 100 * n() / total_px
      ) |>
      mutate(
        top_pct_threshold = pct,
        suit_cutoff = suit_cut,
        total_area_km2 = total_area_km2
      )
    
    out[[i]] <- list(
      lulc = lulc_summary,
      af   = af_summary
    )
  }
  
  out
}

# =========================================================
# 3. MAIN LOOP
# =========================================================

lulc_r <- rast(lulc_path)

all_lulc <- list()
all_af   <- list()

for (scen in names(suit_paths)) {
  
  message("Processing: ", scen)
  
  suit_r <- rast(suit_paths[[scen]])
  df <- extract_suit_lulc(suit_r, lulc_r, scen)
  res <- analyse_thresholds(df)
  
  lulc_rows <- bind_rows(lapply(res, `[[`, "lulc")) |>
    mutate(scenario = scen)
  
  af_rows <- bind_rows(lapply(res, `[[`, "af")) |>
    mutate(scenario = scen)
  
  all_lulc[[scen]] <- lulc_rows
  all_af[[scen]] <- af_rows
}

lulc_combined <- bind_rows(all_lulc)
af_combined   <- bind_rows(all_af)

write.csv(lulc_combined,
          file.path(out_dir, "baseline_vs_ssp_lulc_by_threshold.csv"),
          row.names = FALSE)

write.csv(af_combined,
          file.path(out_dir, "baseline_vs_ssp_af_area_by_threshold.csv"),
          row.names = FALSE)

# =========================================================
# 4. PLOTS
# =========================================================

lulc_colours <- c(
  "Cropland" = "#D4A017",
  "Grassland" = "#90EE90",
  "Shrubland" = "#8B6914",
  "Sparse vegetation" = "#D2B48C",
  "Tree cover" = "#228B22",
  "Water" = "#4682B4",
  "Bare" = "#A9A9A9",
  "Wetland / flooded" = "#20B2AA",
  "Mangrove" = "#006400",
  "Urban" = "#FF4500",
  "No data" = "#E0E0E0"
)

scenario_cols <- c(
  "Baseline" = "#4D4D4D",
  "SSP126"   = "#2166AC",
  "SSP370"   = "#FDAE61",
  "SSP585"   = "#D73027"
)

# ---------------------------------------------------------
# Plot 1. LULC composition in top thresholds
# ---------------------------------------------------------
p_stack <- lulc_combined |>
  mutate(
    top_pct_threshold = factor(
      paste0("Top ", top_pct_threshold, "%"),
      levels = paste0("Top ", sort(top_pct_thresholds), "%")
    )
  ) |>
  ggplot(aes(x = top_pct_threshold, y = pct_of_total,
             fill = lulc_name, alpha = af_appropriate)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = lulc_colours, name = "LULC class") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.45), guide = "none") +
  facet_wrap(~ scenario, ncol = 2) +
  scale_y_continuous(labels = label_percent(scale = 1)) +
  labs(
    title = "LULC composition within top suitability thresholds",
    subtitle = "Baseline compared with future SSP suitability scenarios",
    x = "Suitability threshold",
    y = "% of pixels within threshold",
    caption = "AF-appropriate classes: Cropland, Grassland, Shrubland, Sparse vegetation"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "fig1_baseline_vs_ssp_lulc_composition.png"),
       p_stack, width = 13, height = 8, dpi = 180)

# ---------------------------------------------------------
# Plot 2. AF-appropriate area by threshold
# ---------------------------------------------------------
p_af_area <- af_combined |>
  mutate(
    top_pct_threshold = factor(
      paste0("Top ", top_pct_threshold, "%"),
      levels = paste0("Top ", sort(top_pct_thresholds), "%")
    )
  ) |>
  ggplot(aes(x = top_pct_threshold, y = af_area_km2 / 1e6, fill = scenario)) +
  geom_bar(stat = "identity", position = "dodge", colour = "white") +
  scale_fill_manual(values = scenario_cols, name = "Scenario") +
  scale_y_continuous(labels = label_comma(suffix = "M km²")) +
  labs(
    title = "AF-appropriate land area within top suitability thresholds",
    subtitle = "Baseline vs future SSP scenarios",
    x = "Suitability threshold",
    y = "Area (million km²)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "fig2_baseline_vs_ssp_af_area.png"),
       p_af_area, width = 10, height = 6, dpi = 180)

# ---------------------------------------------------------
# Plot 3. AF-appropriate share of threshold area
# ---------------------------------------------------------
p_af_pct <- af_combined |>
  mutate(
    top_pct_threshold = factor(
      paste0("Top ", top_pct_threshold, "%"),
      levels = paste0("Top ", sort(top_pct_thresholds), "%")
    )
  ) |>
  ggplot(aes(x = scenario, y = pct_af_of_total, fill = scenario)) +
  geom_bar(stat = "identity", colour = "white") +
  geom_text(aes(label = paste0(round(pct_af_of_total, 1), "%")),
            vjust = -0.35, size = 3.4) +
  scale_fill_manual(values = scenario_cols, guide = "none") +
  facet_wrap(~ top_pct_threshold, ncol = 4) +
  scale_y_continuous(limits = c(0, 110), labels = label_percent(scale = 1)) +
  labs(
    title = "AF-appropriate land as % of total pixels within each threshold",
    subtitle = "Comparison of current baseline and future SSP suitability maps",
    x = "Scenario",
    y = "% AF-appropriate"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "fig3_baseline_vs_ssp_af_share.png"),
       p_af_pct, width = 12, height = 6, dpi = 180)

# ---------------------------------------------------------
# Plot 4. Heatmap: AF-appropriate area by LULC class
# ---------------------------------------------------------
p_heatmap <- lulc_combined |>
  filter(af_appropriate) |>
  mutate(
    top_pct_threshold = factor(
      paste0("Top ", top_pct_threshold, "%"),
      levels = paste0("Top ", sort(top_pct_thresholds), "%")
    ),
    lulc_name = factor(
      lulc_name,
      levels = c("Cropland", "Grassland", "Shrubland", "Sparse vegetation")
    )
  ) |>
  ggplot(aes(x = top_pct_threshold, y = scenario, fill = area_km2 / 1000)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = comma(round(area_km2 / 1000))),
            size = 2.8, colour = "black") +
  scale_fill_distiller(
    palette = "YlOrRd",
    direction = 1,
    name = "Area\n('000 km²)",
    labels = label_comma()
  ) +
  facet_wrap(~ lulc_name, ncol = 2) +
  labs(
    title = "AF-appropriate area by LULC class, threshold, and scenario",
    subtitle = "Baseline included alongside SSP scenarios",
    x = "Suitability threshold",
    y = "Scenario"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(file.path(out_dir, "fig4_baseline_vs_ssp_heatmap.png"),
       p_heatmap, width = 12, height = 8, dpi = 180)

