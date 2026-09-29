# =============================================================================
# SECTION 14 — AF3-CLASS INTEGRATION
# Replaces the cropland-only AF transition runs in Comparisons A, B, C with
# the AF3-class (crop+grass+shrub) runs. All downstream plots (fig1-fig5)
# now reflect the AF3-class transition instead of cropland-only.
# =============================================================================
# Requires: invest_output_root, run_meta (with AF3 rows appended), ws_metrics,
# out14, load_watershed(), compare_runs(), summarise_comp(), pct_change(),
# climate_cols, thresh_cols — all defined in the fixed Section 14 script.
# =============================================================================

# -----------------------------------------------------------------------
# Reference run IDs (no-AF baselines) — unchanged, used by all 3 comparisons
# -----------------------------------------------------------------------

ref_baseline <- "baseline_cropland_0.34"
ref_ssp370   <- "ssp370_cropland_0.34"

# =============================================================================
# COMPARISON A — AF3-CLASS EFFECT BY SUITABILITY THRESHOLD (same climate)
# OVERWRITES: previously compared cropland-only AF runs vs no-AF reference.
# Now compares AF3-class (crop+grass+shrub) runs vs the same no-AF reference.
# =============================================================================

message("Block A — AF3-class effect by threshold (crop+grass+shrub)...")

comp_A_pairs <- list(
  list(ref = ref_baseline, scen = "baseline_af3_top10", group = "Baseline"),
  list(ref = ref_baseline, scen = "baseline_af3_top20", group = "Baseline"),
  list(ref = ref_baseline, scen = "baseline_af3_top25", group = "Baseline"),
  list(ref = ref_ssp370,   scen = "ssp370_af3_top10",   group = "SSP370"),
  list(ref = ref_ssp370,   scen = "ssp370_af3_top20",   group = "SSP370"),
  list(ref = ref_ssp370,   scen = "ssp370_af3_top25",   group = "SSP370")
)

comp_A_long <- map_dfr(comp_A_pairs, function(p) {
  comp <- compare_runs(p$ref, p$scen)
  if (is.null(comp)) return(NULL)
  s <- summarise_comp(comp, p$ref, p$scen)
  s$group     <- p$group
  s$threshold <- run_meta$threshold[run_meta$run_id == p$scen]
  s
})

write.csv(comp_A_long, file.path(out14, "A_af3_threshold_effect.csv"), row.names = FALSE)

# =============================================================================
# COMPARISON B — CLIMATE PENALTY (no AF, baseline vs SSP370)
# UNCHANGED: this comparison never touched AF at all — kept identical so the
# climate-only signal remains a stable reference point for Plot 3.
# =============================================================================

message("Block B — Climate penalty (unchanged, no AF)...")

comp_B <- compare_runs(ref_baseline, ref_ssp370)
if (!is.null(comp_B)) {
  comp_B_summary <- summarise_comp(comp_B, ref_baseline, ref_ssp370) |>
    mutate(group = "Climate penalty (no AF)")
  write.csv(comp_B_summary, file.path(out14, "B_climate_penalty.csv"), row.names = FALSE)
} else {
  message("  Block B skipped — one or both watershed shapefiles missing.")
}

# =============================================================================
# COMPARISON C — COMBINED: SSP370 AF3-class runs vs Baseline no-AF
# OVERWRITES: previously used SSP370 cropland-only AF runs.
# Now answers: does AF3-class (crop+grass+shrub) transition under SSP370
# offset the climate penalty relative to today's no-AF baseline?
# =============================================================================

message("Block C — Combined climate + AF3-class effect...")

comp_C_pairs <- list(
  list(ref = ref_baseline, scen = "ssp370_af3_top10"),
  list(ref = ref_baseline, scen = "ssp370_af3_top20"),
  list(ref = ref_baseline, scen = "ssp370_af3_top25")
)

comp_C_long <- map_dfr(comp_C_pairs, function(p) {
  comp <- compare_runs(p$ref, p$scen)
  if (is.null(comp)) return(NULL)
  s <- summarise_comp(comp, p$ref, p$scen)
  s$threshold <- run_meta$threshold[run_meta$run_id == p$scen]
  s
})

write.csv(comp_C_long, file.path(out14, "C_combined_climate_af3_effect.csv"), row.names = FALSE)

# =============================================================================
# PIXEL-LEVEL DIFFERENCE MAPS — usle.tif and sed_export.tif
# OVERWRITES: run_ids updated to point at AF3-class output folders instead
# of the cropland-only AF folders.
# =============================================================================

message("Computing pixel-level difference rasters (AF3-class)...")

diff_pairs <- list(
  list(ref = ref_baseline, scen = "baseline_af3_top10", layer = "usle"),
  list(ref = ref_baseline, scen = "baseline_af3_top20", layer = "usle"),
  list(ref = ref_baseline, scen = "baseline_af3_top25", layer = "usle"),
  list(ref = ref_baseline, scen = "baseline_af3_top10", layer = "sed_export"),
  list(ref = ref_baseline, scen = "baseline_af3_top20", layer = "sed_export"),
  list(ref = ref_baseline, scen = "baseline_af3_top25", layer = "sed_export"),
  list(ref = ref_baseline, scen = ref_ssp370,            layer = "usle"),
  list(ref = ref_baseline, scen = ref_ssp370,            layer = "sed_export"),
  list(ref = ref_baseline, scen = "ssp370_af3_top10",   layer = "usle"),
  list(ref = ref_baseline, scen = "ssp370_af3_top25",   layer = "usle")
)

raster_file <- list(
  usle       = "usle.tif",
  sed_export = "sed_export.tif"
)

diff_dir <- file.path(out14, "difference_rasters_af3")
dir.create(diff_dir, recursive = TRUE, showWarnings = FALSE)

for (dp in diff_pairs) {
  ref_path  <- file.path(invest_output_root, dp$ref,  raster_file[[dp$layer]])
  scen_path <- file.path(invest_output_root, dp$scen, raster_file[[dp$layer]])
  if (!file.exists(ref_path) || !file.exists(scen_path)) {
    message("  Skipping diff (files not found): ", dp$scen, " vs ", dp$ref,
            " [", dp$layer, "]")
    next
  }
  r_ref  <- rast(ref_path)
  r_scen <- rast(scen_path)
  if (!compareGeom(r_ref, r_scen, stopOnError = FALSE)) {
    r_scen <- resample(r_scen, r_ref, method = "bilinear")
  }
  r_diff   <- r_scen - r_ref
  out_name <- paste0("diff_", dp$layer, "_", dp$scen, "_minus_", dp$ref, ".tif")
  out_path <- file.path(diff_dir, out_name)
  writeRaster(r_diff, out_path, datatype = "FLT4S", NAflag = -9999,
              overwrite = TRUE)
  message("  Written: ", out_name)
  rm(r_ref, r_scen, r_diff); gc()
}

# =============================================================================
# PLOTS — regenerated from comp_A_long / comp_B_summary / comp_C_long, now
# reflecting AF3-class (crop+grass+shrub) instead of cropland-only.
# File names suffixed "_af3" so cropland-only cropland outputs (if you kept
# them from a prior run) are not overwritten.
# =============================================================================

message("Generating AF3-class plots...")

thr_levels <- c("Full area", "Top 10%", "Top 20%", "Top 25%")

# ---------------------------------------------------------------------------
# Plot 1 (AF3): sed_export and usle_tot — mean % change by threshold, climate
# ---------------------------------------------------------------------------
p1_af3 <- comp_A_long |>
  filter(metric %in% c("sed_export", "usle_tot")) |>
  mutate(
    threshold = factor(threshold, levels = thr_levels),
    metric_lab = recode(metric,
                        sed_export = "Sediment Export (sed_export)",
                        usle_tot   = "Total Soil Loss (usle_tot)"
    )
  ) |>
  ggplot(aes(x = threshold, y = mean_pct,
             ymin = mean_pct - sd_pct, ymax = mean_pct + sd_pct,
             fill = group, colour = group)) +
  geom_col(position = position_dodge(0.75), width = 0.65, alpha = 0.85) +
  geom_errorbar(position = position_dodge(0.75), width = 0.3, linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = climate_cols, name = "Climate") +
  scale_colour_manual(values = climate_cols, name = "Climate") +
  scale_y_continuous(labels = label_percent(scale = 1, suffix = "%")) +
  facet_wrap(~ metric_lab, scales = "free_y", ncol = 2) +
  labs(
    title = "AF3-class (crop+grass+shrub) transition effect by suitability threshold",
    subtitle = "% change vs no-AF reference run (mean \u00b1 SD across watersheds)",
    x = "Land targeted for AF3-class transition",
    y = "Mean % change in metric",
    caption = "Negative = AF3 reduces erosion/export (beneficial)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(file.path(out14, "fig1_af3_threshold_effect_sed_usle.png"),
       p1_af3, width = 12, height = 6, dpi = 180)

# ---------------------------------------------------------------------------
# Plot 2 (AF3): all metrics, faceted
# ---------------------------------------------------------------------------
p2_af3 <- comp_A_long |>
  mutate(
    threshold = factor(threshold, levels = thr_levels),
    metric = factor(metric, levels = c("usle_tot", "sed_export", "avoid_eros", "avoid_exp"))
  ) |>
  ggplot(aes(x = threshold, y = mean_pct, fill = group)) +
  geom_col(position = position_dodge(0.75), width = 0.65, alpha = 0.85) +
  geom_errorbar(
    aes(ymin = mean_pct - sd_pct, ymax = mean_pct + sd_pct, colour = group),
    position = position_dodge(0.75), width = 0.3, linewidth = 0.5
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = climate_cols, name = "Climate") +
  scale_colour_manual(values = climate_cols, name = "Climate") +
  scale_y_continuous(labels = label_percent(scale = 1)) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(
    title = "AF3-class effect on all SDR metrics by suitability threshold",
    subtitle = "% change vs no-AF reference; mean \u00b1 SD across watersheds",
    x = "Suitability threshold", y = "Mean % change"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(file.path(out14, "fig2_af3_all_metrics_threshold.png"),
       p2_af3, width = 13, height = 9, dpi = 180)

# ---------------------------------------------------------------------------
# Plot 3 (AF3): climate penalty vs AF3 offset on sed_export
# ---------------------------------------------------------------------------
combined_sed_af3 <- bind_rows(
  if (!is.null(comp_B) && exists("comp_B_summary")) {
    comp_B_summary |>
      filter(metric == "sed_export") |>
      mutate(label = "SSP370 \u2014 No AF\n(climate penalty)", colour_group = "Climate penalty")
  },
  comp_C_long |>
    filter(metric == "sed_export") |>
    mutate(label = paste0("SSP370 \u2014 AF3 ", threshold), colour_group = "AF3 offsets climate")
) |>
  mutate(
    label = factor(label, levels = c(
      "SSP370 \u2014 No AF\n(climate penalty)",
      "SSP370 \u2014 AF3 Top 10%",
      "SSP370 \u2014 AF3 Top 20%",
      "SSP370 \u2014 AF3 Top 25%"
    ))
  )

p3_af3 <- ggplot(combined_sed_af3,
                 aes(x = label, y = mean_pct, fill = colour_group,
                     ymin = mean_pct - sd_pct, ymax = mean_pct + sd_pct)) +
  geom_col(width = 0.65, alpha = 0.9) +
  geom_errorbar(width = 0.3, linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(
    values = c("Climate penalty" = "#D73027", "AF3 offsets climate" = "#2166AC"),
    name = NULL
  ) +
  scale_y_continuous(labels = label_percent(scale = 1)) +
  labs(
    title = "Climate penalty vs AF3-class offset on sediment export",
    subtitle = "All relative to Baseline no-AF reference | mean \u00b1 SD across watersheds",
    x = NULL, y = "% change in sed_export vs Baseline no-AF",
    caption = "Negative values = less sediment export than Baseline no-AF"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 15, hjust = 1))

ggsave(file.path(out14, "fig3_af3_climate_vs_offset_sedexport.png"),
       p3_af3, width = 10, height = 6, dpi = 180)

# ---------------------------------------------------------------------------
# Plot 4 (AF3): watershed-level distribution of % change in sed_export
# ---------------------------------------------------------------------------
ws_pct_df_af3 <- map_dfr(comp_A_pairs, function(p) {
  comp <- compare_runs(p$ref, p$scen)
  if (is.null(comp) || !"sed_export_pct" %in% names(comp)) return(NULL)
  tibble(
    group = p$group,
    threshold = run_meta$threshold[run_meta$run_id == p$scen],
    pct_chg = comp$sed_export_pct
  )
}) |>
  mutate(threshold = factor(threshold, levels = thr_levels))

p4_af3 <- ggplot(ws_pct_df_af3,
                 aes(x = threshold, y = pct_chg, fill = group, colour = group)) +
  geom_violin(position = position_dodge(0.8), alpha = 0.35, linewidth = 0.4,
              scale = "width", trim = TRUE) +
  geom_boxplot(position = position_dodge(0.8), width = 0.25,
               outlier.size = 0.8, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = climate_cols, name = "Climate") +
  scale_colour_manual(values = climate_cols, name = "Climate") +
  scale_y_continuous(labels = label_percent(scale = 1)) +
  coord_cartesian(ylim = c(
    max(-100, quantile(ws_pct_df_af3$pct_chg, 0.01, na.rm = TRUE) * 1.2),
    min(100,  quantile(ws_pct_df_af3$pct_chg, 0.99, na.rm = TRUE) * 1.2)
  )) +
  labs(
    title = "Distribution of watershed-level sediment export change (AF3-class)",
    subtitle = "AF3 (crop+grass+shrub) transition vs no-AF reference by suitability threshold",
    x = "Land targeted for AF3-class transition",
    y = "% change in sed_export per watershed",
    caption = "Violin = full distribution | Box = IQR + median | Clipped to 1st-99th %ile"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(out14, "fig4_af3_watershed_sedexport_distribution.png"),
       p4_af3, width = 11, height = 7, dpi = 180)

message("AF3-class integration complete.")
message("  A_af3_threshold_effect.csv, B_climate_penalty.csv, C_combined_climate_af3_effect.csv")
message("  fig1_af3_..., fig2_af3_..., fig3_af3_..., fig4_af3_... written to: ", out14) 

# =============================================================================
# COMPARISON D — AF3-CLASS EFFECT (crop+grass+shrub) vs single-class cropland AF
# Isolates the incremental benefit of adding grassland + shrubland transitions
# on top of the cropland-only AF already assessed in Comparison A
# =============================================================================

message("Block D — AF3-class effect (crop+grass+shrub) vs cropland-only AF and vs no-AF...")

comp_D_pairs <- list(
  # vs no-AF reference (same climate) — total AF3 benefit
  list(ref = "baseline_cropland_0.34", scen = "baseline_af3_top10", group = "Baseline", vs = "noAF"),
  list(ref = "baseline_cropland_0.34", scen = "baseline_af3_top20", group = "Baseline", vs = "noAF"),
  list(ref = "baseline_cropland_0.34", scen = "baseline_af3_top25", group = "Baseline", vs = "noAF"),
  list(ref = "ssp370_cropland_0.34",   scen = "ssp370_af3_top10",   group = "SSP370",   vs = "noAF"),
  list(ref = "ssp370_cropland_0.34",   scen = "ssp370_af3_top20",   group = "SSP370",   vs = "noAF"),
  list(ref = "ssp370_cropland_0.34",   scen = "ssp370_af3_top25",   group = "SSP370",   vs = "noAF"),
  # vs cropland-only AF (same climate, same threshold) — incremental grass+shrub benefit
  list(ref = "baseline_cropland_0.34_to_0.2832_top10", scen = "baseline_af3_top10", group = "Baseline", vs = "croplandAF"),
  list(ref = "baseline_cropland_0.34_to_0.2832_top20", scen = "baseline_af3_top20", group = "Baseline", vs = "croplandAF"),
  list(ref = "baseline_cropland_0.34_to_0.2832_top25", scen = "baseline_af3_top25", group = "Baseline", vs = "croplandAF"),
  list(ref = "ssp370_cropland_0.34_to_0.2832_top10",   scen = "ssp370_af3_top10",   group = "SSP370",   vs = "croplandAF"),
  list(ref = "ssp370_cropland_0.34_to_0.2832_top20",   scen = "ssp370_af3_top20",   group = "SSP370",   vs = "croplandAF"),
  list(ref = "ssp370_cropland_0.34_to_0.2832_top25",   scen = "ssp370_af3_top25",   group = "SSP370",   vs = "croplandAF")
)

comp_D_long <- map_dfr(comp_D_pairs, function(p) {
  comp <- compare_runs(p$ref, p$scen)
  if (is.null(comp)) return(NULL)
  s <- summarise_comp(comp, p$ref, p$scen)
  s$group <- p$group
  s$vs <- p$vs
  s$threshold <- run_meta$threshold[run_meta$run_id == p$scen]
  s
})

write.csv(comp_D_long, file.path(out14, "D_af3class_effect.csv"), row.names = FALSE)


cat("\n========================================\n")
cat("SECTION 14 COMPLETE — InVEST Impact Analysis\n")
cat("Outputs written to:", out14, "\n\n")
cat("CSVs:\n")
cat("  watershed_all_runs_long.csv         — raw watershed data, all runs\n")
cat("  A_af_threshold_effect.csv           — AF vs no-AF by threshold\n")
cat("  B_climate_penalty.csv               — SSP370 vs Baseline, no AF\n")
cat("  C_combined_climate_af_effect.csv    — SSP370 AF runs vs Baseline no-AF\n")
cat("  MASTER_impact_summary.csv           — all comparisons combined\n\n")
cat("Figures:\n")
cat("  fig1_af_threshold_effect_sed_usle.png   — mean % change, sed+usle\n")
cat("  fig2_all_metrics_af_threshold.png       — all 4 metrics by threshold\n")
cat("  fig3_climate_vs_af_offset_sedexport.png — penalty vs offset (sed_export)\n")
cat("  fig4_watershed_sedexport_distribution.png — violin/box distributions\n")
cat("  fig5_marginal_benefit_curve.png         — cost curve: 10→20→25%\n\n")
cat("Difference rasters (difference_rasters/):\n")
cat("  diff_usle_<scen>_minus_<ref>.tif\n")
cat("  diff_sed_export_<scen>_minus_<ref>.tif\n")
cat("========================================\n") 


# ---------------------------------------------------------------------------
# Plot 6: AF3-class total benefit vs cropland-only AF — sed_export & usle_tot
# ---------------------------------------------------------------------------
p6 <- comp_D_long |>
  filter(metric %in% c("sed_export", "usle_tot")) |>
  mutate(
    threshold = factor(threshold, levels = thr_levels),
    metric_lab = recode(metric, sed_export = "Sediment Export", usle_tot = "Total Soil Loss"),
    vs_lab = recode(vs, noAF = "vs No-AF Reference", croplandAF = "vs Cropland-Only AF")
  ) |>
  ggplot(aes(x = threshold, y = mean_pct, fill = group,
             ymin = mean_pct - sd_pct, ymax = mean_pct + sd_pct)) +
  geom_col(position = position_dodge(0.75), width = 0.65, alpha = 0.85) +
  geom_errorbar(position = position_dodge(0.75), width = 0.3, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = climate_cols, name = "Climate") +
  scale_y_continuous(labels = label_percent(scale = 1)) +
  facet_grid(metric_lab ~ vs_lab, scales = "free_y") +
  labs(
    title = "AF3-class (crop+grassland+shrubland) transition effect",
    subtitle = "Total benefit vs no-AF, and incremental benefit vs cropland-only AF",
    x = "Suitability threshold", y = "Mean % change",
    caption = "Negative = AF3-class reduces erosion/export further"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(file.path(out14, "fig6_af3class_effect.png"), p6, width = 12, height = 8, dpi = 180)

cat("\nSection 14 extended with AF3-class comparisons.\n")
cat("  D_af3class_effect.csv — total + incremental benefit of grass/shrub transitions\n")
cat("  fig6_af3class_effect.png\n")



# =============================================================================
# SECTION 15 (REWRITTEN) — SPATIAL MAPS WITH OUTLIER FILTER INTEGRATED
# Replaces BOTH prior Section 15 blocks (the choropleth/bivariate block and
# the separate "AF3-class vs Baseline" block, which duplicated map logic).
#
# Requires from Section 14 (fixed): invest_output_root, ref_baseline,
# run_meta, out14, load_watershed(), pct_change(), compare_runs(),
# compare_runs_robust(), min_baseline_by_metric.
# =============================================================================

library(sf)
library(classInt)
library(scales)

out15 <- file.path(out_dir, "spatial_benefit_maps")
dir.create(out15, recursive = TRUE, showWarnings = FALSE)

map_metric <- "sed_export"   # change to "usle_tot" for total soil loss maps
cutoff_used <- min_baseline_by_metric[[map_metric]]

# -----------------------------------------------------------------------
# Load watershed geometry once (identical polygons across all runs)
# -----------------------------------------------------------------------

ws_geom_path <- file.path(invest_output_root, ref_baseline, "watershed_results_sdr.shp")
stopifnot(file.exists(ws_geom_path))

ws_geom <- st_read(ws_geom_path, quiet = TRUE)
names(ws_geom) <- tolower(names(ws_geom))
if (!"ws_id" %in% names(ws_geom)) ws_geom$ws_id <- seq_len(nrow(ws_geom))
ws_geom <- ws_geom[, "ws_id"]

# -----------------------------------------------------------------------
# All run pairs to map — cropland-only AND AF3-class, both climates,
# all three thresholds. Consolidates what was previously two separate
# pair lists (comp_map1/map2_pairs and af3_map_pairs).
# -----------------------------------------------------------------------

map_pairs <- list(
  list(ref = "baseline_cropland_0.34", scen = "baseline_cropland_0.34_to_0.2832_top10", climate = "Baseline", threshold = "Top 10%", class = "Cropland-only"),
  list(ref = "baseline_cropland_0.34", scen = "baseline_cropland_0.34_to_0.2832_top20", climate = "Baseline", threshold = "Top 20%", class = "Cropland-only"),
  list(ref = "baseline_cropland_0.34", scen = "baseline_cropland_0.34_to_0.2832_top25", climate = "Baseline", threshold = "Top 25%", class = "Cropland-only"),
  list(ref = "ssp370_cropland_0.34",   scen = "ssp370_cropland_0.34_to_0.2832_top10",   climate = "SSP370",   threshold = "Top 10%", class = "Cropland-only"),
  list(ref = "ssp370_cropland_0.34",   scen = "ssp370_cropland_0.34_to_0.2832_top20",   climate = "SSP370",   threshold = "Top 20%", class = "Cropland-only"),
  list(ref = "ssp370_cropland_0.34",   scen = "ssp370_cropland_0.34_to_0.2832_top25",   climate = "SSP370",   threshold = "Top 25%", class = "Cropland-only"),
  list(ref = "baseline_cropland_0.34", scen = "baseline_af3_top10", climate = "Baseline", threshold = "Top 10%", class = "AF3-class"),
  list(ref = "baseline_cropland_0.34", scen = "baseline_af3_top20", climate = "Baseline", threshold = "Top 20%", class = "AF3-class"),
  list(ref = "baseline_cropland_0.34", scen = "baseline_af3_top25", climate = "Baseline", threshold = "Top 25%", class = "AF3-class"),
  list(ref = "ssp370_cropland_0.34",   scen = "ssp370_af3_top10",   climate = "SSP370",   threshold = "Top 10%", class = "AF3-class"),
  list(ref = "ssp370_cropland_0.34",   scen = "ssp370_af3_top20",   climate = "SSP370",   threshold = "Top 20%", class = "AF3-class"),
  list(ref = "ssp370_cropland_0.34",   scen = "ssp370_af3_top25",   climate = "SSP370",   threshold = "Top 25%", class = "AF3-class")
)

# -----------------------------------------------------------------------
# Build one merged sf dataset for ALL pairs at once, tagged with the
# low_baseline flag from compare_runs_robust(). Single source of truth
# replacing the old duplicate build_af3_map_data()/compare_runs() paths.
# -----------------------------------------------------------------------

build_map_row <- function(p) {
  comp <- compare_runs_robust(p$ref, p$scen)
  if (is.null(comp)) {
    message("  Skipping (missing data): ", p$scen, " vs ", p$ref)
    return(NULL)
  }
  pct_col   <- paste0(map_metric, "_pct")
  ref_col   <- paste0(map_metric, "_ref")
  scen_col  <- paste0(map_metric, "_scen")
  delta_col <- paste0(map_metric, "_delta")
  flag_col  <- paste0(map_metric, "_low_baseline")
  if (!all(c(pct_col, ref_col, scen_col, delta_col, flag_col) %in% names(comp))) {
    warning("Metric columns missing for ", p$scen); return(NULL)
  }
  tibble(
    ws_id        = comp$ws_id,
    pct_chg      = comp[[pct_col]],
    delta        = comp[[delta_col]],
    ref_val      = comp[[ref_col]],
    scen_val     = comp[[scen_col]],
    low_baseline = comp[[flag_col]],
    climate      = p$climate,
    threshold    = p$threshold,
    class        = p$class,
    ref_id       = p$ref,
    scen_id      = p$scen
  )
}

map_data_long <- map_dfr(map_pairs, build_map_row) |>
  mutate(
    threshold = factor(threshold, levels = c("Top 10%", "Top 20%", "Top 25%")),
    climate   = factor(climate, levels = c("Baseline", "SSP370")),
    class     = factor(class, levels = c("Cropland-only", "AF3-class")),
    # Category used for map fill: NA out excluded watersheds so they render
    # as a distinct grey "excluded" bucket instead of distorting the scale
    pct_chg_plot = ifelse(low_baseline, NA_real_, pct_chg)
  )

write.csv(map_data_long, file.path(out15, paste0("map_data_", map_metric, "_all_runs.csv")),
          row.names = FALSE)

map_sf_all <- ws_geom |>
  left_join(map_data_long, by = "ws_id") |>
  st_as_sf()

cat("Excluded watersheds per run (low baseline, threshold =", cutoff_used, map_metric, "):\n")
excl_summary <- map_data_long |>
  group_by(class, climate, threshold) |>
  summarise(n_total = n(), n_excluded = sum(low_baseline, na.rm = TRUE),
            pct_excluded = round(100 * n_excluded / n_total, 1), .groups = "drop")
print(excl_summary)
write.csv(excl_summary, file.path(out15, "exclusion_summary_by_run.csv"), row.names = FALSE)

# =============================================================================
# MAP 1 — % change choropleth, AF3-class only, Baseline top25 (headline map)
# Colour limits now computed from FILTERED data only.
# =============================================================================

map1_sf <- map_sf_all |> filter(class == "AF3-class", climate == "Baseline", threshold == "Top 25%")

pct_limit_1 <- max(abs(quantile(map1_sf$pct_chg_plot, c(0.02, 0.98), na.rm = TRUE)))

p_map1 <- ggplot(map1_sf) +
  geom_sf(aes(fill = pct_chg_plot), colour = "white", linewidth = 0.15) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
    limits = c(-pct_limit_1, pct_limit_1), oob = scales::squish,
    na.value = "grey80",
    name = "% change\nsed_export",
    labels = label_percent(scale = 1)
  ) +
  labs(
    title = "Sediment export change \u2014 AF3-class transition (crop+grass+shrub)",
    subtitle = "Baseline climate, top 25% suitability, vs no-AF reference",
    caption = paste0("Grey = excluded (baseline ", map_metric, " < ", cutoff_used, " t/yr, unreliable % denominator) | ",
                     "Colour scale from filtered watersheds only, clipped to 2nd\u201398th pct")
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), legend.position = "right") +
  coord_sf(datum = NA)

ggsave(file.path(out15, "map1_af3_sedexport_change_baseline_top25_filtered.png"),
       p_map1, width = 10, height = 9, dpi = 220)

# =============================================================================
# MAP 2 — Small multiples: % change, AF3-class, all thresholds x climates
# =============================================================================

map2_sf <- map_sf_all |> filter(class == "AF3-class")
pct_limit_2 <- max(abs(quantile(map2_sf$pct_chg_plot, c(0.02, 0.98), na.rm = TRUE)))

p_map2 <- ggplot(map2_sf) +
  geom_sf(aes(fill = pct_chg_plot), colour = "white", linewidth = 0.1) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
    limits = c(-pct_limit_2, pct_limit_2), oob = scales::squish,
    na.value = "grey80",
    name = "% change\nsed_export", labels = label_percent(scale = 1)
  ) +
  facet_grid(climate ~ threshold) +
  labs(
    title = "Sediment export change \u2014 AF3-class transition across thresholds",
    subtitle = "Common colour scale (filtered watersheds only) for direct visual comparison",
    caption = paste0("Grey = excluded (baseline < ", cutoff_used, " t/yr) | See exclusion_summary_by_run.csv for counts per panel")
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  coord_sf(datum = NA)

ggsave(file.path(out15, "map2_af3_sedexport_smallmultiples_filtered.png"),
       p_map2, width = 13, height = 9, dpi = 220)

# =============================================================================
# MAP 2B (NEW) — Same small multiples for Cropland-only, for direct visual
# comparison against Map 2's AF3-class panels using the SAME colour scale
# =============================================================================

map2b_sf <- map_sf_all |> filter(class == "Cropland-only")

p_map2b <- ggplot(map2b_sf) +
  geom_sf(aes(fill = pct_chg_plot), colour = "white", linewidth = 0.1) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
    limits = c(-pct_limit_2, pct_limit_2), oob = scales::squish,   # SAME scale as Map 2
    na.value = "grey80",
    name = "% change\nsed_export", labels = label_percent(scale = 1)
  ) +
  facet_grid(climate ~ threshold) +
  labs(
    title = "Sediment export change \u2014 Cropland-only transition across thresholds",
    subtitle = "Same colour scale as AF3-class map above \u2014 directly comparable magnitude of benefit",
    caption = paste0("Grey = excluded (baseline < ", cutoff_used, " t/yr)")
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  coord_sf(datum = NA)

ggsave(file.path(out15, "map2b_cropland_sedexport_smallmultiples_filtered.png"),
       p_map2b, width = 13, height = 9, dpi = 220)

# =============================================================================
# MAP 3 — Exclusion-summary panel: what fraction of each map is masked
# =============================================================================

p_exclusion <- ggplot(excl_summary, aes(x = threshold, y = pct_excluded, fill = class)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_text(aes(label = paste0(pct_excluded, "%")),
            position = position_dodge(0.7), vjust = -0.4, size = 3.2) +
  facet_wrap(~ climate) +
  scale_fill_manual(values = c("Cropland-only" = "#2166AC", "AF3-class" = "#D73027")) +
  labs(
    title = "Share of watersheds excluded from percent-change maps",
    subtitle = paste0("Excluded where baseline ", map_metric, " < ", cutoff_used, " t/yr"),
    x = "Suitability threshold", y = "% of watersheds excluded", fill = "Transition"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(out15, "map3_exclusion_share_by_run.png"),
       p_exclusion, width = 9, height = 6, dpi = 180)

# =============================================================================
# MAP 4 — Bivariate raster: magnitude (Baseline) x agreement (SSP370)
# Terciles now classified on FILTERED watersheds only; excluded watersheds
# get their own explicit code (0) rather than contaminating tercile breaks
# =============================================================================

message("Building bivariate raster: magnitude x cross-climate agreement (filtered)...")

usle_ref_path <- file.path(invest_output_root, ref_baseline, "usle.tif")
template_r <- rast(usle_ref_path)

magA_sf <- map_sf_all |> filter(class == "AF3-class", climate == "Baseline", threshold == "Top 25%")
magB_sf <- map_sf_all |> filter(class == "AF3-class", climate == "SSP370",  threshold == "Top 25%")

magA_r <- rasterize(vect(magA_sf), template_r, field = "pct_chg_plot")
magB_r <- rasterize(vect(magB_sf), template_r, field = "pct_chg_plot")

classify_tercile <- function(r) {
  vals <- values(r, na.rm = TRUE)
  brks <- quantile(vals, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  brks[1] <- brks[1] - 1e-6
  classify(r, rcl = matrix(c(
    brks[1], brks[2], 1,
    brks[2], brks[3], 2,
    brks[3], brks[4], 3
  ), ncol = 3, byrow = TRUE))
}

classA_r <- classify_tercile(magA_r)
classB_r <- classify_tercile(magB_r)

bivar_r <- (classA_r - 1) * 3 + classB_r
# Watersheds excluded in EITHER dimension get code 0 (excluded), overwriting
# any tercile code assigned from the other dimension's non-NA value
excluded_mask <- is.na(magA_r) | is.na(magB_r)
bivar_r <- ifel(excluded_mask, 0, bivar_r)
names(bivar_r) <- "bivar_class"

bivar_out <- file.path(out15, "bivariate_sedexport_magnitude_agreement_filtered.tif")
writeRaster(bivar_r, bivar_out, datatype = "INT1U", NAflag = 255, overwrite = TRUE)

bivar_lookup <- tibble(
  code = 0:9,
  class_A_magnitude = c("Excluded (low baseline)", rep(c("Low benefit (Baseline)", "Medium benefit (Baseline)", "High benefit (Baseline)"), each = 3)),
  class_B_agreement  = c("Excluded (low baseline)", rep(c("Low benefit (SSP370)", "Medium benefit (SSP370)", "High benefit (SSP370)"), times = 3)),
  suggested_hex = c(
    "#bfbfbf",                          # 0: excluded
    "#e8e8e8", "#b0d5df", "#64acbe",     # low A: low/med/high B
    "#e4acac", "#ad9ea5", "#627f8c",     # med A: low/med/high B
    "#c85a5a", "#985356", "#574249"      # high A: low/med/high B
  )
)
write_csv(bivar_lookup, file.path(out15, "bivariate_lookup_table_filtered.csv"))

p_bivar <- ggplot() +
  geom_raster(data = as.data.frame(bivar_r, xy = TRUE),
              aes(x = x, y = y, fill = factor(bivar_class))) +
  scale_fill_manual(values = setNames(bivar_lookup$suggested_hex, bivar_lookup$code),
                    name = "Bivariate\nclass", guide = "none") +
  coord_equal() +
  labs(
    title = "Bivariate raster: sed_export benefit magnitude \u00d7 cross-climate agreement",
    subtitle = "3x3 classification (filtered terciles) | Grey = excluded low-baseline watersheds",
    caption = "Darker/saturated = high benefit AND high agreement across climates. Import .tif into ArcGIS Pro."
  ) +
  theme_void(base_size = 12)

ggsave(file.path(out15, "map4_bivariate_quickview_filtered.png"),
       p_bivar, width = 9, height = 8, dpi = 200)

message("\nSection 15 (filtered) complete. Outputs in: ", out15)
message("  map1_af3_sedexport_change_baseline_top25_filtered.png")
message("  map2_af3_sedexport_smallmultiples_filtered.png")
message("  map2b_cropland_sedexport_smallmultiples_filtered.png  (NEW \u2014 same scale, direct comparison)")
message("  map3_exclusion_share_by_run.png  (NEW \u2014 shows % excluded per panel)")
message("  map4_bivariate_quickview_filtered.png + bivariate_sedexport_magnitude_agreement_filtered.tif")
message("  exclusion_summary_by_run.csv, map_data_", map_metric, "_all_runs.csv")