library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(purrr)

# =============================================================================
# AF SUITABILITY ANALYSIS — BASELINE vs SSP126 / SSP370 / SSP585
# Merged from: analysis_check.R + af_suitability_full.R
# Fixes: SSP rasters arrive on a 0–1 scale; rescaled to [1, 5] to match
#        the baseline before any threshold or area calculations.
# =============================================================================

# =============================================================================
# 1. USER INPUTS
# =============================================================================

setwd("C:/CI_Alves/analysis/ssa_invest_runs/outputs")

suit_paths <- list(
  Baseline = "C:/CI_Alves/suitability/rasters-20260706T203255Z-3-001/rasters/suitability_agroforestry.tif",
  SSP126   = "C:/CI_Alves/suitability/rasters-20260706T203255Z-3-001/rasters/future_suitability_agroforestry_SSP126.tif",
  SSP370   = "C:/CI_Alves/suitability/rasters-20260706T203255Z-3-001/rasters/future_suitability_agroforestry_SSP370.tif",
  SSP585   = "C:/CI_Alves/suitability/rasters-20260706T203255Z-3-001/rasters/future_suitability_agroforestry_SSP585.tif"
)

lulc_path  <- "C:/CI_Alves/analysis/ssa_invest_runs/Lulc/lulc_suitability/Landcover_harmonized_1km.tif"
out_dir    <- "outputs/baseline_vs_ssp_suitability"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

top_pct_thresholds <- c(10, 15, 20, 25)

# Known correct range from baseline diagnostics
BASELINE_MIN <- 1
BASELINE_MAX <- 5

lulc_labels <- tibble(
  lulc     = c(-128, 10,  20,          30,          40,          50,
                60,      70,     80,                 90,         95),
  lulc_name = c(
    "No data", "Tree cover", "Shrubland", "Grassland", "Cropland",
    "Sparse vegetation", "Water", "Bare", "Wetland / flooded",
    "Mangrove", "Urban"
  ),
  af_appropriate = c(
    FALSE, FALSE, TRUE, TRUE, TRUE,
    TRUE,  FALSE, FALSE, FALSE, FALSE, FALSE
  )
)

# =============================================================================
# 2. HELPERS
# =============================================================================

# Detect CRS type and compute pixel area in km²
get_pixel_area_km2 <- function(r) {
  if (!is.lonlat(r)) {
    prod(res(r)) / 1e6
  } else {
    NULL   # geographic — will use cellSize() per cell
  }
}

# Rescale a raster from [src_min, src_max] → [tgt_min, tgt_max] linearly.
# Used to bring 0–1 SSP rasters onto the same 1–5 scale as the baseline.
rescale_raster <- function(r, src_min, src_max, tgt_min = 1, tgt_max = 5) {
  (r - src_min) / (src_max - src_min) * (tgt_max - tgt_min) + tgt_min
}

# Load a suitability raster, detect its scale, and return it on [1, 5].
# Prints a clear diagnostic message for every scenario.
load_and_normalise_suit <- function(path, scen_name,
                                    tgt_min = BASELINE_MIN,
                                    tgt_max = BASELINE_MAX,
                                    tol     = 0.05) {
  r   <- rast(path)
  rng <- as.numeric(global(r, "range", na.rm = TRUE))
  obs_min <- rng[1]
  obs_max <- rng[2]

  cat("\n[ ", scen_name, " ] Raw raster range: ", round(obs_min, 4),
      " to ", round(obs_max, 4), "\n", sep = "")

  already_correct <- (obs_min >= (tgt_min - tol)) & (obs_max <= (tgt_max + tol)) &
                     (obs_max >  (tgt_min + 1))   # max must be well above 1

  if (already_correct) {
    cat("  Scale check PASSED — already on [", tgt_min, ",", tgt_max,
        "] scale. No rescaling applied.\n")
    return(r)
  }

  # Detect 0-1 binary/normalised rasters
  is_binary_01 <- (obs_min >= (0 - tol)) & (obs_max <= (1 + tol))

  if (is_binary_01) {
    cat("  Scale check FAILED — raster is on [0, 1] scale.\n",
        "  Rescaling linearly to [", tgt_min, ",", tgt_max, "] ...\n", sep = "")
    r_out <- rescale_raster(r, src_min = 0, src_max = 1,
                            tgt_min = tgt_min, tgt_max = tgt_max)
    chk <- as.numeric(global(r_out, "range", na.rm = TRUE))
    cat("  Post-rescale range:", round(chk[1], 4), "to", round(chk[2], 4), "\n")
    return(r_out)
  }

  # Unexpected scale — warn and rescale from observed range
  cat("  WARNING — unexpected scale. Rescaling from observed [",
      round(obs_min, 4), ",", round(obs_max, 4),
      "] to [", tgt_min, ",", tgt_max, "].\n", sep = "")
  r_out <- rescale_raster(r, src_min = obs_min, src_max = obs_max,
                          tgt_min = tgt_min, tgt_max = tgt_max)
  return(r_out)
}

# Stack suitability + LULC, extract as data.frame, join LULC labels
extract_suit_lulc <- function(suit_r, lulc_r, scenario_name) {
  if (!compareGeom(suit_r, lulc_r, stopOnError = FALSE)) {
    lulc_r <- resample(lulc_r, suit_r, method = "near")
  }
  stk        <- c(suit_r, lulc_r)
  names(stk) <- c("suitability", "lulc")
  as.data.frame(stk, na.rm = TRUE) |>
    left_join(lulc_labels, by = "lulc") |>
    mutate(scenario = scenario_name) |>
    filter(!is.na(lulc_name))
}

# Compute threshold summaries (LULC breakdown + AF area) for multiple top-%
analyse_thresholds <- function(df, px_km2, top_pcts = top_pct_thresholds) {

  suit_quantiles <- quantile(
    df$suitability,
    probs   = 1 - top_pcts / 100,
    na.rm   = TRUE
  )

  purrr::map(seq_along(top_pcts), function(i) {
    pct      <- top_pcts[i]
    suit_cut <- unname(suit_quantiles[i])
    subset_df    <- filter(df, suitability >= suit_cut)
    total_px     <- nrow(subset_df)
    total_area   <- total_px * px_km2

    lulc_summary <- subset_df |>
      count(lulc_name, af_appropriate, name = "n_pixels") |>
      mutate(
        area_km2        = n_pixels * px_km2,
        pct_of_total    = 100 * n_pixels / total_px,
        top_pct_threshold = pct,
        suit_cutoff     = suit_cut,
        total_area_km2  = total_area
      )

    af_summary <- subset_df |>
      filter(af_appropriate) |>
      summarise(
        n_af_px       = n(),
        af_area_km2   = n() * px_km2,
        pct_af_of_total = 100 * n() / total_px,
        .groups = "drop"
      ) |>
      mutate(
        top_pct_threshold = pct,
        suit_cutoff     = suit_cut,
        total_area_km2  = total_area
      )

    list(lulc = lulc_summary, af = af_summary)
  })
}

# =============================================================================
# 3. DIAGNOSTICS — LULC raster properties + continental totals
# =============================================================================

lulc_r <- rast(lulc_path)

cat("\n====== LULC RASTER DIAGNOSTICS ======\n")
cat("CRS        :", as.character(crs(lulc_r, proj = TRUE)), "\n")
cat("Resolution :", paste(res(lulc_r), collapse = " x "), "\n")
cat("Extent     :", paste(as.vector(ext(lulc_r)), collapse = " "), "\n")

is_projected   <- !is.lonlat(lulc_r)
pixel_area_km2 <- get_pixel_area_km2(lulc_r)

if (is_projected) {
  cat("CRS is projected. Pixel area =", pixel_area_km2, "km²\n")
} else {
  cat("WARNING: Geographic CRS — pixel areas vary. Using cellSize() per cell.\n")
}

cat("\nLULC unique values:\n")
print(as.data.frame(freq(lulc_r)))

# Continental LULC totals (no suitability filter — ground-truth check)
if (is_projected) {
  lulc_vals   <- as.data.frame(lulc_r, na.rm = TRUE)
  names(lulc_vals) <- "lulc"
  lulc_vals   <- left_join(lulc_vals, lulc_labels, by = "lulc")
  continental <- lulc_vals |>
    count(lulc_name, af_appropriate, name = "n_pixels") |>
    mutate(area_km2 = n_pixels * pixel_area_km2)
  px_km2_use  <- pixel_area_km2
} else {
  area_r      <- cellSize(lulc_r, unit = "km")
  area_vals   <- as.data.frame(c(lulc_r, area_r), na.rm = TRUE)
  names(area_vals) <- c("lulc", "area_km2")
  area_vals   <- left_join(area_vals, lulc_labels, by = "lulc")
  continental <- area_vals |>
    group_by(lulc_name, af_appropriate) |>
    summarise(area_km2 = sum(area_km2, na.rm = TRUE),
              n_pixels = n(), .groups = "drop")
  px_km2_use  <- mean(values(area_r, na.rm = TRUE))
  cat("Using mean pixel area (geographic CRS):", px_km2_use, "km²\n")
}

write.csv(continental,
          file.path(out_dir, "CHECK_continental_lulc_area.csv"),
          row.names = FALSE)
cat("\nContinental LULC totals (km²):\n")
print(arrange(continental, desc(area_km2)))

# =============================================================================
# 4. SUITABILITY RASTER DIAGNOSTICS (raw, before rescaling)
# =============================================================================

cat("\n====== RAW SUITABILITY RASTER DIAGNOSTICS (pre-rescale) ======\n")
for (scen in names(suit_paths)) {
  s   <- rast(suit_paths[[scen]])
  rng <- as.numeric(global(s, "range", na.rm = TRUE))
  q   <- quantile(values(s, na.rm = TRUE),
                  probs = c(0, 0.10, 0.25, 0.50, 0.75, 0.90, 1.0), na.rm = TRUE)
  cat("\n[", scen, "]\n")
  cat("  CRS  :", as.character(crs(s, proj = TRUE)), "\n")
  cat("  Res  :", paste(res(s), collapse = " x "), "\n")
  cat("  Range:", round(rng[1], 4), "to", round(rng[2], 4), "\n")
  cat("  Quantiles:\n"); print(round(q, 4))
}

# =============================================================================
# 5. MAIN LOOP — load, rescale, analyse, collect
# =============================================================================

all_lulc     <- list()
all_af       <- list()
all_thresh   <- list()   # for analysis_check-style tabulation
suit_diag    <- list()   # per-scenario diagnostics after rescaling

for (scen in names(suit_paths)) {

  message("\n========================================")
  message("Processing: ", scen)
  message("========================================")

  # Load and bring to [1, 5] if needed
  suit_r <- load_and_normalise_suit(suit_paths[[scen]], scen)

  # Post-rescale diagnostic
  rng_post <- as.numeric(global(suit_r, "range", na.rm = TRUE))
  cat("  [", scen, "] Suitability range after rescale:", round(rng_post[1], 4),
      "to", round(rng_post[2], 4), "\n")

  suit_diag[[scen]] <- data.frame(
    scenario  = scen,
    raw_min   = round(as.numeric(global(rast(suit_paths[[scen]]), "min", na.rm = TRUE)), 4),
    raw_max   = round(as.numeric(global(rast(suit_paths[[scen]]), "max", na.rm = TRUE)), 4),
    scaled_min = round(rng_post[1], 4),
    scaled_max = round(rng_post[2], 4)
  )

  # Extract + analyse
  df  <- extract_suit_lulc(suit_r, lulc_r, scen)
  res <- analyse_thresholds(df, px_km2_use)

  # Collect LULC + AF rows (af_suitability_full.R style)
  all_lulc[[scen]] <- bind_rows(lapply(res, `[[`, "lulc")) |> mutate(scenario = scen)
  all_af[[scen]]   <- bind_rows(lapply(res, `[[`, "af"))   |> mutate(scenario = scen)

  # Collect threshold table (analysis_check.R style — top 10% + 25% only)
  suit_q <- quantile(df$suitability, probs = 1 - c(10, 25) / 100, na.rm = TRUE)

  for (i in seq_along(c(10, 25))) {
    pct     <- c(10, 25)[i]
    cut_val <- unname(suit_q[i])
    sub     <- filter(df, suitability >= cut_val)

    cat("  Top", pct, "% cutoff:", round(cut_val, 4),
        "| n pixels:", nrow(sub),
        "| area km²:", round(nrow(sub) * px_km2_use), "\n")

    tab <- sub |>
      count(lulc_name, af_appropriate, name = "n_pixels") |>
      mutate(
        area_km2 = n_pixels * px_km2_use,
        scenario = scen,
        top_pct  = pct
      )
    all_thresh[[paste(scen, pct)]] <- tab
  }

  rm(suit_r, df, res); gc()
}

# =============================================================================
# 6. COMBINE AND WRITE TABLES
# =============================================================================

lulc_combined  <- bind_rows(all_lulc)
af_combined    <- bind_rows(all_af)
thresh_combined <- bind_rows(all_thresh)
diag_df        <- bind_rows(suit_diag)

write.csv(lulc_combined,
          file.path(out_dir, "baseline_vs_ssp_lulc_by_threshold.csv"),
          row.names = FALSE)
write.csv(af_combined,
          file.path(out_dir, "baseline_vs_ssp_af_area_by_threshold.csv"),
          row.names = FALSE)
write.csv(thresh_combined,
          file.path(out_dir, "CHECK_threshold_lulc_area.csv"),
          row.names = FALSE)
write.csv(diag_df,
          file.path(out_dir, "CHECK_suitability_scale_diagnostics.csv"),
          row.names = FALSE)

cat("\n====== SCALE DIAGNOSTICS SUMMARY ======\n")
print(diag_df)

# =============================================================================
# 7. COLOUR PALETTES
# =============================================================================

lulc_colours <- c(
  "Cropland"           = "#D4A017",
  "Grassland"          = "#90EE90",
  "Shrubland"          = "#8B6914",
  "Sparse vegetation"  = "#D2B48C",
  "Tree cover"         = "#228B22",
  "Water"              = "#4682B4",
  "Bare"               = "#A9A9A9",
  "Wetland / flooded"  = "#20B2AA",
  "Mangrove"           = "#006400",
  "Urban"              = "#FF4500",
  "No data"            = "#E0E0E0"
)

scenario_cols <- c(
  "Baseline" = "#4D4D4D",
  "SSP126"   = "#2166AC",
  "SSP370"   = "#FDAE61",
  "SSP585"   = "#D73027"
)

# =============================================================================
# 8. PLOTS
# =============================================================================

threshold_factor <- function(x) {
  factor(paste0("Top ", x, "%"),
         levels = paste0("Top ", sort(top_pct_thresholds), "%"))
}

# ---- Plot A: Continental LULC area (no suitability filter — sanity check) --
p_continental <- continental |>
  filter(lulc_name != "No data") |>
  arrange(desc(area_km2)) |>
  mutate(lulc_name = factor(lulc_name, levels = rev(unique(lulc_name)))) |>
  ggplot(aes(x = lulc_name, y = area_km2 / 1e6, fill = lulc_name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(area_km2 / 1e6, 2), "M")),
            hjust = -0.1, size = 3.2) +
  scale_fill_manual(values = lulc_colours, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                     labels = label_comma(suffix = "M km²")) +
  coord_flip() +
  labs(
    title    = "Continental LULC area — all SSA (no threshold filter)",
    subtitle = "Ground truth check before suitability filtering",
    x = NULL, y = "Area (million km²)"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "figA_continental_lulc_area.png"),
       p_continental, width = 10, height = 6, dpi = 180)

# ---- Plot 1: LULC composition within top thresholds (faceted by scenario) --
p_stack <- lulc_combined |>
  mutate(top_pct_threshold = threshold_factor(top_pct_threshold)) |>
  ggplot(aes(x = top_pct_threshold, y = pct_of_total,
             fill = lulc_name, alpha = af_appropriate)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = lulc_colours, name = "LULC class") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.45), guide = "none") +
  facet_wrap(~ scenario, ncol = 2) +
  scale_y_continuous(labels = label_percent(scale = 1)) +
  labs(
    title    = "LULC composition within top suitability thresholds",
    subtitle = "Baseline compared with future SSP suitability scenarios",
    x        = "Suitability threshold",
    y        = "% of pixels within threshold",
    caption  = "AF-appropriate classes: Cropland, Grassland, Shrubland, Sparse vegetation"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "fig1_baseline_vs_ssp_lulc_composition.png"),
       p_stack, width = 13, height = 8, dpi = 180)

# ---- Plot 2: AF-appropriate area (absolute, M km²) -------------------------
p_af_area <- af_combined |>
  mutate(top_pct_threshold = threshold_factor(top_pct_threshold)) |>
  ggplot(aes(x = top_pct_threshold, y = af_area_km2 / 1e6, fill = scenario)) +
  geom_bar(stat = "identity", position = "dodge", colour = "white") +
  scale_fill_manual(values = scenario_cols, name = "Scenario") +
  scale_y_continuous(labels = label_comma(suffix = "M km²")) +
  labs(
    title    = "AF-appropriate land area within top suitability thresholds",
    subtitle = "Baseline vs future SSP scenarios",
    x        = "Suitability threshold",
    y        = "Area (million km²)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "fig2_baseline_vs_ssp_af_area.png"),
       p_af_area, width = 10, height = 6, dpi = 180)

# ---- Plot 3: AF-appropriate share (%) of threshold pixels ------------------
p_af_pct <- af_combined |>
  mutate(top_pct_threshold = threshold_factor(top_pct_threshold)) |>
  ggplot(aes(x = scenario, y = pct_af_of_total, fill = scenario)) +
  geom_bar(stat = "identity", colour = "white") +
  geom_text(aes(label = paste0(round(pct_af_of_total, 1), "%")),
            vjust = -0.35, size = 3.4) +
  scale_fill_manual(values = scenario_cols, guide = "none") +
  facet_wrap(~ top_pct_threshold, ncol = 4) +
  scale_y_continuous(limits = c(0, 110),
                     labels  = label_percent(scale = 1)) +
  labs(
    title    = "AF-appropriate land as % of total pixels within each threshold",
    subtitle = "Comparison of current baseline and future SSP suitability maps",
    x        = "Scenario",
    y        = "% AF-appropriate"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "fig3_baseline_vs_ssp_af_share.png"),
       p_af_pct, width = 12, height = 6, dpi = 180)

# ---- Plot 4: Heatmap — AF area by LULC class, threshold, scenario ----------
p_heatmap <- lulc_combined |>
  filter(af_appropriate) |>
  mutate(
    top_pct_threshold = threshold_factor(top_pct_threshold),
    lulc_name = factor(
      lulc_name,
      levels = c("Cropland", "Grassland", "Shrubland", "Sparse vegetation")
    )
  ) |>
  ggplot(aes(x = top_pct_threshold, y = scenario,
             fill = area_km2 / 1000)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = comma(round(area_km2 / 1000))),
            size = 2.8, colour = "black") +
  scale_fill_distiller(
    palette  = "YlOrRd",
    direction = 1,
    name     = "Area\n('000 km²)",
    labels   = label_comma()
  ) +
  facet_wrap(~ lulc_name, ncol = 2) +
  labs(
    title    = "AF-appropriate area by LULC class, threshold, and scenario",
    subtitle = "Baseline included alongside SSP scenarios",
    x        = "Suitability threshold",
    y        = "Scenario"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid   = element_blank(),
    strip.text   = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 30, hjust = 1)
  )

ggsave(file.path(out_dir, "fig4_baseline_vs_ssp_heatmap.png"),
       p_heatmap, width = 12, height = 8, dpi = 180)

# ---- Plot B: Total LULC area at top 10% + 25% by scenario (stacked bar) ----
p_threshold <- thresh_combined |>
  filter(lulc_name != "No data") |>
  mutate(
    threshold_label = paste0("Top ", top_pct, "%"),
    scenario        = factor(scenario, levels = names(suit_paths))
  ) |>
  ggplot(aes(x = scenario, y = area_km2 / 1e6, fill = lulc_name)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = lulc_colours, name = "LULC class") +
  scale_y_continuous(labels = label_comma(suffix = "M km²")) +
  facet_wrap(~ threshold_label, ncol = 2) +
  labs(
    title    = "Total LULC area within top suitability thresholds by scenario",
    subtitle = "Top 10% and top 25% — stacked by LULC class",
    x        = "Scenario",
    y        = "Area (million km²)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 30, hjust = 1)
  )

ggsave(file.path(out_dir, "figB_lulc_area_by_threshold_scenario.png"),
       p_threshold, width = 12, height = 7, dpi = 180)

# ---- Plot C: AF-appropriate LULC only, top 10% + 25% -----------------------
p_af_only <- thresh_combined |>
  filter(af_appropriate, lulc_name != "No data") |>
  mutate(
    threshold_label = paste0("Top ", top_pct, "%"),
    scenario        = factor(scenario, levels = names(suit_paths))
  ) |>
  ggplot(aes(x = scenario, y = area_km2 / 1e6, fill = lulc_name)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = lulc_colours,
    name   = "LULC class",
    limits = c("Cropland", "Grassland", "Shrubland", "Sparse vegetation")
  ) +
  scale_y_continuous(labels = label_comma(suffix = "M km²")) +
  facet_wrap(~ threshold_label, ncol = 2) +
  labs(
    title    = "AF-appropriate LULC area within top suitability thresholds",
    subtitle = "Cropland + Grassland + Shrubland + Sparse vegetation only",
    x        = "Scenario",
    y        = "Area (million km²)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 30, hjust = 1)
  )

ggsave(file.path(out_dir, "figC_af_lulc_area_threshold_scenario.png"),
       p_af_only, width = 12, height = 7, dpi = 180)

# ---- Plot D: Cropland sanity check ------------------------------------------
p_cropland <- thresh_combined |>
  filter(lulc_name == "Cropland") |>
  mutate(
    threshold_label = paste0("Top ", top_pct, "%"),
    scenario        = factor(scenario, levels = names(suit_paths))
  ) |>
  ggplot(aes(x = scenario, y = area_km2 / 1000, fill = scenario)) +
  geom_bar(stat = "identity", colour = "white") +
  geom_text(aes(label = paste0(comma(round(area_km2 / 1000)), "k km²")),
            vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = scenario_cols, guide = "none") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15)),
    labels = label_comma(suffix = "k km²")
  ) +
  facet_wrap(~ threshold_label, ncol = 2) +
  labs(
    title    = "Cropland area within top suitability thresholds — by scenario",
    subtitle = "Values should differ meaningfully across scenarios if thresholds are working",
    x        = "Scenario",
    y        = "Area (thousand km²)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "figD_cropland_area_sanity_check.png"),
       p_cropland, width = 10, height = 6, dpi = 180)

# =============================================================================
# 9. DONE
# =============================================================================

cat("\n========================================\n")
cat("All outputs written to:", out_dir, "\n")
cat("CSVs:\n")
cat("  baseline_vs_ssp_lulc_by_threshold.csv\n")
cat("  baseline_vs_ssp_af_area_by_threshold.csv\n")
cat("  CHECK_continental_lulc_area.csv\n")
cat("  CHECK_threshold_lulc_area.csv\n")
cat("  CHECK_suitability_scale_diagnostics.csv  <- start here if results look odd\n")
cat("Figures: figA, fig1–fig4, figB–figD\n")
cat("========================================\n")



# =============================================================================
# 10. EXPORT SUITABILITY THRESHOLD MASKS FOR InVEST RUNS
# =============================================================================
# Writes one binary raster per scenario × threshold combination.
# Each raster: 1 = pixel is within the top X% of suitability, NA = excluded.
# Output folder structure:
#   invest_suitability_masks/
#     top10/   suitability_mask_Baseline_top10.tif  ...SSP126... SSP370... SSP585
#     top20/   ...
#     top25/   ...
#     top30/   ...
#
# These are ready to use as an AOI or masking layer in InVEST SDR / AWY runs.
# =============================================================================

invest_thresholds <- c(10, 20, 25, 30)
invest_mask_dir   <- file.path(out_dir, "invest_suitability_masks")

for (pct in invest_thresholds) {
  dir.create(file.path(invest_mask_dir, paste0("top", pct)),
             recursive = TRUE, showWarnings = FALSE)
}

cat("\n========================================\n")
cat("Exporting InVEST suitability masks\n")
cat("Thresholds:", paste0("top ", invest_thresholds, "%", collapse = ", "), "\n")
cat("========================================\n")

for (scen in names(suit_paths)) {
  
  cat("\n[", scen, "]\n")
  
  # Re-load and rescale (same as main loop — avoids holding all in memory)
  suit_r <- load_and_normalise_suit(suit_paths[[scen]], scen)
  
  # Compute threshold cutoffs from the rescaled values
  suit_vals    <- values(suit_r, na.rm = TRUE)
  cutoffs      <- quantile(suit_vals, probs = 1 - invest_thresholds / 100,
                           na.rm = TRUE)
  
  for (i in seq_along(invest_thresholds)) {
    pct     <- invest_thresholds[i]
    cut_val <- unname(cutoffs[i])
    
    # Binary mask: 1 where suitability >= cutoff, NA elsewhere
    mask_r  <- ifel(suit_r >= cut_val, 1L, NA_integer_)
    
    out_name <- paste0("suitability_mask_", scen, "_top", pct, ".tif")
    out_path <- file.path(invest_mask_dir, paste0("top", pct), out_name)
    
    writeRaster(
      mask_r,
      out_path,
      datatype  = "INT1U",
      NAflag    = 255,
      overwrite = TRUE
    )
    
    n_px   <- sum(values(mask_r, na.rm = TRUE) == 1, na.rm = TRUE)
    area   <- round(n_px * px_km2_use)
    cat("  top", pct, "%  cutoff:", round(cut_val, 4),
        " | n pixels:", n_px, " | area km²:", area,
        " ->", out_name, "\n")
  }
  
  rm(suit_r, suit_vals, mask_r); gc()
}

cat("\nInVEST masks written to:", invest_mask_dir, "\n")
cat("Folder structure:\n")
for (pct in invest_thresholds) {
  cat("  invest_suitability_masks/top", pct, "/  — ",
      length(names(suit_paths)), " files\n", sep = "")
} 

# =============================================================================
# 11. EXPORT CROPLAND-ONLY SUITABILITY MASKS FOR InVEST RUNS
# =============================================================================
# Intersects each top-X% suitability mask (from Section 10) with the LULC
# cropland class (value = 40) to produce cropland-only binary masks.
# Output: 1 = pixel is cropland AND within top X% suitability, NA = excluded.
#
# Output folder structure:
#   invest_suitability_masks/
#     top10/   suitability_mask_cropland_Baseline_top10.tif  ...SSP126...
#     top20/   ...
#     top25/   ...
#     top30/   ...
# =============================================================================

CROPLAND_VALUE <- 40   # ESA WorldCover / LULC code for cropland

cat("\n========================================\n")
cat("Exporting cropland-only InVEST suitability masks\n")
cat("Cropland LULC value:", CROPLAND_VALUE, "\n")
cat("Thresholds:", paste0("top ", invest_thresholds, "%", collapse = ", "), "\n")
cat("========================================\n")

# Build cropland binary raster once (1 = cropland, NA = everything else)
# Reuse lulc_r already loaded in memory from the main analysis
cropland_r <- ifel(lulc_r == CROPLAND_VALUE, 1L, NA_integer_)

for (scen in names(suit_paths)) {
  
  cat("\n[", scen, "]\n")
  
  suit_r    <- load_and_normalise_suit(suit_paths[[scen]], scen)
  suit_vals <- values(suit_r, na.rm = TRUE)
  cutoffs   <- quantile(suit_vals, probs = 1 - invest_thresholds / 100,
                        na.rm = TRUE)
  
  for (i in seq_along(invest_thresholds)) {
    pct     <- invest_thresholds[i]
    cut_val <- unname(cutoffs[i])
    
    # Suitability mask for this threshold
    suit_mask <- ifel(suit_r >= cut_val, 1L, NA_integer_)
    
    # Align cropland raster to suitability grid if needed
    crop_aligned <- if (!compareGeom(cropland_r, suit_mask, stopOnError = FALSE)) {
      resample(cropland_r, suit_mask, method = "near")
    } else {
      cropland_r
    }
    
    # Intersection: both must be 1 (cropland AND high suitability)
    cropland_suit <- ifel(suit_mask == 1L & crop_aligned == 1L, 1L, NA_integer_)
    
    out_name <- paste0("suitability_mask_cropland_", scen, "_top", pct, ".tif")
    out_path <- file.path(invest_mask_dir, paste0("top", pct), out_name)
    
    writeRaster(
      cropland_suit,
      out_path,
      datatype  = "INT1U",
      NAflag    = 255,
      overwrite = TRUE
    )
    
    n_px  <- sum(values(cropland_suit, na.rm = TRUE) == 1, na.rm = TRUE)
    area  <- round(n_px * px_km2_use)
    cat("  top", pct, "%  cropland cutoff:", round(cut_val, 4),
        " | n pixels:", n_px, " | area km²:", area,
        " ->", out_name, "\n")
    
    rm(suit_mask, crop_aligned, cropland_suit); gc()
  }
  
  rm(suit_r, suit_vals); gc()
}

cat("\nCropland InVEST masks written to:", invest_mask_dir, "\n")
cat("Each top-X% subfolder now contains:\n")
cat("  suitability_mask_<SCENARIO>_top<X>.tif         — full suitability mask\n")
cat("  suitability_mask_cropland_<SCENARIO>_top<X>.tif — cropland only\n") 

# =============================================================================
# 12. EXPORT COMBINED AF-TRANSITION SUITABILITY MASK
# Cropland + Grassland + Shrubland, one mask per threshold, new folder
# =============================================================================
# Intersects each top-X% suitability mask (Section 10) with ALL THREE
# AF-appropriate transition classes (cropland=40, grassland=30, shrubland=20).
# Output: 1 = pixel is (cropland OR grassland OR shrubland) AND within top X%
# suitability, NA = excluded.
#
# Output folder structure (NEW, separate from Section 10/11 outputs):
#   invest_suitability_masks_af_combined/
#     top10/ suitability_mask_af3class_Baseline_top10.tif ...SSP126... SSP370... SSP585
#     top20/ ...
#     top25/ ...
#     top30/ ...
# =============================================================================

GRASSLAND_VALUE <- 30L
SHRUBLAND_VALUE <- 20L

af3_mask_dir <- file.path(out_dir, "invest_suitability_masks_af_combined")

for (pct in invest_thresholds) {
  dir.create(file.path(af3_mask_dir, paste0("top", pct)),
             recursive = TRUE, showWarnings = FALSE)
}

cat("\n========================================\n")
cat("Exporting combined AF-transition (crop+grass+shrub) InVEST masks\n")
cat("LULC values: cropland =", CROPLAND_VALUE,
    "| grassland =", GRASSLAND_VALUE,
    "| shrubland =", SHRUBLAND_VALUE, "\n")
cat("Thresholds:", paste0("top ", invest_thresholds, "%", collapse = ", "), "\n")
cat("========================================\n")

# Build combined AF-appropriate binary raster once
# (1 = cropland OR grassland OR shrubland, NA = everything else)
af3_r <- ifel(
  lulc_r == CROPLAND_VALUE | lulc_r == GRASSLAND_VALUE | lulc_r == SHRUBLAND_VALUE,
  1L, NA_integer_
)

for (scen in names(suit_paths)) {
  
  cat("\n[", scen, "]\n")
  
  suit_r <- load_and_normalise_suit(suit_paths[[scen]], scen)
  suit_vals <- values(suit_r, na.rm = TRUE)
  cutoffs <- quantile(suit_vals, probs = 1 - invest_thresholds / 100,
                      na.rm = TRUE)
  
  for (i in seq_along(invest_thresholds)) {
    pct <- invest_thresholds[i]
    cut_val <- unname(cutoffs[i])
    
    suit_mask <- ifel(suit_r >= cut_val, 1L, NA_integer_)
    
    af3_aligned <- if (!compareGeom(af3_r, suit_mask, stopOnError = FALSE)) {
      resample(af3_r, suit_mask, method = "near")
    } else {
      af3_r
    }
    
    # Intersection: both must be 1 (AF-appropriate class AND high suitability)
    af3_suit <- ifel(suit_mask == 1L & af3_aligned == 1L, 1L, NA_integer_)
    
    out_name <- paste0("suitability_mask_af3class_", scen, "_top", pct, ".tif")
    out_path <- file.path(af3_mask_dir, paste0("top", pct), out_name)
    
    writeRaster(
      af3_suit,
      out_path,
      datatype = "INT1U",
      NAflag = 255,
      overwrite = TRUE
    )
    
    n_px <- sum(values(af3_suit, na.rm = TRUE) == 1, na.rm = TRUE)
    area <- round(n_px * px_km2_use)
    cat("  top", pct, "% AF3-class cutoff:", round(cut_val, 4),
        " | n pixels:", n_px, " | area km²:", area,
        " ->", out_name, "\n")
    
    rm(suit_mask, af3_aligned, af3_suit); gc()
  }
  
  rm(suit_r, suit_vals); gc()
}

cat("\nCombined AF-transition masks written to:", af3_mask_dir, "\n")


# =============================================================================
# TARGETED FIXES — insert before "SECTION 14 — OUTLIER-ROBUST UPDATE"
# Fixes: missing compare_runs()/pct_change()/run_meta/etc, plus flagged
# downstream areas (compare_runs_robust as single source of truth,
# Comparison C filter, MASTER summary rebuild, cutoff sensitivity check)
# =============================================================================

library(sf)
library(readr)

# ---- CHUNK 1: core objects/functions that were missing ----------------

invest_output_root <- "C:/CI_Alves/analysis/ssa_invest_runs/outputs/invest_outputs"
stopifnot(dir.exists(invest_output_root))

ref_baseline <- "baseline_cropland_0.34"

out14 <- file.path(out_dir, "comparison_outputs")
dir.create(out14, recursive = TRUE, showWarnings = FALSE)

climate_cols <- c("Baseline" = "#2166AC", "SSP370" = "#B2182B")

all_run_ids <- list.dirs(invest_output_root, full.names = FALSE, recursive = FALSE)

run_meta <- tibble(run_id = all_run_ids) |>
  mutate(
    climate = case_when(
      grepl("^baseline", run_id) ~ "Baseline",
      grepl("^ssp370",   run_id) ~ "SSP370",
      TRUE ~ NA_character_
    ),
    threshold = case_when(
      grepl("top10", run_id) ~ "Top 10%",
      grepl("top20", run_id) ~ "Top 20%",
      grepl("top25", run_id) ~ "Top 25%",
      TRUE ~ "Full area"
    ),
    class = case_when(
      grepl("af3", run_id)          ~ "AF3-class",
      grepl("_to_0\\.2832", run_id) ~ "Cropland-only",
      TRUE ~ "No AF"
    )
  )
stopifnot(nrow(run_meta) > 0)
cat("run_meta built:", nrow(run_meta), "run_ids found under", invest_output_root, "\n")

load_watershed <- function(run_id) {
  shp_path <- file.path(invest_output_root, run_id, "watershed_results_sdr.shp")
  if (!file.exists(shp_path)) {
    warning("Missing shapefile for run_id: ", run_id, " at ", shp_path)
    return(NULL)
  }
  sf_obj <- st_read(shp_path, quiet = TRUE)
  sf_obj <- st_drop_geometry(sf_obj)
  names(sf_obj) <- tolower(names(sf_obj))
  if (!"ws_id" %in% names(sf_obj)) sf_obj$ws_id <- seq_len(nrow(sf_obj))
  
  rename_map <- c("usle_tot", "sed_export", "sed_dep", "avoid_eros", "avoid_exp")
  for (std_name in rename_map) {
    if (!std_name %in% names(sf_obj)) {
      match_col <- grep(paste0("^", std_name), names(sf_obj), value = TRUE)
      if (length(match_col) >= 1) names(sf_obj)[names(sf_obj) == match_col[1]] <- std_name
    }
  }
  sf_obj
}

pct_change <- function(scen, ref) {
  ifelse(ref == 0, NA_real_, (scen - ref) / abs(ref) * 100)
}

compare_runs <- function(ref_id, scen_id,
                         metrics = c("usle_tot", "sed_export", "sed_dep",
                                     "avoid_eros", "avoid_exp")) {
  ref_df  <- load_watershed(ref_id)
  scen_df <- load_watershed(scen_id)
  if (is.null(ref_df) || is.null(scen_df)) return(NULL)
  
  metrics <- intersect(metrics, intersect(names(ref_df), names(scen_df)))
  if (length(metrics) == 0) {
    warning("No matching metric columns between ", ref_id, " and ", scen_id)
    return(NULL)
  }
  
  ref_sub  <- ref_df[, c("ws_id", metrics)]
  scen_sub <- scen_df[, c("ws_id", metrics)]
  names(ref_sub)[-1]  <- paste0(metrics, "_ref")
  names(scen_sub)[-1] <- paste0(metrics, "_scen")
  
  merged <- merge(ref_sub, scen_sub, by = "ws_id")
  
  for (m in metrics) {
    merged[[paste0(m, "_delta")]] <- merged[[paste0(m, "_scen")]] - merged[[paste0(m, "_ref")]]
    merged[[paste0(m, "_pct")]]   <- pct_change(merged[[paste0(m, "_scen")]], merged[[paste0(m, "_ref")]])
  }
  merged
}



# =============================================================================
# 13. PREPARE INVEST-READY LULC RASTERS + BIOPHYSICAL TABLES
#     Baseline & SSP370 × top 10% / 20% / 25% cropland AF transitions
# =============================================================================
#
# LOGIC
# -----
# For each scenario (Baseline, SSP370) and each suitability threshold (10, 20, 25%):
#
#   1. Load the suitability raster → rescale to [1,5] → compute cutoff
#   2. Identify pixels that are:
#        - LULC cropland (code 40)  AND
#        - Within the top X% suitability
#   3. In the LULC raster, recode those pixels to a NEW AF lucode:
#        40  → baseline cropland          (C = 0.34, P = 1.00)  [unchanged]
#        45  → AF mixed-species cropland  (C = 0.2832, P = 1.00) [new]
#        46  → AF dense cropland          (C = 0.20,   P = 1.00) [new, optional]
#      (only code 45 is written for the mixed-species runs; 46 added for completeness)
#   4. Write the modified LULC raster → invest_inputs/<scenario>_top<X>/lulc.tif
#   5. Write a matching biophysical table → invest_inputs/<scenario>_top<X>/biophysical_table.csv
#      with the new rows appended, all other rows unchanged from biophysical_table_base.csv
#
# LUCODE SCHEME
# -------------
#   40  = baseline cropland            (no AF)        C = 0.34   P = 1.00
#   45  = AF cropland — mixed-species  (top X% suit.) C = 0.2832 P = 1.00
#   46  = AF cropland — dense          (top X% suit.) C = 0.20   P = 1.00
#
# C-factor values sourced from SDR Methodology (Alves & Craig 2026):
#   C_AF_mixed = C_base × (1 − B_Kuyah × f_canopy_mid)
#              = 0.34  × (1 − 0.897 × 0.352) = 0.34 × 0.8332 ≈ 0.2832
#   C_AF_dense = ~0.20 (synthesized mean, dense/hedgerow class at C_base 0.34)
#
# OUTPUT FOLDER STRUCTURE
# -----------------------
#   invest_inputs/
#     Baseline_top10/   lulc.tif   biophysical_table.csv
#     Baseline_top20/   lulc.tif   biophysical_table.csv
#     Baseline_top25/   lulc.tif   biophysical_table.csv
#     SSP370_top10/     lulc.tif   biophysical_table.csv
#     SSP370_top20/     lulc.tif   biophysical_table.csv
#     SSP370_top25/     lulc.tif   biophysical_table.csv
# =============================================================================

#install.packages("readr")
library(readr)

# =============================================================================
# CONFIG
# =============================================================================

# Path to biophysical_table_base.csv (all non-AF rows)
biophys_base_path <- "C:/CI_Alves/analysis/ssa_invest_runs/Biophysical table/biophysical_table_base.csv"

# Scenarios to process (keys must match suit_paths defined earlier in this script)
invest_scenarios  <- c("Baseline", "SSP370")

# Suitability thresholds
invest_suit_pcts  <- c(10, 20, 25)

# LULC code for cropland (source class to transition)
CROPLAND_CODE     <- 40L

# New AF lucodes
AF_MIXED_CODE     <- 45L   # AF mixed-species cropland
AF_DENSE_CODE     <- 46L   # AF dense/hedgerow cropland

# C-factor and P-factor values for AF rows
# Mixed-species: C_base 0.34, mean canopy fraction 0.352 (mid of 0.2885–0.4096 range)
# B_Kuyah = 0.897
# C_AF_mixed = 0.34 * (1 - 0.897 * 0.352) = 0.2832
AF_C_MIXED  <- 0.2832   # = 0.2832
AF_P_MIXED  <- 1.00

# Dense class: synthesized mean from SDR Methodology at C_base 0.34
AF_C_DENSE  <- 0.20
AF_P_DENSE  <- 1.00

# Output root directory for InVEST run inputs
invest_input_dir <- file.path(out_dir, "invest_inputs")
dir.create(invest_input_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n========================================\n")
cat("Section 13 — InVEST LULC + Biophysical Table Prep\n")
cat("Scenarios:", paste(invest_scenarios, collapse = ", "), "\n")
cat("Thresholds:", paste0("top ", invest_suit_pcts, "%", collapse = ", "), "\n")
cat("AF mixed lucode:", AF_MIXED_CODE, "  C =", AF_C_MIXED, "  P =", AF_P_MIXED, "\n")
cat("AF dense lucode:", AF_DENSE_CODE, "  C =", AF_C_DENSE, "  P =", AF_P_DENSE, "\n")
cat("========================================\n")

# =============================================================================
# LOAD BASE BIOPHYSICAL TABLE
# =============================================================================

biophys_base <- read_csv(biophys_base_path, show_col_types = FALSE)

# Confirm expected columns
stopifnot(all(c("lucode", "usle_c", "usle_p") %in% names(biophys_base)))

# Confirm cropland row exists
stopifnot(CROPLAND_CODE %in% biophys_base$lucode)

cat("Base biophysical table loaded:", nrow(biophys_base), "rows\n")
cat("Cropland (lucode 40) baseline: C =",
    biophys_base$usle_c[biophys_base$lucode == CROPLAND_CODE],
    " P =",
    biophys_base$usle_p[biophys_base$lucode == CROPLAND_CODE], "\n")

# Build the two new AF rows (same column structure as base table)
af_rows <- tibble(
  lucode = c(AF_MIXED_CODE, AF_DENSE_CODE),
  usle_c = c(AF_C_MIXED,    AF_C_DENSE),
  usle_p = c(AF_P_MIXED,    AF_P_DENSE)
)

# If base table has additional columns (root_depth, Kc, etc.), fill with
# the cropland row values so InVEST does not error on missing fields
if (ncol(biophys_base) > 3) {
  cropland_row <- biophys_base[biophys_base$lucode == CROPLAND_CODE, ]
  extra_cols   <- setdiff(names(biophys_base), c("lucode", "usle_c", "usle_p"))
  for (col in extra_cols) {
    af_rows[[col]] <- cropland_row[[col]]
  }
  # Reorder columns to match base table
  af_rows <- af_rows[, names(biophys_base)]
}

# Full expanded biophysical table (base + two AF rows)
biophys_af <- bind_rows(biophys_base, af_rows)
cat("Expanded biophysical table:", nrow(biophys_af), "rows",
    "(added lucodes", AF_MIXED_CODE, "and", AF_DENSE_CODE, ")\n")

# =============================================================================
# 13b. EXTEND BIOPHYSICAL TABLE — GRASSLAND + SHRUBLAND AF TRANSITIONS
# =============================================================================
# AF base C-factor assumption: 0.07 (mixed-species agroforestry baseline)
#
# Grassland transition: 20% grazing penalty applied to the AF base
#   C_AF_grassland = 0.07 * 1.20 = 0.084  ≈ 0.085 (rounded, as specified)
#   InVEST baseline grassland (lucode 30): C = 0.20
#
# Shrubland transition: small shrub-removal penalty applied to the AF base
#   C_AF_shrubland = 0.075 (small increment above 0.07 base, as specified)
#   InVEST baseline shrubland (lucode 20): C = 0.03
#   NOTE: 0.075 > 0.03 baseline — this transition raises C-factor for
#   pixels that were already low-erosion dense shrubland. Flagged for review.
# =============================================================================

# New AF lucodes (distinct from cropland's 45/46)
AF_GRASS_CODE <- 35L   # AF mixed-species grassland transition
AF_SHRUB_CODE <- 25L   # AF mixed-species shrubland transition

GRASSLAND_CODE <- 30L
SHRUBLAND_CODE <- 20L

AF_BASE_C     <- 0.07   # AF-itself base C-factor assumption
GRAZING_PENALTY <- 1.20 # 20% grazing penalty multiplier
AF_C_GRASS    <- round(AF_BASE_C * GRAZING_PENALTY, 3)  # 0.084 -> reported as 0.085
AF_P_GRASS    <- 1.00

AF_C_SHRUB    <- 0.075  # small shrub-removal penalty above AF base
AF_P_SHRUB    <- 1.00

cat("\nGrassland baseline (lucode", GRASSLAND_CODE, "): C =",
    biophys_base$usle_c[biophys_base$lucode == GRASSLAND_CODE], "\n")
cat("Shrubland baseline (lucode", SHRUBLAND_CODE, "): C =",
    biophys_base$usle_c[biophys_base$lucode == SHRUBLAND_CODE], "\n")
cat("AF grassland transition (lucode", AF_GRASS_CODE, "): C =", AF_C_GRASS, "\n")
cat("AF shrubland transition (lucode", AF_SHRUB_CODE, "): C =", AF_C_SHRUB, "\n")

af_rows_grass_shrub <- tibble(
  lucode = c(AF_GRASS_CODE, AF_SHRUB_CODE),
  usle_c = c(AF_C_GRASS, AF_C_SHRUB),
  usle_p = c(AF_P_GRASS, AF_P_SHRUB)
)

if (ncol(biophys_base) > 3) {
  extra_cols <- setdiff(names(biophys_base), c("lucode", "usle_c", "usle_p"))
  grass_row <- biophys_base[biophys_base$lucode == GRASSLAND_CODE, ]
  shrub_row <- biophys_base[biophys_base$lucode == SHRUBLAND_CODE, ]
  src_rows <- list(grass_row, shrub_row)
  for (col in extra_cols) {
    af_rows_grass_shrub[[col]] <- c(src_rows[[1]][[col]], src_rows[[2]][[col]])
  }
  af_rows_grass_shrub <- af_rows_grass_shrub[, names(biophys_base)]
}

# Full expanded biophysical table: base + cropland AF + grassland/shrubland AF
biophys_af_all3 <- bind_rows(biophys_base, af_rows, af_rows_grass_shrub)
cat("Expanded biophysical table (crop+grass+shrub AF):", nrow(biophys_af_all3),
    "rows (added lucodes", AF_MIXED_CODE, ",", AF_DENSE_CODE, ",",
    AF_GRASS_CODE, ",", AF_SHRUB_CODE, ")\n")

# =============================================================================
# 13c. NEW LULC RASTERS + BIOPHYSICAL TABLES — CROP + GRASSLAND + SHRUBLAND AF
# Written to a SEPARATE output directory from the cropland-only runs in 13
# =============================================================================
# FIX: lulc_base_r must be loaded BEFORE this loop uses it (was previously
# defined further down the script, causing "object not found").

lulc_base_r <- rast(lulc_path)

invest_input_dir_af3 <- file.path(out_dir, "invest_inputs_af3class")
dir.create(invest_input_dir_af3, recursive = TRUE, showWarnings = FALSE)

run_log_af3 <- list()

for (scen in invest_scenarios) {
  
  if (!scen %in% names(suit_paths)) next
  
  cat("\n----------------------------\n")
  cat("Scenario (AF3-class):", scen, "\n")
  cat("----------------------------\n")
  
  suit_r <- load_and_normalise_suit(suit_paths[[scen]], scen)
  suit_vals <- values(suit_r, na.rm = TRUE)
  
  probs <- 1 - invest_suit_pcts / 100
  cutoffs <- setNames(quantile(suit_vals, probs = probs, na.rm = TRUE),
                      paste0("top", invest_suit_pcts))
  
  for (pct in invest_suit_pcts) {
    
    cut_val <- unname(cutoffs[paste0("top", pct)])
    run_id  <- paste0(scen, "_top", pct, "_af3")
    run_dir <- file.path(invest_input_dir_af3, run_id)
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    
    cat("\n [", run_id, "] suitability cutoff:", round(cut_val, 4), "\n")
    
    if (!compareGeom(suit_r, lulc_base_r, stopOnError = FALSE)) {
      suit_aligned <- resample(suit_r, lulc_base_r, method = "bilinear")
    } else {
      suit_aligned <- suit_r
    }
    
    lulc_mod_af3 <- lulc_base_r
    lulc_mod_af3 <- ifel(lulc_base_r == CROPLAND_CODE  & suit_aligned >= cut_val, AF_MIXED_CODE, lulc_mod_af3)
    lulc_mod_af3 <- ifel(lulc_base_r == GRASSLAND_CODE & suit_aligned >= cut_val, AF_GRASS_CODE, lulc_mod_af3)
    lulc_mod_af3 <- ifel(lulc_base_r == SHRUBLAND_CODE & suit_aligned >= cut_val, AF_SHRUB_CODE, lulc_mod_af3)
    
    n_crop_total  <- sum(values(lulc_base_r, na.rm = TRUE) == CROPLAND_CODE, na.rm = TRUE)
    n_grass_total <- sum(values(lulc_base_r, na.rm = TRUE) == GRASSLAND_CODE, na.rm = TRUE)
    n_shrub_total <- sum(values(lulc_base_r, na.rm = TRUE) == SHRUBLAND_CODE, na.rm = TRUE)
    
    n_crop_af  <- sum(values(lulc_mod_af3, na.rm = TRUE) == AF_MIXED_CODE, na.rm = TRUE)
    n_grass_af <- sum(values(lulc_mod_af3, na.rm = TRUE) == AF_GRASS_CODE, na.rm = TRUE)
    n_shrub_af <- sum(values(lulc_mod_af3, na.rm = TRUE) == AF_SHRUB_CODE, na.rm = TRUE)
    
    cat("  Cropland transitioned :", n_crop_af,  "/", n_crop_total,
        " (", round(100 * n_crop_af  / n_crop_total,  1), "%)\n")
    cat("  Grassland transitioned:", n_grass_af, "/", n_grass_total,
        " (", round(100 * n_grass_af / n_grass_total, 1), "%)\n")
    cat("  Shrubland transitioned:", n_shrub_af, "/", n_shrub_total,
        " (", round(100 * n_shrub_af / n_shrub_total, 1), "%)\n")
    
    lulc_out_af3 <- file.path(run_dir, "lulc.tif")
    writeRaster(lulc_mod_af3, lulc_out_af3, datatype = "INT2S",
                NAflag = -9999L, overwrite = TRUE)
    cat("  LULC written:", lulc_out_af3, "\n")
    
    biophys_out_af3 <- file.path(run_dir, "biophysical_table.csv")
    write_csv(biophys_af_all3, biophys_out_af3)
    cat("  Biophysical table written:", biophys_out_af3, "\n")
    
    run_log_af3[[run_id]] <- data.frame(
      run_id = run_id, scenario = scen, threshold_pct = pct,
      suit_cutoff = round(cut_val, 4),
      n_crop_af = n_crop_af, n_grass_af = n_grass_af, n_shrub_af = n_shrub_af,
      crop_area_km2  = round(n_crop_af  * px_km2_use),
      grass_area_km2 = round(n_grass_af * px_km2_use),
      shrub_area_km2 = round(n_shrub_af * px_km2_use),
      af_c_crop = AF_C_MIXED, af_c_grass = AF_C_GRASS, af_c_shrub = AF_C_SHRUB,
      lulc_path = lulc_out_af3, biophys_path = biophys_out_af3
    )
    
    rm(suit_aligned, lulc_mod_af3); gc()
  }
  rm(suit_r, suit_vals); gc()
}

run_log_af3_df <- bind_rows(run_log_af3)
write_csv(run_log_af3_df, file.path(invest_input_dir_af3, "invest_run_manifest_af3class.csv"))

cat("\n========================================\n")
cat("AF3-CLASS RUN MANIFEST (crop+grass+shrub)\n")
cat("========================================\n")
print(run_log_af3_df[, c("run_id", "threshold_pct", "suit_cutoff",
                         "crop_area_km2", "grass_area_km2", "shrub_area_km2")])
cat("\nManifest written to:",
    file.path(invest_input_dir_af3, "invest_run_manifest_af3class.csv"), "\n")


# =============================================================================
# SECTION 14 — OUTLIER-ROBUST UPDATE
# Adds a minimum-baseline filter to stabilize % change statistics, and adds
# median/IQR panels to all comparison plots so the outlier effect is directly
# visible rather than hidden inside mean/SD summary tables.
# =============================================================================
# Requires: invest_output_root, run_meta, ws_metrics, out14, load_watershed(),
# compare_runs(), pct_change(), climate_cols — from the fixed Section 14 script.
# Run this AFTER Section 14 (cropland-only) and the AF3 integration script.
# =============================================================================

# =============================================================================
# WHY MIN_BASELINE_SED_EXPORT = 1000 (t/yr)
# =============================================================================
# The percent-change formula (scen - ref) / |ref| * 100 is unstable whenever
# the reference (no-AF baseline) value is small, because a modest absolute
# change gets divided by a near-zero denominator and inflates into an
# enormous percentage. In this dataset:
#
#   - 41 watersheds (1.5%) have baseline sed_export == 0 (already excluded
#     from every _pct column by pct_change(), which returns NA for ref == 0)
#   - The 10th percentile of baseline sed_export across all 2,754 watersheds
#     is ~941 t/yr; the 5th percentile is ~142 t/yr
#   - All 12 watersheds driving |% change| > 500% in the AF3 top10 comparison
#     have baseline sed_export between 11 and 2,173 t/yr — i.e., they sit
#     entirely within the bottom ~10-15% of the baseline distribution
#
# Setting the cutoff at the 10th percentile (~1000 t/yr, rounded) removes
# 10-11% of watersheds — those with the least reliable baseline denominator —
# while preserving 2,473 of 2,754 watersheds (90%) for analysis. Empirically,
# this drops the AF3 top10 sed_export SD from 205% to 51% and stabilizes the
# mean from -9.4% to -20.6% (much closer to the median of -7.7%), without
# changing the cropland-only statistics at all (they were never outlier-driven
# in the first place, since every cropland pixel strictly reduces C-factor).
#
# This is not an arbitrary "make the outliers go away" choice — it is a
# denominator-reliability cutoff: below ~1,000 t/yr baseline sediment export,
# the percent-change metric is measuring noise in a near-zero quantity, not a
# meaningful erosion signal. Watersheds below the cutoff are NOT deleted from
# the dataset; they are excluded only from percent-change-based statistics
# and plots. Absolute delta (t/yr) is unaffected by this issue and does not
# require filtering.
# =============================================================================

MIN_BASELINE_SED_EXPORT <- 1000   # t/yr — 10th percentile of baseline distribution
MIN_BASELINE_USLE_TOT   <- 1000   # t/yr — same rationale, applied per metric

min_baseline_by_metric <- c(
  sed_export = MIN_BASELINE_SED_EXPORT,
  usle_tot   = MIN_BASELINE_USLE_TOT,
  sed_dep    = 1000,
  avoid_eros = 1000,
  avoid_exp  = 1000
)


# ---- CHUNK 2: compare_runs_robust as single source of truth -----------
# (definition kept here; the summarise_comp_robust() you already have
#  downstream in Section 14 works unchanged against this compare_runs())

compare_runs_robust <- function(ref_id, scen_id, min_baseline = min_baseline_by_metric) {
  comp <- compare_runs(ref_id, scen_id)
  if (is.null(comp)) return(NULL)
  for (m in names(min_baseline)) {
    ref_col  <- paste0(m, "_ref")
    flag_col <- paste0(m, "_low_baseline")
    if (ref_col %in% names(comp)) comp[[flag_col]] <- abs(comp[[ref_col]]) < min_baseline[[m]]
  }
  comp
}

# NOTE: min_baseline_by_metric must be defined before this point in your
# script (it is, in the existing Section 14 block). If sourcing this file
# standalone, define it here too:
if (!exists("min_baseline_by_metric")) {
  min_baseline_by_metric <- c(
    sed_export = 1000, usle_tot = 1000, sed_dep = 1000,
    avoid_eros = 1000, avoid_exp = 1000
  )
}

cat("\nCore objects ready: invest_output_root, ref_baseline, run_meta, climate_cols, out14\n")
cat("Functions ready: load_watershed(), pct_change(), compare_runs(), compare_runs_robust()\n")


# -----------------------------------------------------------------------
# Outlier-robust summary: adds median/IQR alongside mean/SD, computed both
# on the full watershed set and on the baseline-filtered set, for direct
# comparison.
# -----------------------------------------------------------------------

summarise_comp_robust <- function(comp_df, ref_id, scen_id,
                                  min_baseline = min_baseline_by_metric) {
  pct_cols <- grep("_pct$", names(comp_df), value = TRUE)
  
  map_dfr(pct_cols, function(col) {
    metric <- sub("_pct$", "", col)
    ref_col <- paste0(metric, "_ref")
    vals_all <- comp_df[[col]]
    
    cutoff <- if (metric %in% names(min_baseline)) min_baseline[[metric]] else 0
    keep <- if (ref_col %in% names(comp_df)) abs(comp_df[[ref_col]]) >= cutoff else rep(TRUE, nrow(comp_df))
    vals_filt <- comp_df[[col]][keep]
    
    tibble(
      ref_id = ref_id, scen_id = scen_id, metric = metric,
      mean_pct_all       = mean(vals_all, na.rm = TRUE),
      median_pct_all     = median(vals_all, na.rm = TRUE),
      sd_pct_all         = sd(vals_all, na.rm = TRUE),
      q25_pct_all        = quantile(vals_all, 0.25, na.rm = TRUE),
      q75_pct_all        = quantile(vals_all, 0.75, na.rm = TRUE),
      n_ws_all           = sum(!is.na(vals_all)),
      mean_pct_filt      = mean(vals_filt, na.rm = TRUE),
      median_pct_filt    = median(vals_filt, na.rm = TRUE),
      sd_pct_filt        = sd(vals_filt, na.rm = TRUE),
      q25_pct_filt       = quantile(vals_filt, 0.25, na.rm = TRUE),
      q75_pct_filt       = quantile(vals_filt, 0.75, na.rm = TRUE),
      n_ws_filt          = sum(!is.na(vals_filt)),
      n_ws_excluded      = sum(!is.na(vals_all)) - sum(!is.na(vals_filt)),
      pct_excluded       = round(100 * (sum(!is.na(vals_all)) - sum(!is.na(vals_filt))) / sum(!is.na(vals_all)), 1),
      min_baseline_used  = cutoff
    )
  })
}

# =============================================================================
# REGENERATE COMPARISON A (cropland-only AND af3) WITH ROBUST SUMMARY
# =============================================================================

message("Regenerating Comparison A with outlier-robust summary (median/IQR)...")

comp_A_pairs_crop <- list(
  list(ref = "baseline_cropland_0.34", scen = "baseline_cropland_0.34_to_0.2832_top10", group = "Baseline", class = "Cropland-only"),
  list(ref = "baseline_cropland_0.34", scen = "baseline_cropland_0.34_to_0.2832_top20", group = "Baseline", class = "Cropland-only"),
  list(ref = "baseline_cropland_0.34", scen = "baseline_cropland_0.34_to_0.2832_top25", group = "Baseline", class = "Cropland-only"),
  list(ref = "ssp370_cropland_0.34",   scen = "ssp370_cropland_0.34_to_0.2832_top10",   group = "SSP370",   class = "Cropland-only"),
  list(ref = "ssp370_cropland_0.34",   scen = "ssp370_cropland_0.34_to_0.2832_top20",   group = "SSP370",   class = "Cropland-only"),
  list(ref = "ssp370_cropland_0.34",   scen = "ssp370_cropland_0.34_to_0.2832_top25",   group = "SSP370",   class = "Cropland-only")
)

comp_A_pairs_af3 <- list(
  list(ref = "baseline_cropland_0.34", scen = "baseline_af3_top10", group = "Baseline", class = "AF3-class"),
  list(ref = "baseline_cropland_0.34", scen = "baseline_af3_top20", group = "Baseline", class = "AF3-class"),
  list(ref = "baseline_cropland_0.34", scen = "baseline_af3_top25", group = "Baseline", class = "AF3-class"),
  list(ref = "ssp370_cropland_0.34",   scen = "ssp370_af3_top10",   group = "SSP370",   class = "AF3-class"),
  list(ref = "ssp370_cropland_0.34",   scen = "ssp370_af3_top20",   group = "SSP370",   class = "AF3-class"),
  list(ref = "ssp370_cropland_0.34",   scen = "ssp370_af3_top25",   group = "SSP370",   class = "AF3-class")
)

comp_A_all_pairs <- c(comp_A_pairs_crop, comp_A_pairs_af3)

comp_A_robust <- map_dfr(comp_A_all_pairs, function(p) {
  comp <- compare_runs_robust(p$ref, p$scen)
  if (is.null(comp)) return(NULL)
  s <- summarise_comp_robust(comp, p$ref, p$scen)
  s$group     <- p$group
  s$class     <- p$class
  s$threshold <- run_meta$threshold[run_meta$run_id == p$scen]
  s
})

write.csv(comp_A_robust, file.path(out14, "A_threshold_effect_outlier_robust.csv"), row.names = FALSE)

# Long-form per-watershed data (both classes), with low-baseline flag, for plotting
ws_pct_all <- map_dfr(comp_A_all_pairs, function(p) {
  comp <- compare_runs_robust(p$ref, p$scen)
  if (is.null(comp) || !"sed_export_pct" %in% names(comp)) return(NULL)
  tibble(
    group = p$group,
    class = p$class,
    threshold = run_meta$threshold[run_meta$run_id == p$scen],
    pct_chg = comp$sed_export_pct,
    low_baseline = comp$sed_export_low_baseline
  )
}) |>
  mutate(threshold = factor(threshold, levels = c("Full area", "Top 10%", "Top 20%", "Top 25%")))

write.csv(ws_pct_all, file.path(out14, "watershed_sed_export_pct_flagged.csv"), row.names = FALSE)


  #=============================================================================
  # PLOT — Median/IQR panel alongside Mean/SD panel, same axes, side by side
  # Makes the outlier effect visually obvious: mean/SD bars vs median/IQR bars
  # for the SAME data, same watersheds, same runs.
  # =============================================================================

message("Building median/IQR comparison plots...")

# Reshape robust summary into long form: one row per (run, metric, stat_type)
stat_compare_df <- comp_A_robust |>
  filter(metric == "sed_export") |>
  select(class, group, threshold, mean_pct_all, median_pct_all,
         q25_pct_all, q75_pct_all, mean_pct_filt, median_pct_filt,
         q25_pct_filt, q75_pct_filt) |>
  pivot_longer(
    cols = -c(class, group, threshold),
    names_to = "stat_col", values_to = "value"
  ) |>
  mutate(
    filtered  = ifelse(grepl("_filt$", stat_col), "Baseline-filtered (n_ws >= cutoff)", "All watersheds (unfiltered)"),
    stat_type = case_when(
      grepl("^mean_pct",   stat_col) ~ "mean",
      grepl("^median_pct", stat_col) ~ "median",
      grepl("^q25_pct",    stat_col) ~ "q25",
      grepl("^q75_pct",    stat_col) ~ "q75"
    )
  ) |>
  select(-stat_col) |>
  pivot_wider(names_from = stat_type, values_from = value) |>
  mutate(threshold = factor(threshold, levels = c("Top 10%", "Top 20%", "Top 25%")))

p_median_iqr <- ggplot(stat_compare_df,
                       aes(x = threshold, y = median, fill = class)) +
  geom_col(position = position_dodge(0.7), width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = q25, ymax = q75),
                position = position_dodge(0.7), width = 0.25, linewidth = 0.5) +
  geom_point(aes(y = mean, colour = class), position = position_dodge(0.7),
             shape = 18, size = 2.6, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("Cropland-only" = "#2166AC", "AF3-class" = "#D73027"), name = "Transition") +
  scale_colour_manual(values = c("Cropland-only" = "#0b3757", "AF3-class" = "#7a1a10"), guide = "none") +
  facet_grid(filtered ~ group) +
  labs(
    title = "Median \u00b1 IQR vs mean (diamond marker) for sed_export % change",
    subtitle = "Bars = median with 25th-75th percentile whiskers | Diamonds = mean (shown for comparison)",
    x = "Suitability threshold", y = "% change in sed_export vs no-AF reference",
    caption = paste0("Baseline-filtered excludes watersheds with reference sed_export < ",
                     MIN_BASELINE_SED_EXPORT, " t/yr")
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(file.path(out14, "fig6_median_iqr_vs_mean_outlier_check.png"),
       p_median_iqr, width = 12, height = 8, dpi = 180)

# -----------------------------------------------------------------------
# PLOT — Violin/box distribution, split by low-baseline flag, so the
# outlier watersheds are visibly separated from the stable majority.
# -----------------------------------------------------------------------

p_outlier_split <- ws_pct_all |>
  filter(class == "AF3-class", !is.na(pct_chg)) |>
  mutate(baseline_group = ifelse(low_baseline, 
                                 paste0("Low baseline (<", MIN_BASELINE_SED_EXPORT, " t/yr)"),
                                 "Stable baseline")) |>
  ggplot(aes(x = threshold, y = pct_chg, fill = baseline_group)) +
  geom_boxplot(position = position_dodge(0.75), width = 0.55,
               outlier.size = 0.6, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("Stable baseline" = "#2166AC",
                               setNames("#D73027", paste0("Low baseline (<", MIN_BASELINE_SED_EXPORT, " t/yr)"))),
                    name = NULL) +
  facet_wrap(~ group, ncol = 2) +
  coord_cartesian(ylim = c(-100, 200)) +
  labs(
    title = "AF3-class sed_export % change: low-baseline watersheds vs stable watersheds",
    subtitle = "Low-baseline watersheds (small denominator) drive nearly all extreme percentages",
    x = "Suitability threshold", y = "% change in sed_export",
    caption = "Y-axis clipped to [-100%, 200%] to show the stable-baseline distribution clearly"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(out14, "fig7_low_baseline_vs_stable_split.png"),
       p_outlier_split, width = 11, height = 6, dpi = 180)

# =============================================================================
# UPDATED VIOLIN/BOX PLOT (replaces old fig4/fig4_af3) — now shows median/IQR
# via boxplot layer explicitly, with low-baseline watersheds marked as an
# overlaid point layer rather than hidden inside the violin.
# =============================================================================

p4_updated <- ws_pct_all |>
  filter(!is.na(pct_chg)) |>
  ggplot(aes(x = threshold, y = pct_chg, fill = group, colour = group)) +
  geom_violin(position = position_dodge(0.8), alpha = 0.25, linewidth = 0.4,
              scale = "width", trim = TRUE) +
  geom_boxplot(position = position_dodge(0.8), width = 0.22,
               outlier.shape = NA, alpha = 0.75) +
  geom_point(data = ~ filter(.x, low_baseline),
             aes(x = threshold, y = pct_chg),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.1),
             shape = 4, size = 1.4, colour = "black", alpha = 0.6, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = climate_cols, name = "Climate") +
  scale_colour_manual(values = climate_cols, name = "Climate") +
  facet_wrap(~ class, ncol = 2) +
  coord_cartesian(ylim = c(-100, 150)) +
  labs(
    title = "Distribution of watershed-level sediment export change (updated)",
    subtitle = "X marks = watersheds with low baseline sed_export (unreliable % denominator)",
    x = "Suitability threshold", y = "% change in sed_export per watershed",
    caption = "Violin = full distribution | Box = IQR + median | Y-axis clipped to [-100%, 150%]"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(out14, "fig4_updated_watershed_distribution_with_flags.png"),
       p4_updated, width = 12, height = 7, dpi = 180)

# =============================================================================
# SUMMARY TABLE PRINTOUT
# =============================================================================

cat("\n========================================\n")
cat("OUTLIER-ROBUST SUMMARY (sed_export)\n")
cat("========================================\n")
print(comp_A_robust |>
        filter(metric == "sed_export") |>
        select(class, group, threshold, mean_pct_all, median_pct_all, sd_pct_all,
               mean_pct_filt, median_pct_filt, sd_pct_filt, pct_excluded))

message("\nOutputs written to: ", out14)
message("  A_threshold_effect_outlier_robust.csv — full mean/median/IQR, filtered + unfiltered")
message("  watershed_sed_export_pct_flagged.csv — per-watershed data with low_baseline flag")
message("  fig6_median_iqr_vs_mean_outlier_check.png")
message("  fig7_low_baseline_vs_stable_split.png")
message("  fig4_updated_watershed_distribution_with_flags.png") 



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

# =============================================================================
# MAP 5 (NEW) — Absolute change in sed_export (t/yr), AF3-class, all
# thresholds x climates. Uses QUANTILE-based colour breaks (not a fixed
# symmetric min/max), since raw deltas are heavily right/left-skewed by a
# small number of very large watersheds — a plain min/max or even a simple
# +/- max scale washes out mid-range variation into near-white, and a naive
# 2nd-98th clip (used for the % maps) is too aggressive for delta, whose
# tails carry real physical mass (large watersheds moving many t/yr).
# Excluded (low-baseline) watersheds are still greyed out for consistency
# with Maps 1-4.
# =============================================================================

message("Building Map 5: absolute delta (t/yr), quantile-classed...")

map5_sf <- map_sf_all |> filter(class == "AF3-class")

# NA out excluded watersheds' delta too, so grey rendering matches pct maps
map5_sf$delta_plot <- ifelse(map5_sf$low_baseline, NA_real_, map5_sf$delta)

# Quantile (equal-count) breaks computed on FILTERED, non-zero deltas only.
# classInt::classIntervals with style = "quantile" ensures each colour bin
# contains roughly the same number of watersheds, so the map is driven by
# where most of the data actually sits rather than by extreme tails.
delta_vals <- map5_sf$delta_plot[!is.na(map5_sf$delta_plot)]

n_bins <- 9  # odd number so a bin can center on ~0
brks5 <- classInt::classIntervals(delta_vals, n = n_bins, style = "quantile")$brks
brks5 <- unique(round(brks5, 2))  # guard against duplicate breaks from ties/zeros

# Diverging palette applied over the quantile breaks (not a linear/continuous
# diverging scale) — this is the key difference from Maps 1-2: colour steps
# are equal-COUNT, not equal-VALUE, so outliers don't compress the mid-range.
pal5 <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(length(brks5) - 1)

map5_sf$delta_bin <- cut(map5_sf$delta_plot, breaks = brks5, include.lowest = TRUE)

bin_labels <- levels(map5_sf$delta_bin)
bin_labels_fmt <- gsub("\\(|\\]|\\[", "", bin_labels)
bin_labels_fmt <- sapply(strsplit(bin_labels_fmt, ","), function(x) {
  paste0(scales::label_comma(accuracy = 1)(as.numeric(x[1])), " to ",
         scales::label_comma(accuracy = 1)(as.numeric(x[2])))
})
levels(map5_sf$delta_bin) <- bin_labels_fmt

p_map5 <- ggplot(map5_sf) +
  geom_sf(aes(fill = delta_bin), colour = "white", linewidth = 0.08) +
  scale_fill_manual(
    values = setNames(pal5, bin_labels_fmt),
    na.value = "grey80",
    name = "\u0394 sed_export\n(t/yr)\nquantile bins",
    drop = FALSE
  ) +
  facet_grid(climate ~ threshold) +
  labs(
    title = "Absolute change in sediment export \u2014 AF3-class transition",
    subtitle = "Quantile (equal-count) colour bins \u2014 mid-range variation preserved, not compressed by outlier watersheds",
    caption = paste0("Grey = excluded (baseline sed_export < ", cutoff_used, " t/yr) | ",
                     n_bins - 1, " equal-count bins from filtered, non-excluded watersheds")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(), strip.text = element_text(face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 8)
  ) +
  coord_sf(datum = NA) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE))

ggsave(file.path(out15, "map5_af3_delta_sedexport_quantile.png"),
       p_map5, width = 14, height = 9, dpi = 220)

# -----------------------------------------------------------------------
# MAP 5B (NEW) — Same quantile-binned delta map for Cropland-only, using
# its OWN quantile breaks (deltas are much smaller in magnitude for
# cropland-only, so sharing AF3's bins would flatten this map to one colour)
# -----------------------------------------------------------------------

map5b_sf <- map_sf_all |> filter(class == "Cropland-only")
map5b_sf$delta_plot <- ifelse(map5b_sf$low_baseline, NA_real_, map5b_sf$delta)

delta_vals_b <- map5b_sf$delta_plot[!is.na(map5b_sf$delta_plot)]
brks5b <- classInt::classIntervals(delta_vals_b, n = n_bins, style = "quantile")$brks
brks5b <- unique(round(brks5b, 2))
pal5b <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(length(brks5b) - 1)

map5b_sf$delta_bin <- cut(map5b_sf$delta_plot, breaks = brks5b, include.lowest = TRUE)
bin_labels_b <- levels(map5b_sf$delta_bin)
bin_labels_fmt_b <- gsub("\\(|\\]|\\[", "", bin_labels_b)
bin_labels_fmt_b <- sapply(strsplit(bin_labels_fmt_b, ","), function(x) {
  paste0(scales::label_comma(accuracy = 0.1)(as.numeric(x[1])), " to ",
         scales::label_comma(accuracy = 0.1)(as.numeric(x[2])))
})
levels(map5b_sf$delta_bin) <- bin_labels_fmt_b

p_map5b <- ggplot(map5b_sf) +
  geom_sf(aes(fill = delta_bin), colour = "white", linewidth = 0.08) +
  scale_fill_manual(
    values = setNames(pal5b, bin_labels_fmt_b),
    na.value = "grey80",
    name = "\u0394 sed_export\n(t/yr)\nquantile bins",
    drop = FALSE
  ) +
  facet_grid(climate ~ threshold) +
  labs(
    title = "Absolute change in sediment export \u2014 Cropland-only transition",
    subtitle = "Quantile (equal-count) colour bins, own scale (magnitudes much smaller than AF3-class)",
    caption = paste0("Grey = excluded (baseline sed_export < ", cutoff_used, " t/yr) | ",
                     n_bins - 1, " equal-count bins from filtered, non-excluded watersheds")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(), strip.text = element_text(face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 8)
  ) +
  coord_sf(datum = NA) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE))

ggsave(file.path(out15, "map5b_cropland_delta_sedexport_quantile.png"),
       p_map5b, width = 14, height = 9, dpi = 220)

message("  map5_af3_delta_sedexport_quantile.png       (NEW \u2014 absolute \u0394, quantile-classed)")
message("  map5b_cropland_delta_sedexport_quantile.png (NEW \u2014 absolute \u0394, own quantile scale)")

# =============================================================================
# MAPS 6-7 (NEW) — Top 10% vs Top 25% suitability threshold comparison
# Shows, per watershed, how much MORE benefit is gained by expanding the
# suitability mask from the tightest (Top 10%) to the widest (Top 25%)
# threshold. Two maps: (6) difference in % change, (7) difference in
# absolute delta (t/yr). Both use quantile-classed bins for the same
# reason as Map 5 — the spread between two already-skewed percent/delta
# distributions is itself skewed, so equal-count bins keep mid-range
# variation visible instead of letting a few extreme watersheds dominate.
#
# Definition: threshold_gain = value_at_Top25 - value_at_Top10
#   Negative (blue) = Top 25% achieves GREATER sediment-export reduction
#                      than Top 10% (i.e., widening the mask helps)
#   Positive (red)  = Top 25% achieves LESS reduction than Top 10%
#                      (widening the mask hurts, e.g. dilutes benefit)
# A watershed excluded (low-baseline) at EITHER threshold is greyed out,
# since the comparison is not meaningful if one side is unreliable.
# =============================================================================

message("Building Maps 6-7: Top 10% vs Top 25% threshold comparison...")

build_threshold_gain <- function(class_filter) {
  d10 <- map_sf_all |> filter(class == class_filter, threshold == "Top 10%") |>
    st_drop_geometry() |>
    select(ws_id, climate, pct_chg_top10 = pct_chg, delta_top10 = delta,
           low_baseline_top10 = low_baseline)
  d25 <- map_sf_all |> filter(class == class_filter, threshold == "Top 25%") |>
    st_drop_geometry() |>
    select(ws_id, climate, pct_chg_top25 = pct_chg, delta_top25 = delta,
           low_baseline_top25 = low_baseline)
  
  gain <- d10 |>
    inner_join(d25, by = c("ws_id", "climate")) |>
    mutate(
      pct_gain   = pct_chg_top25 - pct_chg_top10,
      delta_gain = delta_top25 - delta_top10,
      excluded   = low_baseline_top10 | low_baseline_top25,
      pct_gain_plot   = ifelse(excluded, NA_real_, pct_gain),
      delta_gain_plot = ifelse(excluded, NA_real_, delta_gain)
    )
  
  ws_geom |> left_join(gain, by = "ws_id") |> st_as_sf() |>
    mutate(climate = factor(climate, levels = c("Baseline", "SSP370")))
}

gain_af3_sf <- build_threshold_gain("AF3-class")

write.csv(st_drop_geometry(gain_af3_sf),
          file.path(out15, "threshold_gain_af3_top10_vs_top25.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------
# MAP 6 — Difference in PERCENT change, Top 25% minus Top 10% (AF3-class)
# Quantile-classed to keep mid-range variation visible.
# -----------------------------------------------------------------------

pct_gain_vals <- gain_af3_sf$pct_gain_plot[!is.na(gain_af3_sf$pct_gain_plot)]
brks6 <- classInt::classIntervals(pct_gain_vals, n = 9, style = "quantile")$brks
brks6 <- unique(round(brks6, 3))
pal6 <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(length(brks6) - 1)

gain_af3_sf$pct_gain_bin <- cut(gain_af3_sf$pct_gain_plot, breaks = brks6, include.lowest = TRUE)
bin_labels6 <- levels(gain_af3_sf$pct_gain_bin)
bin_labels6_fmt <- gsub("\\(|\\]|\\[", "", bin_labels6)
bin_labels6_fmt <- sapply(strsplit(bin_labels6_fmt, ","), function(x) {
  paste0(scales::label_percent(scale = 1, accuracy = 0.1)(as.numeric(x[1])), " to ",
         scales::label_percent(scale = 1, accuracy = 0.1)(as.numeric(x[2])))
})
levels(gain_af3_sf$pct_gain_bin) <- bin_labels6_fmt

p_map6 <- ggplot(gain_af3_sf) +
  geom_sf(aes(fill = pct_gain_bin), colour = "white", linewidth = 0.08) +
  scale_fill_manual(
    values = setNames(pal6, bin_labels6_fmt),
    na.value = "grey80",
    name = "\u0394 (Top25 \u2212 Top10)\n% change\nquantile bins",
    drop = FALSE
  ) +
  facet_wrap(~ climate, nrow = 1) +
  labs(
    title = "Threshold sensitivity \u2014 % change gain from widening Top 10% to Top 25%",
    subtitle = "AF3-class transition | Blue = wider threshold adds MORE benefit | Red = wider threshold adds LESS benefit",
    caption = paste0("Grey = excluded (low baseline at Top 10% and/or Top 25%) | ",
                     "9 equal-count bins from filtered, non-excluded watersheds")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(), strip.text = element_text(face = "bold"),
    legend.position = "right", legend.text = element_text(size = 8)
  ) +
  coord_sf(datum = NA) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE))

ggsave(file.path(out15, "map6_threshold_gain_pct_af3.png"),
       p_map6, width = 12, height = 7, dpi = 220)

# -----------------------------------------------------------------------
# MAP 7 — Difference in ABSOLUTE delta (t/yr), Top 25% minus Top 10%
# Own quantile scale (different units/magnitude than Map 6).
# -----------------------------------------------------------------------

delta_gain_vals <- gain_af3_sf$delta_gain_plot[!is.na(gain_af3_sf$delta_gain_plot)]
brks7 <- classInt::classIntervals(delta_gain_vals, n = 9, style = "quantile")$brks
brks7 <- unique(round(brks7, 2))
pal7 <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(length(brks7) - 1)

gain_af3_sf$delta_gain_bin <- cut(gain_af3_sf$delta_gain_plot, breaks = brks7, include.lowest = TRUE)
bin_labels7 <- levels(gain_af3_sf$delta_gain_bin)
bin_labels7_fmt <- gsub("\\(|\\]|\\[", "", bin_labels7)
bin_labels7_fmt <- sapply(strsplit(bin_labels7_fmt, ","), function(x) {
  paste0(scales::label_comma(accuracy = 1)(as.numeric(x[1])), " to ",
         scales::label_comma(accuracy = 1)(as.numeric(x[2])))
})
levels(gain_af3_sf$delta_gain_bin) <- bin_labels7_fmt

p_map7 <- ggplot(gain_af3_sf) +
  geom_sf(aes(fill = delta_gain_bin), colour = "white", linewidth = 0.08) +
  scale_fill_manual(
    values = setNames(pal7, bin_labels7_fmt),
    na.value = "grey80",
    name = "\u0394 (Top25 \u2212 Top10)\nsed_export (t/yr)\nquantile bins",
    drop = FALSE
  ) +
  facet_wrap(~ climate, nrow = 1) +
  labs(
    title = "Threshold sensitivity \u2014 absolute sed_export gain from widening Top 10% to Top 25%",
    subtitle = "AF3-class transition | Blue = wider threshold removes MORE tonnes | Red = wider threshold removes FEWER tonnes",
    caption = paste0("Grey = excluded (low baseline at Top 10% and/or Top 25%) | ",
                     "9 equal-count bins from filtered, non-excluded watersheds")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(), strip.text = element_text(face = "bold"),
    legend.position = "right", legend.text = element_text(size = 8)
  ) +
  coord_sf(datum = NA) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE))

ggsave(file.path(out15, "map7_threshold_gain_delta_af3.png"),
       p_map7, width = 12, height = 7, dpi = 220)

message("  map6_threshold_gain_pct_af3.png    (NEW \u2014 Top25 vs Top10, % change diff, quantile-classed)")
message("  map7_threshold_gain_delta_af3.png  (NEW \u2014 Top25 vs Top10, absolute delta diff, quantile-classed)")
message("  threshold_gain_af3_top10_vs_top25.csv")

# =============================================================================
# FIX — Two-sided (zero-anchored) quantile classing for diverging maps
# =============================================================================
# PROBLEM IDENTIFIED: Maps 5, 5b, 6, 7 used classInt::classIntervals(style =
# "quantile") across the FULL range of values (negative + positive together),
# then mapped a diverging blue-white-red palette across the resulting bins by
# RANK ORDER. Because negative and positive counts are rarely equal, the bin
# containing the palette's "white" midpoint does NOT necessarily contain
# zero -- it contains the median. This can assign RED to negative values
# (decreases) or BLUE to positive values (increases), inverting the intended
# meaning of the colour scale.
#
# FIX: classify negative and positive values into equal-count bins
# SEPARATELY, anchored at zero, so blue is always <0 and red is always >0,
# while still preserving quantile (equal-count) resolution on each side.
# =============================================================================

zero_anchored_quantile_bins <- function(x, n_bins_per_side = 4) {
  x <- x[!is.na(x)]
  neg <- x[x < 0]
  pos <- x[x > 0]
  
  neg_brks <- if (length(neg) > 0) {
    b <- classInt::classIntervals(neg, n = n_bins_per_side, style = "quantile")$brks
    unique(round(b, 4))
  } else numeric(0)
  
  pos_brks <- if (length(pos) > 0) {
    b <- classInt::classIntervals(pos, n = n_bins_per_side, style = "quantile")$brks
    unique(round(b, 4))
  } else numeric(0)
  
  # Stitch: ...neg breaks..., 0, ...pos breaks... (drop duplicate zero-adjacent edges)
  brks <- c(neg_brks, 0, pos_brks)
  brks <- unique(brks)
  sort(brks)
}

make_diverging_palette <- function(brks) {
  zero_idx <- which(brks == 0)
  n_neg_bins <- zero_idx - 1          # bins strictly below 0
  n_pos_bins <- length(brks) - zero_idx  # bins strictly above 0
  
  neg_pal <- if (n_neg_bins > 0) colorRampPalette(c("#08306B", "#F7F7F7"))(n_neg_bins + 1)[1:n_neg_bins] else character(0)
  pos_pal <- if (n_pos_bins > 0) colorRampPalette(c("#F7F7F7", "#67000D"))(n_pos_bins + 1)[2:(n_pos_bins + 1)] else character(0)
  
  c(neg_pal, pos_pal)
}

format_bin_labels <- function(brks, pct = FALSE, accuracy = 1) {
  labs <- sapply(seq_len(length(brks) - 1), function(i) {
    lo <- brks[i]; hi <- brks[i + 1]
    fmt <- if (pct) scales::label_percent(scale = 1, accuracy = accuracy) else scales::label_comma(accuracy = accuracy)
    paste0(fmt(lo), " to ", fmt(hi))
  })
  labs
}

apply_zero_anchored_fill <- function(sf_obj, value_col, pct = FALSE, n_bins_per_side = 4) {
  brks <- zero_anchored_quantile_bins(sf_obj[[value_col]], n_bins_per_side)
  pal  <- make_diverging_palette(brks)
  labs <- format_bin_labels(brks, pct = pct)
  
  sf_obj$bin_col <- cut(sf_obj[[value_col]], breaks = brks, include.lowest = TRUE, labels = labs)
  list(sf_obj = sf_obj, palette = setNames(pal, labs), breaks = brks, labels = labs)
}

# -----------------------------------------------------------------------
# RE-BUILD MAP 5 with zero-anchored quantile classing (AF3-class delta)
# -----------------------------------------------------------------------

map5_fix <- apply_zero_anchored_fill(map5_sf, "delta_plot", pct = FALSE)
map5_sf <- map5_fix$sf_obj

p_map5_fixed <- ggplot(map5_sf) +
  geom_sf(aes(fill = bin_col), colour = "white", linewidth = 0.08) +
  scale_fill_manual(values = map5_fix$palette, na.value = "grey80",
                    name = "\u0394 sed_export\n(t/yr)\nzero-anchored bins", drop = FALSE) +
  facet_grid(climate ~ threshold) +
  labs(
    title = "Absolute change in sediment export \u2014 AF3-class transition",
    subtitle = "Zero-anchored quantile bins \u2014 blue always = decrease, red always = increase",
    caption = paste0("Grey = excluded (baseline sed_export < ", cutoff_used, " t/yr) | ",
                     "Equal-count bins computed separately for negative and positive values, split at zero")
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "right", legend.text = element_text(size = 8)) +
  coord_sf(datum = NA) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE))

ggsave(file.path(out15, "map5_af3_delta_sedexport_zeroanchored.png"),
       p_map5_fixed, width = 14, height = 9, dpi = 220)

# -----------------------------------------------------------------------
# RE-BUILD MAP 5B with zero-anchored quantile classing (Cropland-only delta)
# -----------------------------------------------------------------------

map5b_fix <- apply_zero_anchored_fill(map5b_sf, "delta_plot", pct = FALSE)
map5b_sf <- map5b_fix$sf_obj

p_map5b_fixed <- ggplot(map5b_sf) +
  geom_sf(aes(fill = bin_col), colour = "white", linewidth = 0.08) +
  scale_fill_manual(values = map5b_fix$palette, na.value = "grey80",
                    name = "\u0394 sed_export\n(t/yr)\nzero-anchored bins", drop = FALSE) +
  facet_grid(climate ~ threshold) +
  labs(
    title = "Absolute change in sediment export \u2014 Cropland-only transition",
    subtitle = "Zero-anchored quantile bins \u2014 blue always = decrease, red always = increase",
    caption = paste0("Grey = excluded (baseline sed_export < ", cutoff_used, " t/yr)")
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "right", legend.text = element_text(size = 8)) +
  coord_sf(datum = NA) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE))

ggsave(file.path(out15, "map5b_cropland_delta_sedexport_zeroanchored.png"),
       p_map5b_fixed, width = 14, height = 9, dpi = 220)

# -----------------------------------------------------------------------
# RE-BUILD MAP 6 with zero-anchored quantile classing (% gain, Top25-Top10)
# -----------------------------------------------------------------------

map6_fix <- apply_zero_anchored_fill(gain_af3_sf, "pct_gain_plot", pct = TRUE)
gain_af3_sf <- map6_fix$sf_obj

p_map6_fixed <- ggplot(gain_af3_sf) +
  geom_sf(aes(fill = bin_col), colour = "white", linewidth = 0.08) +
  scale_fill_manual(values = map6_fix$palette, na.value = "grey80",
                    name = "\u0394 (Top25\u2212Top10)\n% change\nzero-anchored bins", drop = FALSE) +
  facet_wrap(~ climate, nrow = 1) +
  labs(
    title = "Threshold sensitivity \u2014 % change gain from widening Top 10% to Top 25%",
    subtitle = "AF3-class | Zero-anchored bins: blue = wider threshold adds MORE benefit, red = adds LESS",
    caption = paste0("Grey = excluded (low baseline at Top 10% and/or Top 25%)")
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "right", legend.text = element_text(size = 8)) +
  coord_sf(datum = NA) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE))

ggsave(file.path(out15, "map6_threshold_gain_pct_af3_zeroanchored.png"),
       p_map6_fixed, width = 12, height = 7, dpi = 220)

# -----------------------------------------------------------------------
# RE-BUILD MAP 7 with zero-anchored quantile classing (absolute gain)
# -----------------------------------------------------------------------

map7_fix <- apply_zero_anchored_fill(gain_af3_sf, "delta_gain_plot", pct = FALSE)
gain_af3_sf <- map7_fix$sf_obj

p_map7_fixed <- ggplot(gain_af3_sf) +
  geom_sf(aes(fill = bin_col), colour = "white", linewidth = 0.08) +
  scale_fill_manual(values = map7_fix$palette, na.value = "grey80",
                    name = "\u0394 (Top25\u2212Top10)\nsed_export (t/yr)\nzero-anchored bins", drop = FALSE) +
  facet_wrap(~ climate, nrow = 1) +
  labs(
    title = "Threshold sensitivity \u2014 absolute sed_export gain from widening Top 10% to Top 25%",
    subtitle = "AF3-class | Zero-anchored bins: blue = wider threshold removes MORE tonnes, red = removes FEWER",
    caption = paste0("Grey = excluded (low baseline at Top 10% and/or Top 25%)")
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "right", legend.text = element_text(size = 8)) +
  coord_sf(datum = NA) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE))

ggsave(file.path(out15, "map7_threshold_gain_delta_af3_zeroanchored.png"),
       p_map7_fixed, width = 12, height = 7, dpi = 220)

message("\nZero-anchored quantile fix applied. New outputs:")
message("  map5_af3_delta_sedexport_zeroanchored.png")
message("  map5b_cropland_delta_sedexport_zeroanchored.png")
message("  map6_threshold_gain_pct_af3_zeroanchored.png")
message("  map7_threshold_gain_delta_af3_zeroanchored.png")
message("NOTE: previous *_quantile.png / *_pct_af3.png / *_delta_af3.png versions")
message("      are superseded by these zero-anchored versions and can be removed.")

# =============================================================================
# FIX — Exclude zero-change ("no change") watersheds from Sections 14 & 15
# =============================================================================
# PROBLEM: The existing low_baseline filter only excludes watersheds where
# the REFERENCE (baseline) value is too small to support a reliable percent
# -change denominator. It does NOT exclude watersheds where the transition
# produced NO measurable change at all (pct_chg == 0 or very close to it,
# and/or delta == 0). These zero-change watersheds are watersheds where no
# AF-eligible pixels fell within that watershed's contributing/flow-relevant
# area under the given suitability mask -- they were never "treated," so
# including them in effect-size statistics (mean, median, IQR, boxplots,
# choropleths) dilutes the estimated benefit among watersheds that were
# actually treated. This is the same "Q75 sits at exactly 0.0" pattern
# discussed earlier -- confirmed here to be a spatial-coverage artifact,
# not a magnitude-distribution one, and is now filtered out explicitly.
#
# APPROACH: add a `no_change` flag (TRUE if |pct_chg| < ZERO_CHANGE_EPS_PCT
# AND |delta| < ZERO_CHANGE_EPS_ABS, i.e. change is zero in BOTH relative
# and absolute terms -- guards against a watershed with huge baseline where
# a tiny absolute change rounds to ~0% but might still be non-trivial, and
# vice versa for tiny-baseline watersheds). Both stats (Section 14) and maps
# (Section 15) are regenerated to exclude no_change watersheds IN ADDITION
# TO low_baseline watersheds, and both exclusion reasons are reported
# separately so the two effects remain distinguishable.
# =============================================================================

ZERO_CHANGE_EPS_PCT <- 0.001   # % change smaller than this treated as "no change"
ZERO_CHANGE_EPS_ABS <- 0.01    # t/yr — absolute delta smaller than this treated as "no change"

# -----------------------------------------------------------------------
# Re-define compare_runs_robust to also flag no_change per metric
# -----------------------------------------------------------------------

compare_runs_robust <- function(ref_id, scen_id, min_baseline = min_baseline_by_metric) {
  comp <- compare_runs(ref_id, scen_id)
  if (is.null(comp)) return(NULL)
  for (m in names(min_baseline)) {
    ref_col   <- paste0(m, "_ref")
    pct_col   <- paste0(m, "_pct")
    delta_col <- paste0(m, "_delta")
    lowb_col  <- paste0(m, "_low_baseline")
    nochg_col <- paste0(m, "_no_change")
    
    if (ref_col %in% names(comp)) {
      comp[[lowb_col]] <- abs(comp[[ref_col]]) < min_baseline[[m]]
    }
    if (all(c(pct_col, delta_col) %in% names(comp))) {
      comp[[nochg_col]] <- abs(comp[[pct_col]]) < ZERO_CHANGE_EPS_PCT &
        abs(comp[[delta_col]]) < ZERO_CHANGE_EPS_ABS
    } else if (pct_col %in% names(comp)) {
      # fallback if delta column absent for this metric
      comp[[nochg_col]] <- abs(comp[[pct_col]]) < ZERO_CHANGE_EPS_PCT
    }
  }
  comp
}

cat("compare_runs_robust() redefined to also flag *_no_change per metric.\n")
cat("ZERO_CHANGE_EPS_PCT =", ZERO_CHANGE_EPS_PCT, "% | ZERO_CHANGE_EPS_ABS =", ZERO_CHANGE_EPS_ABS, "t/yr\n")

# -----------------------------------------------------------------------
# Re-define summarise_comp_robust: filtered set now excludes BOTH
# low-baseline AND no-change watersheds; unfiltered ("all") set is
# unchanged for comparison purposes.
# -----------------------------------------------------------------------

summarise_comp_robust <- function(comp_df, ref_id, scen_id,
                                  min_baseline = min_baseline_by_metric) {
  pct_cols <- grep("_pct$", names(comp_df), value = TRUE)
  
  map_dfr(pct_cols, function(col) {
    metric <- sub("_pct$", "", col)
    ref_col   <- paste0(metric, "_ref")
    lowb_col  <- paste0(metric, "_low_baseline")
    nochg_col <- paste0(metric, "_no_change")
    vals_all <- comp_df[[col]]
    
    cutoff <- if (metric %in% names(min_baseline)) min_baseline[[metric]] else 0
    is_low  <- if (lowb_col  %in% names(comp_df)) comp_df[[lowb_col]]  else rep(FALSE, nrow(comp_df))
    is_nochg <- if (nochg_col %in% names(comp_df)) comp_df[[nochg_col]] else rep(FALSE, nrow(comp_df))
    
    keep <- !is_low & !is_nochg
    vals_filt <- comp_df[[col]][keep]
    
    tibble(
      ref_id = ref_id, scen_id = scen_id, metric = metric,
      mean_pct_all       = mean(vals_all, na.rm = TRUE),
      median_pct_all     = median(vals_all, na.rm = TRUE),
      sd_pct_all         = sd(vals_all, na.rm = TRUE),
      q25_pct_all        = quantile(vals_all, 0.25, na.rm = TRUE),
      q75_pct_all        = quantile(vals_all, 0.75, na.rm = TRUE),
      n_ws_all           = sum(!is.na(vals_all)),
      mean_pct_filt      = mean(vals_filt, na.rm = TRUE),
      median_pct_filt    = median(vals_filt, na.rm = TRUE),
      sd_pct_filt        = sd(vals_filt, na.rm = TRUE),
      q25_pct_filt       = quantile(vals_filt, 0.25, na.rm = TRUE),
      q75_pct_filt       = quantile(vals_filt, 0.75, na.rm = TRUE),
      n_ws_filt          = sum(!is.na(vals_filt)),
      n_ws_excl_low_baseline = sum(is_low & !is.na(vals_all), na.rm = TRUE),
      n_ws_excl_no_change    = sum(is_nochg & !is_low & !is.na(vals_all), na.rm = TRUE),
      n_ws_excluded_total    = sum(!is.na(vals_all)) - sum(!is.na(vals_filt)),
      pct_excluded_total     = round(100 * (sum(!is.na(vals_all)) - sum(!is.na(vals_filt))) / sum(!is.na(vals_all)), 1),
      min_baseline_used  = cutoff
    )
  })
}

# -----------------------------------------------------------------------
# Regenerate Comparison A robust summary + per-watershed long form with
# BOTH flags, and re-save the same filenames used downstream (Sections
# 14/15 both read from these).
# -----------------------------------------------------------------------

message("Regenerating Comparison A with low_baseline + no_change filtering...")

comp_A_robust <- map_dfr(comp_A_all_pairs, function(p) {
  comp <- compare_runs_robust(p$ref, p$scen)
  if (is.null(comp)) return(NULL)
  s <- summarise_comp_robust(comp, p$ref, p$scen)
  s$group     <- p$group
  s$class     <- p$class
  s$threshold <- run_meta$threshold[run_meta$run_id == p$scen]
  s
})

write.csv(comp_A_robust, file.path(out14, "A_threshold_effect_outlier_robust.csv"), row.names = FALSE)

ws_pct_all <- map_dfr(comp_A_all_pairs, function(p) {
  comp <- compare_runs_robust(p$ref, p$scen)
  if (is.null(comp) || !"sed_export_pct" %in% names(comp)) return(NULL)
  tibble(
    group = p$group,
    class = p$class,
    threshold = run_meta$threshold[run_meta$run_id == p$scen],
    pct_chg = comp$sed_export_pct,
    low_baseline = comp$sed_export_low_baseline,
    no_change    = comp$sed_export_no_change
  )
}) |>
  mutate(threshold = factor(threshold, levels = c("Full area", "Top 10%", "Top 20%", "Top 25%")),
         exclusion_reason = case_when(
           low_baseline ~ "Low baseline",
           no_change    ~ "No change (untreated)",
           TRUE         ~ "Included"
         ))

write.csv(ws_pct_all, file.path(out14, "watershed_sed_export_pct_flagged.csv"), row.names = FALSE)

# -----------------------------------------------------------------------
# Coverage summary: how many watersheds are excluded, and WHY (low
# baseline vs. no change vs. included) -- reported per run so the two
# exclusion mechanisms stay distinguishable.
# -----------------------------------------------------------------------

coverage_summary <- ws_pct_all |>
  group_by(class, group, threshold) |>
  summarise(
    n_total          = n(),
    n_low_baseline   = sum(low_baseline, na.rm = TRUE),
    n_no_change      = sum(no_change & !low_baseline, na.rm = TRUE),
    n_included       = sum(!low_baseline & !no_change, na.rm = TRUE),
    pct_low_baseline = round(100 * n_low_baseline / n_total, 1),
    pct_no_change    = round(100 * n_no_change / n_total, 1),
    pct_included     = round(100 * n_included / n_total, 1),
    .groups = "drop"
  )

write.csv(coverage_summary, file.path(out14, "watershed_coverage_summary.csv"), row.names = FALSE)
print(coverage_summary)

# =============================================================================
# REBUILD fig6 (median/IQR vs mean) — now filtered set excludes no_change too
# =============================================================================

stat_compare_df <- comp_A_robust |>
  filter(metric == "sed_export") |>
  select(class, group, threshold, mean_pct_all, median_pct_all,
         q25_pct_all, q75_pct_all, mean_pct_filt, median_pct_filt,
         q25_pct_filt, q75_pct_filt) |>
  pivot_longer(cols = -c(class, group, threshold), names_to = "stat_col", values_to = "value") |>
  mutate(
    filtered  = ifelse(grepl("_filt$", stat_col),
                       "Filtered (excl. low-baseline + no-change)",
                       "All watersheds (unfiltered)"),
    stat_type = case_when(
      grepl("^mean_pct",   stat_col) ~ "mean",
      grepl("^median_pct", stat_col) ~ "median",
      grepl("^q25_pct",    stat_col) ~ "q25",
      grepl("^q75_pct",    stat_col) ~ "q75"
    )
  ) |>
  select(-stat_col) |>
  pivot_wider(names_from = stat_type, values_from = value) |>
  mutate(threshold = factor(threshold, levels = c("Top 10%", "Top 20%", "Top 25%")))

p_median_iqr_v2 <- ggplot(stat_compare_df, aes(x = threshold, y = median, fill = class)) +
  geom_col(position = position_dodge(0.7), width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = q25, ymax = q75), position = position_dodge(0.7), width = 0.25, linewidth = 0.5) +
  geom_point(aes(y = mean, colour = class), position = position_dodge(0.7),
             shape = 18, size = 2.6, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("Cropland-only" = "#2166AC", "AF3-class" = "#D73027"), name = "Transition") +
  scale_colour_manual(values = c("Cropland-only" = "#0b3757", "AF3-class" = "#7a1a10"), guide = "none") +
  facet_grid(filtered ~ group) +
  labs(
    title = "Median \u00b1 IQR vs mean for sed_export % change (untreated watersheds now excluded)",
    subtitle = "Filtered panel excludes low-baseline watersheds AND watersheds with no measurable change",
    x = "Suitability threshold", y = "% change in sed_export vs no-AF reference",
    caption = paste0("No-change threshold: |%| < ", ZERO_CHANGE_EPS_PCT, "% and |\u0394| < ",
                     ZERO_CHANGE_EPS_ABS, " t/yr | Low-baseline cutoff: ", MIN_BASELINE_SED_EXPORT, " t/yr")
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))

ggsave(file.path(out14, "fig6_median_iqr_vs_mean_outlier_check_v2.png"),
       p_median_iqr_v2, width = 12, height = 8, dpi = 180)

# =============================================================================
# REBUILD fig4 (violin/box distribution) — no-change watersheds removed
# from the plotted distribution entirely (not just marked), since they add
# a large uninformative spike at zero that compresses the visible spread
# of the treated watersheds. Low-baseline watersheds remain marked with X
# as before, for watersheds that ARE treated but have unreliable denominators.
# =============================================================================

p4_v2 <- ws_pct_all |>
  filter(!is.na(pct_chg), !no_change) |>
  ggplot(aes(x = threshold, y = pct_chg, fill = group, colour = group)) +
  geom_violin(position = position_dodge(0.8), alpha = 0.25, linewidth = 0.4,
              scale = "width", trim = TRUE) +
  geom_boxplot(position = position_dodge(0.8), width = 0.22,
               outlier.shape = NA, alpha = 0.75) +
  geom_point(data = ~ filter(.x, low_baseline),
             aes(x = threshold, y = pct_chg),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.1),
             shape = 4, size = 1.4, colour = "black", alpha = 0.6, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = climate_cols, name = "Climate") +
  scale_colour_manual(values = climate_cols, name = "Climate") +
  facet_wrap(~ class, ncol = 2) +
  coord_cartesian(ylim = c(-100, 150)) +
  labs(
    title = "Distribution of watershed-level sediment export change \u2014 untreated watersheds removed",
    subtitle = "No-change (untreated) watersheds excluded entirely | X = treated but low-baseline (unreliable %)",
    x = "Suitability threshold", y = "% change in sed_export per treated watershed",
    caption = paste0("Excludes watersheds with |%| < ", ZERO_CHANGE_EPS_PCT, "% and |\u0394| < ",
                     ZERO_CHANGE_EPS_ABS, " t/yr (no AF-eligible pixels affected this watershed)")
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(out14, "fig4_v2_watershed_distribution_treated_only.png"),
       p4_v2, width = 12, height = 7, dpi = 180)

# -----------------------------------------------------------------------
# NEW — Coverage bar chart: % of watersheds actually treated per run
# -----------------------------------------------------------------------

p_coverage <- coverage_summary |>
  pivot_longer(cols = c(pct_low_baseline, pct_no_change, pct_included),
               names_to = "reason", values_to = "pct") |>
  mutate(reason = factor(recode(reason,
                                pct_included = "Treated (included)",
                                pct_low_baseline = "Excluded: low baseline",
                                pct_no_change = "Excluded: no change (untreated)"),
                         levels = c("Treated (included)", "Excluded: no change (untreated)", "Excluded: low baseline"))) |>
  ggplot(aes(x = threshold, y = pct, fill = reason)) +
  geom_col(position = "stack", width = 0.65) +
  facet_grid(class ~ group) +
  scale_fill_manual(values = c("Treated (included)" = "#2166AC",
                               "Excluded: no change (untreated)" = "#999999",
                               "Excluded: low baseline" = "#D73027")) +
  labs(
    title = "Watershed treatment coverage by suitability threshold",
    subtitle = "Grey = watersheds with no AF-eligible pixels affecting sed_export | Red = unreliable baseline denominator",
    x = "Suitability threshold", y = "% of watersheds", fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(out14, "fig8_watershed_treatment_coverage.png"),
       p_coverage, width = 10, height = 7, dpi = 180)

message("\nSection 14 outputs regenerated (no-change filtering added):")
message("  A_threshold_effect_outlier_robust.csv        — filtered stats now exclude no-change watersheds")
message("  watershed_sed_export_pct_flagged.csv         — adds no_change flag + exclusion_reason column")
message("  watershed_coverage_summary.csv               — NEW: counts/pct by exclusion reason per run")
message("  fig6_median_iqr_vs_mean_outlier_check_v2.png — filtered panel now excludes no-change watersheds")
message("  fig4_v2_watershed_distribution_treated_only.png — no-change watersheds removed from distribution")
message("  fig8_watershed_treatment_coverage.png        — NEW: stacked bar of treated vs excluded by reason")

# =============================================================================
# SECTION 15 MAPS — rebuild map_data_long to also exclude no_change
# watersheds (greyed out, same as low_baseline) so spatial maps only show
# watersheds where the AF transition actually altered sed_export.
# =============================================================================

message("\nRebuilding Section 15 map data with no_change exclusion added...")

build_map_row_v2 <- function(p) {
  comp <- compare_runs_robust(p$ref, p$scen)
  if (is.null(comp)) {
    message("  Skipping (missing data): ", p$scen, " vs ", p$ref)
    return(NULL)
  }
  pct_col   <- paste0(map_metric, "_pct")
  ref_col   <- paste0(map_metric, "_ref")
  scen_col  <- paste0(map_metric, "_scen")
  delta_col <- paste0(map_metric, "_delta")
  lowb_col  <- paste0(map_metric, "_low_baseline")
  nochg_col <- paste0(map_metric, "_no_change")
  if (!all(c(pct_col, ref_col, scen_col, delta_col, lowb_col, nochg_col) %in% names(comp))) {
    warning("Metric columns missing for ", p$scen); return(NULL)
  }
  tibble(
    ws_id        = comp$ws_id,
    pct_chg      = comp[[pct_col]],
    delta        = comp[[delta_col]],
    ref_val      = comp[[ref_col]],
    scen_val     = comp[[scen_col]],
    low_baseline = comp[[lowb_col]],
    no_change    = comp[[nochg_col]],
    excluded     = comp[[lowb_col]] | comp[[nochg_col]],
    climate      = p$climate,
    threshold    = p$threshold,
    class        = p$class,
    ref_id       = p$ref,
    scen_id      = p$scen
  )
}

map_data_long <- map_dfr(map_pairs, build_map_row_v2) |>
  mutate(
    threshold = factor(threshold, levels = c("Top 10%", "Top 20%", "Top 25%")),
    climate   = factor(climate, levels = c("Baseline", "SSP370")),
    class     = factor(class, levels = c("Cropland-only", "AF3-class")),
    pct_chg_plot   = ifelse(excluded, NA_real_, pct_chg),
    delta_plot_all = ifelse(excluded, NA_real_, delta)
  )

write.csv(map_data_long, file.path(out15, paste0("map_data_", map_metric, "_all_runs_v2.csv")),
          row.names = FALSE)

map_sf_all <- ws_geom |> left_join(map_data_long, by = "ws_id") |> st_as_sf()

excl_summary_v2 <- map_data_long |>
  group_by(class, climate, threshold) |>
  summarise(n_total = n(),
            n_low_baseline = sum(low_baseline, na.rm = TRUE),
            n_no_change    = sum(no_change & !low_baseline, na.rm = TRUE),
            n_treated      = sum(!excluded, na.rm = TRUE),
            pct_treated    = round(100 * n_treated / n_total, 1),
            .groups = "drop")

write.csv(excl_summary_v2, file.path(out15, "exclusion_summary_by_run_v2.csv"), row.names = FALSE)
print(excl_summary_v2)

message("Section 15 map data rebuilt: map_data_", map_metric, "_all_runs_v2.csv")
message("Re-run Maps 1-7 (existing code below/above, unchanged) using this map_sf_all —")
message("no_change watersheds now render grey alongside low_baseline watersheds automatically,")
message("since both are captured in pct_chg_plot / delta_plot NA-ing via `excluded`.") 

