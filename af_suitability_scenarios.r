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
# AF SUITABILITY + SCENARIO C-FACTOR INTEGRATION
# Builds on: af_suitability_merged-5.r
# This script ADDS three new sections (16, 17, 18) and STREAMLINES
# Sections 6-11 of the prior script into a single parameterized loop.
# Sections 1-5, 12-15 (diagnostics, masks, comparison/robust stats, maps)
# are UNCHANGED and should be run from af_suitability_merged-5.r first,
# or sourced before this script if you want the full pipeline in one run.
# =============================================================================
#
# SECTIONS IN THIS FILE
# ----------------------
#   SECTION 6-11 (STREAMLINED) - generalized threshold/LULC/AF-area summary
#     loop, replacing the scenario-specific duplicated blocks in the merged
#     script. Produces the same CSV/plot outputs as before but driven by one
#     function so future scenario additions do not require new copy-pasted
#     blocks.
#   SECTION 16 - SCENARIO C-FACTOR DEFINITIONS (conservative/central/optimistic)
#     Defines baseline and AF usle_c values for cropland, grassland, shrubland
#     under each of the 3 scenarios, taken directly from the finalized
#     C-factor summary tab of ES-Quantification-Methods-and-Data-Requirements.
#   SECTION 17 - INVEST-READY LAYERS FOR 3 SCENARIOS x 3 LAND USES x 3 SUIT %
#     For each scenario (conservative/central/optimistic), builds:
#       - suitability masks at 10/20/25% (reused from Section 10 masks)
#       - LULC raster with cropland/grassland/shrubland AF pixels recoded
#       - biophysical_table.csv with scenario-specific usle_c values
#   SECTION 18 - DEGRADED LANDS ANALYSIS
#     Uses band 10 of the degradation-trend raster (-1 = degraded, 0 = neutral,
#     1 = improving) recoded to a binary degraded mask (0 = degraded,
#     1 = neutral/improving). Combines this with the degraded-land baseline
#     and AF C-factor values (conservative/central/optimistic) to build
#     InVEST-ready LULC rasters, biophysical tables, and suitability-masked
#     layers for degraded-land runs.
#
# =============================================================================

library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(purrr)
library(readr)

# =============================================================================
# ASSUMED OBJECTS FROM af_suitability_merged-5.r (Sections 1-5)
# If running this file standalone, uncomment and set these paths.
# =============================================================================
setwd("C:/CI_Alves/analysis/ssa_invest_runs/outputs")
suit_paths <- list(
  Baseline = "C:/CI_Alves/suitability/rasters-20260706T203255Z-3-001/rasters/suitability_agroforestry.tif",
   SSP126   = "C:/CI_Alves/suitability/rasters-20260706T203255Z-3-001/rasters/future_suitability_agroforestry_SSP126.tif",
   SSP370   = "C:/CI_Alves/suitability/rasters-20260706T203255Z-3-001/rasters/future_suitability_agroforestry_SSP370.tif",
   SSP585   = "C:/CI_Alves/suitability/rasters-20260706T203255Z-3-001/rasters/future_suitability_agroforestry_SSP585.tif"
 )
 lulc_path <- "C:/CI_Alves/analysis/ssa_invest_runs/Lulc/lulc_suitability/Landcover_harmonized_1km.tif"
 out_dir <- "outputs/baseline_vs_ssp_suitability"
 lulc_r <- rast(lulc_path)
# px_km2_use <- ... (from Section 3 diagnostics)
#
# load_and_normalise_suit() must already be defined (Section 2 helper).

stopifnot(exists("suit_paths"), exists("lulc_path"), exists("out_dir"),
          exists("lulc_r"), exists("px_km2_use"),
          exists("load_and_normalise_suit"))

# LULC codes used throughout (unchanged from merged script)
CROPLAND_CODE   <- 40L
GRASSLAND_CODE  <- 30L
SHRUBLAND_CODE  <- 20L

invest_suit_pcts <- c(10, 20, 25)   # suitability thresholds used everywhere below

# =============================================================================
# SECTION 6-11 (STREAMLINED) — GENERALIZED THRESHOLD / LULC / AF-AREA SUMMARY
# =============================================================================
# WHY THIS REPLACES THE OLD SECTIONS 6-11
# ----------------------------------------
# The original merged script computed continental LULC totals, per-scenario
# threshold summaries, cropland-only masks, and AF3-class (crop+grass+shrub)
# masks in four separate, heavily duplicated code blocks (Sections 6, 10, 11,
# 12). Each block re-loaded and re-rescaled the same suitability rasters.
# This section collapses that into ONE reusable function, run once per mask
# "flavor" (full AF3-class, cropland-only, grassland-only, shrubland-only),
# which is also what Section 17 needs for the scenario-specific runs below.
# All outputs (CSVs, mask rasters) are written to the SAME paths as before,
# so nothing downstream (Sections 12-15) breaks.
# =============================================================================

invest_thresholds <- c(10, 20, 25, 30)
invest_mask_dir    <- file.path(out_dir, "invest_suitability_masks")
af3_mask_dir       <- file.path(out_dir, "invest_suitability_masks_af_combined")

for (d in c(invest_mask_dir, af3_mask_dir)) {
  for (pct in invest_thresholds) {
    dir.create(file.path(d, paste0("top", pct)), recursive = TRUE, showWarnings = FALSE)
  }
}

# Build one binary "class of interest" raster: 1 = pixel matches lulc_codes, NA otherwise
build_class_mask <- function(lulc_r, lulc_codes) {
  ifel(lulc_r %in% lulc_codes, 1L, NA_integer_)
}

# Core reusable function: for one scenario (Baseline/SSP126/...), write the
# full suitability mask + a class-specific mask (cropland only, af3-class,
# grassland only, shrubland only - whatever lulc_codes specifies) at every
# threshold in `thresholds`. Returns a manifest data.frame of what was written.
export_suitability_masks <- function(scen, suit_paths, lulc_r, px_km2_use,
                                      thresholds = invest_thresholds,
                                      class_lulc_codes = NULL,
                                      class_label = "full",
                                      out_root = invest_mask_dir) {

  suit_r    <- load_and_normalise_suit(suit_paths[[scen]], scen)
  suit_vals <- values(suit_r, na.rm = TRUE)
  cutoffs   <- quantile(suit_vals, probs = 1 - thresholds / 100, na.rm = TRUE)

  class_r <- if (!is.null(class_lulc_codes)) build_class_mask(lulc_r, class_lulc_codes) else NULL

  manifest <- list()

  for (i in seq_along(thresholds)) {
    pct     <- thresholds[i]
    cut_val <- unname(cutoffs[i])
    suit_mask <- ifel(suit_r >= cut_val, 1L, NA_integer_)

    if (!is.null(class_r)) {
      class_aligned <- if (!compareGeom(class_r, suit_mask, stopOnError = FALSE)) {
        resample(class_r, suit_mask, method = "near")
      } else class_r
      out_mask <- ifel(suit_mask == 1L & class_aligned == 1L, 1L, NA_integer_)
      out_name <- paste0("suitability_mask_", class_label, "_", scen, "_top", pct, ".tif")
    } else {
      out_mask <- suit_mask
      out_name <- paste0("suitability_mask_", scen, "_top", pct, ".tif")
    }

    out_path <- file.path(out_root, paste0("top", pct), out_name)
    writeRaster(out_mask, out_path, datatype = "INT1U", NAflag = 255, overwrite = TRUE)

    n_px <- sum(values(out_mask, na.rm = TRUE) == 1, na.rm = TRUE)
    manifest[[paste(scen, class_label, pct)]] <- data.frame(
      scenario = scen, class_label = class_label, top_pct = pct,
      suit_cutoff = round(cut_val, 4), n_pixels = n_px,
      area_km2 = round(n_px * px_km2_use), out_path = out_path
    )
    rm(suit_mask, out_mask); gc()
  }
  rm(suit_r, suit_vals); gc()
  bind_rows(manifest)
}

# Run for every scenario x mask flavor (full suitability, cropland-only,
# AF3-class combined). This reproduces the outputs of old Sections 10-12.
mask_flavors <- list(
  full      = list(codes = NULL,                                       label = NULL),
  cropland  = list(codes = CROPLAND_CODE,                               label = "cropland"),
  af3class  = list(codes = c(CROPLAND_CODE, GRASSLAND_CODE, SHRUBLAND_CODE), label = "af3class")
)

streamlined_manifest <- map_dfr(names(suit_paths), function(scen) {
  map_dfr(names(mask_flavors), function(flavor) {
    fl <- mask_flavors[[flavor]]
    lbl <- if (is.null(fl$label)) "full" else fl$label
    out_root <- if (flavor == "af3class") af3_mask_dir else invest_mask_dir
    export_suitability_masks(scen, suit_paths, lulc_r, px_km2_use,
                              thresholds = invest_thresholds,
                              class_lulc_codes = fl$codes,
                              class_label = lbl,
                              out_root = out_root)
  })
})

write_csv(streamlined_manifest, file.path(out_dir, "CHECK_streamlined_mask_manifest.csv"))
cat("\n[Sections 6-11 streamlined] Mask manifest written:",
    file.path(out_dir, "CHECK_streamlined_mask_manifest.csv"), "\n")
cat("Rows:", nrow(streamlined_manifest), "(scenarios x flavors x thresholds)\n")

# =============================================================================
# SECTION 16 — FINAL NON-DEGRADED C-FACTOR DEFINITIONS (FIXED BASELINE, P = 1)
# SECTION 17 — INVEST-READY NON-DEGRADED SUITABILITY SCENARIOS (P = 1)
# =============================================================================
# Source values for C-factors come from the workbook-derived final values used
# in the prior rewrite, with the additional standardization that usle_p = 1.00
# for BOTH baseline and AF rows to isolate C-factor effects and maintain
# symmetry across comparisons.
#
# Per current instruction, Section 20 is NOT rewritten here.
# =============================================================================

library(terra)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(tibble)

stopifnot(exists("suit_paths"), exists("lulc_path"), exists("out_dir"),
          exists("lulc_r"), exists("px_km2_use"),
          exists("load_and_normalise_suit"), exists("af3_mask_dir"))

CROPLAND_CODE  <- 40L
GRASSLAND_CODE <- 30L
SHRUBLAND_CODE <- 20L
AF_CROP_CODE   <- 45L
AF_GRASS_CODE  <- 35L
AF_SHRUB_CODE  <- 25L
DEFAULT_P      <- 1.00
invest_suit_pcts         <- c(10, 20, 25)
invest_scenarios_climate <- c("Baseline", "SSP370")

# =============================================================================
# SECTION 16 — SCENARIO C-FACTOR TABLE (FIXED BASELINE, ALL P = 1)
# =============================================================================

baseline_fixed_tbl <- tribble(
  ~land_use,   ~lucode_base,    ~c_base,
  "cropland",  CROPLAND_CODE,   0.2129,
  "grassland", GRASSLAND_CODE,  0.0698,
  "shrubland", SHRUBLAND_CODE,  0.0693
)

af_scenario_tbl <- tribble(
  ~scenario,      ~land_use,   ~lucode_af,    ~c_af,
  "conservative", "cropland",  AF_CROP_CODE,  0.1971,
  "central",      "cropland",  AF_CROP_CODE,  0.1816,
  "optimistic",   "cropland",  AF_CROP_CODE,  0.0999,
  "conservative", "grassland", AF_GRASS_CODE, 0.5000,
  "central",      "grassland", AF_GRASS_CODE, 0.0510,
  "optimistic",   "grassland", AF_GRASS_CODE, 0.0155,
  "conservative", "shrubland", AF_SHRUB_CODE, 0.3600,
  "central",      "shrubland", AF_SHRUB_CODE, 0.0655,
  "optimistic",   "shrubland", AF_SHRUB_CODE, 0.0135
)

scenario_c_factors <- af_scenario_tbl |>
  left_join(baseline_fixed_tbl, by = "land_use") |>
  mutate(p_base = DEFAULT_P, p_af = DEFAULT_P) |>
  select(scenario, land_use, lucode_base, c_base, p_base, lucode_af, c_af, p_af)

stopifnot(
  nrow(scenario_c_factors) == 9,
  !anyDuplicated(scenario_c_factors[, c("scenario", "land_use")]),
  all(scenario_c_factors$lucode_base %in% c(SHRUBLAND_CODE, GRASSLAND_CODE, CROPLAND_CODE)),
  all(scenario_c_factors$lucode_af   %in% c(AF_SHRUB_CODE, AF_GRASS_CODE, AF_CROP_CODE)),
  all(scenario_c_factors$c_base >= 0 & scenario_c_factors$c_base <= 1),
  all(scenario_c_factors$c_af   >= 0 & scenario_c_factors$c_af   <= 1),
  all(scenario_c_factors$p_base == 1),
  all(scenario_c_factors$p_af == 1)
)

baseline_check <- scenario_c_factors |>
  group_by(land_use) |>
  summarise(n_distinct_base = n_distinct(c_base), .groups = "drop")
stopifnot(all(baseline_check$n_distinct_base == 1))

write_csv(scenario_c_factors,
          file.path(out_dir, "scenario_c_factors_nondegraded_FINAL_P1.csv"))

scenario_c_factor_summary <- scenario_c_factors |>
  arrange(match(land_use, c("cropland", "grassland", "shrubland")),
          match(scenario, c("conservative", "central", "optimistic"))) |>
  mutate(delta_c = c_af - c_base,
         pct_change_vs_base = 100 * (c_af - c_base) / c_base)

write_csv(scenario_c_factor_summary,
          file.path(out_dir, "SUMMARY_scenario_c_factors_nondegraded_FINAL_P1.csv"))

cat("\n[Section 16] Fixed-baseline scenario C-factor table written (P = 1 for baseline and AF).\n")
print(scenario_c_factor_summary)

# =============================================================================
# SECTION 17 — INVEST-READY LAYERS FOR NON-DEGRADED SCENARIOS
# =============================================================================

invest_input_dir_scen <- file.path(out_dir, "invest_inputs_scenarios_nondegraded_FINAL_P1")
dir.create(invest_input_dir_scen, recursive = TRUE, showWarnings = FALSE)

lulc_base_r <- if (exists("lulc_base_r")) lulc_base_r else rast(lulc_path)

canonical_standard_biophys <- tribble(
  ~lucode, ~usle_c, ~usle_p,
  -128L, 0.0000000000, 0.00,
  0L, 0.0000000000, 1.00,
  10L, 0.0050000000, 1.00,
  20L, 0.0693000000, 1.00,
  30L, 0.0698000000, 1.00,
  40L, 0.2129000000, 1.00,
  50L, 0.5000000000, 1.00,
  60L, 1.0000000000, 1.00,
  70L, 0.0000000000, 1.00,
  80L, 0.0000000000, 1.00,
  90L, 0.1000000000, 1.00,
  95L, 0.0100000000, 1.00,
  100L, 0.1500000000, 1.00,
  127L, 0.1500000000, 1.00
)

canon_check <- canonical_standard_biophys |>
  filter(lucode %in% c(SHRUBLAND_CODE, GRASSLAND_CODE, CROPLAND_CODE)) |>
  arrange(lucode)
expect_check <- baseline_fixed_tbl |>
  arrange(lucode_base) |>
  pull(c_base)
stopifnot(isTRUE(all.equal(canon_check$usle_c, expect_check, tolerance = 1e-9)))
stopifnot(all(canonical_standard_biophys$usle_p[canonical_standard_biophys$lucode %in% c(20L,30L,40L)] == 1))

build_scenario_biophys <- function(canonical_tbl, scen_c_tbl) {
  stopifnot(nrow(scen_c_tbl) == 3)
  
  af_rows <- scen_c_tbl |>
    transmute(lucode = as.integer(lucode_af), usle_c = c_af, usle_p = DEFAULT_P)
  
  out <- bind_rows(canonical_tbl, af_rows) |>
    distinct(lucode, .keep_all = TRUE) |>
    arrange(lucode)
  
  stopifnot(!anyDuplicated(out$lucode))
  stopifnot(all(c(SHRUBLAND_CODE, GRASSLAND_CODE, CROPLAND_CODE,
                  AF_SHRUB_CODE, AF_GRASS_CODE, AF_CROP_CODE) %in% out$lucode))
  stopifnot(all(out$usle_c >= 0 & out$usle_c <= 1))
  stopifnot(all(out$usle_p[out$lucode %in% c(20L,30L,40L,25L,35L,45L)] == 1))
  
  base_rows_out <- out |> filter(lucode %in% c(SHRUBLAND_CODE, GRASSLAND_CODE, CROPLAND_CODE)) |> arrange(lucode)
  base_rows_in  <- canonical_tbl |> filter(lucode %in% c(SHRUBLAND_CODE, GRASSLAND_CODE, CROPLAND_CODE)) |> arrange(lucode)
  stopifnot(isTRUE(all.equal(base_rows_out$usle_c, base_rows_in$usle_c, tolerance = 0)))
  
  out
}

align_mask_to_lulc <- function(mask_r, lulc_r) {
  if (!compareGeom(mask_r, lulc_r, stopOnError = FALSE)) {
    mask_r <- resample(mask_r, lulc_r, method = "near")
  }
  vals <- unique(na.omit(values(mask_r)))
  stopifnot(all(vals %in% 1L))
  mask_r
}

run_log_scenarios <- list()
change_log <- list()

for (scenario_name in unique(scenario_c_factors$scenario)) {
  scen_c_tbl   <- scenario_c_factors |> filter(scenario == scenario_name)
  biophys_scen <- build_scenario_biophys(canonical_standard_biophys, scen_c_tbl)
  
  for (climate in invest_scenarios_climate) {
    if (!climate %in% names(suit_paths)) next
    
    for (pct in invest_suit_pcts) {
      run_id  <- paste0(scenario_name, "_", climate, "_top", pct)
      run_dir <- file.path(invest_input_dir_scen, run_id)
      dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
      
      mask_path <- file.path(af3_mask_dir, paste0("top", pct),
                             paste0("suitability_mask_af3class_", climate, "_top", pct, ".tif"))
      stopifnot(file.exists(mask_path))
      af3_mask_r <- rast(mask_path)
      af3_mask_aligned <- align_mask_to_lulc(af3_mask_r, lulc_base_r)
      
      lulc_mod <- lulc_base_r
      for (i in seq_len(nrow(scen_c_tbl))) {
        row <- scen_c_tbl[i, ]
        lulc_mod <- ifel(lulc_base_r == row$lucode_base & af3_mask_aligned == 1L,
                         row$lucode_af, lulc_mod)
      }
      
      changed_idx <- which(!is.na(values(af3_mask_aligned)))
      if (length(changed_idx) > 0) {
        orig_vals <- values(lulc_base_r)[changed_idx]
        new_vals  <- values(lulc_mod)[changed_idx]
        valid_pairs <- rbind(
          c(CROPLAND_CODE,  AF_CROP_CODE),
          c(GRASSLAND_CODE, AF_GRASS_CODE),
          c(SHRUBLAND_CODE, AF_SHRUB_CODE)
        )
        ok <- rep(FALSE, length(orig_vals))
        for (p in seq_len(nrow(valid_pairs))) {
          ok <- ok | (orig_vals == valid_pairs[p, 1] & new_vals == valid_pairs[p, 2])
        }
        ok <- ok | (orig_vals == new_vals)
        stopifnot(all(ok, na.rm = TRUE))
      }
      
      lulc_out <- file.path(run_dir, "lulc.tif")
      writeRaster(lulc_mod, lulc_out, datatype = "INT2S", NAflag = -9999L, overwrite = TRUE)
      
      biophys_out <- file.path(run_dir, "biophysical_table.csv")
      write_csv(biophys_scen, biophys_out)
      
      n_by_class <- sapply(scen_c_tbl$lucode_af, function(lc) sum(values(lulc_mod, na.rm = TRUE) == lc, na.rm = TRUE))
      names(n_by_class) <- scen_c_tbl$land_use
      
      run_log_scenarios[[run_id]] <- tibble(
        run_id = run_id, scenario = scenario_name, climate = climate, top_pct = pct,
        crop_area_km2  = round(n_by_class["cropland"]  * px_km2_use),
        grass_area_km2 = round(n_by_class["grassland"] * px_km2_use),
        shrub_area_km2 = round(n_by_class["shrubland"] * px_km2_use),
        c_base_crop  = scen_c_tbl$c_base[scen_c_tbl$land_use == "cropland"],
        c_af_crop    = scen_c_tbl$c_af[scen_c_tbl$land_use == "cropland"],
        c_base_grass = scen_c_tbl$c_base[scen_c_tbl$land_use == "grassland"],
        c_af_grass   = scen_c_tbl$c_af[scen_c_tbl$land_use == "grassland"],
        c_base_shrub = scen_c_tbl$c_base[scen_c_tbl$land_use == "shrubland"],
        c_af_shrub   = scen_c_tbl$c_af[scen_c_tbl$land_use == "shrubland"],
        p_base = 1,
        p_af = 1,
        lulc_path = lulc_out,
        biophys_path = biophys_out
      )
      
      change_log[[run_id]] <- tibble(
        run_id = run_id, climate = climate, top_pct = pct,
        land_use = scen_c_tbl$land_use,
        lucode_base = scen_c_tbl$lucode_base,
        lucode_af = scen_c_tbl$lucode_af,
        c_base = scen_c_tbl$c_base,
        c_af = scen_c_tbl$c_af,
        p_base = 1,
        p_af = 1,
        delta_c = scen_c_tbl$c_af - scen_c_tbl$c_base,
        n_transitioned_pixels = as.numeric(n_by_class[scen_c_tbl$land_use])
      )
      
      rm(af3_mask_r, af3_mask_aligned, lulc_mod); gc()
    }
  }
}

run_log_scenarios_df <- bind_rows(run_log_scenarios)
change_log_df <- bind_rows(change_log)

write_csv(run_log_scenarios_df,
          file.path(invest_input_dir_scen, "invest_run_manifest_scenarios_nondegraded_FINAL_P1.csv"))
write_csv(change_log_df,
          file.path(invest_input_dir_scen, "SUMMARY_scenario_cfactor_changes_by_run_P1.csv"))

stopifnot(all(file.exists(run_log_scenarios_df$lulc_path)))
stopifnot(all(file.exists(run_log_scenarios_df$biophys_path)))
stopifnot(nrow(run_log_scenarios_df) == length(unique(scenario_c_factors$scenario)) *
            length(intersect(invest_scenarios_climate, names(suit_paths))) *
            length(invest_suit_pcts))

cat("\n[Section 17] Non-degraded scenario InVEST-input layers written (P = 1 for baseline and AF).\n")
cat("Run manifest:", file.path(invest_input_dir_scen, "invest_run_manifest_scenarios_nondegraded_FINAL_P1.csv"), "\n")



# =============================================================================
# SECTION 17B — NO-AF REFERENCE RUNS (BASELINE LULC, NO TRANSITION, P = 1)
# =============================================================================
# Builds the reference "no-AF" InVEST input pair (lulc.tif + biophysical_table.csv)
# for each climate scenario (Baseline, SSP370). These runs use the UNCHANGED
# LULC raster (no suitability-based cropland/grassland/shrubland -> AF
# recoding) and a biophysical table that reflects the fixed, current baseline
# C-factors from Section 16 with usle_p = 1.00 for all rows (matching the
# treatment-run standardization from Section 16/17).
#
# Written into the SAME invest_input_dir_scen root used for the AF-transition
# scenario runs so that Section 20 can find them via nondeg_reference_by_climate.
#
# Output run_id naming matches what Section 20 expects:
#   "baseline_noAF_reference"
#   "ssp370_noAF_reference"
# =============================================================================

stopifnot(exists("out_dir"), exists("lulc_path"), exists("canonical_standard_biophys"),
          exists("invest_input_dir_scen"), exists("px_km2_use"))

noaf_reference_run_ids <- c(Baseline = "baseline_noAF_reference", SSP370 = "ssp370_noAF_reference")
noaf_reference_climates <- names(noaf_reference_run_ids)

lulc_base_r <- if (exists("lulc_base_r")) lulc_base_r else rast(lulc_path)

# ---- Biophysical table for the no-AF reference: identical structure to the
# canonical table used for scenario runs, but with NO AF lucodes (25/35/45)
# present at all — because no transition occurs, these codes never appear in
# the no-AF LULC raster. usle_p = 1.00 throughout, matching the P standardization. ----
noaf_biophys <- canonical_standard_biophys |>
  filter(!lucode %in% c(25L, 35L, 45L)) |>
  arrange(lucode)

stopifnot(
  all(noaf_biophys$usle_p == 1 | noaf_biophys$lucode == -128L),
  !anyDuplicated(noaf_biophys$lucode),
  all(c(20L, 30L, 40L) %in% noaf_biophys$lucode)
)

# Confirm the baseline C-factors embedded here match Section 16's fixed
# baseline values exactly (no drift between scenario runs and the reference).
noaf_check <- noaf_biophys |> filter(lucode %in% c(20L, 30L, 40L)) |> arrange(lucode)
expected_check <- baseline_fixed_tbl |> arrange(lucode_base) |> pull(c_base)
stopifnot(isTRUE(all.equal(noaf_check$usle_c, expected_check, tolerance = 1e-9)))

noaf_run_log <- list()

for (climate in noaf_reference_climates) {
  run_id  <- noaf_reference_run_ids[[climate]]
  run_dir <- file.path(invest_input_dir_scen, run_id)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  
  # LULC is UNCHANGED — no suitability mask applied, no recoding to AF lucodes.
  # We still write a fresh copy so every InVEST run_dir is self-contained and
  # the raster is guaranteed aligned/typed identically to the treatment runs.
  lulc_out <- file.path(run_dir, "lulc.tif")
  writeRaster(lulc_base_r, lulc_out, datatype = "INT2S", NAflag = -9999L, overwrite = TRUE)
  
  biophys_out <- file.path(run_dir, "biophysical_table.csv")
  write_csv(noaf_biophys, biophys_out)
  
  n_crop  <- sum(values(lulc_base_r, na.rm = TRUE) == CROPLAND_CODE, na.rm = TRUE)
  n_grass <- sum(values(lulc_base_r, na.rm = TRUE) == GRASSLAND_CODE, na.rm = TRUE)
  n_shrub <- sum(values(lulc_base_r, na.rm = TRUE) == SHRUBLAND_CODE, na.rm = TRUE)
  
  noaf_run_log[[run_id]] <- tibble(
    run_id = run_id, climate = climate,
    crop_area_km2  = round(n_crop  * px_km2_use),
    grass_area_km2 = round(n_grass * px_km2_use),
    shrub_area_km2 = round(n_shrub * px_km2_use),
    c_cropland  = baseline_fixed_tbl$c_base[baseline_fixed_tbl$land_use == "cropland"],
    c_grassland = baseline_fixed_tbl$c_base[baseline_fixed_tbl$land_use == "grassland"],
    c_shrubland = baseline_fixed_tbl$c_base[baseline_fixed_tbl$land_use == "shrubland"],
    p_all = 1,
    lulc_path = lulc_out,
    biophys_path = biophys_out
  )
  
  cat("\n[Section 17B]", run_id, "— no-AF reference written:\n")
  cat("  LULC:     ", lulc_out, "\n")
  cat("  Biophys:  ", biophys_out, "\n")
  cat("  Cropland area (km2): ", round(n_crop * px_km2_use), "\n")
  cat("  Grassland area (km2):", round(n_grass * px_km2_use), "\n")
  cat("  Shrubland area (km2):", round(n_shrub * px_km2_use), "\n")
}

noaf_run_log_df <- bind_rows(noaf_run_log)
write_csv(noaf_run_log_df, file.path(invest_input_dir_scen, "invest_run_manifest_noAF_reference_P1.csv"))

stopifnot(all(file.exists(noaf_run_log_df$lulc_path)))
stopifnot(all(file.exists(noaf_run_log_df$biophys_path)))
stopifnot(nrow(noaf_run_log_df) == length(noaf_reference_climates))

cat("\n[Section 17B COMPLETE] No-AF reference runs written for:", paste(noaf_reference_climates, collapse = ", "), "\n")
cat("Manifest:", file.path(invest_input_dir_scen, "invest_run_manifest_noAF_reference_P1.csv"), "\n")
cat("Run IDs match Section 20 nondeg_reference_by_climate expectations:\n")
print(noaf_reference_run_ids) 



# =============================================================================
# SECTION 18 — DEGRADED LANDS ANALYSIS
# =============================================================================
# INPUT: a multi-band degradation raster where BAND 10 encodes a trend
# classification: -1 = degrading, 0 = neutral/stable, 1 = improving.
# This is recoded to a BINARY degraded mask:
#     0 = degraded  (was -1)
#     1 = neutral or improving (was 0 or 1)
# This binary mask is then intersected with LULC (cropland/grassland/
# shrubland) and the top-X% suitability mask to identify degraded-land
# pixels eligible for AF-transition runs, using the DEGRADED baseline and
# AF C-factor values (conservative/central/optimistic) rather than the
# non-degraded values used in Section 17.
#
# DEGRADED C-FACTOR VALUES (from the finalized C-factor summary tab)
# ---------------------------------------------------------------
#   Cropland:  baseline 0.45 / 0.50 / 0.55 (cons/central/opt);  AF 0.375 / 0.31 / 0.25
#   Grassland: baseline 0.35 / 0.55 / 0.82;                     AF 0.30  / 0.318 / 0.17
#   Shrubland: baseline 0.15 / 0.20 / 0.49;                     AF 0.13  / 0.11  / 0.05
#
# New lucodes for degraded-land AF transitions are offset by +100 from the
# non-degraded AF lucodes used in Section 17, so degraded and non-degraded AF
# pixels can be distinguished in the same biophysical table if ever merged:
#   145 = degraded cropland  -> AF   (vs 45 non-degraded)
#   135 = degraded grassland -> AF   (vs 35 non-degraded)
#   125 = degraded shrubland -> AF   (vs 25 non-degraded)
# Degraded baseline (no-AF) lucodes are similarly offset by +100 relative to
# the standard LULC codes so a degraded cropland pixel that remains
# untransitioned can carry the WORSE degraded C-factor instead of the
# standard 0.34 cropland value:
#   140 = degraded cropland  (no AF)
#   130 = degraded grassland (no AF)
#   120 = degraded shrubland (no AF)
# =============================================================================

degraded_band_path <- "C:/CI_Alves/analysis/ssa_invest_runs/Degraded/band_10_sdg_degraded.tif"
DEGRADED_BAND_INDEX <- 10   # band 10 of the multi-band raster

stopifnot(file.exists(degraded_band_path))
degrad_stack <- rast(degraded_band_path)

degrad_trend_r <- degrad_stack
names(degrad_trend_r) <- "degradation_trend"

degrad_trend_r <- degrad_stack[[DEGRADED_BAND_INDEX]]

degraded_binary_r <- classify(
  degrad_trend_r,
  rcl = matrix(c(-1, 0,
                 0, 1,
                 1, 1), ncol = 2, byrow = TRUE)
)

degraded_only_r <- ifel(degraded_binary_r == 0, 1L, NA_integer_)

degraded_aligned <- if (!compareGeom(degraded_only_r, lulc_base_r, stopOnError = FALSE)) {
  resample(degraded_only_r, lulc_base_r, method = "near")
} else degraded_only_r

suit_aligned <- if (!compareGeom(suit_mask_r, lulc_base_r, stopOnError = FALSE)) {
  resample(suit_mask_r, lulc_base_r, method = "near")
} else suit_mask_r

eligible_r <- ifel(degraded_aligned == 1L & suit_aligned == 1L, 1L, NA_integer_)
residual_deg_r <- ifel(degraded_aligned == 1L & is.na(eligible_r), 1L, NA_integer_)




cat("\n========================================\n")
cat("Section 18 - Degraded Lands Analysis\n")
cat("========================================\n")
cat("Band", DEGRADED_BAND_INDEX, "unique values (raw):\n")
print(as.data.frame(freq(degrad_trend_r)))

# Recode: -1 (degraded) -> 0 ; 0 or 1 (neutral/improving) -> 1
degraded_binary_r <- classify(
  degrad_trend_r,
  rcl = matrix(c(-1, 0,
                  0, 1,
                  1, 1), ncol = 2, byrow = TRUE)
)
names(degraded_binary_r) <- "degraded_binary"

cat("\nRecoded binary degraded mask (0 = degraded, 1 = neutral/improving):\n")
print(as.data.frame(freq(degraded_binary_r)))

degraded_mask_out_dir <- file.path(out_dir, "degraded_lands_masks")
dir.create(degraded_mask_out_dir, recursive = TRUE, showWarnings = FALSE)
writeRaster(degraded_binary_r,
            file.path(degraded_mask_out_dir, "degraded_binary_mask.tif"),
            datatype = "INT1U", NAflag = 255, overwrite = TRUE)

# Isolate ONLY the degraded pixels (value == 0) as the analysis mask
degraded_only_r <- ifel(degraded_binary_r == 0, 1L, NA_integer_)
writeRaster(degraded_only_r,
            file.path(degraded_mask_out_dir, "degraded_only_mask.tif"),
            datatype = "INT1U", NAflag = 255, overwrite = TRUE)

n_degraded_px <- sum(values(degraded_only_r, na.rm = TRUE) == 1, na.rm = TRUE)
cat("\nTotal degraded pixels (band", DEGRADED_BAND_INDEX, "== -1):", n_degraded_px,
    "| area km2:", round(n_degraded_px * px_km2_use), "\n")

# ---- Degraded-land C-factor scenario table ---------------------------------

scenario_c_factors_degraded <- tribble(
  ~scenario,      ~land_use,   ~lucode_base_std, ~lucode_base_deg, ~c_base_deg, ~lucode_af_deg, ~c_af_deg,
  "conservative", "cropland",   CROPLAND_CODE,     140L,             0.45,       145L,           0.375,
  "central",      "cropland",   CROPLAND_CODE,     140L,             0.50,       145L,           0.31,
  "optimistic",   "cropland",   CROPLAND_CODE,     140L,             0.55,       145L,           0.25,

  "conservative", "grassland",  GRASSLAND_CODE,    130L,             0.35,       135L,           0.30,
  "central",      "grassland",  GRASSLAND_CODE,    130L,             0.55,       135L,           0.318,
  "optimistic",   "grassland",  GRASSLAND_CODE,    130L,             0.82,       135L,           0.17,

  "conservative", "shrubland",  SHRUBLAND_CODE,    120L,             0.15,       125L,           0.13,
  "central",      "shrubland",  SHRUBLAND_CODE,    120L,             0.20,       125L,           0.11,
  "optimistic",   "shrubland",  SHRUBLAND_CODE,    120L,             0.49,       125L,           0.05
)

write_csv(scenario_c_factors_degraded, file.path(out_dir, "scenario_c_factors_degraded.csv"))
cat("\nDegraded-land scenario C-factor table written:",
    file.path(out_dir, "scenario_c_factors_degraded.csv"), "\n")
print(scenario_c_factors_degraded)

# Builds a degraded-land biophysical table: base table PLUS new degraded
# baseline rows (140/130/120) at c_base_deg, PLUS new degraded AF rows
# (145/135/125) at c_af_deg. Original standard lucodes (40/30/20) are left
# untouched so non-degraded pixels keep their normal C-factor.
build_degraded_biophys <- function(base_tbl, deg_c_tbl) {
  extra_cols <- setdiff(names(base_tbl), c("lucode", "usle_c", "usle_p"))
  new_rows <- list()
  for (i in seq_len(nrow(deg_c_tbl))) {
    row <- deg_c_tbl[i, ]
    std_row <- base_tbl[base_tbl$lucode == row$lucode_base_std, ]

    base_deg_row <- tibble(lucode = row$lucode_base_deg, usle_c = row$c_base_deg, usle_p = DEFAULT_P)
    af_deg_row   <- tibble(lucode = row$lucode_af_deg,   usle_c = row$c_af_deg,   usle_p = DEFAULT_P)
    if (length(extra_cols) > 0) {
      for (col in extra_cols) {
        base_deg_row[[col]] <- std_row[[col]]
        af_deg_row[[col]]   <- std_row[[col]]
      }
    }
    new_rows[[length(new_rows) + 1]] <- base_deg_row[, names(base_tbl)]
    new_rows[[length(new_rows) + 1]] <- af_deg_row[, names(base_tbl)]
  }
  bind_rows(base_tbl, bind_rows(new_rows))
}

invest_input_dir_degraded <- file.path(out_dir, "invest_inputs_degraded")
dir.create(invest_input_dir_degraded, recursive = TRUE, showWarnings = FALSE)

run_log_degraded <- list()

for (scenario_name in unique(scenario_c_factors_degraded$scenario)) {

  deg_c_tbl <- scenario_c_factors_degraded |> filter(scenario == scenario_name)
  biophys_deg <- build_degraded_biophys(biophys_base, deg_c_tbl)

  for (climate in invest_scenarios_climate) {
    if (!climate %in% names(suit_paths)) next

    for (pct in invest_suit_pcts) {

      run_id  <- paste0("degraded_", scenario_name, "_", climate, "_top", pct)
      run_dir <- file.path(invest_input_dir_degraded, run_id)
      dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

      cat("\n----------------------------\n")
      cat("Degraded-land scenario run:", run_id, "\n")
      cat("----------------------------\n")

      # Reuse af3class suitability mask (Sections 6-11 streamlined)
      mask_path <- file.path(af3_mask_dir, paste0("top", pct),
                              paste0("suitability_mask_af3class_", climate, "_top", pct, ".tif"))
      stopifnot(file.exists(mask_path))
      suit_mask_r <- rast(mask_path)

      # Align degraded-only mask and suitability mask to the LULC grid
      degraded_aligned <- if (!compareGeom(degraded_only_r, lulc_base_r, stopOnError = FALSE)) {
        resample(degraded_only_r, lulc_base_r, method = "near")
      } else degraded_only_r

      suit_aligned <- if (!compareGeom(suit_mask_r, lulc_base_r, stopOnError = FALSE)) {
        resample(suit_mask_r, lulc_base_r, method = "near")
      } else suit_mask_r

      # Combined eligibility mask: degraded AND within top-X% suitability
      eligible_r <- ifel(degraded_aligned == 1L & suit_aligned == 1L, 1L, NA_integer_)

      lulc_mod_deg <- lulc_base_r
      for (i in seq_len(nrow(deg_c_tbl))) {
        row <- deg_c_tbl[i, ]
        # Step 1: recode ALL degraded pixels of this class to the degraded
        # BASELINE lucode (whether or not they get AF), so degraded soil-loss
        # rates are captured even where AF is not applied within this run.
        is_class_degraded <- (lulc_base_r == row$lucode_base_std) & (degraded_aligned == 1L)
        lulc_mod_deg <- ifel(is_class_degraded, row$lucode_base_deg, lulc_mod_deg)

        # Step 2: within degraded pixels that are ALSO in the top-X%
        # suitability mask, recode further to the degraded AF lucode.
        is_class_eligible <- (lulc_base_r == row$lucode_base_std) & (eligible_r == 1L)
        lulc_mod_deg <- ifel(is_class_eligible, row$lucode_af_deg, lulc_mod_deg)
      }

      lulc_out_deg <- file.path(run_dir, "lulc.tif")
      writeRaster(lulc_mod_deg, lulc_out_deg, datatype = "INT2S", NAflag = -9999L, overwrite = TRUE)

      biophys_out_deg <- file.path(run_dir, "biophysical_table.csv")
      write_csv(biophys_deg, biophys_out_deg)

      n_deg_base <- sapply(deg_c_tbl$lucode_base_deg, function(lc) {
        sum(values(lulc_mod_deg, na.rm = TRUE) == lc, na.rm = TRUE)
      })
      n_deg_af <- sapply(deg_c_tbl$lucode_af_deg, function(lc) {
        sum(values(lulc_mod_deg, na.rm = TRUE) == lc, na.rm = TRUE)
      })
      names(n_deg_base) <- deg_c_tbl$land_use
      names(n_deg_af)   <- deg_c_tbl$land_use

      cat(" Degraded (no-AF) pixels : crop =", n_deg_base["cropland"],
          "| grass =", n_deg_base["grassland"], "| shrub =", n_deg_base["shrubland"], "\n")
      cat(" Degraded -> AF pixels   : crop =", n_deg_af["cropland"],
          "| grass =", n_deg_af["grassland"], "| shrub =", n_deg_af["shrubland"], "\n")
      cat(" LULC written:", lulc_out_deg, "\n")
      cat(" Biophysical table written:", biophys_out_deg, "\n")

      run_log_degraded[[run_id]] <- data.frame(
        run_id = run_id, scenario = scenario_name, climate = climate, top_pct = pct,
        deg_base_crop_km2  = round(n_deg_base["cropland"]  * px_km2_use),
        deg_base_grass_km2 = round(n_deg_base["grassland"] * px_km2_use),
        deg_base_shrub_km2 = round(n_deg_base["shrubland"] * px_km2_use),
        deg_af_crop_km2  = round(n_deg_af["cropland"]  * px_km2_use),
        deg_af_grass_km2 = round(n_deg_af["grassland"] * px_km2_use),
        deg_af_shrub_km2 = round(n_deg_af["shrubland"] * px_km2_use),
        c_base_deg_crop  = deg_c_tbl$c_base_deg[deg_c_tbl$land_use == "cropland"],
        c_af_deg_crop    = deg_c_tbl$c_af_deg[deg_c_tbl$land_use == "cropland"],
        c_base_deg_grass = deg_c_tbl$c_base_deg[deg_c_tbl$land_use == "grassland"],
        c_af_deg_grass   = deg_c_tbl$c_af_deg[deg_c_tbl$land_use == "grassland"],
        c_base_deg_shrub = deg_c_tbl$c_base_deg[deg_c_tbl$land_use == "shrubland"],
        c_af_deg_shrub   = deg_c_tbl$c_af_deg[deg_c_tbl$land_use == "shrubland"],
        lulc_path = lulc_out_deg, biophys_path = biophys_out_deg
      )

      rm(degraded_aligned, suit_aligned, eligible_r, lulc_mod_deg, suit_mask_r); gc()
    }
  }
}

run_log_degraded_df <- bind_rows(run_log_degraded)
write_csv(run_log_degraded_df, file.path(invest_input_dir_degraded, "invest_run_manifest_degraded.csv"))

cat("\n========================================\n")
cat("DEGRADED-LAND RUN MANIFEST\n")
cat("========================================\n")
print(run_log_degraded_df[, c("run_id", "top_pct",
                              "deg_base_crop_km2", "deg_af_crop_km2",
                              "deg_base_grass_km2", "deg_af_grass_km2",
                              "deg_base_shrub_km2", "deg_af_shrub_km2")])
cat("\nManifest written to:",
    file.path(invest_input_dir_degraded, "invest_run_manifest_degraded.csv"), "\n")
cat("Total degraded-land runs:", nrow(run_log_degraded_df),
    "(3 scenarios x", length(invest_scenarios_climate), "climates x",
    length(invest_suit_pcts), "thresholds)\n")

# =============================================================================
# SECTION 19 — COMBINED RUN INDEX (ALL RUNS: NON-DEGRADED + DEGRADED)
# =============================================================================
# Convenience index tying together every InVEST-ready run directory produced
# by Section 17 (non-degraded, 3 scenarios) and Section 18 (degraded, 3
# scenarios), for use as the master lookup when actually invoking InVEST SDR
# in a batch loop, and later when re-running the Section 12-15 comparison /
# mapping code against these new run_ids.
# =============================================================================

all_runs_index <- bind_rows(
  run_log_scenarios_df |> mutate(land_condition = "non_degraded"),
  run_log_degraded_df  |> mutate(land_condition = "degraded")
)

write_csv(all_runs_index, file.path(out_dir, "MASTER_scenario_run_index.csv"))

cat("\n========================================\n")
cat("MASTER RUN INDEX WRITTEN\n")
cat("========================================\n")
cat("File:", file.path(out_dir, "MASTER_scenario_run_index.csv"), "\n")
cat("Total runs indexed:", nrow(all_runs_index), "\n")
cat("  Non-degraded:", sum(all_runs_index$land_condition == "non_degraded"), "\n")
cat("  Degraded    :", sum(all_runs_index$land_condition == "degraded"), "\n")
cat("\nEach row's lulc_path + biophys_path are ready to feed directly into\n")
cat("InVEST SDR model_run() calls (Python natcap.invest.sdr.sdr module) as\n")
cat("the lulc_path and biophysical_table_path arguments respectively.\n")

# =============================================================================
# SECTION 20 — NON-DEGRADED SCENARIO-BAND INVEST OUTPUT ANALYSIS
#   (rewritten to match Sections 16-17 fixed-baseline, P = 1 run structure,
#    and updated to consume the SECTION 17B no-AF reference runs)
# =============================================================================
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(sf)
library(scales)
library(readr)


project_root <- "C:/CI_Alves/analysis/ssa_invest_runs" 
out_dir <- file.path( project_root, "outputs", "baseline_vs_ssp_suitability" ) 
out_20 <- file.path( out_dir, "SECTION20_nondegraded_scenario_bands_P1" ) 
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE) 
dir.create(out_20, recursive = TRUE, showWarnings = FALSE) 
stopifnot( dir.exists(out_dir), dir.exists(out_20), file.access(out_20, mode = 2) == 0 ) 
cat("\n[Section 20 output directory]\n", normalizePath(out_20, winslash = "/"), "\n", sep = "") 


if (!exists("min_baseline_by_metric")) {
  min_baseline_by_metric <- c(sed_export = 1000, usle_tot = 1000, sed_dep = 1000,
                              avoid_eros = 1000, avoid_exp = 1000)
}
ZERO_CHANGE_EPS_PCT <- if (exists("ZERO_CHANGE_EPS_PCT")) ZERO_CHANGE_EPS_PCT else 0.5
ZERO_CHANGE_EPS_ABS <- if (exists("ZERO_CHANGE_EPS_ABS")) ZERO_CHANGE_EPS_ABS else 5
OUTLIER_IQR_MULT     <- if (exists("OUTLIER_IQR_MULT")) OUTLIER_IQR_MULT else 3
OUTLIER_DROP_MAX_SHARE_PCT <- if (exists("OUTLIER_DROP_MAX_SHARE_PCT")) OUTLIER_DROP_MAX_SHARE_PCT else 5

if (!exists("pct_change")) {
  pct_change <- function(scen, ref) ifelse(ref == 0, NA_real_, (scen - ref) / abs(ref) * 100)
}

nondeg_output_root <- if (exists("invest_input_dir_scen")) invest_input_dir_scen else
  "C:/CI_Alves/analysis/ssa_invest_runs/Outputs/outputs/baseline_vs_ssp_suitability/invest_inputs_scenarios_nondegraded_FINAL_P1"
stopifnot("nondeg_output_root does not exist — check the path." = dir.exists(nondeg_output_root))

out_20 <- file.path(out_dir, "SECTION20_nondegraded_scenario_bands_P1")
dir.create(out_20, recursive = TRUE, showWarnings = FALSE)

nondeg_scenario_bands <- c("conservative", "central", "optimistic")
nondeg_climates       <- c("Baseline", "SSP370")
nondeg_thresholds     <- c(10, 20, 25)

# The no-AF reference run per climate is produced by SECTION 17B (unchanged
# LULC, no suitability-based transition, fixed baseline C-factors, P = 1 for
# all rows) and written into the SAME invest_input_dir_scen root as the
# scenario-band runs. Run IDs below match Section 17B's noaf_reference_run_ids
# exactly — do not rename independently in either section.
nondeg_reference_by_climate <- if (exists("noaf_reference_run_ids")) {
  setNames(as.character(noaf_reference_run_ids), names(noaf_reference_run_ids))
} else {
  c("Baseline" = "baseline_noAF_reference", "SSP370" = "ssp370_noAF_reference")
}

# ---- Pre-flight check: confirm Section 17B reference INPUT runs exist on
# disk (lulc.tif + biophysical_table.csv) before any comparison is attempted.
# Catches the "forgot to run Section 17B first" failure mode early. ----
for (clim in names(nondeg_reference_by_climate)) {
  ref_run_id  <- nondeg_reference_by_climate[[clim]]
  ref_run_dir <- file.path(nondeg_output_root, ref_run_id)
  ref_lulc    <- file.path(ref_run_dir, "lulc.tif")
  ref_biophys <- file.path(ref_run_dir, "biophysical_table.csv")
  if (!file.exists(ref_lulc) || !file.exists(ref_biophys)) {
    stop("No-AF reference inputs missing for climate '", clim, "' at: ", ref_run_dir,
         ". Run SECTION 17B (no-AF reference run generator) before Section 20.")
  }
  cat("[Pre-flight] No-AF reference inputs found for", clim, "->", ref_run_id, "\n")
}
cat("[Pre-flight] NOTE: checks above confirm INPUT files exist. The watershed_results_sdr.shp\n")
cat("             OUTPUT files are only produced once InVEST SDR has actually been RUN on these\n")
cat("             reference run_dirs (in addition to the scenario-band run_dirs from Section 17).\n")

cat("\n=====================================================\n")
cat("SECTION 20 — Non-Degraded Scenario-Band Analysis (P = 1)\n")
cat("=====================================================\n")

# =============================================================================
# CHUNK 1 — helper functions
# =============================================================================

# NOTE: load_watershed_20() looks for watershed_results_sdr.shp under
# nondeg_output_root/<run_id>/. For the no-AF reference run_ids, these output
# shapefiles are produced by running InVEST SDR on the lulc.tif +
# biophysical_table.csv written by SECTION 17B into that same run_id folder.
load_watershed_20 <- function(run_id) {
  candidate_roots <- c(nondeg_output_root, if (exists("invest_output_root")) invest_output_root else NULL)
  for (root in candidate_roots) {
    shp_path <- file.path(root, run_id, "watershed_results_sdr.shp")
    if (file.exists(shp_path)) {
      sf_obj <- st_read(shp_path, quiet = TRUE)
      sf_obj <- st_drop_geometry(sf_obj)
      names(sf_obj) <- tolower(names(sf_obj))
      if (!"wsid" %in% names(sf_obj)) sf_obj$wsid <- seq_len(nrow(sf_obj))
      return(sf_obj)
    }
  }
  warning("Missing shapefile for run_id: ", run_id, " — checked: ",
          paste(file.path(candidate_roots, run_id), collapse = " | "))
  NULL
}

compare_runs_20 <- function(ref_id, scen_id, metrics = c("usle_tot", "sed_export", "sed_dep", "avoid_eros", "avoid_exp")) {
  ref_df  <- load_watershed_20(ref_id)
  scen_df <- load_watershed_20(scen_id)
  if (is.null(ref_df) || is.null(scen_df)) return(NULL)
  
  metrics <- intersect(metrics, intersect(names(ref_df), names(scen_df)))
  if (length(metrics) == 0) {
    warning("No matching metric columns between ", ref_id, " and ", scen_id)
    return(NULL)
  }
  
  ref_sub  <- ref_df[, c("wsid", metrics)]
  scen_sub <- scen_df[, c("wsid", metrics)]
  names(ref_sub)[-1]  <- paste0(metrics, "_ref")
  names(scen_sub)[-1] <- paste0(metrics, "_scen")
  
  merged <- merge(ref_sub, scen_sub, by = "wsid")
  for (m in metrics) {
    merged[[paste0(m, "_delta")]] <- merged[[paste0(m, "_scen")]] - merged[[paste0(m, "_ref")]]
    merged[[paste0(m, "_pct")]]   <- pct_change(merged[[paste0(m, "_scen")]], merged[[paste0(m, "_ref")]])
  }
  merged
}

# =============================================================================
# CHUNK 2 — VALIDATION / OUTLIER SCREEN
# =============================================================================

flag_outliers_20 <- function(df, group_cols = c("climate", "threshold", "band"),
                             value_col = "pct_chg",
                             min_baseline = 1000, ref_col = "ref_val",
                             delta_col = "delta", eps_abs = ZERO_CHANGE_EPS_ABS,
                             eps_pct = ZERO_CHANGE_EPS_PCT,
                             iqr_mult = OUTLIER_IQR_MULT) {
  
  df <- df |>
    mutate(
      flag_low_baseline = abs(.data[[ref_col]]) < min_baseline,
      flag_no_change    = abs(.data[[delta_col]]) < eps_abs &
        (is.na(.data[[value_col]]) | abs(.data[[value_col]]) < eps_pct)
    )
  
  df <- df |>
    group_by(across(all_of(group_cols))) |>
    mutate(
      .q1 = quantile(.data[[value_col]], 0.25, na.rm = TRUE),
      .q3 = quantile(.data[[value_col]], 0.75, na.rm = TRUE),
      .iqr = .q3 - .q1,
      .lo_fence = .q1 - iqr_mult * .iqr,
      .hi_fence = .q3 + iqr_mult * .iqr,
      flag_statistical_outlier = !is.na(.data[[value_col]]) &
        (.data[[value_col]] < .lo_fence | .data[[value_col]] > .hi_fence)
    ) |>
    ungroup() |>
    select(-.q1, -.q3, -.iqr, -.lo_fence, -.hi_fence)
  
  df |>
    mutate(
      exclusion_reason = case_when(
        flag_low_baseline ~ "Low baseline",
        flag_no_change    ~ "No change (untreated)",
        flag_statistical_outlier ~ "Statistical outlier (Tukey fence)",
        TRUE ~ "Included"
      )
    )
}

report_outlier_shares_20 <- function(flagged_df, group_cols = c("climate", "threshold", "band")) {
  flagged_df |>
    group_by(across(all_of(group_cols))) |>
    summarise(
      n_total = n(),
      n_low_baseline = sum(flag_low_baseline, na.rm = TRUE),
      n_no_change    = sum(flag_no_change & !flag_low_baseline, na.rm = TRUE),
      n_stat_outlier = sum(flag_statistical_outlier & !flag_low_baseline & !flag_no_change, na.rm = TRUE),
      n_included     = sum(exclusion_reason == "Included"),
      pct_low_baseline = round(100 * n_low_baseline / n_total, 1),
      pct_no_change    = round(100 * n_no_change / n_total, 1),
      pct_stat_outlier = round(100 * n_stat_outlier / n_total, 1),
      pct_included     = round(100 * n_included / n_total, 1),
      .groups = "drop"
    ) |>
    mutate(
      stat_outlier_auto_droppable = pct_stat_outlier <= OUTLIER_DROP_MAX_SHARE_PCT,
      recommendation = case_when(
        pct_stat_outlier > OUTLIER_DROP_MAX_SHARE_PCT ~
          paste0("REVIEW: statistical outliers are ", pct_stat_outlier,
                 "% of this group (> ", OUTLIER_DROP_MAX_SHARE_PCT,
                 "% threshold) — do not auto-drop, investigate cause."),
        pct_stat_outlier > 0 ~
          paste0("OK to drop: statistical outliers are only ", pct_stat_outlier,
                 "% of this group."),
        TRUE ~ "No statistical outliers detected."
      )
    )
}

# =============================================================================
# CHUNK 3 — build long-form per-watershed table
# =============================================================================

find_run_id <- function(band, climate, pct, ..., candidates) {
  climate_pat <- tolower(climate)
  pct_pat     <- paste0("top", pct)
  band_pat    <- paste0("^", tolower(band), "_")
  
  candidates_lc <- tolower(candidates)
  matches <- candidates[
    grepl(band_pat, candidates_lc) &
      grepl(climate_pat, candidates_lc) &
      grepl(pct_pat, candidates_lc)
  ]
  
  if (length(matches) == 0) return(NA_character_)
  if (length(matches) > 1) {
    warning("Multiple matches for band=", band,
            " climate=", climate,
            " pct=", pct,
            ": ", paste(matches, collapse = ", "),
            " — using first match.")
  }
  matches[1]
}


all_run_ids_20 <- list.dirs(nondeg_output_root, full.names = FALSE, recursive = FALSE)
cat("\n[Validation] run_ids found under nondeg_output_root:\n")
print(all_run_ids_20)

nondeg_registry <- expand.grid(
  band = nondeg_scenario_bands,
  climate = nondeg_climates,
  pct = nondeg_thresholds,
  stringsAsFactors = FALSE
)

nondeg_registry$run_id <- vapply(
  seq_len(nrow(nondeg_registry)),
  function(i) find_run_id(
    band = nondeg_registry$band[i],
    climate = nondeg_registry$climate[i],
    pct = nondeg_registry$pct[i],
    candidates = all_run_ids_20
  ),
  character(1)
)

cat("\n[Validation] nondeg_registry after matching:\n")
print(nondeg_registry)


n_unmatched <- sum(is.na(nondeg_registry$run_id))
if (n_unmatched > 0) {
  warning(n_unmatched, " of ", nrow(nondeg_registry),
          " band/climate/threshold combinations had NO matching run_id.")
}
stopifnot("No non-degraded scenario-band runs matched — check naming convention." =
            n_unmatched < nrow(nondeg_registry))

nondeg_run_pairs <- purrr::pmap(
  nondeg_registry[!is.na(nondeg_registry$run_id), ],
  function(band, climate, pct, run_id) {
    list(ref = nondeg_reference_by_climate[[climate]], scen = run_id,
         band = band, climate = climate, threshold = paste0("Top ", pct))
  }
)
cat("\n[Validation]", length(nondeg_run_pairs), "of", nrow(nondeg_registry),
    "possible scenario-band runs matched.\n")

build_nondeg_row <- function(p) {
  comp <- compare_runs_20(p$ref, p$scen)
  if (is.null(comp)) {
    message("Skipping missing data: ", p$scen, " vs ", p$ref)
    return(NULL)
  }
  metrics <- c("usle_tot", "sed_export", "sed_dep", "avoid_eros", "avoid_exp")
  metrics <- intersect(metrics, gsub("_pct$", "", grep("_pct$", names(comp), value = TRUE)))
  
  purrr::map_dfr(metrics, function(m) {
    tibble(
      wsid       = comp$wsid,
      metric     = m,
      pct_chg    = comp[[paste0(m, "_pct")]],
      delta      = comp[[paste0(m, "_delta")]],
      ref_val    = comp[[paste0(m, "_ref")]],
      scen_val   = comp[[paste0(m, "_scen")]],
      band       = p$band, climate = p$climate, threshold = p$threshold,
      ref_id = p$ref, scen_id = p$scen
    )
  })
}

nondeg_long_raw <- purrr::map_dfr(nondeg_run_pairs, build_nondeg_row)
stopifnot("No data returned for ANY non-degraded scenario-band run — check reference run_ids and shapefile paths." =
            nrow(nondeg_long_raw) > 0)

nondeg_long_raw <- nondeg_long_raw |>
  mutate(threshold = factor(threshold, levels = paste0("Top ", nondeg_thresholds)),
         climate   = factor(climate, levels = nondeg_climates),
         band      = factor(band, levels = nondeg_scenario_bands))


# =============================================================================
# SECTION 20A — ABSOLUTE SEDIMENT-EXPORT CHANGE ANALYSIS (P = 1)
# =============================================================================
# Insert after SECTION 20 has created `nondeg_long_raw`.
# This analysis deliberately does NOT use percentage-change filtering:
#   delta = scenario sediment export - no-AF reference sediment export
# Negative delta = reduced sediment export (beneficial).
# Positive delta = increased sediment export.
# Units inherit the InVEST watershed `sed_export` field (normally t/yr).
#
# Scientific handling:
#   * Totals and summaries retain every finite watershed delta, including zeros
#     and extreme values, so the sediment mass total is not biased.
#   * Tukey outlier flags are diagnostic only; no observation is removed.
#   * Boxplot/map display limits may clip colours or the visible axis, but the
#     CSV outputs and calculated statistics remain unclipped.
# =============================================================================



required_abs_objects <- c(
  "nondeg_long_raw", "out_20", "nondeg_scenario_bands",
  "nondeg_climates", "nondeg_thresholds", "nondeg_output_root",
  "nondeg_reference_by_climate"
)
missing_abs_objects <- required_abs_objects[
  !vapply(required_abs_objects, exists, logical(1), inherits = TRUE)
]
if (length(missing_abs_objects) > 0) {
  stop(
    "SECTION 20A requires these SECTION 20 objects: ",
    paste(missing_abs_objects, collapse = ", ")
  )
}

dir.create(out_20, recursive = TRUE, showWarnings = FALSE)

abs_band_levels <- as.character(nondeg_scenario_bands)
abs_climate_levels <- as.character(nondeg_climates)
abs_threshold_levels <- paste0("Top ", nondeg_thresholds)
ABS_NEAR_ZERO_TYR <- if (exists("ZERO_CHANGE_EPS_ABS")) ZERO_CHANGE_EPS_ABS else 5
ABS_IQR_MULT <- if (exists("OUTLIER_IQR_MULT")) OUTLIER_IQR_MULT else 3


# Accept either naming convention used in different Section 20 versions.
abs_source <- nondeg_long_raw

if (!"refval" %in% names(abs_source) && "ref_val" %in% names(abs_source)) {
  abs_source <- abs_source |> rename(refval = ref_val)
}

if (!"scenval" %in% names(abs_source) && "scen_val" %in% names(abs_source)) {
  abs_source <- abs_source |> rename(scenval = scen_val)
}

if (!"refid" %in% names(abs_source) && "ref_id" %in% names(abs_source)) {
  abs_source <- abs_source |> rename(refid = ref_id)
}

if (!"scenid" %in% names(abs_source) && "scen_id" %in% names(abs_source)) {
  abs_source <- abs_source |> rename(scenid = scen_id)
}



required_abs_cols <- c(
  "wsid", "metric", "climate", "threshold", "band",
  "refval", "scenval", "delta", "refid", "scenid"
)
missing_abs_cols <- setdiff(required_abs_cols, names(abs_source))
if (length(missing_abs_cols) > 0) {
  stop("Missing required columns in nondeg_long_raw: ",
       paste(missing_abs_cols, collapse = ", "))
}

abs_sed_all <- abs_source |>
  filter(metric == "sed_export") |>
  transmute(
    wsid,
    metric = "sed_export",
    climate = factor(as.character(climate), levels = abs_climate_levels),
    threshold = factor(as.character(threshold), levels = abs_threshold_levels),
    band = factor(as.character(band), levels = abs_band_levels),
    refval = as.numeric(refval),
    scenval = as.numeric(scenval),
    delta = as.numeric(delta),
    refid,
    scenid
  ) |>
  filter(is.finite(refval), is.finite(scenval), is.finite(delta))

if (nrow(abs_sed_all) == 0) {
  stop("No finite sed_export rows were available for absolute-change analysis.")
}

# Key uniqueness is essential before summing watershed amounts.
abs_duplicate_keys <- abs_sed_all |>
  count(wsid, climate, threshold, band, name = "n") |>
  filter(n > 1)

if (nrow(abs_duplicate_keys) > 0) {
  write_csv(
    abs_duplicate_keys,
    file.path(out_20, "VALIDATION_absolute_sedexport_duplicate_keys.csv")
  )
  stop(
    "Duplicate watershed/climate/threshold/band rows detected. ",
    "Do not sum until VALIDATION_absolute_sedexport_duplicate_keys.csv is resolved."
  )
}

# Diagnostic outlier flags only. Nothing is removed from totals or summaries.
abs_sed_flagged <- abs_sed_all |>
  group_by(climate, threshold, band) |>
  mutate(
    abs_q1_group = quantile(delta, 0.25, na.rm = TRUE),
    abs_q3_group = quantile(delta, 0.75, na.rm = TRUE),
    abs_iqr_group = abs_q3_group - abs_q1_group,
    abs_lo_fence = abs_q1_group - ABS_IQR_MULT * abs_iqr_group,
    abs_hi_fence = abs_q3_group + ABS_IQR_MULT * abs_iqr_group,
    flag_abs_outlier = if_else(
      is.finite(abs_iqr_group) & abs_iqr_group > 0,
      delta < abs_lo_fence | delta > abs_hi_fence,
      FALSE
    ),
    flag_near_zero_abs = abs(delta) <= ABS_NEAR_ZERO_TYR
  ) |>
  ungroup()

# -----------------------------------------------------------------------------
# A. Per-scenario absolute summaries and study-area totals
# -----------------------------------------------------------------------------
abs_scenario_summary <- abs_sed_flagged |>
  group_by(climate, threshold, band) |>
  summarise(
    n_watersheds = n(),
    n_reduced = sum(delta < -ABS_NEAR_ZERO_TYR),
    n_increased = sum(delta > ABS_NEAR_ZERO_TYR),
    n_near_zero = sum(flag_near_zero_abs),
    n_abs_outlier_flagged = sum(flag_abs_outlier),
    total_reference_tyr = sum(refval),
    total_scenario_tyr = sum(scenval),
    total_delta_tyr = sum(delta),
    net_reduction_tyr = -sum(delta),
    gross_reduction_tyr = sum(-delta[delta < 0]),
    gross_increase_tyr = sum(delta[delta > 0]),
    mean_delta_tyr = mean(delta),
    sd_delta_tyr = sd(delta),
    min_delta_tyr = min(delta),
    q05_delta_tyr = quantile(delta, 0.05),
    q25_delta_tyr = quantile(delta, 0.25),
    median_delta_tyr = median(delta),
    q75_delta_tyr = quantile(delta, 0.75),
    q95_delta_tyr = quantile(delta, 0.95),
    max_delta_tyr = max(delta),
    .groups = "drop"
  ) |>
  mutate(
    delta_from_aggregate_totals_tyr = total_scenario_tyr - total_reference_tyr,
    mass_balance_error_tyr = total_delta_tyr - delta_from_aggregate_totals_tyr,
    pct_watersheds_reduced = 100 * n_reduced / n_watersheds,
    pct_watersheds_increased = 100 * n_increased / n_watersheds,
    pct_watersheds_near_zero = 100 * n_near_zero / n_watersheds,
    pct_abs_outlier_flagged = 100 * n_abs_outlier_flagged / n_watersheds
  ) |>
  arrange(climate, threshold, band)

write_csv(
  abs_scenario_summary,
  file.path(out_20, "ABS_sedexport_per_scenario_summary.csv")
)

# Scenario-level total range for each climate/threshold.
abs_total_band_range <- abs_scenario_summary |>
  select(climate, threshold, band, total_delta_tyr, net_reduction_tyr) |>
  pivot_wider(
    names_from = band,
    values_from = c(total_delta_tyr, net_reduction_tyr),
    names_glue = "{.value}_{band}"
  ) |>
  rowwise() |>
  mutate(
    total_delta_range_low_tyr = min(c_across(starts_with("total_delta_tyr_")), na.rm = TRUE),
    total_delta_range_high_tyr = max(c_across(starts_with("total_delta_tyr_")), na.rm = TRUE),
    total_delta_range_width_tyr = total_delta_range_high_tyr - total_delta_range_low_tyr,
    net_reduction_range_low_tyr = min(c_across(starts_with("net_reduction_tyr_")), na.rm = TRUE),
    net_reduction_range_high_tyr = max(c_across(starts_with("net_reduction_tyr_")), na.rm = TRUE),
    net_reduction_range_width_tyr = net_reduction_range_high_tyr - net_reduction_range_low_tyr
  ) |>
  ungroup()

write_csv(
  abs_total_band_range,
  file.path(out_20, "ABS_sedexport_total_scenario_band_range.csv")
)

# -----------------------------------------------------------------------------
# B. Watershed scenario range width and pairwise differences
# -----------------------------------------------------------------------------
abs_watershed_wide <- abs_sed_all |>
  select(wsid, climate, threshold, band, delta) |>
  pivot_wider(names_from = band, values_from = delta)

# Ensure all expected scenario columns exist even if a run was absent.
for (bnd in abs_band_levels) {
  if (!bnd %in% names(abs_watershed_wide)) abs_watershed_wide[[bnd]] <- NA_real_
}

abs_watershed_range <- abs_watershed_wide |>
  rowwise() |>
  mutate(
    n_valid_bands = sum(is.finite(c_across(all_of(abs_band_levels)))),
    abs_range_low_tyr = if (n_valid_bands >= 2) {
      min(c_across(all_of(abs_band_levels)), na.rm = TRUE)
    } else NA_real_,
    abs_range_high_tyr = if (n_valid_bands >= 2) {
      max(c_across(all_of(abs_band_levels)), na.rm = TRUE)
    } else NA_real_,
    abs_range_width_tyr = if (n_valid_bands >= 2) {
      abs_range_high_tyr - abs_range_low_tyr
    } else NA_real_
  ) |>
  ungroup()

write_csv(
  abs_watershed_range,
  file.path(out_20, "ABS_sedexport_watershed_scenario_range_wide.csv")
)

abs_range_summary <- abs_watershed_range |>
  filter(is.finite(abs_range_width_tyr)) |>
  group_by(climate, threshold) |>
  summarise(
    n_watersheds = n(),
    mean_range_width_tyr = mean(abs_range_width_tyr),
    median_range_width_tyr = median(abs_range_width_tyr),
    q25_range_width_tyr = quantile(abs_range_width_tyr, 0.25),
    q75_range_width_tyr = quantile(abs_range_width_tyr, 0.75),
    q95_range_width_tyr = quantile(abs_range_width_tyr, 0.95),
    max_range_width_tyr = max(abs_range_width_tyr),
    .groups = "drop"
  )

write_csv(
  abs_range_summary,
  file.path(out_20, "ABS_sedexport_watershed_range_summary.csv")
)

abs_pairwise <- abs_watershed_wide |>
  mutate(
    optimistic_minus_conservative_tyr = optimistic - conservative,
    optimistic_minus_central_tyr = optimistic - central,
    central_minus_conservative_tyr = central - conservative
  ) |>
  select(
    wsid, climate, threshold,
    optimistic_minus_conservative_tyr,
    optimistic_minus_central_tyr,
    central_minus_conservative_tyr
  ) |>
  pivot_longer(
    cols = ends_with("_tyr"),
    names_to = "scenario_pair",
    values_to = "delta_difference_tyr"
  ) |>
  mutate(
    scenario_pair = recode(
      scenario_pair,
      optimistic_minus_conservative_tyr = "optimistic minus conservative",
      optimistic_minus_central_tyr = "optimistic minus central",
      central_minus_conservative_tyr = "central minus conservative"
    )
  )

write_csv(
  abs_pairwise,
  file.path(out_20, "ABS_sedexport_pairwise_scenario_differences.csv")
)

# Full finite watershed table: retain for detailed follow-up validation.
write_csv(
  abs_sed_flagged,
  file.path(out_20, "ABS_sedexport_watershed_all_flagged.csv")
)

# -----------------------------------------------------------------------------
# C. Boxplots and scenario-total composite
# -----------------------------------------------------------------------------
abs_box_limits <- unname(quantile(abs_sed_all$delta, c(0.01, 0.99), na.rm = TRUE))
if (!all(is.finite(abs_box_limits)) || diff(abs_box_limits) <= 0) {
  abs_box_limits <- range(abs_sed_all$delta, finite = TRUE)
}
if (diff(abs_box_limits) <= 0) abs_box_limits <- abs_box_limits + c(-1, 1)

abs_palette <- c(
  conservative = "#D95F02",
  central = "#1B9E77",
  optimistic = "#7570B3"
)

p_abs_box_all <- ggplot(abs_sed_all, aes(x = band, y = delta, fill = band)) +
  geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.35) +
  geom_boxplot(width = 0.68, outlier.shape = NA, alpha = 0.85) +
  facet_grid(climate ~ threshold) +
  coord_cartesian(ylim = abs_box_limits) +
  scale_fill_manual(values = abs_palette, drop = FALSE) +
  scale_y_continuous(labels = label_number(big.mark = ",", accuracy = 1)) +
  labs(
    title = "Absolute sediment-export change by AF scenario",
    subtitle = "Watershed distributions by climate and suitability threshold (P = 1)",
    x = NULL,
    y = "Scenario minus no-AF sediment export (t/yr)",
    caption = paste0(
      "Negative = reduced export (beneficial). Boxplot axis zoomed to the 1st–99th percentiles; ",
      "all finite values remain in CSV summaries and totals."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

ggsave(
  file.path(out_20, "figure20_absolute_boxplots_all_scenarios.png"),
  p_abs_box_all, width = 14, height = 8, dpi = 220
)

p_abs_box_baseline <- abs_sed_all |>
  filter(climate == "Baseline") |>
  ggplot(aes(x = threshold, y = delta, fill = band)) +
  geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.35) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  coord_cartesian(ylim = abs_box_limits) +
  scale_fill_manual(values = abs_palette, drop = FALSE) +
  scale_y_continuous(labels = label_number(big.mark = ",", accuracy = 1)) +
  labs(
    title = "Absolute sediment-export change — Baseline climate",
    subtitle = "Scenario distributions across top-suitability thresholds (P = 1)",
    x = "Suitability threshold",
    y = "Scenario minus no-AF sediment export (t/yr)",
    fill = "Scenario",
    caption = "Negative = reduced export. Visible axis uses the common 1st–99th percentile zoom."
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave(
  file.path(out_20, "figure20_absolute_boxplots_baseline.png"),
  p_abs_box_baseline, width = 11, height = 7, dpi = 220
)

p_abs_total_composite <- ggplot(
  abs_scenario_summary,
  aes(x = threshold, y = net_reduction_tyr, colour = band, group = band)
) +
  geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.35) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.6) +
  facet_wrap(~climate, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = abs_palette, drop = FALSE) +
  scale_y_continuous(labels = label_number(big.mark = ",", accuracy = 1)) +
  labs(
    title = "Total sediment-export reduction by scenario",
    subtitle = "Sum across all finite watershed changes; positive values indicate a net reduction",
    x = "Suitability threshold",
    y = "Net reduction in sediment export (t/yr)",
    colour = "Scenario",
    caption = "Net reduction = no-AF total minus AF-scenario total; no outlier removal or visual clipping."
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave(
  file.path(out_20, "figure20_absolute_total_scenario_composite.png"),
  p_abs_total_composite, width = 12, height = 6.5, dpi = 220
)

# -----------------------------------------------------------------------------
# D. Absolute-change maps: scenario composite, range width, pairwise differences
# -----------------------------------------------------------------------------
abs_map_threshold <- "Top 25"
abs_map_geometry_found <- FALSE
abs_map_n_geometry <- NA_integer_
abs_map_n_unmatched <- NA_integer_

wsgeom_abs_path <- file.path(
  nondeg_output_root,
  nondeg_reference_by_climate[["Baseline"]],
  "watershed_results_sdr.shp"
)
if (!file.exists(wsgeom_abs_path) && exists("invest_output_root")) {
  alt_abs_path <- file.path(
    invest_output_root,
    nondeg_reference_by_climate[["Baseline"]],
    "watershed_results_sdr.shp"
  )
  if (file.exists(alt_abs_path)) wsgeom_abs_path <- alt_abs_path
}

if (file.exists(wsgeom_abs_path)) {
  wsgeom_abs <- st_read(wsgeom_abs_path, quiet = TRUE)
  names(wsgeom_abs) <- tolower(names(wsgeom_abs))
  if (!"wsid" %in% names(wsgeom_abs)) wsgeom_abs$wsid <- seq_len(nrow(wsgeom_abs))
  wsgeom_abs <- wsgeom_abs["wsid"]
  
  abs_map_geometry_found <- TRUE
  abs_map_n_geometry <- nrow(wsgeom_abs)
  
  abs_map_data <- abs_sed_all |>
    filter(as.character(threshold) == abs_map_threshold) |>
    select(wsid, climate, threshold, band, delta, refval, scenval)
  
  abs_map_n_unmatched <- sum(!abs_map_data$wsid %in% wsgeom_abs$wsid)
  map_abs_sf <- wsgeom_abs |> left_join(abs_map_data, by = "wsid")
  
  abs_map_values <- map_abs_sf$delta[is.finite(map_abs_sf$delta)]
  abs_map_limit <- unname(quantile(abs(abs_map_values), 0.98, na.rm = TRUE))
  if (!is.finite(abs_map_limit) || abs_map_limit <= 0) {
    abs_map_limit <- max(abs(abs_map_values), na.rm = TRUE)
  }
  if (!is.finite(abs_map_limit) || abs_map_limit <= 0) abs_map_limit <- 1
  
  p_abs_map_composite <- ggplot(map_abs_sf) +
    geom_sf(aes(fill = delta), colour = "white", linewidth = 0.06) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0,
      limits = c(-abs_map_limit, abs_map_limit),
      oob = squish,
      na.value = "grey82",
      name = expression(Delta~"sediment export"~"(t/yr)")
    ) +
    facet_grid(climate ~ band) +
    labs(
      title = "Absolute sediment-export change — scenario composite",
      subtitle = paste0(abs_map_threshold, "% suitability; common colour scale across all panels (P = 1)"),
      caption = "Blue = reduced export; red = increased export; grey = unavailable. Colours clipped at the 98th percentile of |change|."
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text = element_blank(), axis.ticks = element_blank(),
      panel.grid = element_blank(), strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    coord_sf(datum = NA)
  
  ggsave(
    file.path(out_20, "map20_absolute_scenario_composite_top25.png"),
    p_abs_map_composite, width = 15, height = 9, dpi = 220
  )
  
  abs_range_map_data <- abs_watershed_range |>
    filter(as.character(threshold) == abs_map_threshold) |>
    select(wsid, climate, threshold, abs_range_low_tyr,
           abs_range_high_tyr, abs_range_width_tyr, n_valid_bands)
  
  map_abs_range_sf <- wsgeom_abs |> left_join(abs_range_map_data, by = "wsid")
  abs_range_values <- map_abs_range_sf$abs_range_width_tyr[
    is.finite(map_abs_range_sf$abs_range_width_tyr)
  ]
  abs_range_limit <- unname(quantile(abs_range_values, 0.98, na.rm = TRUE))
  if (!is.finite(abs_range_limit) || abs_range_limit <= 0) {
    abs_range_limit <- max(abs_range_values, na.rm = TRUE)
  }
  if (!is.finite(abs_range_limit) || abs_range_limit <= 0) abs_range_limit <- 1
  
  p_abs_range <- ggplot(map_abs_range_sf) +
    geom_sf(aes(fill = abs_range_width_tyr), colour = "white", linewidth = 0.08) +
    scale_fill_distiller(
      palette = "YlOrRd", direction = 1,
      limits = c(0, abs_range_limit), oob = squish,
      na.value = "grey82", name = "Scenario range\nwidth (t/yr)"
    ) +
    facet_wrap(~climate, nrow = 1) +
    labs(
      title = "Absolute-change scenario range width",
      subtitle = paste0(abs_map_threshold, "% suitability; maximum minus minimum scenario delta per watershed"),
      caption = "Larger values indicate greater sensitivity to the selected C-factor scenario. Colours clipped at the 98th percentile."
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text = element_blank(), axis.ticks = element_blank(),
      panel.grid = element_blank(), strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    coord_sf(datum = NA)
  
  ggsave(
    file.path(out_20, "map20_absolute_range_width_top25.png"),
    p_abs_range, width = 13, height = 7, dpi = 220
  )
  
  abs_pair_map_data <- abs_pairwise |>
    filter(as.character(threshold) == abs_map_threshold)
  map_abs_pair_sf <- wsgeom_abs |> left_join(abs_pair_map_data, by = "wsid")
  
  abs_pair_values <- map_abs_pair_sf$delta_difference_tyr[
    is.finite(map_abs_pair_sf$delta_difference_tyr)
  ]
  abs_pair_limit <- unname(quantile(abs(abs_pair_values), 0.98, na.rm = TRUE))
  if (!is.finite(abs_pair_limit) || abs_pair_limit <= 0) {
    abs_pair_limit <- max(abs(abs_pair_values), na.rm = TRUE)
  }
  if (!is.finite(abs_pair_limit) || abs_pair_limit <= 0) abs_pair_limit <- 1
  
  p_abs_pair <- ggplot(map_abs_pair_sf) +
    geom_sf(aes(fill = delta_difference_tyr), colour = "white", linewidth = 0.05) +
    scale_fill_distiller(
      palette = "PuOr", direction = 1,
      limits = c(-abs_pair_limit, abs_pair_limit), oob = squish,
      na.value = "grey82", name = "Difference in\nabsolute change (t/yr)"
    ) +
    facet_grid(climate ~ scenario_pair) +
    labs(
      title = "Pairwise scenario differences in absolute sediment-export change",
      subtitle = paste0(abs_map_threshold, "% suitability; larger magnitude means stronger scenario sensitivity"),
      caption = "Each value is the first named scenario delta minus the second. Colours clipped at the 98th percentile of |difference|."
    ) +
    theme_minimal(base_size = 10.5) +
    theme(
      axis.text = element_blank(), axis.ticks = element_blank(),
      panel.grid = element_blank(), strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    coord_sf(datum = NA)
  
  ggsave(
    file.path(out_20, "map20_absolute_pairwise_differences_top25.png"),
    p_abs_pair, width = 15, height = 9, dpi = 220
  )
  
  write_csv(
    st_drop_geometry(map_abs_sf),
    file.path(out_20, "ABS_mapdata_scenario_composite_top25.csv")
  )
  write_csv(
    st_drop_geometry(map_abs_range_sf),
    file.path(out_20, "ABS_mapdata_range_width_top25.csv")
  )
  write_csv(
    st_drop_geometry(map_abs_pair_sf),
    file.path(out_20, "ABS_mapdata_pairwise_differences_top25.csv")
  )
} else {
  message("Absolute-change maps skipped; watershed geometry not found at: ",
          wsgeom_abs_path)
}

# -----------------------------------------------------------------------------
# E. Compact validation bundle — upload this CSV back into the thread
# -----------------------------------------------------------------------------
abs_group_completeness <- abs_scenario_summary |>
  group_by(climate, threshold) |>
  summarise(
    n_scenario_bands_found = n_distinct(band),
    all_scenario_bands_present =
      n_scenario_bands_found == length(abs_band_levels),
    reference_total_range_across_bands_tyr =
      max(total_reference_tyr) - min(total_reference_tyr),
    .groups = "drop"
  )

abs_reference_wsid_check <- abs_sed_all |>
  group_by(wsid, climate, threshold) |>
  summarise(
    reference_spread_tyr = max(refval) - min(refval),
    .groups = "drop"
  ) |>
  group_by(climate, threshold) |>
  summarise(
    n_reference_mismatch_watersheds = sum(reference_spread_tyr > 1e-8),
    max_reference_spread_tyr = max(reference_spread_tyr),
    .groups = "drop"
  )

abs_validation_bundle <- abs_scenario_summary |>
  left_join(abs_group_completeness, by = c("climate", "threshold")) |>
  left_join(abs_reference_wsid_check, by = c("climate", "threshold")) |>
  mutate(
    duplicate_key_count = nrow(abs_duplicate_keys),
    geometry_found = abs_map_geometry_found,
    n_geometry_features = abs_map_n_geometry,
    n_map_rows_without_geometry = abs_map_n_unmatched,
    mass_balance_tolerance_tyr = pmax(
      1e-6,
      1e-10 * pmax(abs(total_reference_tyr), abs(total_scenario_tyr), 1)
    ),
    mass_balance_pass =
      abs(mass_balance_error_tyr) <= mass_balance_tolerance_tyr,
    reference_consistency_pass =
      n_reference_mismatch_watersheds == 0 &
      reference_total_range_across_bands_tyr <= mass_balance_tolerance_tyr,
    validation_status = if_else(
      duplicate_key_count == 0 &
        all_scenario_bands_present &
        mass_balance_pass &
        reference_consistency_pass,
      "PASS",
      "CHECK"
    ),
    delta_definition = "scenario minus no-AF reference",
    sign_interpretation = "negative delta = reduced export; positive net_reduction = benefit",
    totals_filtering = "all finite watershed values retained; outliers diagnostic only"
  ) |>
  arrange(climate, threshold, band)

write_csv(
  abs_validation_bundle,
  file.path(out_20, "VALIDATION_absolute_sedexport_bundle.csv")
)

cat("\nSECTION 20A COMPLETE — absolute sediment-export analysis\n")
cat("Primary validation file to upload:\n  ",
    file.path(out_20, "VALIDATION_absolute_sedexport_bundle.csv"), "\n", sep = "")
cat("Detailed watershed file for deeper validation:\n  ",
    file.path(out_20, "ABS_sedexport_watershed_all_flagged.csv"), "\n", sep = "")
cat("Validation status counts:\n")
print(abs_validation_bundle |> count(validation_status))
cat("Absolute scenario totals:\n")
print(
  abs_scenario_summary |>
    select(climate, threshold, band, n_watersheds,
           total_delta_tyr, net_reduction_tyr, median_delta_tyr)
) 



nondeg_long <- flag_outliers_20(nondeg_long_raw, group_cols = c("metric", "climate", "threshold", "band"))

outlier_share_report <- report_outlier_shares_20(nondeg_long, group_cols = c("metric", "climate", "threshold", "band"))
write_csv(outlier_share_report, file.path(out_20, "VALIDATION_outlier_share_report.csv"))

cat("\n[Validation] Outlier / exclusion share report (sed_export only, printed):\n")
print(outlier_share_report |> filter(metric == "sed_export"))

groups_needing_review <- outlier_share_report |> filter(!stat_outlier_auto_droppable, pct_stat_outlier > 0)
if (nrow(groups_needing_review) > 0) {
  cat("\n*** WARNING: the following groups have a LARGE share of statistical outliers",
      "(>", OUTLIER_DROP_MAX_SHARE_PCT, "%) — these are NOT auto-dropped. Review before interpreting. ***\n")
  print(groups_needing_review |> select(metric, climate, threshold, band, pct_stat_outlier, recommendation))
}

nondeg_long <- nondeg_long |>
  left_join(
    outlier_share_report |> select(metric, climate, threshold, band, stat_outlier_auto_droppable),
    by = c("metric", "climate", "threshold", "band")
  ) |>
  mutate(
    drop_for_filtered = case_when(
      exclusion_reason == "Low baseline" ~ TRUE,
      exclusion_reason == "No change (untreated)" ~ TRUE,
      exclusion_reason == "Statistical outlier (Tukey fence)" & stat_outlier_auto_droppable ~ TRUE,
      TRUE ~ FALSE
    )
  )

write_csv(nondeg_long, file.path(out_20, "nondegraded_watershed_ALL_flagged.csv"))
nondeg_filtered <- nondeg_long |> filter(!drop_for_filtered)
write_csv(nondeg_filtered, file.path(out_20, "nondegraded_watershed_FILTERED.csv"))




cat("\n[Validation] Rows total:", nrow(nondeg_long),
    " | Rows retained after filtering:", nrow(nondeg_filtered),
    " (", round(100 * nrow(nondeg_filtered) / nrow(nondeg_long), 1), "% retained)\n")




# =============================================================================
# CHUNK 4 — range/summary tables
# =============================================================================

nondeg_scenario_summary <- nondeg_filtered |>
  group_by(metric, climate, threshold, band) |>
  summarise(
    n_ws = n(),
    mean_pct   = mean(pct_chg, na.rm = TRUE),
    median_pct = median(pct_chg, na.rm = TRUE),
    sd_pct     = sd(pct_chg, na.rm = TRUE),
    q25_pct    = quantile(pct_chg, 0.25, na.rm = TRUE),
    q75_pct    = quantile(pct_chg, 0.75, na.rm = TRUE),
    mean_delta   = mean(delta, na.rm = TRUE),
    median_delta = median(delta, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(nondeg_scenario_summary, file.path(out_20, "PERSCENARIO_summary.csv"))

cat("\n[Chunk 4a] Per-scenario summary (sed_export):\n")
print(nondeg_scenario_summary |> filter(metric == "sed_export"))

nondeg_band_range <- nondeg_scenario_summary |>
  select(metric, climate, threshold, band, median_pct, median_delta) |>
  pivot_wider(names_from = band, values_from = c(median_pct, median_delta),
              names_glue = "{.value}_{band}") |>
  mutate(
    range_pct_low   = pmin(median_pct_conservative, median_pct_central, median_pct_optimistic, na.rm = TRUE),
    range_pct_high  = pmax(median_pct_conservative, median_pct_central, median_pct_optimistic, na.rm = TRUE),
    range_pct_width = range_pct_high - range_pct_low
  )
write_csv(nondeg_band_range, file.path(out_20, "BAND_range_wide.csv"))

cat("\n[Chunk 4b] Conservative-central-optimistic band range (sed_export):\n")
print(nondeg_band_range |> filter(metric == "sed_export"))

# =============================================================================
# CHUNK 5 — plots
# =============================================================================

climate_cols <- c(Baseline = "#4393C3", SSP370 = "#D6604D")
scenario_cols <- c(conservative = "#4393C3", central = "#F4A582", optimistic = "#D6604D")

nondeg_plot_df <- nondeg_scenario_summary |>
  filter(metric == "sed_export") |>
  select(climate, threshold, band, median_pct) |>
  pivot_wider(names_from = band, values_from = median_pct)

p_range_band <- ggplot(nondeg_plot_df, aes(x = threshold, group = climate)) +
  geom_ribbon(aes(ymin = conservative, ymax = optimistic, fill = climate), alpha = 0.25) +
  geom_line(aes(y = central, colour = climate), linewidth = 1) +
  geom_point(aes(y = central, colour = climate), size = 2.4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = climate_cols, name = "Climate") +
  scale_colour_manual(values = climate_cols, name = "Climate") +
  labs(title = "Sediment export change range — non-degraded AF3-class transition (P = 1)",
       subtitle = "Shaded band = conservative-to-optimistic C-factor scenario range. Line = central estimate.",
       x = "Suitability threshold", y = "% change in sed_export vs. no-AF reference",
       caption = "Median values across FILTERED watersheds (outliers/low-baseline/no-change excluded per validation report)") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(out_20, "fig20a_sedexport_range_band.png"), p_range_band, width = 10, height = 7, dpi = 180)

p_band_box <- nondeg_filtered |>
  filter(metric == "sed_export") |>
  ggplot(aes(x = threshold, y = pct_chg, fill = band)) +
  geom_boxplot(position = position_dodge(0.75), width = 0.6, outlier.size = 0.5, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = scenario_cols, name = "Scenario band") +
  facet_wrap(~climate, ncol = 2) +
  coord_cartesian(ylim = c(-125, 100)) +
  labs(title = "Watershed-level sed_export change by scenario band (P = 1, filtered)",
       subtitle = "Filtered set excludes low-baseline, untreated, and statistical-outlier watersheds",
       x = "Suitability threshold", y = "% change in sed_export",
       caption = "Y-axis clipped to [-125, 100] for readability") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(out_20, "fig20b_band_boxplot_distribution.png"), p_band_box, width = 11, height = 6, dpi = 180)

p_perscenario_bar <- nondeg_scenario_summary |>
  filter(metric == "sed_export", climate == "Baseline") |>
  ggplot(aes(x = band, y = median_pct, fill = band)) +
  geom_col(width = 0.6, alpha = 0.9) +
  geom_errorbar(aes(ymin = q25_pct, ymax = q75_pct), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = scenario_cols, guide = "none") +
  facet_wrap(~threshold, nrow = 1) +
  labs(title = "Individual scenario summary — sed_export % change (Baseline climate, P = 1)",
       subtitle = "Bar = median; error bars = IQR (25th-75th percentile), filtered watersheds",
       x = "Scenario band", y = "% change in sed_export") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
ggsave(file.path(out_20, "fig20c_perscenario_bar_summary.png"), p_perscenario_bar, width = 10, height = 5, dpi = 180)

cat("\n[Chunk 5] Plots written: fig20a (range band), fig20b (boxplot), fig20c (per-scenario bars)\n")

# =============================================================================
# CHUNK 6 — spatial maps
# =============================================================================

wsgeom_path <- file.path(nondeg_output_root, nondeg_reference_by_climate[["Baseline"]], "watershed_results_sdr.shp")
if (!file.exists(wsgeom_path) && exists("invest_output_root")) {
  alt <- file.path(invest_output_root, nondeg_reference_by_climate[["Baseline"]], "watershed_results_sdr.shp")
  if (file.exists(alt)) wsgeom_path <- alt
}

if (file.exists(wsgeom_path)) {
  wsgeom <- st_read(wsgeom_path, quiet = TRUE)
  names(wsgeom) <- tolower(names(wsgeom))
  if (!"wsid" %in% names(wsgeom)) wsgeom$wsid <- seq_len(nrow(wsgeom))
  wsgeom <- wsgeom["wsid"]
  
  map_data_scen <- nondeg_long |>
    filter(metric == "sed_export", climate == "Baseline", threshold == "Top 25") |>
    select(wsid, band, pct_chg, delta, exclusion_reason)
  
  map_sf_scen <- wsgeom |> left_join(map_data_scen, by = "wsid")
  pct_limit_scen <- max(abs(quantile(map_sf_scen$pct_chg, c(0.02, 0.98), na.rm = TRUE)))
  
  for (bnd in nondeg_scenario_bands) {
    df_b <- map_sf_scen |> filter(band == bnd)
    p_ind <- ggplot(df_b) +
      geom_sf(aes(fill = pct_chg), colour = "white", linewidth = 0.1) +
      scale_fill_distiller(palette = "RdBu", direction = 1,
                           limits = c(-pct_limit_scen, pct_limit_scen), oob = scales::squish,
                           name = "% change\nsed_export", na.value = "grey85") +
      labs(title = paste0("Sediment export change — ", bnd, " scenario (P = 1)"),
           subtitle = "Baseline climate, top 25% suitability, vs. no-AF reference",
           caption = "Blue = reduced sediment export (beneficial); Red = increased export. Grey = excluded.") +
      theme_minimal(base_size = 12) +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            panel.grid = element_blank(), legend.position = "right") +
      coord_sf(datum = NA)
    ggsave(file.path(out_20, paste0("map20_individual_", bnd, "_sedexport.png")),
           p_ind, width = 9, height = 8, dpi = 200)
  }
  
  p_composite <- ggplot(map_sf_scen) +
    geom_sf(aes(fill = pct_chg), colour = "white", linewidth = 0.08) +
    scale_fill_distiller(palette = "RdBu", direction = 1,
                         limits = c(-pct_limit_scen, pct_limit_scen), oob = scales::squish,
                         name = "% change\nsed_export", na.value = "grey85") +
    facet_wrap(~band, nrow = 1) +
    labs(title = "Sediment export change — scenario composite comparison (P = 1)",
         subtitle = "Baseline climate, top 25% suitability. Common colour scale across panels for direct comparison.",
         caption = "Blue = reduced export (beneficial); Red = increased export. Grey = excluded (see validation report).") +
    theme_minimal(base_size = 12) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), strip.text = element_text(face = "bold"),
          legend.position = "bottom") +
    coord_sf(datum = NA)
  ggsave(file.path(out_20, "map20_composite_all_scenarios.png"), p_composite, width = 15, height = 6, dpi = 200)
  
  map_wide_scen <- map_data_scen |>
    select(wsid, band, pct_chg) |>
    pivot_wider(names_from = band, values_from = pct_chg)
  
  diff_pairs <- list(
    c("optimistic", "conservative"),
    c("optimistic", "central"),
    c("central", "conservative")
  )
  
  diff_maps_list <- list()
  for (pr in diff_pairs) {
    hi <- pr[1]; lo <- pr[2]
    diff_col <- paste0(hi, "_minus_", lo)
    map_wide_scen[[diff_col]] <- map_wide_scen[[hi]] - map_wide_scen[[lo]]
    diff_maps_list[[diff_col]] <- map_wide_scen |> select(wsid, value = all_of(diff_col)) |>
      mutate(pair = paste0(hi, " minus ", lo))
  }
  diff_long <- bind_rows(diff_maps_list)
  diff_sf <- wsgeom |> left_join(diff_long, by = "wsid")
  
  diff_limit <- max(abs(quantile(diff_sf$value, c(0.02, 0.98), na.rm = TRUE)))
  p_diff <- ggplot(diff_sf) +
    geom_sf(aes(fill = value), colour = "white", linewidth = 0.08) +
    scale_fill_distiller(palette = "PuOr", direction = 1,
                         limits = c(-diff_limit, diff_limit), oob = scales::squish,
                         name = "Scenario\ndifference\n(pct pts)", na.value = "grey85") +
    facet_wrap(~pair, nrow = 1) +
    labs(title = "Where does scenario choice matter most? Pairwise scenario differences (P = 1)",
         subtitle = "Baseline climate, top 25% suitability. Difference in % change in sed_export between scenarios.",
         caption = "Larger magnitude = greater sensitivity of local outcome to which C-factor scenario is assumed.") +
    theme_minimal(base_size = 12) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), strip.text = element_text(face = "bold"),
          legend.position = "bottom") +
    coord_sf(datum = NA)
  ggsave(file.path(out_20, "map20_pairwise_scenario_differences.png"), p_diff, width = 15, height = 6, dpi = 200)
  
  write_csv(st_drop_geometry(diff_sf), file.path(out_20, "mapdata_pairwise_scenario_differences.csv"))
  
  
  ws_band_wide <- nondeg_long |>
    filter(
      metric == "sed_export",
      climate == "Baseline",
      threshold == "Top 25"
    ) |>
    select(wsid, band, pct_chg, exclusion_reason) |>
    pivot_wider(names_from = band, values_from = pct_chg) |>
    rowwise() |>
    mutate(
      n_valid_bands = sum(is.finite(c(conservative, central, optimistic))),
      range_width = if_else(
        n_valid_bands >= 2,
        max(c(conservative, central, optimistic), na.rm = TRUE) -
          min(c(conservative, central, optimistic), na.rm = TRUE),
        NA_real_
      )
    ) |>
    ungroup()
  
  
  map_range_sf <- wsgeom |> left_join(ws_band_wide, by = "wsid")
  range_limit <- quantile(map_range_sf$range_width, 0.98, na.rm = TRUE)
  
  p_map_range <- ggplot(map_range_sf) +
    geom_sf(aes(fill = range_width), colour = "white", linewidth = 0.15) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1,
                         limits = c(0, range_limit), oob = scales::squish, na.value = "grey80",
                         name = "Range width\n(pct pts)") +
    labs(title = "Uncertainty range width — sed_export change (Baseline, top 25, P = 1)",
         subtitle = "Width of conservative-to-optimistic scenario band per watershed",
         caption = "Grey = excluded (low baseline / untreated / stat. outlier). Colour clipped at 98th percentile.") +
    theme_minimal(base_size = 12) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), legend.position = "right") +
    coord_sf(datum = NA)
  ggsave(file.path(out_20, "map20d_range_width_choropleth.png"), p_map_range, width = 10, height = 9, dpi = 220)
  
  p_map_central <- ggplot(map_range_sf) +
    geom_sf(aes(fill = central), colour = "white", linewidth = 0.15) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
                         na.value = "grey80", name = "% change\n(central)") +
    labs(title = "Sediment export change — central estimate (Baseline, top 25, P = 1)",
         subtitle = "Non-degraded AF3-class transition, central C-factor scenario",
         caption = "Grey = excluded (low baseline / untreated / stat. outlier)") +
    theme_minimal(base_size = 12) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), legend.position = "right") +
    coord_sf(datum = NA)
  ggsave(file.path(out_20, "map20e_central_estimate_choropleth.png"), p_map_central, width = 10, height = 9, dpi = 220)
  
  
  
  # ---- MAP: absolute sediment-export change, central scenario ----------------
  # delta is scenario sediment export minus reference sediment export.
  # Negative values = less sediment export under AF (beneficial).
  # Positive values = more sediment export under AF.
  
  ws_delta_abs <- nondeg_long |>
    filter(
      metric == "sed_export",
      climate == "Baseline",
      threshold == "Top 25",
      band == "central"
    ) |>
    select(wsid, delta, ref_val, scen_val, exclusion_reason)
  
  map_delta_abs_sf <- wsgeom |>
    left_join(ws_delta_abs, by = "wsid")
  
  # Use the 98th percentile of absolute valid changes for a symmetric,
  # zero-centred colour scale. This prevents a few extreme watersheds
  # from flattening the colour contrast across the map.
  delta_vals <- map_delta_abs_sf$delta[
    is.finite(map_delta_abs_sf$delta)
  ]
  
  stopifnot(length(delta_vals) > 0)
  
  delta_lim <- unname(quantile(abs(delta_vals), 0.98, na.rm = TRUE))
  
  # Failsafe in case all valid changes are exactly zero.
  if (!is.finite(delta_lim) || delta_lim == 0) {
    delta_lim <- max(abs(delta_vals), na.rm = TRUE)
  }
  if (!is.finite(delta_lim) || delta_lim == 0) {
    delta_lim <- 1
  }
  
  p_map_delta_abs <- ggplot(map_delta_abs_sf) +
    geom_sf(aes(fill = delta), colour = "white", linewidth = 0.15) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-delta_lim, delta_lim),
      oob = scales::squish,
      na.value = "grey80",
      name = expression(Delta~"sediment export"~"(t yr"^{-1}*")")
    ) +
    labs(
      title = "Sediment export amount change — central estimate (Baseline, top 25, P = 1)",
      subtitle = "Absolute watershed change relative to the no-AF reference",
      caption = paste0(
        "Blue = reduced sediment export (beneficial); red = increased export. ",
        "Grey = excluded, untreated, or unavailable. Colours clipped at the 98th percentile of |change|."
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    coord_sf(datum = NA)
  
  ggsave(
    file.path(out_20, "map20f_central_sedexport_absolute_change.png"),
    p_map_delta_abs,
    width = 10,
    height = 9,
    dpi = 220
  )
  
  # Export the exact map values and their reference/scenario amounts.
  write_csv(
    st_drop_geometry(map_delta_abs_sf),
    file.path(out_20, "mapdata_central_sedexport_absolute_change_baseline_top25.csv")
  )
  
  
  write_csv(st_drop_geometry(map_range_sf), file.path(out_20, "mapdata_range_width_baseline_top25.csv"))
  
} else {
  message("Skipping spatial maps — watershed geometry not found at: ", wsgeom_path)
}

cat("\n=====================================================\n")
cat("SECTION 20 COMPLETE — outputs written to:", out_20, "\n")
cat("=====================================================\n") 


