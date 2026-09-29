# ================================================================
# InVEST SDR Watershed Results: Multi-Scenario Sensitivity Analysis
# ================================================================
# Blocks:
#   1. AF vs Baseline comparisons (within climate)
#   2. Climate comparisons (baseline vs SSP3-7.0, same C)
#   3. Cross comparisons (baseline vs SSP3-7.0 + AF)
# ================================================================

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# ----------------------------------------------------------------
# CONFIGURATION — update paths to your InVEST output directories
# ----------------------------------------------------------------

runs_parent <- "C:/CI_Alves/ssa_invest_runs/outputs"  # UPDATE

# Baseline (current erosivity) run folders
# Naming convention suggested: scenario_climate_c[value]
run_dirs <- list(
  # Baseline climate runs
  base_034        = file.path(runs_parent, "baseline_0.34"),
  base_034_af     = file.path(runs_parent, "baseline_0.34_to_0.2832"),
  base_045        = file.path(runs_parent, "baseline_0.45"),
  base_045_af     = file.path(runs_parent, "baseline_0.45_to_0.3748"),
  base_034_dense  = file.path(runs_parent, "baseline_hedge_0.34_to_0.2"),   # 0.34 -> 0.20
  base_045_dense  = file.path(runs_parent, "baseline_hedge_0.45_to_0.2643"),   # 0.45 -> 0.2643
  
  # SSP3-7.0 erosivity runs (same C structure, future R)
  ssp_034         = file.path(runs_parent, "ssp370_0.34"),
  ssp_034_af      = file.path(runs_parent, "ssp370_0.34_to_0.2832"),
  ssp_045         = file.path(runs_parent, "ssp370_0.45"),
  ssp_045_af      = file.path(runs_parent, "ssp370_0.45_to_0.3748"),
  ssp_034_dense   = file.path(runs_parent, "ssp370_hedge_0.34_to_0.2"),
  ssp_045_dense   = file.path(runs_parent, "ssp370_hedge_0.45_to_0.2643")
)

# ----------------------------------------------------------------
# DIAGNOSTIC — run this block on its own to see what's actually
# in your run folders before touching any other code
# ----------------------------------------------------------------

# 1. Print every path in run_dirs and whether it exists
message("=== run_dirs path check ===")
for (k in names(run_dirs)) {
  folder_exists <- dir.exists(run_dirs[[k]])
  shp_path      <- file.path(run_dirs[[k]], "watershed_results_sdr.shp")
  shp_exists    <- file.exists(shp_path)
  cat(sprintf("%-20s | folder: %-5s | shp: %-5s | %s\n",
              k,
              folder_exists,
              shp_exists,
              run_dirs[[k]]))
}

# 2. List ALL .shp files found anywhere under runs_parent
# This tells you the real filename and subfolder structure InVEST used
message("\n=== All .shp files found under runs_parent ===")
all_shps <- list.files(runs_parent, pattern = "\\.shp$",
                       recursive = TRUE, full.names = TRUE)
if (length(all_shps) == 0) {
  message("No .shp files found. Check that runs_parent is correct: ", runs_parent)
} else {
  print(all_shps)
}

# 3. Also show the top-level subfolders so you can match to run_dirs keys
message("\n=== Subfolders directly inside runs_parent ===")
print(list.dirs(runs_parent, recursive = FALSE))


# Scenario metadata table — used for clean labeling throughout
scenario_meta <- tribble(
  ~run_key,          ~climate,    ~c_base, ~c_af,   ~af_type,         ~label,
  "base_034",        "baseline",  0.34,    NA,       "none",           "Base 0.34",
  "base_034_af",     "baseline",  0.34,    0.2832,   "mixed_AF",       "Base 0.34 + AF (0.2832)",
  "base_045",        "baseline",  0.45,    NA,       "none",           "Base 0.45",
  "base_045_af",     "baseline",  0.45,    0.3748,   "mixed_AF",       "Base 0.45 + AF (0.3748)",
  "base_034_dense",  "baseline",  0.34,    0.20,     "dense_AF",       "Base 0.34 + Dense AF (0.20)",
  "base_045_dense",  "baseline",  0.45,    0.2643,   "dense_AF",       "Base 0.45 + Dense AF (0.2643)",
  "ssp_034",         "ssp370",    0.34,    NA,       "none",           "SSP370 0.34",
  "ssp_034_af",      "ssp370",    0.34,    0.2832,   "mixed_AF",       "SSP370 0.34 + AF (0.2832)",
  "ssp_045",         "ssp370",    0.45,    NA,       "none",           "SSP370 0.45",
  "ssp_045_af",      "ssp370",    0.45,    0.3748,   "mixed_AF",       "SSP370 0.45 + AF (0.3748)",
  "ssp_034_dense",   "ssp370",    0.34,    0.20,     "dense_AF",       "SSP370 0.34 + Dense AF (0.20)",
  "ssp_045_dense",   "ssp370",    0.45,    0.2643,   "dense_AF",       "SSP370 0.45 + Dense AF (0.2643)"
)

# Watershed metrics to extract from watershed_results_sdr.shp
ws_metrics <- c("sed_export", "usle_tot", "sed_dep", "avoid_exp", "avoid_eros")

# Raster layers available in each InVEST output folder
raster_layers <- c("sed_export.tif", "sed_deposition.tif", "rkls.tif",
                   "usle.tif", "avoided_export.tif", "avoided_erosion.tif")

# Output directory
out_dir <- file.path(runs_parent, "sensitivity_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ================================================================
# UTILITY FUNCTIONS
# ================================================================

# Load watershed_results_sdr.shp and return as data frame with ws_id
load_ws <- function(run_key) {
  run_dir <- run_dirs[[run_key]]
  shp <- file.path(run_dir, "watershed_results_sdr.shp")
  if (!file.exists(shp)) {
    warning("Missing: ", shp)
    return(NULL)
  }
  sf_obj <- st_read(shp, quiet = TRUE)
  df <- st_drop_geometry(sf_obj)
  if (!"ws_id" %in% names(df)) df$ws_id <- seq_len(nrow(df))
  df$run_key <- run_key
  df
}

# Compute % change: (scenario - reference) / |reference| * 100
pct_change <- function(scen, ref) {
  ifelse(ref == 0 | is.na(ref), NA_real_, (scen - ref) / abs(ref) * 100)
}

# Arc elasticity: proportional change in output / proportional change in C
arc_elasticity <- function(c_vals, metric_means) {
  if (length(c_vals) < 2) return(NA_real_)
  c1 <- c_vals[1]; c2 <- c_vals[length(c_vals)]
  m1 <- metric_means[1]; m2 <- metric_means[length(metric_means)]
  if (c1 == 0 || m1 == 0) return(NA_real_)
  ((m2 - m1) / ((m1 + m2) / 2)) / ((c2 - c1) / ((c1 + c2) / 2))
}

# Build pairwise comparison: reference run vs. scenario run
# Returns per-watershed delta and % change for all ws_metrics
compare_runs <- function(ref_key, scen_key) {
  df_ref  <- load_ws(ref_key)
  df_scen <- load_ws(scen_key)
  if (is.null(df_ref) || is.null(df_scen)) return(NULL)
  if (nrow(df_ref) != nrow(df_scen)) {
    warning("Row mismatch: ", ref_key, " vs ", scen_key)
    return(NULL)
  }
  
  valid_metrics <- intersect(ws_metrics, intersect(names(df_ref), names(df_scen)))
  if (length(valid_metrics) == 0) {
    warning("No matching metrics between ", ref_key, " and ", scen_key)
    return(NULL)
  }
  
  out <- data.frame(ws_id = df_ref$ws_id,
                    ref_key = ref_key,
                    scen_key = scen_key)
  for (m in valid_metrics) {
    out[[paste0(m, "_ref")]]   <- df_ref[[m]]
    out[[paste0(m, "_scen")]]  <- df_scen[[m]]
    out[[paste0(m, "_delta")]] <- df_scen[[m]] - df_ref[[m]]
    out[[paste0(m, "_pct")]]   <- pct_change(df_scen[[m]], df_ref[[m]])
  }
  out
}

# Summarise a comparison data frame to mean/sd/median % change per metric
summarise_comparison <- function(comp_df, ref_key, scen_key) {
  pct_cols <- grep("_pct$", names(comp_df), value = TRUE)
  map_dfr(pct_cols, function(col) {
    vals <- comp_df[[col]]
    tibble(
      ref_key      = ref_key,
      scen_key     = scen_key,
      ref_label    = scenario_meta$label[scenario_meta$run_key == ref_key],
      scen_label   = scenario_meta$label[scenario_meta$run_key == scen_key],
      metric       = sub("_pct$", "", col),
      mean_pct_chg = mean(vals, na.rm = TRUE),
      median_pct_chg = median(vals, na.rm = TRUE),
      sd_pct_chg   = sd(vals, na.rm = TRUE),
      n_watersheds = sum(!is.na(vals))
    )
  })
}


# ----------------------------------------------------------------
# REPLACE the broken column-check block with this
# ----------------------------------------------------------------

message("=== Column check: inspecting first available run ===")

# Find the first run folder that actually has the shapefile
first_valid_key <- NULL
for (k in names(run_dirs)) {
  shp <- file.path(run_dirs[[k]], "watershed_results_sdr.shp")
  if (file.exists(shp)) {
    first_valid_key <- k
    break
  }
}

if (is.null(first_valid_key)) {
  stop("No watershed_results_sdr.shp found in any run_dirs folder. ",
       "Check your run_dirs paths.")
}

message("Using run '", first_valid_key, "' for column inspection.")
first_df <- load_ws(first_valid_key)

if (!is.null(first_df)) {
  message("Columns found in watershed_results_sdr.shp:")
  print(names(first_df))
  
  # Restrict ws_metrics to only columns that actually exist
  ws_metrics <- intersect(ws_metrics, names(first_df))
  message("ws_metrics after filtering to available columns: ",
          paste(ws_metrics, collapse = ", "))
  
  if (length(ws_metrics) == 0) {
    stop("None of the requested ws_metrics were found. ",
         "Check column names above and update ws_metrics manually.")
  }
} else {
  stop("load_ws() returned NULL for '", first_valid_key, "'. ",
       "Check that the .shp and .dbf files are both present.")
} 




# ================================================================
# BLOCK 1: AF vs BASELINE comparisons (within climate forcing)
# ================================================================
# Pairs:
#   baseline 0.34 vs 0.34+AF(0.2832)
#   baseline 0.45 vs 0.45+AF(0.3748)
#   baseline 0.34 vs 0.34+dense(0.20)
#   baseline 0.45 vs 0.45+dense(0.2643)
#   SSP 0.34 vs SSP 0.34+AF(0.2832)
#   SSP 0.45 vs SSP 0.45+AF(0.3748)
#   SSP 0.34 vs SSP 0.34+dense(0.20)
#   SSP 0.45 vs SSP 0.45+dense(0.2643)
# ================================================================

af_pairs <- list(
  # Baseline climate AF pairs
  list(ref = "base_034",       scen = "base_034_af",    group = "baseline_mixed_AF"),
  list(ref = "base_045",       scen = "base_045_af",    group = "baseline_mixed_AF"),
  list(ref = "base_034",       scen = "base_034_dense", group = "baseline_dense_AF"),
  list(ref = "base_045",       scen = "base_045_dense", group = "baseline_dense_AF"),
  # SSP3-7.0 AF pairs
  list(ref = "ssp_034",        scen = "ssp_034_af",     group = "ssp370_mixed_AF"),
  list(ref = "ssp_045",        scen = "ssp_045_af",     group = "ssp370_mixed_AF"),
  list(ref = "ssp_034",        scen = "ssp_034_dense",  group = "ssp370_dense_AF"),
  list(ref = "ssp_045",        scen = "ssp_045_dense",  group = "ssp370_dense_AF")
)

message("=== BLOCK 1: AF vs Baseline comparisons ===")

block1_long <- map_dfr(af_pairs, function(p) {
  comp <- compare_runs(p$ref, p$scen)
  if (is.null(comp)) return(NULL)
  s <- summarise_comparison(comp, p$ref, p$scen)
  s$group <- p$group
  s
})

write.csv(block1_long,
          file.path(out_dir, "block1_AF_vs_baseline_summary.csv"),
          row.names = FALSE)

# --- Block 1 plot: % change in sed_export by pair and C configuration ---
p1 <- block1_long %>%
  filter(metric == "sed_export") %>%
  mutate(pair_label = paste0(ref_label, "\nvs\n", scen_label)) %>%
  ggplot(aes(x = pair_label, y = mean_pct_chg,
             ymin = mean_pct_chg - sd_pct_chg,
             ymax = mean_pct_chg + sd_pct_chg,
             fill = group)) +
  geom_col(position = "dodge") +
  geom_errorbar(width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Sediment Export Change: AF vs Baseline (Mean ± SD across watersheds)",
       x = NULL, y = "% Change in sed_export",
       fill = "AF group") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(out_dir, "block1_sed_export_AF_vs_baseline.png"),
       p1, width = 14, height = 6, dpi = 150)

# --- Block 1 plot: all metrics, AF vs baseline ---
p1b <- block1_long %>%
  mutate(pair_label = paste0(sub(".*label = '(.*)'.*", "\\1", ref_label),
                             "\nvs ", scen_label)) %>%
  ggplot(aes(x = scen_label, y = mean_pct_chg, fill = group)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = mean_pct_chg - sd_pct_chg,
                    ymax = mean_pct_chg + sd_pct_chg), width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~metric, scales = "free_y") +
  labs(title = "All Metrics: AF vs Baseline comparisons",
       x = NULL, y = "Mean % change", fill = "Group") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave(file.path(out_dir, "block1_all_metrics_AF_vs_baseline.png"),
       p1b, width = 16, height = 10, dpi = 150)

message("Block 1 complete. Files written to: ", out_dir)

# ================================================================
# BLOCK 2: Climate comparisons (baseline vs SSP3-7.0, same C config)
# ================================================================
# Pairs:
#   baseline 0.34 vs SSP 0.34    — erosivity effect on current land
#   baseline 0.45 vs SSP 0.45    — erosivity effect on degraded land
#   baseline 0.34+AF vs SSP 0.34+AF  — erosivity effect under mixed AF
#   baseline 0.45+AF vs SSP 0.45+AF
#   baseline 0.34+dense vs SSP 0.34+dense
#   baseline 0.45+dense vs SSP 0.45+dense
# These isolate the pure climate change signal at each land config.
# ================================================================

climate_pairs <- list(
  list(ref = "base_034",       scen = "ssp_034",       group = "no_AF",    c_base = 0.34),
  list(ref = "base_045",       scen = "ssp_045",       group = "no_AF",    c_base = 0.45),
  list(ref = "base_034_af",    scen = "ssp_034_af",    group = "mixed_AF", c_base = 0.34),
  list(ref = "base_045_af",    scen = "ssp_045_af",    group = "mixed_AF", c_base = 0.45),
  list(ref = "base_034_dense", scen = "ssp_034_dense", group = "dense_AF", c_base = 0.34),
  list(ref = "base_045_dense", scen = "ssp_045_dense", group = "dense_AF", c_base = 0.45)
)

message("=== BLOCK 2: Climate (baseline vs SSP3-7.0) comparisons ===")

block2_long <- map_dfr(climate_pairs, function(p) {
  comp <- compare_runs(p$ref, p$scen)
  if (is.null(comp)) return(NULL)
  s <- summarise_comparison(comp, p$ref, p$scen)
  s$group    <- p$group
  s$c_base   <- p$c_base
  s
})

write.csv(block2_long,
          file.path(out_dir, "block2_climate_effect_summary.csv"),
          row.names = FALSE)

# --- Block 2 plot: climate change signal by AF configuration ---
p2 <- block2_long %>%
  filter(metric == "sed_export") %>%
  ggplot(aes(x = group, y = mean_pct_chg,
             fill = factor(c_base),
             ymin = mean_pct_chg - sd_pct_chg,
             ymax = mean_pct_chg + sd_pct_chg)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(position = position_dodge(0.8), width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Climate Change Effect on Sediment Export (Baseline → SSP3-7.0)",
       subtitle = "Isolates erosivity change at constant land cover / C-factor",
       x = "AF configuration", y = "Mean % change in sed_export",
       fill = "C base") +
  theme_bw(base_size = 11)

ggsave(file.path(out_dir, "block2_climate_effect_sed_export.png"),
       p2, width = 9, height = 6, dpi = 150)

# --- Block 2 plot: all metrics ---
p2b <- block2_long %>%
  ggplot(aes(x = group, y = mean_pct_chg, fill = factor(c_base))) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_pct_chg - sd_pct_chg,
                    ymax = mean_pct_chg + sd_pct_chg),
                position = position_dodge(0.8), width = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~metric, scales = "free_y") +
  labs(title = "Climate Effect on All Metrics (Baseline → SSP3-7.0)",
       x = "AF configuration", y = "Mean % change", fill = "C base") +
  theme_bw(base_size = 9)

ggsave(file.path(out_dir, "block2_all_metrics_climate_effect.png"),
       p2b, width = 14, height = 9, dpi = 150)

message("Block 2 complete.")

# ================================================================
# BLOCK 3: Cross comparisons (climate × AF sensitivity analysis)
# ================================================================
# This is your full sensitivity decomposition. It compares:
#   A. baseline 0.34  vs SSP 0.34+AF(0.2832)  — combined effect
#   B. baseline 0.34  vs SSP 0.45+AF(0.3748)  — worst-case mitigation
#   C. baseline 0.34  vs SSP 0.34+dense(0.20) — climate + dense AF
#   D. baseline 0.34  vs SSP 0.45+dense(0.2643)
#   E. baseline 0.45  vs SSP 0.34 (degraded lands improved by AF)
#   ...and the arc elasticity of outputs to the C-factor gradient,
#      computed separately within baseline climate and SSP climate.
# ================================================================

cross_pairs <- list(
  # Combined climate + mixed AF (key policy-relevant comparisons)
  list(ref = "base_034", scen = "ssp_034_af",    label = "Base_034 → SSP+AF(0.2832)"),
  list(ref = "base_034", scen = "ssp_045_af",    label = "Base_034 → SSP_045+AF(0.3748)"),
  list(ref = "base_034", scen = "ssp_034_dense", label = "Base_034 → SSP+Dense(0.20)"),
  list(ref = "base_034", scen = "ssp_045_dense", label = "Base_034 → SSP_045+Dense(0.2643)"),
  list(ref = "base_045", scen = "ssp_034_af",    label = "Base_045 → SSP_034+AF(0.2832)"),
  list(ref = "base_045", scen = "ssp_045_af",    label = "Base_045 → SSP_045+AF(0.3748)"),
  list(ref = "base_045", scen = "ssp_034_dense", label = "Base_045 → SSP_034+Dense(0.20)"),
  list(ref = "base_045", scen = "ssp_045_dense", label = "Base_045 → SSP_045+Dense(0.2643)")
)

message("=== BLOCK 3: Cross comparisons (climate x AF) ===")

block3_long <- map_dfr(cross_pairs, function(p) {
  comp <- compare_runs(p$ref, p$scen)
  if (is.null(comp)) return(NULL)
  s <- summarise_comparison(comp, p$ref, p$scen)
  s$cross_label <- p$label
  s
})

write.csv(block3_long,
          file.path(out_dir, "block3_cross_climate_AF_summary.csv"),
          row.names = FALSE)

# --- Block 3 plot: combined climate + AF effect ---
p3 <- block3_long %>%
  filter(metric == "sed_export") %>%
  ggplot(aes(x = reorder(cross_label, mean_pct_chg), y = mean_pct_chg,
             fill = mean_pct_chg < 0)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_pct_chg - sd_pct_chg,
                    ymax = mean_pct_chg + sd_pct_chg), width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato"),
                    labels = c("TRUE" = "Reduction", "FALSE" = "Increase"),
                    name = NULL) +
  coord_flip() +
  labs(title = "Sediment Export: Combined Climate × AF Effect",
       subtitle = "Negative = AF offsets or reverses climate-driven increase",
       x = NULL, y = "Mean % change in sed_export") +
  theme_bw(base_size = 10)

ggsave(file.path(out_dir, "block3_cross_sed_export.png"),
       p3, width = 11, height = 7, dpi = 150)

# ================================================================
# BLOCK 3b: Arc Elasticity of SDR outputs to C-factor
# Computed separately for baseline and SSP3-7.0 climates
# ================================================================

message("=== BLOCK 3b: Arc elasticity across C-factor gradient ===")

# For each climate, order the C-factor runs from highest to lowest
# (high C = less AF cover = more erosion, so decreasing C = adding AF)

c_gradient <- list(
  baseline = list(
    c_vals    = c(0.45,    0.3748,  0.34,    0.2832,  0.2643,  0.20),
    run_keys  = c("base_045", "base_045_af", "base_034",
                  "base_034_af", "base_045_dense", "base_034_dense")
  ),
  ssp370 = list(
    c_vals    = c(0.45,    0.3748,  0.34,    0.2832,  0.2643,  0.20),
    run_keys  = c("ssp_045", "ssp_045_af", "ssp_034",
                  "ssp_034_af", "ssp_045_dense", "ssp_034_dense")
  )
)

elasticity_rows <- list()

for (climate_name in names(c_gradient)) {
  cg      <- c_gradient[[climate_name]]
  c_vals  <- cg$c_vals
  run_keys <- cg$run_keys
  
  # Load all watershed dfs for this climate gradient
  dfs <- map(run_keys, load_ws)
  valid <- !map_lgl(dfs, is.null)
  
  if (sum(valid) < 2) {
    warning("Not enough valid runs for elasticity in ", climate_name)
    next
  }
  
  c_vals_valid  <- c_vals[valid]
  dfs_valid     <- dfs[valid]
  
  # For each metric, compute mean across watersheds at each C level
  valid_metrics <- intersect(ws_metrics, names(dfs_valid[[1]]))
  
  for (m in valid_metrics) {
    means_at_c <- map_dbl(dfs_valid, ~ mean(.x[[m]], na.rm = TRUE))
    
    elast <- arc_elasticity(c_vals_valid, means_at_c)
    interp <- case_when(
      is.na(elast)      ~ "insufficient data",
      abs(elast) > 1    ~ "Elastic (output sensitive to C)",
      abs(elast) == 1   ~ "Unit elastic",
      abs(elast) < 1    ~ "Inelastic (output buffered from C)"
    )
    elasticity_rows[[paste(climate_name, m, sep = "_")]] <- tibble(
      climate        = climate_name,
      metric         = m,
      arc_elasticity = round(elast, 4),
      interpretation = interp
    )
  }
}

elasticity_df <- bind_rows(elasticity_rows)
write.csv(elasticity_df,
          file.path(out_dir, "block3b_arc_elasticity.csv"),
          row.names = FALSE)

# --- Elasticity plot: baseline vs SSP side by side ---
p_elast <- elasticity_df %>%
  filter(!is.na(arc_elasticity)) %>%
  ggplot(aes(x = reorder(metric, abs(arc_elasticity)),
             y = arc_elasticity,
             fill = climate)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 0, colour = "black") +
  coord_flip() +
  scale_fill_manual(values = c("baseline" = "steelblue", "ssp370" = "tomato")) +
  labs(title = "Arc Elasticity of SDR Metrics to C-Factor",
       subtitle = "Dashed lines at ±1 (unit elasticity threshold)",
       x = "Metric", y = "Arc elasticity (ΔOutput/ΔC)",
       fill = "Climate") +
  theme_bw(base_size = 11)

ggsave(file.path(out_dir, "block3b_elasticity_baseline_vs_ssp.png"),
       p_elast, width = 8, height = 5, dpi = 150)

# ================================================================
# BLOCK 3c: Rank stability — do watersheds rank consistently
# as AF intensity increases, within each climate?
# ================================================================

message("=== BLOCK 3c: Rank stability across C-factor gradient ===")

rank_rows <- list()

for (climate_name in names(c_gradient)) {
  cg       <- c_gradient[[climate_name]]
  c_vals   <- cg$c_vals
  run_keys <- cg$run_keys
  dfs      <- map(run_keys, load_ws)
  valid    <- !map_lgl(dfs, is.null)
  
  if (sum(valid) < 2) next
  
  c_vals_valid <- c_vals[valid]
  dfs_valid    <- dfs[valid]
  
  valid_metrics <- intersect(ws_metrics, names(dfs_valid[[1]]))
  
  for (m in valid_metrics) {
    val_mat <- do.call(cbind, map(dfs_valid, ~ .x[[m]]))
    n_scen  <- ncol(val_mat)
    
    for (j in seq_len(n_scen - 1)) {
      for (k in (j + 1):n_scen) {
        rho <- cor(val_mat[, j], val_mat[, k],
                   method = "spearman", use = "complete.obs")
        rank_rows[[paste(climate_name, m, j, k, sep = "_")]] <- tibble(
          climate     = climate_name,
          metric      = m,
          c_from      = c_vals_valid[j],
          c_to        = c_vals_valid[k],
          spearman_rho = round(rho, 4),
          stable      = rho >= 0.9
        )
      }
    }
  }
}

rank_df <- bind_rows(rank_rows)
write.csv(rank_df,
          file.path(out_dir, "block3c_rank_stability.csv"),
          row.names = FALSE)

# --- Rank stability heatmap: sed_export, both climates ---
for (clim in c("baseline", "ssp370")) {
  sub <- rank_df %>% filter(climate == clim, metric == "sed_export")
  if (nrow(sub) == 0) next
  
  # Make symmetric
  sym <- bind_rows(
    sub,
    sub %>% rename(c_from = c_to, c_to = c_from)
  ) %>%
    bind_rows(tibble(climate = clim, metric = "sed_export",
                     c_from = unique(c(sub$c_from, sub$c_to)),
                     c_to   = unique(c(sub$c_from, sub$c_to)),
                     spearman_rho = 1, stable = TRUE))
  
  pheat <- ggplot(sym, aes(x = factor(c_from), y = factor(c_to),
                           fill = spearman_rho)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = round(spearman_rho, 2)), size = 3.5) +
    scale_fill_gradient2(low = "tomato", mid = "lightyellow", high = "steelblue",
                         midpoint = 0.9, limits = c(0, 1), name = "Spearman ρ") +
    labs(title = paste0("Rank Stability: sed_export — ", clim),
         subtitle = "ρ ≥ 0.9 = stable watershed rankings across C-factor scenarios",
         x = "C-factor", y = "C-factor") +
    theme_bw(base_size = 11)
  
  ggsave(file.path(out_dir, paste0("block3c_rankheatmap_sedexport_", clim, ".png")),
         pheat, width = 7, height = 6, dpi = 150)
}

message("Block 3c complete.")

# ================================================================
# SUMMARY TABLE: All comparisons, all metrics
# ================================================================

message("=== Writing master summary table ===")

master_summary <- bind_rows(
  block1_long %>% mutate(block = "Block1_AF_vs_Baseline"),
  block2_long %>% mutate(block = "Block2_Climate_Effect"),
  block3_long %>% mutate(block = "Block3_Cross_Climate_AF")
) %>%
  select(block, ref_key, scen_key, ref_label, scen_label, metric,
         mean_pct_chg, median_pct_chg, sd_pct_chg, n_watersheds,
         everything())

write.csv(master_summary,
          file.path(out_dir, "MASTER_sensitivity_summary.csv"),
          row.names = FALSE)

# Print console summary
message("\n========================================")
message("SENSITIVITY ANALYSIS COMPLETE")
message("========================================")
message("Outputs written to: ", out_dir)
message("")
message("Files produced:")
message("  Block 1: AF vs Baseline (same climate)")
message("    block1_AF_vs_baseline_summary.csv")
message("    block1_sed_export_AF_vs_baseline.png")
message("    block1_all_metrics_AF_vs_baseline.png")
message("")
message("  Block 2: Climate effect (baseline vs SSP3-7.0, same C)")
message("    block2_climate_effect_summary.csv")
message("    block2_climate_effect_sed_export.png")
message("    block2_all_metrics_climate_effect.png")
message("")
message("  Block 3: Cross climate x AF comparisons")
message("    block3_cross_climate_AF_summary.csv")
message("    block3_cross_sed_export.png")
message("")
message("  Block 3b: Arc elasticity")
message("    block3b_arc_elasticity.csv")
message("    block3b_elasticity_baseline_vs_ssp.png")
message("")
message("  Block 3c: Rank stability")
message("    block3c_rank_stability.csv")
message("    block3c_rankheatmap_sedexport_baseline.png")
message("    block3c_rankheatmap_sedexport_ssp370.png")
message("")
message("  Master summary: MASTER_sensitivity_summary.csv")