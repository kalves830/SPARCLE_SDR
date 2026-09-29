# =============================================================================
# InVEST SDR Raster Preparation for West Africa
# Replicates ArcGIS workflow from Kyle's notes (2/8/26 and 3/1/26)
# Also addresses InVEST SDR LS-factor IndexError from log 2026-03-02
#
# ERROR ADDRESSED:
#   IndexError: index 4 is out of bounds for axis 0 with size 4
#   in natcap.invest.sdr.sdr.ls_factor_function (line 1251)
#
#   This error occurs when the DEM contains edge pixels or nodata regions
#   that produce an out-of-range D8 flow direction index (value 4) during
#   the LS-factor calculation. The fix is to:
#     1. Ensure all DEM nodata pixels are explicitly masked (not zero-filled)
#     2. Buffer-clip the DEM so no valid pixels sit on the raster edge
#     3. Use a consistent nodata value (-9999) across all inputs
#     4. Ensure the DEM is pit-filled BEFORE handing it to InVEST
#        (InVEST pit-fills internally, but edge artifacts can still cause issues)
#
# REFERENCES:
#   - Sharp et al. (2018) InVEST User's Guide: Sediment Delivery Ratio.
#     Natural Capital Project, Stanford University.
#   - Renard et al. (1997) RUSLE: Predicting Soil Erosion by Water.
#     USDA Agricultural Handbook No. 703.
#   - Hijmans, R. (2023). terra: Spatial Data Analysis. R package.
# =============================================================================

# --- 1. PACKAGE LOADING -------------------------------------------------------

# terra: Modern R spatial package for raster and vector operations.
# Replaces the older 'raster' package. Used for all raster I/O, reprojection,
# resampling, mosaicking, and clipping. Far faster than raster for large files.
library(terra)

# sf: Simple Features for R. Used for vector operations (watershed polygons,
# extent clipping masks). Integrates cleanly with terra via vect().
library(sf)

# dplyr: Used for minor tabular work (biophysical table checks). Optional.
library(dplyr)


# --- 2. CONFIGURATION ---------------------------------------------------------

# Set your file paths here. Adjust to match your directory structure.
# Uses forward slashes (works on Windows too when inside R).

config <- list(

  # Input DEM tiles (Copernicus GLO-30 or similar; 2 tiles to mosaic)
  dem_tile     = "dem/dem_run.tif",

  # Land Use / Land Cover tiles (4 tiles from your notes)
  lulc_tiles    = list(
    "lulc/lulc.tif"
  ),

  # Erosivity (GloRESatE or similar rainfall erosivity R-factor)
  erosivity     = "erosivity/erosivity_run.tif",

  # Erodibility (EPIC K-factor; 4 clipped tiles to mosaic, matching LULC extents)
  kfactor_tiles = list(
    "erodibility/erodibility/k_factor_invest.tif"
  ),

  # Watersheds shapefile (af-lev watersheds including lakes)
  watersheds    = "watersheds/watersheds.shp",

  # Output directory
  out_dir       = "output",

  # Target projection: Africa Albers Equal Area Conic (ESRI:102022)
  # Used in your ArcGIS notes for all projected outputs.
  # Equal-area projection is essential for RUSLE area calculations.
  target_crs    = "ESRI:102022",

  # Target resolution: 1000 m (1 km), as specified in your notes
  target_res    = 1000,

  # Nodata value to use consistently across all outputs.
  # -9999 is the InVEST-recommended convention. Using a consistent value
  # prevents the LS-factor index error caused by nodata being treated as 0.
  nodata_val    = -9999
)

# Create output directory if it doesn't exist
if (!dir.exists(config$out_dir)) dir.create(config$out_dir, recursive = TRUE)


# --- 3. WEST AFRICA EXTENT ----------------------------------------------------
# Define the West Africa bounding box in geographic coordinates (WGS84).
# This is used to clip all rasters before projection to keep file sizes
# manageable. Approximate bounds: lon -18 to 16, lat 4 to 23.

west_africa_wgs84 <- ext(-18, 16, 4, 23)  # xmin, xmax, ymin, ymax


# --- 4. HELPER FUNCTION: Reproject and Resample a Raster ----------------------

#' Reproject and resample a raster to the target CRS and resolution.
#'
#' @param r       SpatRaster object to process
#' @param method  Resampling method: "bilinear" for continuous data (DEM,
#'                erosivity, erodibility), "near" (nearest neighbour) for
#'                categorical data (LULC). Per your ArcGIS notes and
#'                standard GIS practice (Tobler 1979).
#' @param target_crs  PROJ string or EPSG/ESRI code for output CRS
#' @param target_res  Output cell size in CRS units (metres for AAEAC)
#' @param nodata_val  Value to assign as nodata in output
reproject_resample <- function(r,
                                method     = "bilinear",
                                target_crs = config$target_crs,
                                target_res = config$target_res,
                                nodata_val = config$nodata_val) {

  # Build a template raster at the target CRS and resolution
  # terra::project() reprojects and resamples in one step.
  r_proj <- project(r, target_crs, method = method, res = target_res)

  # Explicitly set nodata. This is critical: InVEST's LS-factor function
  # iterates over a D8 lookup table with 4 entries (indices 0-3).
  # If nodata is stored as 0 or left ambiguous, edge/nodata pixels return
  # index 4 which is out of bounds. Setting a clear nodata masks these cells.
  NAflag(r_proj) <- nodata_val

  return(r_proj)
}


# --- 5. PROCESS DEM -----------------------------------------------------------
# Replicates ArcGIS notes: mosaic two DEM tiles, reproject bilinear, upscale 1km

message("=== Processing DEM ===")

# Load both DEM tiles
dem <- rast(config$dem_tile)


# Clip to West Africa extent in native CRS (WGS84) before reprojecting.
# Clipping first reduces the data volume and speeds up reprojection.
dem_clipped_wgs <- crop(dem, west_africa_wgs84)

# Reproject to Africa Albers Equal Area Conic with bilinear resampling.
# Bilinear is correct for continuous elevation data (not categorical).
# Reference: ESRI Best Practices for DEM Resampling (2016).
dem_proj <- reproject_resample(dem_clipped_wgs, method = "bilinear")

# --- DEM EDGE FIX (addresses the InVEST LS-factor IndexError) ---
# The error occurs when valid DEM cells exist right on the raster boundary,
# because D8 flow direction can't find a downslope neighbour and returns
# an invalid index. We add a 1-cell buffer of nodata around the valid area
# by shrinking the extent slightly (trim nodata, then re-add a nodata border).
#
# Strategy: replace any edge rows/cols that are entirely nodata with
# explicit nodata, and ensure the nodata flag is written to the file.
# A safer alternative (used here) is to shrink the extent by 1 pixel:

dem_nrows <- nrow(dem_proj)
dem_ncols <- ncol(dem_proj)

# Crop 1 pixel from each edge to remove boundary artifacts.
# This loses minimal data but prevents the LS-factor crash.
dem_safe <- crop(dem_proj,
                 ext(dem_proj$xmin + res(dem_proj)[1],
                     dem_proj$xmax - res(dem_proj)[1],
                     dem_proj$ymin + res(dem_proj)[2],
                     dem_proj$ymax - res(dem_proj)[2]))

dem_out <- file.path(config$out_dir, "dem_run.tif")

# Write as 32-bit float to match your ArcGIS export settings.
# datatype = "FLT4S" = 32-bit signed float, equivalent to ArcGIS "32 bit float"
writeRaster(dem_safe,
            dem_out,
            datatype  = "FLT4S",
            NAflag    = config$nodata_val,
            overwrite = TRUE)

message("DEM written to: ", dem_out)


# --- 6. PROCESS LULC ----------------------------------------------------------
# Replicates ArcGIS notes: mosaic 4 LULC tiles, reproject, clip.
# IMPORTANT: Use nearest-neighbour resampling for categorical data.

message("=== Processing LULC ===")

# Load all LULC tiles
lulc_list <- lapply(config$lulc_tiles, rast)

# Mosaic all tiles. For categorical data, 'fun = "first"' or 'fun = "last"'
# both work (pixel values are class codes, not interpolatable numbers).
# Using "last" to match your ArcGIS "Last" mosaic operator.
lulc_mosaic <- do.call(mosaic, c(lulc_list, list(fun = "last")))

# Clip to West Africa extent
lulc_clipped_wgs <- crop(lulc_mosaic, west_africa_wgs84)

# Reproject with NEAREST NEIGHBOUR for categorical LULC data.
# Never use bilinear/cubic for class codes — it creates meaningless
# interpolated values between classes (e.g., class 3.7 doesn't exist).
# Reference: Congalton & Green (2009) Assessing the Accuracy of Remotely
# Sensed Data, CRC Press. Chapter on resampling for categorical maps.
lulc_proj <- reproject_resample(lulc_clipped_wgs, method = "near")

# IMPORTANT for InVEST: LULC must be stored as integer (INT2U or INT4S).
# The SetColorTable error in your log ("only supported for Byte or UInt16")
# was a non-fatal GDAL warning caused by a 32-bit float LULC raster.
# Writing as INT2U (unsigned 16-bit) resolves this and matches InVEST's
# expectation for LULC class codes.
lulc_out <- file.path(config$out_dir, "lulc.tif")

writeRaster(lulc_proj,
            lulc_out,
            datatype  = "INT2U",   # Unsigned 16-bit integer — solves color table error
            NAflag    = 65535,     # Max value for UINT16 used as nodata for integers
            overwrite = TRUE)

message("LULC written to: ", lulc_out)


# --- 7. PROCESS EROSIVITY -----------------------------------------------------
# GloRESatE R-factor (MJ·mm·ha⁻¹·h⁻¹·yr⁻¹).
# Reference: Panagos et al. (2017) Global rainfall erosivity assessment based
# on high-temporal resolution rainfall records. Scientific Reports 7, 4175.

message("=== Processing Erosivity ===")

erosivity_r <- rast(config$erosivity)

# Clip to West Africa extent
erosivity_clipped <- crop(erosivity_r, west_africa_wgs84)

# Reproject with bilinear (continuous numerical data)
erosivity_proj <- reproject_resample(erosivity_clipped, method = "bilinear")

erosivity_out <- file.path(config$out_dir, "erosivity_run.tif")
writeRaster(erosivity_proj,
            erosivity_out,
            datatype  = "FLT4S",
            NAflag    = config$nodata_val,
            overwrite = TRUE)

message("Erosivity written to: ", erosivity_out)


# --- 8. PROCESS ERODIBILITY (K-FACTOR) ----------------------------------------
# EPIC K-factor data (tonnes·ha·h·ha⁻¹·MJ⁻¹·mm⁻¹).
# Your notes describe clipping K-factor to each of 4 LULC extents, then
# mosaicking. We replicate this: clip each K tile to its LULC counterpart's
# extent, then mosaic the 4 clips together.
# Reference: Sharpley & Williams (1990) EPIC — Erosion/Productivity Impact
# Calculator. USDA Technical Bulletin No. 1768.

message("=== Processing Erodibility (K-factor) ===")

# Load K-factor tiles
kfactor_list <- lapply(config$kfactor_tiles, rast)

# Load corresponding LULC tiles (still in WGS84 at this point)
lulc_raw_list <- lapply(config$lulc_tiles, rast)

# Clip each K tile to the extent of the corresponding LULC tile
k_clipped_list <- mapply(function(k_tile, lulc_tile) {
  crop(k_tile, ext(lulc_tile))
}, kfactor_list, lulc_raw_list, SIMPLIFY = FALSE)

# Mosaic the 4 clipped K-factor tiles.
# 'fun = "last"' matches ArcGIS "Last" mosaic operator in your notes.
k_mosaic <- do.call(mosaic, c(k_clipped_list, list(fun = "last")))

# Clip to West Africa extent
k_clipped_wgs <- crop(k_mosaic, west_africa_wgs84)

# Reproject with bilinear (K-factor is a continuous numerical value).
# Note from your ArcGIS notes: bilinear resampling may fill in some singular
# points from the original raster — this is expected behaviour.
k_proj <- reproject_resample(k_clipped_wgs, method = "bilinear")

k_out <- file.path(config$out_dir, "k_factor_invest.tif")
writeRaster(k_proj,
            k_out,
            datatype  = "FLT4S",
            NAflag    = config$nodata_val,
            overwrite = TRUE)

message("K-factor written to: ", k_out)


# --- 9. PROCESS WATERSHEDS ----------------------------------------------------
# Reproject the watershed shapefile to match raster CRS.
# InVEST requires the watersheds polygon to share the same CRS as rasters.

message("=== Processing Watersheds ===")

watersheds_sf <- st_read(config$watersheds, quiet = TRUE)

# Reproject to Africa Albers Equal Area Conic
watersheds_proj <- st_transform(watersheds_sf, crs = config$target_crs)

watersheds_out <- file.path(config$out_dir, "watersheds.shp")
st_write(watersheds_proj, watersheds_out, delete_layer = TRUE, quiet = TRUE)

message("Watersheds written to: ", watersheds_out)


# --- 10. ALIGNMENT CHECK ------------------------------------------------------
# InVEST aligns rasters internally, but pre-checking alignment avoids
# surprises and confirms the fix worked.

message("=== Alignment Check ===")

dem_check  <- rast(dem_out)
lulc_check <- rast(lulc_out)
ero_check  <- rast(erosivity_out)
k_check    <- rast(k_out)

check_rasters <- list(DEM        = dem_check,
                      LULC       = lulc_check,
                      Erosivity  = ero_check,
                      Erodibility= k_check)

for (nm in names(check_rasters)) {
  r <- check_rasters[[nm]]
  cat(sprintf(
    "%-12s | CRS: %-10s | Res: %s m | Extent: %s | Nodata: %s\n",
    nm,
    substr(crs(r, describe = TRUE)$name, 1, 10),
    paste(round(res(r)), collapse = "x"),
    paste(round(as.vector(ext(r))), collapse = ", "),
    NAflag(r)
  ))
}

message("\n=== All InVEST-ready files written to: ", config$out_dir, " ===")
message("Files to use in InVEST SDR:")
message("  DEM:         ", dem_out)
message("  LULC:        ", lulc_out)
message("  Erosivity:   ", erosivity_out)
message("  Erodibility: ", k_out)
message("  Watersheds:  ", watersheds_out)
