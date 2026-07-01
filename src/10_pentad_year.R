# =========================================================================
# PENTAD (5-DAY) MULTI-PANEL SUMMARY: SINGLE-YEAR ABSOLUTE + ANOMALY GRIDS
# =========================================================================
# Per-year companion to 06_plot_climo.R's 5x5 climatology grid. For a
# chosen year, generates two 5x5 panels covering the same 25 five-day
# blocks (June 10 - Oct 12):
#   1. Absolute 5-day means: TCWV fill, Z500/Q700 contours, one "H" ridge
#      marker per block - same visual language as the climatology grid,
#      for direct side-by-side comparison (e.g. 2020 vs. 2021).
#   2. True 5-day anomalies vs. the 1991-2025 pixel baseline: TCWV anomaly
#      fill, Z500 anomaly contours - shows exactly when/where the season
#      ran wet, dry, or displaced.
#
# Assumes working directory = project root (monsoonTracker/), consistent
# with scripts 01-09.
# =========================================================================

library(terra)
library(ggplot2)
library(tidyterra)
library(maps)

# -------------------------------------------------------------------------
# 0. HELPER: 25 x 5-DAY BLOCK DEFINITIONS (June 10 - Oct 12)
# -------------------------------------------------------------------------
build_pentad_blocks <- function(start_day = 10) {
  blocks <- lapply(0:24, function(i) (start_day + i * 5):(start_day + i * 5 + 4))
  base_date <- as.Date("2021-06-01")  # non-leap reference year, for labels only
  start_dates <- base_date + sapply(blocks, min) - 1
  end_dates   <- base_date + sapply(blocks, max) - 1
  list(
    blocks = blocks,
    labels = sprintf("%s - %s", format(start_dates, "%d %b"), format(end_dates, "%d %b"))
  )
}

# -------------------------------------------------------------------------
# 1. CORE FUNCTION
# -------------------------------------------------------------------------
generate_pentad_year <- function(target_year,
                                  slice_dir = "data/processed/4slices",
                                  climo_dir = "data/processed/4slices/climo",
                                  out_dir   = "NAM_tracker/output/pentad") {

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  months_to_process <- c("06", "07", "08", "09", "10")
  pentad      <- build_pentad_blocks()
  blocks_5day <- pentad$blocks
  layer_names <- pentad$labels

  cat(sprintf("Building pentad grids for %d...\n", target_year))

  # -----------------------------------------------------------------------
  # 1a. LOAD THE FULL 153-DAY SEASON: TCWV, Z500, Q700
  # -----------------------------------------------------------------------
  tcwv_list <- list(); z500_list <- list(); q700_list <- list()

  for (m in months_to_process) {
    nc_single <- list.files(slice_dir, pattern = sprintf("single_%d_%s\\.nc$", target_year, m), full.names = TRUE)[1]
    nc_press  <- list.files(slice_dir, pattern = sprintf("pressure_%d_%s\\.nc$", target_year, m), full.names = TRUE)[1]
    if (is.na(nc_single) || is.na(nc_press)) {
      stop(sprintf("Missing processed data for %d-%s in %s", target_year, m, slice_dir))
    }

    r_single <- rast(nc_single)
    r_press  <- rast(nc_press)
    days_in_month <- nlyr(r_press) / 8   # 8 stacked variables per month block (q850..z250)

    tcwv_list[[m]] <- r_single
    q700_list[[m]] <- r_press[[(1:days_in_month) + (1 * days_in_month)]] * 1000     # block 2 = q700
    z500_list[[m]] <- r_press[[(1:days_in_month) + (6 * days_in_month)]] / 9.80665  # block 7 = z500
  }

  season_tcwv <- rast(tcwv_list)
  season_q700 <- rast(q700_list)
  season_z500 <- rast(z500_list)

  # -----------------------------------------------------------------------
  # 1b. LOAD THE 1991-2025 PIXEL BASELINE (RAW DAILY CLIMATOLOGY)
  # -----------------------------------------------------------------------
  climo_tcwv_path <- file.path(climo_dir, "climatology_tcwv_raw.tif")
  climo_z500_path <- file.path(climo_dir, "climatology_z500_raw.tif")
  if (!file.exists(climo_tcwv_path) || !file.exists(climo_z500_path)) {
    stop("Missing raw climatology baselines. Run 03A/03B first.")
  }

  climo_tcwv <- rast(climo_tcwv_path)
  climo_z500 <- rast(climo_z500_path)

  # -----------------------------------------------------------------------
  # 2. COLLAPSE TO 25 FIVE-DAY BLOCK MEANS (ACTUAL YEAR + CLIMATOLOGY)
  # -----------------------------------------------------------------------
  tcwv_actual <- rast(lapply(blocks_5day, function(idx) mean(season_tcwv[[idx]], na.rm = TRUE)))
  z500_actual <- rast(lapply(blocks_5day, function(idx) mean(season_z500[[idx]], na.rm = TRUE)))
  q700_actual <- rast(lapply(blocks_5day, function(idx) mean(season_q700[[idx]], na.rm = TRUE)))

  tcwv_climo_block <- rast(lapply(blocks_5day, function(idx) mean(climo_tcwv[[idx]], na.rm = TRUE)))
  z500_climo_block <- rast(lapply(blocks_5day, function(idx) mean(climo_z500[[idx]], na.rm = TRUE)))

  names(tcwv_actual) <- names(z500_actual) <- names(q700_actual) <- layer_names
  names(tcwv_climo_block) <- names(z500_climo_block) <- layer_names

  tcwv_anom <- tcwv_actual - tcwv_climo_block
  z500_anom <- z500_actual - z500_climo_block
  names(tcwv_anom) <- names(z500_anom) <- layer_names

  # Spatial smoothing for contouring (matches 06_plot_climo.R)
  spatial_smooth_mat <- matrix(1 / 25, nrow = 5, ncol = 5)
  q700_smooth             <- focal(q700_actual, w = spatial_smooth_mat, fun = mean, na.rm = TRUE)
  tcwv_smooth_for_contour <- focal(tcwv_actual, w = spatial_smooth_mat, fun = mean, na.rm = TRUE)

  # -----------------------------------------------------------------------
  # 3. TRACK THE "H" RIDGE PER 5-DAY BLOCK (ABSOLUTE GRID ONLY)
  # -----------------------------------------------------------------------
  tracker_smooth_mat <- matrix(1 / 225, nrow = 15, ncol = 15)
  visual_smooth_mat  <- matrix(1 / 25, nrow = 5, ncol = 5)
  monsoon_box <- ext(-114, -95, 15, 45)

  h_lon <- numeric(nlyr(z500_actual)); h_lat <- numeric(nlyr(z500_actual))

  for (i in 1:nlyr(z500_actual)) {
    z_search <- crop(z500_actual[[i]], monsoon_box)
    z_macro <- focal(z_search, w = tracker_smooth_mat, fun = mean, na.rm = TRUE)
    macro_max <- global(z_macro, "max", na.rm = TRUE)[1, 1]
    macro_core <- ifel(z_macro >= (macro_max - 10), 1, NA)

    blobs <- patches(macro_core); blob_sizes <- freq(blobs)
    if (nrow(blob_sizes) > 0) {
      largest_blob_id <- blob_sizes$value[which.max(blob_sizes$count)]
      primary_macro_core <- ifel(blobs == largest_blob_id, 1, NA)
      z_visual <- focal(z_search, w = visual_smooth_mat, fun = mean, na.rm = TRUE)
      z_isolated <- mask(z_visual, primary_macro_core)
      true_max <- global(z_isolated, "max", na.rm = TRUE)[1, 1]
      final_core <- ifel(z_isolated >= (true_max - 2), 1, NA)
      centroid <- colMeans(crds(as.points(final_core)))
      h_lon[i] <- centroid[1]; h_lat[i] <- centroid[2]
    } else {
      h_lon[i] <- NA; h_lat[i] <- NA
    }
  }

  h_pts <- data.frame(lon = h_lon, lat = h_lat, lyr = factor(layer_names, levels = layer_names))

  # -----------------------------------------------------------------------
  # 4. PLOT 1: ABSOLUTE 5x5 GRID
  # -----------------------------------------------------------------------
  feature_colors <- c(
    "Primary Ridge (H)" = "darkred",
    "TCWV > 25 kg/m2"   = "orange",
    "Q700 > 6 g/kg"     = "darkblue",
    "Z500 Heights"      = "gray30"
  )
  obs_pal <- colorRampPalette(c(
    "#4a2c11", "#bc966c", "#f4efe6", "#31a354", "#008080", "#064cb5", "#5106b5"
  ))(100)

  world_map <- map_data("world"); state_map <- map_data("state")
  b <- ext(tcwv_actual)

  p_absolute <- ggplot() +
    geom_spatraster(data = tcwv_actual) +
    scale_fill_gradientn(colors = obs_pal, name = "Mean TCWV\n(kg/m2)", na.value = "transparent",
                          limits = c(5, 65), oob = scales::squish) +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.5) +
    geom_polygon(data = state_map, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.3) +
    geom_spatraster_contour(data = z500_actual, aes(color = "Z500 Heights"), breaks = seq(5500, 6000, by = 20), linewidth = 0.3) +
    geom_spatraster_contour(data = q700_smooth, aes(color = "Q700 > 6 g/kg"), breaks = 6, linewidth = 0.6, linetype = "dashed") +
    geom_spatraster_contour(data = tcwv_smooth_for_contour, aes(color = "TCWV > 25 kg/m2"), breaks = 25, linewidth = 0.8) +
    geom_text(data = h_pts, aes(x = lon, y = lat, label = "H"), color = "blue1", fontface = "bold", size = 5) +
    scale_color_manual(name = "Features", values = feature_colors) +
    coord_sf(xlim = c(b[1], b[2]), ylim = c(b[3], b[4]), expand = FALSE) +
    facet_wrap(~ factor(lyr, levels = layer_names), ncol = 5) +
    labs(title = sprintf("%d Monsoon Season: 5-Day Evolution (Observed)", target_year), x = "", y = "") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 22, hjust = 0.5, margin = margin(b = 20)),
      strip.text = element_text(face = "bold", size = 9),
      strip.background = element_rect(fill = "gray90", color = NA),
      axis.text = element_blank(), axis.ticks = element_blank(), panel.spacing = unit(0.4, "lines"),
      legend.position = "right", legend.key.height = unit(1.3, "cm")
    )

  out_absolute <- file.path(out_dir, sprintf("pentad_absolute_%d.png", target_year))
  ggsave(out_absolute, plot = p_absolute, width = 16, height = 12, dpi = 300, bg = "white")
  cat(sprintf("Successfully saved: %s\n", out_absolute))

  # -----------------------------------------------------------------------
  # 5. PLOT 2: ANOMALY 5x5 GRID (vs. 1991-2025 PIXEL BASELINE)
  # -----------------------------------------------------------------------
  anom_pal <- c("#543005", "#8c510a", "#bf812d", "#f5f5f5", "#35978f", "#01665e", "#003c30")
  feature_colors_anom <- c(
    "Z500 Anom (+)" = "hotpink",
    "Z500 Anom (-)" = "steelblue"
  )

  p_anomaly <- ggplot() +
    geom_spatraster(data = tcwv_anom) +
    scale_fill_gradientn(colors = anom_pal, name = "TCWV Anom\n(kg/m2)", limits = c(-15, 15),
                          oob = scales::squish, na.value = "transparent") +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.5) +
    geom_polygon(data = state_map, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.3) +
    geom_spatraster_contour(data = z500_anom, aes(color = "Z500 Anom (-)"), breaks = seq(-300, -10, by = 10), linewidth = 0.4, linetype = "dashed") +
    geom_spatraster_contour(data = z500_anom, aes(color = "Z500 Anom (+)"), breaks = seq(10, 300, by = 10), linewidth = 0.4, linetype = "solid") +
    scale_color_manual(name = "Features", values = feature_colors_anom) +
    coord_sf(xlim = c(b[1], b[2]), ylim = c(b[3], b[4]), expand = FALSE) +
    facet_wrap(~ factor(lyr, levels = layer_names), ncol = 5) +
    labs(title = sprintf("%d Monsoon Season: 5-Day Evolution (Anomaly vs. 1991-2025)", target_year), x = "", y = "") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 22, hjust = 0.5, margin = margin(b = 20)),
      strip.text = element_text(face = "bold", size = 9),
      strip.background = element_rect(fill = "gray90", color = NA),
      axis.text = element_blank(), axis.ticks = element_blank(), panel.spacing = unit(0.4, "lines"),
      legend.position = "right", legend.key.height = unit(1.3, "cm")
    )

  out_anomaly <- file.path(out_dir, sprintf("pentad_anomaly_%d.png", target_year))
  ggsave(out_anomaly, plot = p_anomaly, width = 16, height = 12, dpi = 300, bg = "white")
  cat(sprintf("Successfully saved: %s\n", out_anomaly))

  invisible(list(absolute = p_absolute, anomaly = p_anomaly))
}

# =========================================================================
# EXECUTION
# =========================================================================

# --- Single year ---
#generate_pentad_year(target_year = 2020)

#--- Batch: every year with processed data on hand ---
for (yr in 1991:2025) {
  tryCatch({
    generate_pentad_year(target_year = yr)
  }, error = function(e) {
    cat(sprintf("Skipped %d due to error: %s\n", yr, e$message))
  })
}

cat("\n========================================\n")
cat("PENTAD GRID GENERATION COMPLETE!\n")
cat("========================================\n")
