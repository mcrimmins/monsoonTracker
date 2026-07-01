# =========================================================================
# FULL-SEASON MACRO SUMMARY (Mean State + Ridge Evolution, Anomaly Footprint)
# =========================================================================
# Two-panel "season report card" for a single year:
#   Panel A - Mean State & Ridge Evolution: June-Oct mean TCWV/Z500, with
#             the full-season "H" ridge trajectory plotted on top (monthly
#             waypoints labeled Jun-Oct) so the season's *movement*, not
#             just its average, is visible in one image.
#   Panel B - Seasonal Anomaly Footprint: net TCWV/Z500 anomaly vs. the
#             1991-2025 climatology, showing where the season's moisture
#             anomaly actually sat (Sonora core, displaced into NM, or
#             largely absent).
#
# Assumes working directory = project root (monsoonTracker/), consistent
# with scripts 01-08.
# =========================================================================

library(terra)
library(ggplot2)
library(tidyterra)
library(maps)
library(patchwork)

# -------------------------------------------------------------------------
# 1. CORE FUNCTION
# -------------------------------------------------------------------------
generate_seasonal_summary <- function(target_year,
                                       slice_dir = "data/processed/4slices",
                                       anom_dir  = "data/processed/4slices/anomalies",
                                       out_dir   = "NAM_tracker/output/seasonal_summary") {

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  months_to_process <- c("06", "07", "08", "09", "10")

  cat(sprintf("Building full-season macro summary for %d...\n", target_year))

  # -----------------------------------------------------------------------
  # 1a. LOAD THE FULL 153-DAY SEASON: TCWV & Z500
  # -----------------------------------------------------------------------
  tcwv_list <- list()
  z500_list <- list()

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
    # Block 7 of 8 is z500 (see var_metadata in 04_calculate_daily_anomalies.R)
    z500_list[[m]] <- r_press[[(1:days_in_month) + (6 * days_in_month)]] / 9.80665
  }

  season_tcwv <- rast(tcwv_list)   # 153 layers
  season_z500 <- rast(z500_list)   # 153 layers

  date_vals <- seq(as.Date(sprintf("%d-06-01", target_year)),
                    as.Date(sprintf("%d-10-31", target_year)), by = "day")

  mean_tcwv <- mean(season_tcwv, na.rm = TRUE)
  mean_z500 <- mean(season_z500, na.rm = TRUE)

  # -----------------------------------------------------------------------
  # 1b. LOAD THE PRECOMPUTED DAILY ANOMALY STACKS & COLLAPSE TO A SEASON MEAN
  # -----------------------------------------------------------------------
  anom_tcwv_path <- file.path(anom_dir, sprintf("anomaly_tcwv_%d.tif", target_year))
  anom_z500_path <- file.path(anom_dir, sprintf("anomaly_z500_%d.tif", target_year))
  if (!file.exists(anom_tcwv_path) || !file.exists(anom_z500_path)) {
    stop(sprintf("Missing anomaly stacks for %d. Run 04_calculate_daily_anomalies.R first.", target_year))
  }

  mean_tcwv_anom <- mean(rast(anom_tcwv_path), na.rm = TRUE)
  mean_z500_anom <- mean(rast(anom_z500_path), na.rm = TRUE)

  # -----------------------------------------------------------------------
  # 2. TRACK THE "H" RIDGE EVERY DAY OF THE SEASON
  # -----------------------------------------------------------------------
  # Same two-pass macro/micro focal-smoothing tracker used in 05B and
  # 06_plot_climo.R, run once per day across the full 153-day season.
  visual_smooth_mat  <- matrix(1 / 25, nrow = 5, ncol = 5)
  tracker_smooth_mat <- matrix(1 / 225, nrow = 15, ncol = 15)
  monsoon_box <- ext(-114, -95, 15, 45)

  h_lon <- numeric(nlyr(season_z500))
  h_lat <- numeric(nlyr(season_z500))

  for (i in 1:nlyr(season_z500)) {
    z_search <- crop(season_z500[[i]], monsoon_box)

    # PASS 1: MACRO TRACKER
    z_macro <- focal(z_search, w = tracker_smooth_mat, fun = mean, na.rm = TRUE)
    macro_max <- global(z_macro, "max", na.rm = TRUE)[1, 1]
    macro_core <- ifel(z_macro >= (macro_max - 10), 1, NA)

    blobs <- patches(macro_core)
    blob_sizes <- freq(blobs)

    if (nrow(blob_sizes) > 0) {
      largest_blob_id <- blob_sizes$value[which.max(blob_sizes$count)]
      primary_macro_core <- ifel(blobs == largest_blob_id, 1, NA)

      # PASS 2: MICRO PINPOINT
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

  trajectory_df <- data.frame(date = date_vals, lon = h_lon, lat = h_lat)
  trajectory_df <- trajectory_df[stats::complete.cases(trajectory_df), ]

  # Monthly waypoints (1st of each month) so direction/progression along
  # the path can be read without needing a continuous color gradient.
  waypoint_dates <- as.Date(sprintf("%d-%s-01", target_year, months_to_process))
  waypoints_df <- trajectory_df[trajectory_df$date %in% waypoint_dates, ]
  waypoints_df$label <- format(waypoints_df$date, "%b")

  # -----------------------------------------------------------------------
  # 3. PLOT
  # -----------------------------------------------------------------------
  # Reuse the project's standard TCWV color ramps for visual consistency
  # with the daily synoptic dashboard (05B), Hovmoller (08), and climo grids.
  obs_pal <- colorRampPalette(c(
    "#4a2c11", "#bc966c", "#f4efe6", "#31a354", "#008080", "#064cb5", "#5106b5"
  ))(100)

  anom_pal <- c("#543005", "#8c510a", "#bf812d", "#f5f5f5", "#35978f", "#01665e", "#003c30")

  feature_colors <- c(
    "Z500 Heights"          = "gray30",
    "Ridge Path (Jun-Oct)"  = "darkred",
    "Z500 Anom (+)"         = "hotpink",
    "Z500 Anom (-)"         = "steelblue"
  )

  world_map <- map_data("world")
  state_map <- map_data("state")
  b <- ext(mean_tcwv)

  # --- PANEL A: MEAN STATE + RIDGE EVOLUTION ---
  p_mean <- ggplot() +
    geom_spatraster(data = mean_tcwv) +
    scale_fill_gradientn(colors = obs_pal, name = "Mean TCWV\n(kg/m2)",
                          na.value = "transparent", limits = c(5, 65), oob = scales::squish) +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.5) +
    geom_polygon(data = state_map, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.3) +
    geom_spatraster_contour(data = mean_z500, aes(color = "Z500 Heights"),
                             breaks = seq(5500, 6000, by = 20), linewidth = 0.3) +
    geom_path(data = trajectory_df, aes(x = lon, y = lat, color = "Ridge Path (Jun-Oct)"),
              linewidth = 0.7, alpha = 0.85) +
    geom_point(data = waypoints_df, aes(x = lon, y = lat), color = "darkred", size = 2.2) +
    geom_text(data = waypoints_df, aes(x = lon, y = lat, label = label),
              color = "darkred", fontface = "bold", size = 3, vjust = -0.9) +
    geom_point(data = trajectory_df[1, ], aes(x = lon, y = lat),
               shape = 21, fill = "white", color = "darkred", size = 3, stroke = 1.1) +
    geom_point(data = trajectory_df[nrow(trajectory_df), ], aes(x = lon, y = lat),
               shape = 22, fill = "darkred", color = "darkred", size = 3) +
    scale_color_manual(name = "Features", values = feature_colors, drop = FALSE) +
    coord_sf(xlim = c(b[1], b[2]), ylim = c(b[3], b[4]), expand = FALSE) +
    labs(title = "Mean State & Ridge Evolution", x = "", y = "") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text = element_text(color = "black", size = 10)
    )

  # --- PANEL B: SEASONAL ANOMALY FOOTPRINT ---
  p_anom <- ggplot() +
    geom_spatraster(data = mean_tcwv_anom) +
    scale_fill_gradientn(colors = anom_pal, name = "TCWV Anom\n(kg/m2)",
                          limits = c(-15, 15), oob = scales::squish, na.value = "transparent") +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.5) +
    geom_polygon(data = state_map, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.3) +
    geom_spatraster_contour(data = mean_z500_anom, aes(color = "Z500 Anom (-)"),
                             breaks = seq(-300, -10, by = 10), linewidth = 0.5, linetype = "dashed") +
    geom_spatraster_contour(data = mean_z500_anom, aes(color = "Z500 Anom (+)"),
                             breaks = seq(10, 300, by = 10), linewidth = 0.5, linetype = "solid") +
    scale_color_manual(name = "Features", values = feature_colors, drop = FALSE) +
    coord_sf(xlim = c(b[1], b[2]), ylim = c(b[3], b[4]), expand = FALSE) +
    labs(title = "Seasonal Anomaly Footprint", x = "", y = "") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.x = element_text(color = "black", size = 10),
      axis.text.y = element_blank(), axis.ticks.y = element_blank()
    )

  footnote_text <- paste0(
    "Data Source: ERA5 Atmospheric Reanalysis | Climatological Base Period: 1991-2025\n",
    "Ridge path = daily 500 hPa height maximum, tracked Jun 1 - Oct 31; ",
    "circle = season start, square = season end, labeled dots = 1st of month.\n",
    "Produced by: Climate Science Applications Program, University of Arizona"
  )

  combo <- p_mean + p_anom +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = sprintf("%d Monsoon Season: Full-Season Macro Summary", target_year),
      caption = footnote_text,
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 15)),
        plot.caption = element_text(
          size = 8.5, color = "grey35", hjust = 0, face = "italic",
          lineheight = 1.2, margin = margin(t = 15, l = 10)
        )
      )
    )

  out_file <- file.path(out_dir, sprintf("seasonal_summary_%d.png", target_year))
  ggsave(out_file, plot = combo, width = 16, height = 7.5, dpi = 150, bg = "white")
  cat(sprintf("Successfully saved: %s\n", out_file))

  invisible(combo)
}

# =========================================================================
# EXECUTION
# =========================================================================

# --- Single year ---
generate_seasonal_summary(target_year = 2020)

# --- Batch: every year with processed data + anomalies on hand ---
# for (yr in 1991:2025) {
#   tryCatch({
#     generate_seasonal_summary(target_year = yr)
#   }, error = function(e) {
#     cat(sprintf("Skipped %d due to error: %s\n", yr, e$message))
#   })
# }

cat("\n========================================\n")
cat("SEASONAL SUMMARY GENERATION COMPLETE!\n")
cat("========================================\n")
