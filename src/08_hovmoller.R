# =========================================================================
# HOVMÖLLER DIAGRAM GENERATOR (Date x Latitude, Zonal-Mean TCWV)
# =========================================================================
# Collapses a single monsoon season (June 1 - Oct 31) into one dense 2D
# image: the zonal (longitudinal) mean of TCWV across a chosen corridor,
# plotted as Date (x) by Latitude (y). Reveals the northward moisture
# "surge" up the Gulf of California, and how it compares to the
# 1991-2025 climatological baseline.
#
# Assumes working directory = project root (monsoonTracker/), consistent
# with scripts 01-05 and 07.
# =========================================================================

library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(maps)

# -------------------------------------------------------------------------
# 0. HELPER: ZONAL (LONGITUDINAL) MEAN, ONE VALUE PER LATITUDE ROW/LAYER
# -------------------------------------------------------------------------
# Returns a nrow(r) x nlyr(r) matrix - the across-longitude mean for every
# latitude row, for every layer. Built directly from values() (which terra
# always returns in row-major cell order: row 1 left-to-right, then row 2,
# etc.) instead of relying on terra::aggregate()'s fact = c(h, v) argument,
# which aggregates ROWS first and is easy to get backwards.
zonal_mean_by_row <- function(r) {
  n_row <- nrow(r); n_col <- ncol(r); n_lyr <- nlyr(r)
  v <- values(r)
  out <- matrix(NA_real_, nrow = n_row, ncol = n_lyr)
  for (L in 1:n_lyr) {
    out[, L] <- rowMeans(matrix(v[, L], nrow = n_row, ncol = n_col, byrow = TRUE), na.rm = TRUE)
  }
  out
}

# -------------------------------------------------------------------------
# 1. CORE FUNCTION
# -------------------------------------------------------------------------
generate_hovmoller <- function(target_year,
                                lon_range = c(-115, -105),   # Gulf of California corridor
                                lat_range = c(15, 45),       # Full monsoon domain
                                slice_dir = "data/processed/4slices",
                                climo_dir = "data/processed/4slices/climo",
                                out_dir   = "NAM_tracker/output/hovmoller") {

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  months_to_process <- c("06", "07", "08", "09", "10")
  corridor_box <- ext(lon_range[1], lon_range[2], lat_range[1], lat_range[2])

  cat(sprintf("Building Hovmoller diagram for %d (corridor %d°W-%d°W)...\n",
              target_year, abs(lon_range[1]), abs(lon_range[2])))

  # -----------------------------------------------------------------------
  # 1a. PULL & ZONALLY AVERAGE THE TARGET YEAR'S DAILY TCWV
  # -----------------------------------------------------------------------
  month_zonal_list <- list()
  lat_vals <- NULL

  for (m in months_to_process) {
    nc_file <- list.files(slice_dir, pattern = sprintf("single_%d_%s\\.nc$", target_year, m),
                           full.names = TRUE)[1]
    if (is.na(nc_file) || !file.exists(nc_file)) {
      stop(sprintf("Missing single-level TCWV file for %d-%s", target_year, m))
    }

    r_month <- crop(rast(nc_file), corridor_box)

    # Collapse every longitude column into a zonal mean, leaving each
    # latitude row (and each daily layer) intact.
    month_zonal_list[[m]] <- zonal_mean_by_row(r_month)

    if (is.null(lat_vals)) lat_vals <- yFromRow(r_month, 1:nrow(r_month))
  }

  obs_mat <- do.call(cbind, month_zonal_list)   # n_lat_rows x 153

  # -----------------------------------------------------------------------
  # 1b. PULL & ZONALLY AVERAGE THE 153-DAY CLIMATOLOGICAL BASELINE
  # -----------------------------------------------------------------------
  climo_path <- file.path(climo_dir, "climatology_tcwv_raw.tif")
  if (!file.exists(climo_path)) stop("Missing TCWV climatology baseline: ", climo_path)

  climo_full <- crop(rast(climo_path), corridor_box)

  if (nrow(climo_full) != length(lat_vals)) {
    stop(sprintf(
      "Grid mismatch: season data has %d latitude rows but climatology has %d. Check that both share the same ERA5 grid/resolution.",
      length(lat_vals), nrow(climo_full)
    ))
  }

  climo_mat <- zonal_mean_by_row(climo_full)   # n_lat_rows x 153

  # -----------------------------------------------------------------------
  # 2. ASSEMBLE THE DATE x LATITUDE MATRICES
  # -----------------------------------------------------------------------
  date_vals <- seq(as.Date(sprintf("%d-06-01", target_year)),
                    as.Date(sprintf("%d-10-31", target_year)), by = "day")

  anom_mat <- obs_mat - climo_mat

  colnames(obs_mat)  <- as.character(date_vals)
  colnames(anom_mat) <- as.character(date_vals)

  obs_df <- as.data.frame(obs_mat) %>%
    mutate(lat = lat_vals) %>%
    tidyr::pivot_longer(cols = -lat, names_to = "date", values_to = "tcwv") %>%
    mutate(date = as.Date(date))

  anom_df <- as.data.frame(anom_mat) %>%
    mutate(lat = lat_vals) %>%
    tidyr::pivot_longer(cols = -lat, names_to = "date", values_to = "tcwv_anom") %>%
    mutate(date = as.Date(date))

  # Climatological zonal-mean TCWV, for overlaying a "normal" 25 kg/m2
  # threshold line on the observed panel for direct comparison.
  climo_mat_named <- climo_mat
  colnames(climo_mat_named) <- as.character(date_vals)

  climo_df <- as.data.frame(climo_mat_named) %>%
    mutate(lat = lat_vals) %>%
    tidyr::pivot_longer(cols = -lat, names_to = "date", values_to = "tcwv_climo") %>%
    mutate(date = as.Date(date))

  # -----------------------------------------------------------------------
  # 3. PLOT: OBSERVED ZONAL TCWV (TOP) + ANOMALY VS. CLIMATOLOGY (BOTTOM)
  # -----------------------------------------------------------------------
  # Reuse the project's standard TCWV color ramps for visual consistency
  # with the daily synoptic dashboard (05B) and climatology grids (06/06B).
  obs_pal <- colorRampPalette(c(
    "#4a2c11", "#bc966c", "#f4efe6", "#31a354", "#008080", "#064cb5", "#5106b5"
  ))(100)

  anom_pal <- c("#543005", "#8c510a", "#bf812d", "#f5f5f5", "#35978f", "#01665e", "#003c30")

  # -----------------------------------------------------------------------
  # 3a. LOCATOR INSET: WHERE IS THIS CORRIDOR ON THE MAP?
  # -----------------------------------------------------------------------
  # Small reference map of the full project domain with the zonal-average
  # corridor (lon_range x lat_range) highlighted, so the Hovmoller's
  # latitude axis can be read against real geography at a glance.
  world_map <- map_data("world")
  state_map <- map_data("state")

  lat_breaks <- seq(15, 45, by = 5)

  p_locator <- ggplot() +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
                 fill = "grey85", color = "grey45", linewidth = 0.15) +
    geom_polygon(data = state_map, aes(x = long, y = lat, group = group),
                 fill = NA, color = "grey45", linewidth = 0.1) +
    geom_hline(yintercept = lat_breaks, color = "grey40", linewidth = 0.15, linetype = "dashed") +
    geom_rect(aes(xmin = lon_range[1], xmax = lon_range[2],
                  ymin = lat_range[1], ymax = lat_range[2]),
              fill = "red", alpha = 0.25, color = "red", linewidth = 0.5) +
    scale_y_continuous(breaks = lat_breaks, position = "right", expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    coord_sf(xlim = c(-125, -95), ylim = c(15, 45), expand = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "aliceblue", color = "black", linewidth = 0.4),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.4),
      axis.text.y = element_text(size = 5, color = "grey20", margin = margin(l = 1)),
      axis.ticks.y = element_line(linewidth = 0.2, color = "grey20"),
      axis.ticks.length.y = unit(1, "pt"),
      plot.margin = margin(0, 0, 0, 0)
    )

  obs_label   <- sprintf("Observed (%d)", target_year)
  climo_label <- "1991-2025 Mean"

  p_obs <- ggplot(obs_df, aes(x = date, y = lat, fill = tcwv)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradientn(colors = obs_pal, name = "Zonal TCWV\n(kg/m2)",
                          limits = c(5, 65), oob = scales::squish) +
    geom_contour(aes(z = tcwv, color = obs_label, linetype = obs_label), breaks = 25, linewidth = 0.45) +
    geom_contour(data = climo_df, aes(x = date, y = lat, z = tcwv_climo, color = climo_label, linetype = climo_label),
                 breaks = 25, linewidth = 0.5, inherit.aes = FALSE) +
    scale_color_manual(name = "25 kg/m² Threshold",
                        values = setNames(c("black", "deeppink"), c(obs_label, climo_label))) +
    scale_linetype_manual(name = "25 kg/m² Threshold",
                           values = setNames(c("solid", "dashed"), c(obs_label, climo_label))) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b 1", expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      title = sprintf("%d Monsoon Hovmöller: Zonal-Mean TCWV (%d°W–%d°W corridor)",
                       target_year, abs(lon_range[1]), abs(lon_range[2])),
      subtitle = "25 kg/m² monsoon moisture threshold: solid = observed, dashed = 1991-2025 climatological mean",
      x = NULL, y = "Latitude (°N)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30")
    )

  # Drop the locator map into the top-right corner of the observed panel
  p_obs <- p_obs + inset_element(p_locator, left = 0.87, bottom = 0.76, right = 0.995, top = 0.985)

  p_anom <- ggplot(anom_df, aes(x = date, y = lat, fill = tcwv_anom)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradientn(colors = anom_pal, name = "TCWV Anom\n(kg/m2)",
                          limits = c(-15, 15), oob = scales::squish, na.value = "transparent") +
    scale_x_date(date_breaks = "1 month", date_labels = "%b 1", expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = "Departure from 1991–2025 Climatology", x = "Date", y = "Latitude (°N)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  footnote_text <- paste0(
    "Data Source: ERA5 Atmospheric Reanalysis | Climatological Base Period: 1991-2025\n",
    "Zonal mean computed by collapsing all longitudes within the corridor at each latitude/day.\n",
    "Produced by: Climate Science Applications Program, University of Arizona"
  )

  hov_combo <- p_obs / p_anom +
    plot_layout(guides = "collect") +
    plot_annotation(
      caption = footnote_text,
      theme = theme(
        plot.caption = element_text(size = 8, color = "grey35", hjust = 0,
                                     face = "italic", lineheight = 1.2,
                                     margin = margin(t = 10))
      )
    )

  out_file <- file.path(out_dir, sprintf("hovmoller_%d.png", target_year))
  ggsave(out_file, plot = hov_combo, width = 12, height = 10, dpi = 200, bg = "white")
  cat(sprintf("Successfully saved: %s\n", out_file))

  invisible(hov_combo)
}

# =========================================================================
# EXECUTION
# =========================================================================

# --- Single year ---
generate_hovmoller(target_year = 2021)

#--- Batch: every year with processed data on hand ---
for (yr in 1991:2025) {
  tryCatch({
    generate_hovmoller(target_year = yr)
  }, error = function(e) {
    cat(sprintf("Skipped %d due to error: %s\n", yr, e$message))
  })
}

cat("\n========================================\n")
cat("HOVMÖLLER GENERATION COMPLETE!\n")
cat("========================================\n")
