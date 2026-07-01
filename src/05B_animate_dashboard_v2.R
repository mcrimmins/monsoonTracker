# =========================================================================
# SYNOPTIC DASHBOARD ANIMATION GENERATOR (V4)
# =========================================================================

library(terra)
library(ggplot2)
library(tidyterra) 
library(maps)
library(patchwork)

# -------------------------------------------------------------------------
# 1. DEFINE THE CORE FRAME GENERATION FUNCTION
# -------------------------------------------------------------------------
generate_synoptic_frame <- function(test_year, test_month, test_day, 
                                    slice_dir = "data/processed/4slices",
                                    climo_dir = "data/processed/4slices/climo", # Path to your climo folder
                                    out_dir = "output/frames",
                                    show_trajectory = TRUE, 
                                    tail_length_days = 3) {
  
  # Ensure output directory exists
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # --- GLOBAL COLOR SCALES ---
  feature_colors <- c(
    "Ridge Trajectory (H)" = "darkred",
    "TCWV > 25 kg/m2"   = "orange",
    "Q700 > 6 g/kg"     = "darkblue",
    "Z500 Heights"      = "gray30",
    "Z500 Anom (+)"     = "hotpink",    
    "Z500 Anom (-)"     = "steelblue"    
  )
  
  # --- DYNAMIC TITLE & DATE HANDLING ---
  map_date <- as.Date(sprintf("%d-%s-%02d", test_year, test_month, test_day))
  formatted_date <- format(map_date, "%d %b %Y")
  dashboard_title <- sprintf("NAM Daily Atmospheric Analysis\n%s", formatted_date)
  
  cat(sprintf("Rendering frame for: %s...\n", formatted_date))
  
  # --- LOAD DAILY OBSERVED DATA ---
  nc_single <- list.files(slice_dir, pattern = sprintf("single_%d_%s\\.nc$", test_year, test_month), full.names = TRUE)[1]
  r_single  <- rast(nc_single)
  tcwv_day  <- r_single[[test_day]] 
  
  nc_press <- list.files(slice_dir, pattern = sprintf("pressure_%d_%s\\.nc$", test_year, test_month), full.names = TRUE)[1]
  r_press  <- rast(nc_press)
  days_in_month <- nlyr(r_press) / 8
  
  q700_day <- r_press[[test_day + (1 * days_in_month)]] * 1000 
  z500_day <- r_press[[test_day + (6 * days_in_month)]] / 9.80665 
  
  # -------------------------------------------------------------------------
  # 2. CALCULATE TRUE PIXEL-BY-PIXEL CLIMATOLOGICAL ANOMALIES
  # -------------------------------------------------------------------------
  # Load the 153-day climatology stacks
  tcwv_climo <- rast(file.path(climo_dir, "climatology_tcwv_raw.tif"))
  z500_climo <- rast(file.path(climo_dir, "climatology_z500_raw.tif"))
  
  # Calculate index relative to June 1st (June 1 = Layer 1, Oct 31 = Layer 153)
  start_of_monsoon <- as.Date(sprintf("%d-06-01", test_year))
  day_idx <- as.integer(map_date - start_of_monsoon) + 1
  
  # Extract the historical baseline grid for this specific calendar day
  tcwv_baseline <- tcwv_climo[[day_idx]]
  z500_baseline <- z500_climo[[day_idx]]
  
  # True grid spatial subtraction
  tcwv_anom <- tcwv_day - tcwv_baseline
  z500_anom <- z500_day - z500_baseline
  
  # -------------------------------------------------------------------------
  # 3. ADVANCED "H" TRACKING
  # -------------------------------------------------------------------------
  if (show_trajectory && tail_length_days > 0) {
    start_tail <- max(1, test_day - tail_length_days)
  } else {
    start_tail <- test_day 
  }
  tail_days <- start_tail:test_day
  
  h_lon <- numeric(length(tail_days))
  h_lat <- numeric(length(tail_days))
  
  visual_smooth_mat  <- matrix(1/25, nrow = 5, ncol = 5)    
  tracker_smooth_mat <- matrix(1/225, nrow = 15, ncol = 15) 
  
  monsoon_box <- ext(-125, -95, 15, 45)
  
  for (i in seq_along(tail_days)) {
    d <- tail_days[i]
    z_d <- r_press[[d + (6 * days_in_month)]] / 9.80665
    z_search <- crop(z_d, monsoon_box)
    
    # PASS 1: MACRO TRACKER 
    z_macro <- focal(z_search, w = tracker_smooth_mat, fun = mean, na.rm = TRUE)
    macro_max <- global(z_macro, "max", na.rm = TRUE)[1, 1]
    macro_core <- ifel(z_macro >= (macro_max - 10), 1, NA)
    
    blobs <- patches(macro_core)
    blob_sizes <- freq(blobs)
    largest_blob_id <- blob_sizes$value[which.max(blob_sizes$count)]
    primary_macro_core <- ifel(blobs == largest_blob_id, 1, NA)
    
    # PASS 2: MICRO PINPOINT 
    z_visual <- focal(z_search, w = visual_smooth_mat, fun = mean, na.rm = TRUE)
    z_isolated <- mask(z_visual, primary_macro_core)
    
    true_max <- global(z_isolated, "max", na.rm = TRUE)[1, 1]
    final_core <- ifel(z_isolated >= (true_max - 2), 1, NA)
    
    core_pts <- crds(as.points(final_core))
    centroid <- colMeans(core_pts)
    
    h_lon[i] <- centroid[1]
    h_lat[i] <- centroid[2]
  }
  
  traj_df <- data.frame(lon = h_lon, lat = h_lat)
  
  # -------------------------------------------------------------------------
  # 4. SPATIAL SMOOTHING FOR CONTOURS
  # -------------------------------------------------------------------------
  q700_smooth             <- focal(q700_day, w = visual_smooth_mat, fun = mean, na.rm = TRUE)
  z500_smooth             <- focal(z500_day, w = visual_smooth_mat, fun = mean, na.rm = TRUE)
  z500_anom_smooth        <- focal(z500_anom, w = visual_smooth_mat, fun = mean, na.rm = TRUE)
  tcwv_smooth_for_contour <- focal(tcwv_day, w = visual_smooth_mat, fun = mean, na.rm = TRUE)
  
  # Base Maps
  world_map <- map_data("world")
  state_map <- map_data("state")
  b <- ext(tcwv_day) 
  #obs_pal <- colorRampPalette(c("wheat", "lightgreen", "deepskyblue3", "darkblue"))(100)
  # Paste this inside your generate_synoptic_frame function around line 121
  obs_pal <- colorRampPalette(c(
    "#4a2c11",  # 5-12 kg/m2:  Severe Desert Dryness (Deep Brown)
    "#bc966c",  # 12-20 kg/m2: Modified Continental Dry Air (Tan)
    "#f4efe6",  # 20-28 kg/m2: The Transition Zone / Sharp Gradient (Cream)
    "#31a354",  # 28-38 kg/m2: Monsoon Periphery / First Surges (Vibrant Green)
    "#008080",  # 38-46 kg/m2: Deep Low-level Moisture (Teal)
    "#064cb5",  # 46-55 kg/m2: Core Monsoon / Gulf Surges (Rich Blue)
    "#5106b5"   # 55-65 kg/m2: Tropical Storms / Extreme PWAT (Deep Purple)
  ))(100)
  
  # -------------------------------------------------------------------------
  # 5. BUILD PLOT A: OBSERVED STATE
  # -------------------------------------------------------------------------
  p_obs <- ggplot() +
    geom_spatraster(data = tcwv_day) +
    scale_fill_gradientn(colors = obs_pal, name = "Observed TCWV\n(kg/m2)", 
                         na.value = "transparent", limits = c(5, 65), oob = scales::squish) +
    
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.5) +
    geom_polygon(data = state_map, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.3) +
    
    geom_spatraster_contour(data = z500_smooth, aes(color = "Z500 Heights"), breaks = seq(5500, 6000, by = 20), linewidth = 0.3) +
    geom_spatraster_contour(data = q700_smooth, aes(color = "Q700 > 6 g/kg"), breaks = 6, linewidth = 0.6, linetype = "dashed") +
    geom_spatraster_contour(data = tcwv_smooth_for_contour, aes(color = "TCWV > 25 kg/m2"), breaks = 25, linewidth = 0.8, linetype = "solid") 
  
  if (show_trajectory && nrow(traj_df) > 1) {
    p_obs <- p_obs + 
      geom_path(data = traj_df, aes(x = lon, y = lat, color = "Ridge Trajectory (H)"), linewidth = 0.8) +
      geom_point(data = traj_df[-nrow(traj_df), ], aes(x = lon, y = lat, color = "Ridge Trajectory (H)"), size = 2)
  }
  
  p_obs <- p_obs + 
    annotate("text", x = traj_df$lon[nrow(traj_df)], y = traj_df$lat[nrow(traj_df)], label = "H", color = "blue1", fontface = "bold", size = 8) +
    scale_color_manual(name = "Contours & Features", values = feature_colors, drop = FALSE) +
    coord_sf(xlim = c(b[1], b[2]), ylim = c(b[3], b[4]), expand = FALSE) + 
    labs(title = "Observed", x = "", y = "") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5), 
      axis.text = element_text(color = "black", size = 10)
    )
  
  # -------------------------------------------------------------------------
  # 6. BUILD PLOT B: CLIMATOLOGICAL ANOMALY
  # -------------------------------------------------------------------------
  p_anom <- ggplot() +
    geom_spatraster(data = tcwv_anom) +
    # scale_fill_gradient2(low = "saddlebrown", mid = "white", high = "forestgreen", midpoint = 0, 
    #                      name = "TCWV Anom\n(kg/m2)", limits = c(-15, 15), oob = scales::squish, na.value = "transparent") +

    scale_fill_gradientn(
      colors = c(
        "#543005",  # Extreme Dry Anomaly (-15 kg/m2)
        "#8c510a",  # Moderate Dry Anomaly
        "#bf812d",  # Light Dry Anomaly
        "#f5f5f5",  # Neutral / Climatological Normal (0 kg/m2)
        "#35978f",  # Light Wet Anomaly
        "#01665e",  # Moderate Wet Anomaly
        "#003c30"   # Extreme Wet Anomaly (+15 kg/m2)
      ),
      name = "TCWV Anom\n(kg/m2)", 
      limits = c(-15, 15), 
      oob = scales::squish, 
      na.value = "transparent"
    ) +
    
    
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.5) +
    geom_polygon(data = state_map, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.3) +
    
    # These will now highlight true anomalous centers (e.g., cut-off lows, anomalous high blocks)
    geom_spatraster_contour(data = z500_anom_smooth, aes(color = "Z500 Anom (-)"), breaks = seq(-300, -10, by = 10), linewidth = 0.6, linetype = "dashed") +
    geom_spatraster_contour(data = z500_anom_smooth, aes(color = "Z500 Anom (+)"), breaks = seq(10, 300, by = 10), linewidth = 0.6, linetype = "solid") +
    
    geom_spatraster_contour_text(data = z500_anom_smooth, breaks = seq(-300, 300, by = 10), 
                                 color = "grey20", size = 3, fontface = "bold") +
    
    scale_color_manual(name = "Contours & Features", values = feature_colors, drop = FALSE) +
    coord_sf(xlim = c(b[1], b[2]), ylim = c(b[3], b[4]), expand = FALSE) + 
    labs(title = "Anomaly", x = "", y = "") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.x = element_text(color = "black", size = 10),
      axis.text.y = element_blank(), axis.ticks.y = element_blank() 
    )
  
  # -------------------------------------------------------------------------
  # 7. COMBINE, ANNOTATE, AND SAVE
  # -------------------------------------------------------------------------
  out_file <- file.path(out_dir, sprintf("dashboard_%d%s%02d.png", test_year, test_month, test_day))
  
  # Construct a clean, technical footnote string
  footnote_text <- paste0(
    "Data Source: ERA5 Atmospheric Reanalysis | Climatological Base Period: 1991-2025\n",
    "Methodology: Anomalies computed from 6,12,18,00Z hourly values; ",
    "Contours processed with a 5x5 focal mean for synoptic-scale smoothing.\n",
    "Produced by: Climate Science Applications Program, University of Arizona"
  )
  
  dashboard <- p_obs + p_anom +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = dashboard_title,
      caption = footnote_text,
      theme = theme(
        # Main Title Styling
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 15)),
        
        # Footnote Styling: Small, muted, left-aligned, italicized
        plot.caption = element_text(
          size = 8.5, 
          color = "grey35", 
          hjust = 0, 
          face = "italic", 
          lineheight = 1.2, 
          margin = margin(t = 15, l = 10)
        )
      )
    )
  
  ggsave(out_file, plot = dashboard, width = 16, height = 7.5, dpi = 150, bg = "white")
  cat(sprintf("Successfully saved: %s\n", out_file))
}

# =========================================================================
# EXECUTION: NESTED LOOPS FOR MULTIPLE YEARS & MONTHS
# =========================================================================

target_years   <- 2023:2023 
target_months  <- c("06", "07", "08", "09", "10") 
#target_months  <- c("06")
base_slice_dir <- "data/processed/4slices"
base_out_dir   <- "NAM_tracker/output/frames"

for (y in target_years) {
  year_out_dir <- file.path(base_out_dir, as.character(y))
  
  for (m in target_months) {
    nc_file <- list.files(base_slice_dir, pattern = sprintf("single_%d_%s\\.nc$", y, m), full.names = TRUE)[1]
    
    if (is.na(nc_file) || !file.exists(nc_file)) {
      cat(sprintf("Data not found for %d-%s. Skipping to next month...\n", y, m))
      next
    }
    
    temp_r <- terra::rast(nc_file)
    num_days <- terra::nlyr(temp_r)
    
    cat(sprintf("\n--- Starting %d-%s (%d days) ---\n", y, m, num_days))
    
    for (d in 1:num_days) {
      tryCatch({
        generate_synoptic_frame(
          test_year = y, 
          test_month = m, 
          test_day = d, 
          slice_dir = base_slice_dir,
          out_dir = year_out_dir,
          show_trajectory = TRUE,
          tail_length_days = 4
        )
      }, error = function(e) {
        cat(sprintf("Skipped Day %02d due to error: %s\n", d, e$message))
      })
    }
  }
}

cat("\n========================================\n")
cat("BATCH PROCESSING COMPLETE!\n")
cat("========================================\n")


# =========================================================================
# SYNOPTIC DASHBOARD: ANIMATOR SCRIPT
# =========================================================================
# library(gifski)
# library(av)
# 
# base_frame_dir  <- "output/frames"
# output_anim_dir <- "output/animations"
# output_format   <- "mp4" 
# gif_delay       <- 0.25 
# mp4_fps         <- 4    
# 
# if (!dir.exists(output_anim_dir)) {
#   dir.create(output_anim_dir, recursive = TRUE)
# }
# 
# year_dirs <- list.dirs(base_frame_dir, full.names = TRUE, recursive = FALSE)
# 
# if (length(year_dirs) == 0) {
#   stop("No year directories found in the specified base_frame_dir.")
# }
# 
# for (y_dir in year_dirs) {
#   year_name <- basename(y_dir)
#   png_files <- list.files(y_dir, pattern = "\\.png$", full.names = TRUE)
#   
#   if (length(png_files) == 0) {
#     cat(sprintf("No PNG frames found for %s. Skipping...\n", year_name))
#     next
#   }
#   
#   cat(sprintf("\nAnimating %d frames for the year %s...\n", length(png_files), year_name))
#   
#   if (output_format == "mp4") {
#     out_file <- file.path(output_anim_dir, sprintf("synoptic_dashboard_%s.mp4", year_name))
#     av_encode_video(
#       input = png_files, 
#       output = out_file, 
#       framerate = mp4_fps,
#       vfilter = "format=yuv420p"
#     )
#   } else if (output_format == "gif") {
#     out_file <- file.path(output_anim_dir, sprintf("synoptic_dashboard_%s.gif", year_name))
#     gifski(
#       png_files = png_files,
#       gif_file = out_file,
#       width = 2400,   
#       height = 1050,  
#       delay = gif_delay,    
#       progress = TRUE
#     )
#   }
#   
#   cat(sprintf("Successfully saved: %s\n", out_file))
# }
# 
# cat("\n========================================\n")
# cat("ALL ANIMATIONS COMPLETED!\n")
# cat("========================================\n")