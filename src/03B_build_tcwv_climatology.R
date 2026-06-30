# src/03_build_tcwv_climatology.R
library(terra)
library(dplyr)
library(stringr)

# -------------------------------------------------------------------------
# 1. SET UP PATHS & CONSTANTS
# -------------------------------------------------------------------------
proc_dir <- "data/processed/4slices"
months_to_process <- c("06", "07", "08", "09", "10") # June to October

# Target output files
climo_raw_path <- "data/processed/climatology_tcwv_raw.tif"
climo_smoothed_path <- "data/processed/climatology_tcwv_5day_smoothed.tif"

# -------------------------------------------------------------------------
# 2. PASS 1: COMPUTE MULTI-YEAR MEAN FOR EACH CALENDAR DAY
# -------------------------------------------------------------------------
cat("Starting Pass 1: Collapsing 35 years into daily means...\n")

monthly_means_list <- list()

for (m in months_to_process) {
  m_files <- list.files(proc_dir, pattern = paste0("single_.*_", m, "\\.nc$"), full.names = TRUE)
  
  cat(sprintf("Processing month %s (%d years found)...\n", m, length(m_files)))
  
  r_list <- lapply(m_files, rast)
  num_days <- nlyr(r_list[[1]]) 
  
  # Initialize empty raster stack for this month
  month_daily_avg <- rast(r_list[[1]], nlyrs = num_days)
  
  for (d in 1:num_days) {
    # FIX: Use rast() instead of sds() to stack the 35 layers cleanly
    day_stack <- rast(lapply(r_list, function(x) x[[d]]))
    
    # This now correctly computes the pixel-by-pixel average across years
    month_daily_avg[[d]] <- mean(day_stack, na.rm = TRUE)
  }
  
  monthly_means_list[[m]] <- month_daily_avg
}

# Combine all months into one continuous 153-layer raw climatology stack
climo_raw <- rast(monthly_means_list)

# FIX: Set layer names explicitly on the final combined stack so they stick
month_names <- c("June", "July", "August", "September", "October")
days_per_month <- sapply(monthly_means_list, nlyr)
all_layer_names <- unlist(lapply(1:5, function(i) {
  sprintf("%s_%02d", month_names[i], 1:days_per_month[i])
}))
names(climo_raw) <- all_layer_names

# Save the raw daily means checkpoint
writeRaster(climo_raw, climo_raw_path, overwrite = TRUE)

# -------------------------------------------------------------------------
# 3. PASS 2: APPLY 5-DAY ROLLING SMOOTH ALONG THE TIME AXIS
# -------------------------------------------------------------------------
cat("\nStarting Pass 2: Applying 5-day rolling mean to smooth daily noise...\n")

climo_smoothed <- climo_raw
n_days <- nlyr(climo_raw)

for (i in 1:n_days) {
  window_idx <- (i - 2):(i + 2)
  window_idx <- window_idx[window_idx >= 1 & window_idx <= n_days]
  
  climo_smoothed[[i]] <- mean(climo_raw[[window_idx]], na.rm = TRUE)
}

# Keep the clean names on the smoothed dataset as well
names(climo_smoothed) <- all_layer_names

writeRaster(climo_smoothed, climo_smoothed_path, overwrite = TRUE)
cat("Success! 5-day smoothed climatology saved to:", climo_smoothed_path, "\n")

##############################################################################



###### testing output --------################################################
# src/10_inspect_climatology.R
# library(terra)
# library(rnaturalearth)
# 
# # -------------------------------------------------------------------------
# # 1. LOAD AND SANITY-CHECK THE OUTPUT FILES
# # -------------------------------------------------------------------------
# climo_raw_path <- "data/processed/climatology_tcwv_raw.tif"
# climo_smoothed_path <- "data/processed/climatology_tcwv_5day_smoothed.tif"
# 
# if (!file.exists(climo_raw_path) || !file.exists(climo_smoothed_path)) {
#   stop("Error: One or both climatology files are missing. Check your previous script output!")
# }
# 
# cat("Loading climatology rasters...\n")
# climo_raw <- rast(climo_raw_path)
# climo_smooth <- rast(climo_smoothed_path)
# 
# # Print structural diagnostics to console
# cat("\n=== DIAGNOSTIC REPORT ===\n")
# cat(sprintf("Raw Climatology Layers:      %d (Expected: 153)\n", nlyr(climo_raw)))
# cat(sprintf("Smoothed Climatology Layers: %d (Expected: 153)\n", nlyr(climo_smooth)))
# cat(sprintf("Spatial Resolution:          %.3f x %.3f degrees\n", res(climo_raw)[1], res(climo_raw)[2]))
# cat("Coordinate Reference System: ", crs(climo_raw, describe = TRUE)$name, "\n")
# 
# # Check value ranges to ensure data isn't corrupted/empty
# raw_range <- minmax(climo_raw)
# cat(sprintf("TCWV Value Range (Raw):      %.2f to %.2f kg/m²\n", min(raw_range[1,]), max(raw_range[2,])))
# cat("=========================\n\n")
# 
# # -------------------------------------------------------------------------
# # 2. FETCH GEOGRAPHIC BOUNDARIES
# # -------------------------------------------------------------------------
# cat("Downloading vector boundaries for mapping...\n")
# countries_sf <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
# countries_vect <- vect(countries_sf)
# 
# states_sf <- ne_states(country = c("United States of America", "Mexico"), returnclass = "sf")
# states_vect <- vect(states_sf)
# 
# # -------------------------------------------------------------------------
# # 3. GRAPHICAL COMPARISON: RAW VS. SMOOTHED (July 15th)
# # -------------------------------------------------------------------------
# # July 15th is Day 45 of our June-October timeline (30 days of June + 15 days of July)
# target_day_idx <- 45 
# target_day_name <- names(climo_raw)[target_day_idx]
# 
# cat(sprintf("Plotting spatial comparison for: %s (Layer Index: %d)...\n", target_day_name, target_day_idx))
# 
# # Set up side-by-side plotting area
# par(mfrow = c(1, 2), mar = c(3, 3, 4, 5))
# 
# # Define a smooth, intuitive moisture color palette (Yellow to Green to Deep Blue)
# moisture_pal <- colorRampPalette(c("#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8", "#253494"))(100)
# 
# # Target a consistent color bar scale based on typical mid-summer PW values
# z_limits <- c(10, 55) 
# 
# # --- Plot 1: Raw 35-Year Mean ---
# plot(climo_raw[[target_day_idx]], 
#      main = sprintf("Raw Climatology Mean\n%s", target_day_name),
#      col = moisture_pal,
#      type = "continuous",
#      range = z_limits,
#      plg = list(title = "kg/m²"))
# 
# # Overlay boundaries
# plot(states_vect, add = TRUE, border = "gray30", lwd = 0.5, lty = "dotted")
# plot(countries_vect, add = TRUE, border = "black", lwd = 1.2)
# 
# # --- Plot 2: 5-Day Smoothed Climatology ---
# plot(climo_smooth[[target_day_idx]], 
#      main = sprintf("5-Day Smoothed Climatology\n%s", target_day_name),
#      col = moisture_pal,
#      type = "continuous",
#      range = z_limits,
#      plg = list(title = "kg/m²"))
# 
# # Overlay boundaries
# plot(states_vect, add = TRUE, border = "gray30", lwd = 0.5, lty = "dotted")
# plot(countries_vect, add = TRUE, border = "black", lwd = 1.2)
# 
# # Reset plotting parameters
# par(mfrow = c(1, 1))
# cat("Inspection complete! Check your RStudio Plots pane.\n")


