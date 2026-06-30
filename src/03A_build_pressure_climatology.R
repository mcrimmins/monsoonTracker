# src/11_build_pressure_climatology.R
library(terra)
library(dplyr)
library(stringr)

# -------------------------------------------------------------------------
# 1. SETUP PATHS & VARIABLE DEFINITIONS
# -------------------------------------------------------------------------
proc_dir <- "data/processed/4slices"
months_to_process <- c("06", "07", "08", "09", "10") # June to October
month_names <- c("June", "July", "August", "September", "October")

# Define our 8 target variables and their corresponding block positions (1-8)
var_metadata <- tribble(
  ~var_id,  ~var_type, ~block_pos, ~unit_scale,
  "q850",   "q",        1,          1000,       # Convert to g/kg
  "q700",   "q",        2,          1000,
  "q500",   "q",        3,          1000,
  "q250",   "q",        4,          1000,
  "z850",   "z",        5,          1/9.80665,  # Convert to geopotential meters
  "z700",   "z",        6,          1/9.80665,
  "z500",   "z",        7,          1/9.80665,
  "z250",   "z",        8,          1/9.80665
)

# -------------------------------------------------------------------------
# 2. MASTER VARIABLE LOOP
# -------------------------------------------------------------------------
for (v in 1:nrow(var_metadata)) {
  target_var <- var_metadata$var_id[v]
  v_type     <- var_metadata$var_type[v]
  b_pos      <- var_metadata$block_pos[v]
  scale_factor <- var_metadata$unit_scale[v]
  
  cat(sprintf("\n==================================================\n"))
  cat(sprintf("PROCESSING VARIABLE: %s\n", toupper(target_var)))
  cat(sprintf("==================================================\n"))
  
  # Master list to hold the mean raster for each of the 5 months
  monthly_means_list <- list()
  
  # --- PASS 1: Collapse 35 Years per Month ---
  for (m in 1:length(months_to_process)) {
    m_num <- months_to_process[m]
    m_files <- list.files(proc_dir, pattern = paste0("pressure_.*_", m_num, "\\.nc$"), full.names = TRUE)
    
    cat(sprintf(" -> Collapsing 35 years for %s...\n", month_names[m]))
    
    r_list <- lapply(m_files, rast)
    
    # Dynamically determine the number of days in this specific month
    # Total layers / 8 variables = number of days
    num_days <- nlyr(r_list[[1]]) / 8
    
    # Initialize empty raster stack for this month's averages
    month_daily_avg <- rast(r_list[[1]][[1]], nlyrs = num_days)
    
    for (d in 1:num_days) {
      # Calculate the exact layer index for this variable on day 'd'
      # Example: For block 5 (z850), day 1 is layer (4 * num_days) + 1
      target_layer_idx <- ((b_pos - 1) * num_days) + d
      
      # Stack day 'd' across all 35 years
      day_stack <- rast(lapply(r_list, function(x) x[[target_layer_idx]]))
      
      # Compute the pixel mean and apply unit conversion scaling
      month_daily_avg[[d]] <- mean(day_stack, na.rm = TRUE) * scale_factor
    }
    
    monthly_means_list[[m_num]] <- month_daily_avg
  }
  
  # Combine months into the continuous 153-layer raw climatology stack
  climo_raw <- rast(monthly_means_list)
  
  # Generate and assign clean layer names (e.g., "June_01")
  days_per_month <- sapply(monthly_means_list, nlyr)
  all_layer_names <- unlist(lapply(1:5, function(i) {
    sprintf("%s_%02d", month_names[i], 1:days_per_month[i])
  }))
  names(climo_raw) <- all_layer_names
  
  # Save raw daily means checkpoint
  raw_out_path <- sprintf("data/processed/climatology_%s_raw.tif", target_var)
  writeRaster(climo_raw, raw_out_path, overwrite = TRUE)
  cat(sprintf(" -> Saved raw daily checkpoint to: %s\n", raw_out_path))
  
  # --- PASS 2: Apply 5-Day Rolling Smooth ---
  cat(" -> Applying 5-day rolling temporal smooth...\n")
  climo_smoothed <- climo_raw
  n_days <- nlyr(climo_raw)
  
  for (i in 1:n_days) {
    window_idx <- (i - 2):(i + 2)
    window_idx <- window_idx[window_idx >= 1 & window_idx <= n_days]
    
    climo_smoothed[[i]] <- mean(climo_raw[[window_idx]], na.rm = TRUE)
  }
  
  names(climo_smoothed) <- all_layer_names
  
  # Save final smoothed baseline
  smooth_out_path <- sprintf("data/processed/climatology_%s_5day_smoothed.tif", target_var)
  writeRaster(climo_smoothed, smooth_out_path, overwrite = TRUE)
  cat(sprintf(" -> Saved smoothed baseline to: %s\n", smooth_out_path))
}

cat("\n==================================================\n")
cat("SUCCESS: All 8 pressure-level baselines completed!\n")
cat("==================================================\n")

###############################################







###############################################
# src/12_inspect_pressure_climo.R
# library(terra)
# library(rnaturalearth)
# 
# # -------------------------------------------------------------------------
# # 1. LOAD TARGET FILES
# # -------------------------------------------------------------------------
# # Define paths to the smoothed 5-day climatologies we want to test
# z500_path <- "data/processed/climatology_z500_5day_smoothed.tif"
# q850_path <- "data/processed/climatology_q850_5day_smoothed.tif"
# q700_path <- "data/processed/climatology_q700_5day_smoothed.tif"
# 
# if (!all(file.exists(c(z500_path, q850_path, q700_path)))) {
#   stop("Error: One or more target smoothed files are missing. Check directory!")
# }
# 
# cat("Loading pressure climatology rasters...\n")
# z500 <- rast(z500_path)
# q850 <- rast(q850_path)
# q700 <- rast(q700_path)
# 
# # Extract Day 45 (July 15th)
# target_idx <- 90
# target_name <- names(z500)[target_idx]
# 
# z500_day <- z500[[target_idx]]
# q850_day <- q850[[target_idx]]
# q700_day <- q700[[target_idx]]
# 
# # -------------------------------------------------------------------------
# # 2. RUN STATISTICAL DIAGNOSTICS (Sanity checking the math)
# # -------------------------------------------------------------------------
# cat("\n=== DIAGNOSTIC REPORT:", target_name, "===\n")
# 
# # Extract min/max to ensure unit conversions worked (z in meters, q in g/kg)
# z500_range <- minmax(z500_day)
# q850_range <- minmax(q850_day)
# q700_range <- minmax(q700_day)
# 
# cat(sprintf("500 mb Heights (z):   %.0f to %.0f meters\n", z500_range[1,], z500_range[2,]))
# cat(sprintf("850 mb Moisture (q):  %.2f to %.2f g/kg\n", q850_range[1,], q850_range[2,]))
# cat(sprintf("700 mb Moisture (q):  %.2f to %.2f g/kg\n", q700_range[1,], q700_range[2,]))
# 
# # Calculate the exact coordinates of the highest 500 mb value (The Subtropical Ridge)
# max_cell_idx <- which.max(values(z500_day))
# high_coords <- xyFromCell(z500_day, max_cell_idx)
# 
# cat(sprintf("Subtropical High (H): Anchored at Lon %.2f, Lat %.2f\n", high_coords[1], high_coords[2]))
# cat("======================================\n\n")
# 
# # -------------------------------------------------------------------------
# # 3. FETCH GEOGRAPHIC BOUNDARIES
# # -------------------------------------------------------------------------
# cat("Downloading vector boundaries for mapping...\n")
# countries_sf <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
# countries_vect <- vect(countries_sf)
# 
# states_sf <- ne_states(country = c("United States of America", "Mexico"), returnclass = "sf")
# states_vect <- vect(states_sf)
# 
# # -------------------------------------------------------------------------
# # 4. PLOT PROTOTYPE SYNOPTIC MAP
# # -------------------------------------------------------------------------
# cat("Rendering prototype map... Check your RStudio Plots pane.\n")
# 
# # Create a heat palette for the 500 mb heights (steering flow)
# height_pal <- colorRampPalette(c("#ffffb2", "#fecc5c", "#fd8d3c", "#f03b20", "#bd0026"))(100)
# 
# # Set map layout and margins
# par(mar = c(3, 3, 4, 5))
# 
# # 1. Base Map: 500 mb Heights as a color grid
# plot(z500_day, 
#      main = sprintf("Monsoon Prototype Map: %s\n500mb Heights & Moisture Boundaries", target_name),
#      col = height_pal,
#      type = "continuous",
#      plg = list(title = "Meters"))
# 
# # 2. Add Geographic Boundaries
# plot(states_vect, add = TRUE, border = "gray40", lwd = 0.5, lty = "dotted")
# plot(countries_vect, add = TRUE, border = "black", lwd = 1.2)
# 
# # 3. Add 500 mb Height Contours (Draw lines every ~20 meters for steering flow context)
# contour(z500_day, levels = seq(5700, 6000, by = 20), col = "gray20", lwd = 0.5, add = TRUE)
# 
# # 4. Add the Tropical Moisture Isolines
# #   - 850 mb boundary (Shallow Surge): Thick Solid Red Line at 10 g/kg
# #   - 700 mb boundary (Terrain Clearer): Thick Dashed Blue Line at 6 g/kg
# contour(q850_day, levels = 10, col = "red", lwd = 3, add = TRUE, drawlabels = FALSE)
# contour(q700_day, levels = 6, col = "blue", lwd = 3, lty = 2, add = TRUE, drawlabels = FALSE)
# 
# # 5. Add the Subtropical High "H" Marker
# text(high_coords[1], high_coords[2], labels = "H", col = "blue", font = 2, cex = 1.8)
# 
# # Add a custom legend for the isolines
# legend("bottomleft", 
#        legend = c("850mb Moisture (10 g/kg)", "700mb Moisture (6 g/kg)", "Subtropical Ridge Center"),
#        col = c("red", "blue", "blue"), 
#        lwd = c(3, 3, NA), 
#        lty = c(1, 2, NA),
#        pch = c(NA, NA, "H"),
#        pt.cex = 1.5,
#        bg = "white",
#        cex = 0.8)
# 
# # Reset plotting parameters
# par(mfrow = c(1, 1))




