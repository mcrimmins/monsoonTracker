# -------------------------------------------------------------------------
# 1. SETUP PARAMETERS & LOAD DATA
# -------------------------------------------------------------------------
library(terra)
library(ggplot2)   # NEW: Replaces base R plotting
library(tidyterra) # NEW: Allows ggplot to read terra rasters natively
library(maps)      # Required for map borders

slice_dir  <- "data/processed/4slices"
test_year  <- 2020
test_month <- "08"
test_day   <- 14 

# --- NEW VISUAL CONTROLS ---
show_trajectory  <- TRUE   # Set to FALSE to turn off the tail completely
tail_length_days <- 4      # How many days back to draw the tail

# --- DYNAMIC TITLE GENERATION ---
# Converts the raw inputs into a formal Date object to format cleanly
map_date <- as.Date(sprintf("%d-%s-%02d", test_year, test_month, test_day))
formatted_date <- format(map_date, "%d %b %Y") # e.g., "14 Aug 2020"

# Use \n to force the date onto the second line
map_title <- sprintf("NAM Synoptic Overview: 500hPa Heights & Moisture\n%s", formatted_date)

cat(sprintf("Rendering...\n%s\n", map_title))

nc_single <- list.files(slice_dir, pattern = sprintf("single_%d_%s\\.nc$", test_year, test_month), full.names = TRUE)[1]
r_single  <- rast(nc_single)
tcwv_day  <- r_single[[test_day]] 

nc_press <- list.files(slice_dir, pattern = sprintf("pressure_%d_%s\\.nc$", test_year, test_month), full.names = TRUE)[1]
r_press  <- rast(nc_press)
days_in_month <- nlyr(r_press) / 8

q850_day <- r_press[[test_day]] * 1000                     
q700_day <- r_press[[test_day + (1 * days_in_month)]] * 1000 
z500_day <- r_press[[test_day + (6 * days_in_month)]] / 9.80665 

# -------------------------------------------------------------------------
# 2. ADVANCED "H" TRACKING: The "Two-Pass" Macro-to-Micro Algorithm
# -------------------------------------------------------------------------
# Adjust tail tracking based on user toggles to save processing time
if (show_trajectory && tail_length_days > 0) {
  start_tail <- max(1, test_day - tail_length_days)
} else {
  start_tail <- test_day # Only compute the current day if the tail is off
}
tail_days <- start_tail:test_day

h_lon <- numeric(length(tail_days))
h_lat <- numeric(length(tail_days))

visual_smooth_mat <- matrix(1/25, nrow = 5, ncol = 5)    
tracker_smooth_mat <- matrix(1/225, nrow = 15, ncol = 15) 

# Open eastern boundary to track ridge progression
monsoon_box <- ext(-114, -95, 15, 45)

for (i in seq_along(tail_days)) {
  d <- tail_days[i]
  z_d <- r_press[[d + (6 * days_in_month)]] / 9.80665
  z_search <- crop(z_d, monsoon_box)
  
  # PASS 1: THE MACRO TRACKER 
  z_macro <- focal(z_search, w = tracker_smooth_mat, fun = mean, na.rm = TRUE)
  macro_max <- global(z_macro, "max", na.rm = TRUE)[1, 1]
  macro_core <- ifel(z_macro >= (macro_max - 10), 1, NA)
  
  blobs <- patches(macro_core)
  blob_sizes <- freq(blobs)
  largest_blob_id <- blob_sizes$value[which.max(blob_sizes$count)]
  primary_macro_core <- ifel(blobs == largest_blob_id, 1, NA)
  
  # PASS 2: THE MICRO PINPOINT 
  z_visual <- focal(z_search, w = visual_smooth_mat, fun = mean, na.rm = TRUE)
  z_isolated <- mask(z_visual, primary_macro_core)
  
  true_max <- global(z_isolated, "max", na.rm = TRUE)[1, 1]
  final_core <- ifel(z_isolated >= (true_max - 2), 1, NA)
  
  core_pts <- crds(as.points(final_core))
  centroid <- colMeans(core_pts)
  
  h_lon[i] <- centroid[1]
  h_lat[i] <- centroid[2]
}

# -------------------------------------------------------------------------
# 3. SPATIAL SMOOTHING FOR PLOTTING
# -------------------------------------------------------------------------
q850_smooth <- focal(q850_day, w = visual_smooth_mat, fun = mean, na.rm = TRUE)
q700_smooth <- focal(q700_day, w = visual_smooth_mat, fun = mean, na.rm = TRUE)
z500_smooth <- focal(z500_day, w = visual_smooth_mat, fun = mean, na.rm = TRUE)

# -------------------------------------------------------------------------
# 4. BUILD THE VISUAL FRAME (ggplot2 Implementation)
# -------------------------------------------------------------------------
out_file <- sprintf("frame_%d%s%02d.png", test_year, test_month, test_day)
tcwv_pal <- colorRampPalette(c("wheat", "lightgreen", "deepskyblue3", "darkblue"))(100)

# Convert trajectory points into a dataframe for ggplot
traj_df <- data.frame(lon = h_lon, lat = h_lat)
b <- ext(tcwv_day) # Get exact coordinate boundaries

# 1. Extract modern map layers to avoid the deprecated borders() warning
world_map <- map_data("world")
state_map <- map_data("state")

# 2. Initialize Plot with Raster, Borders, and Contours
p <- ggplot() +
  # Raster
  geom_spatraster(data = tcwv_day) +
  scale_fill_gradientn(
    colors = tcwv_pal,
    name = "TCWV\n(kg/m2)",
    na.value = "transparent"
  ) +
  
  # Modern Map Borders
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = NA, color = "black", linewidth = 0.8) +
  geom_polygon(data = state_map, aes(x = long, y = lat, group = group), 
               fill = NA, color = "black", linewidth = 0.5) +
  
  # THE FIX: Correct tidyterra contour functions
  geom_spatraster_contour(data = z500_smooth, aes(color = "Z500 Heights"), 
                          breaks = seq(5500, 6000, by = 20), linewidth = 0.4) +
  
  geom_spatraster_contour(data = q850_smooth, aes(color = "Q850 > 10 g/kg"), 
                          breaks = 10, linewidth = 1.2) +
  
  geom_spatraster_contour(data = q700_smooth, aes(color = "Q700 > 6 g/kg"), 
                          breaks = 6, linewidth = 1.2, linetype = "dashed")

# 3. Conditionally add trajectory tail and current High
if (show_trajectory && nrow(traj_df) > 1) {
  p <- p + 
    geom_path(data = traj_df, aes(x = lon, y = lat, color = "Primary Ridge (H)"), linewidth = 1) +
    geom_point(data = traj_df[-nrow(traj_df), ], aes(x = lon, y = lat, color = "Primary Ridge (H)"), size = 2)
}

# Always plot the current 'H' at the end of the trajectory list
p <- p + 
  annotate("text", x = traj_df$lon[nrow(traj_df)], y = traj_df$lat[nrow(traj_df)], 
           label = "H", color = "red", fontface = "bold", size = 8)

# 4. Apply formatting, legends, and absolute clipping walls
p <- p + 
  scale_color_manual(
    name = NULL, 
    values = c("Q850 > 10 g/kg" = "red", "Q700 > 6 g/kg" = "blue", 
               "Z500 Heights" = "gray30", "Primary Ridge (H)" = "darkred")
  ) +
  # The Absolute Bounding Box
  coord_sf(xlim = c(b[1], b[2]), ylim = c(b[3], b[4]), expand = FALSE) + 
  
  labs(title = map_title, x = "", y = "") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 10)),
    axis.text = element_text(color = "black", size = 10),
    
    # Bottom Left Legend Placement
    legend.position.inside = c(0.02, 0.02),
    legend.justification.inside = c(0, 0),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.margin = margin(5, 5, 5, 5)
  ) +
  guides(
    fill = guide_colorbar(order = 1, barheight = 15, title.vjust = 1), # TCWV Guide (Stays Right)
    color = guide_legend(order = 2, position = "inside")               # Contour Guide (Goes Inside)
  )

# 5. Save the plot
ggsave(out_file, plot = p, width = 9, height = 7, dpi = 150, bg = "white")

cat(sprintf("Successfully saved: %s\n", out_file))