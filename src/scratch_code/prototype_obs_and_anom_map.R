# =========================================================================
# SYNOPTIC DASHBOARD: FULL OBSERVED FEATURES VS. ANOMALY (V2)
# =========================================================================

# -------------------------------------------------------------------------
# 1. SETUP PARAMETERS & LOAD DATA
# -------------------------------------------------------------------------
library(terra)
library(ggplot2)
library(tidyterra) # Provides geom_spatraster_contour and _text
library(maps)
library(patchwork)

slice_dir  <- "data/processed/4slices"
test_year  <- 1991
test_month <- "08"
test_day   <- 30 

# --- VISUAL CONTROLS ---
show_trajectory  <- TRUE   
tail_length_days <- 4      

# --- GLOBAL COLOR SCALES (Fixes Duplicated Legends) ---
# Defining this once guarantees patchwork merges the legends perfectly
feature_colors <- c(
  "Primary Ridge (H)" = "blue",
  "Q700 > 6 g/kg"     = "orange",
  "Q850 > 10 g/kg"    = "red1",  # Anchor line (pops against all backgrounds)
  "Z500 Heights"      = "gray30",
  "Z500 Anom (+)"     = "hotpink",    # Standard synoptic warm/ridge color
  "Z500 Anom (-)"     = "steelblue"    # Standard synoptic cool/trough color
)

# --- DYNAMIC TITLE ---
map_date <- as.Date(sprintf("%d-%s-%02d", test_year, test_month, test_day))
formatted_date <- format(map_date, "%d %b %Y")
dashboard_title <- sprintf("Synoptic Overview & Anomalies\n%s", formatted_date)

cat(sprintf("Rendering...\n%s\n", dashboard_title))

# --- LOAD DATA ---
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
# 2. CALCULATE ANOMALIES (Replace with actual climo logic)
# -------------------------------------------------------------------------
tcwv_anom <- tcwv_day - global(tcwv_day, "mean", na.rm = TRUE)[1, 1]
z500_anom <- z500_day - global(z500_day, "mean", na.rm = TRUE)[1, 1]

# -------------------------------------------------------------------------
# 3. ADVANCED "H" TRACKING: The "Two-Pass" Macro-to-Micro Algorithm
# -------------------------------------------------------------------------
if (show_trajectory && tail_length_days > 0) {
  start_tail <- max(1, test_day - tail_length_days)
} else {
  start_tail <- test_day 
}
tail_days <- start_tail:test_day

h_lon <- numeric(length(tail_days))
h_lat <- numeric(length(tail_days))

visual_smooth_mat <- matrix(1/25, nrow = 5, ncol = 5)    
tracker_smooth_mat <- matrix(1/225, nrow = 15, ncol = 15) 

monsoon_box <- ext(-114, -95, 15, 45)

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
# 4. SPATIAL SMOOTHING FOR PLOTTING
# -------------------------------------------------------------------------
q850_smooth <- focal(q850_day, w = visual_smooth_mat, fun = mean, na.rm = TRUE)
q700_smooth <- focal(q700_day, w = visual_smooth_mat, fun = mean, na.rm = TRUE)
z500_smooth <- focal(z500_day, w = visual_smooth_mat, fun = mean, na.rm = TRUE)
z500_anom_smooth <- focal(z500_anom, w = visual_smooth_mat, fun = mean, na.rm = TRUE)

# -------------------------------------------------------------------------
# 5. BUILD PLOT A: OBSERVED (Full Features)
# -------------------------------------------------------------------------
b <- ext(tcwv_day) 
world_map <- map_data("world")
state_map <- map_data("state")

obs_pal <- colorRampPalette(c("wheat", "lightgreen", "deepskyblue3", "darkblue"))(100)

p_obs <- ggplot() +
  geom_spatraster(data = tcwv_day) +
  scale_fill_gradientn(colors = obs_pal, name = "Observed TCWV\n(kg/m2)", 
                       na.value = "transparent", limits = c(5, 65), oob = scales::squish) +
  
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.8) +
  geom_polygon(data = state_map, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.5) +
  
  geom_spatraster_contour(data = z500_smooth, aes(color = "Z500 Heights"), breaks = seq(5500, 6000, by = 20), linewidth = 0.4) +
  geom_spatraster_contour(data = q700_smooth, aes(color = "Q700 > 6 g/kg"), breaks = 6, linewidth = 1, linetype = "dashed") +
  # Anchor Line is now Black
  geom_spatraster_contour(data = q850_smooth, aes(color = "Q850 > 10 g/kg"), breaks = 10, linewidth = 1) 

# Add Trajectory Tail
if (show_trajectory && nrow(traj_df) > 1) {
  p_obs <- p_obs + 
    geom_path(data = traj_df, aes(x = lon, y = lat, color = "Primary Ridge (H)"), linewidth = 1) +
    geom_point(data = traj_df[-nrow(traj_df), ], aes(x = lon, y = lat, color = "Primary Ridge (H)"), size = 2)
}

# Add Current 'H' and Final Touches
p_obs <- p_obs + 
  annotate("text", x = traj_df$lon[nrow(traj_df)], y = traj_df$lat[nrow(traj_df)], label = "H", color = "blue2", fontface = "bold", size = 8) +
  # Use the global feature vector to force legend merging
  scale_color_manual(name = "Contours & Features", values = feature_colors, drop = FALSE) +
  coord_sf(xlim = c(b[1], b[2]), ylim = c(b[3], b[4]), expand = FALSE) + 
  labs(title = "Observed", x = "", y = "") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5), axis.text = element_text(color = "black", size = 10))

# -------------------------------------------------------------------------
# 6. BUILD PLOT B: ANOMALY
# -------------------------------------------------------------------------
p_anom <- ggplot() +
  geom_spatraster(data = tcwv_anom) +
  scale_fill_gradient2(low = "saddlebrown", mid = "white", high = "forestgreen", midpoint = 0, 
                       name = "TCWV Anom\n(kg/m2)", limits = c(-15, 15), oob = scales::squish, na.value = "transparent") +
  
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.8) +
  geom_polygon(data = state_map, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.5) +
  
  # Z500 Anomalies (Now Red/Blue)
  geom_spatraster_contour(data = z500_anom_smooth, aes(color = "Z500 Anom (-)"), breaks = seq(-300, -20, by = 20), linewidth = 0.8, linetype = "dashed") +
  geom_spatraster_contour(data = z500_anom_smooth, aes(color = "Z500 Anom (+)"), breaks = seq(20, 300, by = 20), linewidth = 0.8, linetype = "solid") +
  
  # Add inline contour text labels for anomalies
  geom_spatraster_contour_text(data = z500_anom_smooth, breaks = seq(-300, 300, by = 20), 
                               color = "grey20", size = 3, fontface = "bold") +
  
  # Q850 Anchor Line (Matches Observed exactly)
  geom_spatraster_contour(data = q850_smooth, aes(color = "Q850 > 10 g/kg"), breaks = 10, linewidth = 1, linetype = "solid") +
  
  # Use the global feature vector to force legend merging
  scale_color_manual(name = "Contours & Features", values = feature_colors, drop = FALSE) +
  coord_sf(xlim = c(b[1], b[2]), ylim = c(b[3], b[4]), expand = FALSE) + 
  labs(title = "Anomaly", x = "", y = "") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_blank(), axis.ticks.y = element_blank() # Strip redundant Y-axis
  )

# -------------------------------------------------------------------------
# 7. COMBINE AND SAVE (PATCHWORK)
# -------------------------------------------------------------------------
out_file <- sprintf("dashboard_%d%s%02d.png", test_year, test_month, test_day)

# Stitch together and generate master legend
dashboard <- p_obs + p_anom +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = dashboard_title,
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 15)))
  )

ggsave(out_file, plot = dashboard, width = 16, height = 7, dpi = 150, bg = "white")
cat(sprintf("Successfully saved dashboard: %s\n", out_file))