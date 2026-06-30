# =========================================================================
# CLIMATOLOGY GRID: 5x5 MONSOON EVOLUTION (153 Daily Layers)
# =========================================================================

library(terra)
library(ggplot2)
library(tidyterra)
library(maps)

# -------------------------------------------------------------------------
# 1. SETUP PARAMETERS & 5-DAY AGGREGATION
# -------------------------------------------------------------------------
climo_dir <- "data/processed/4slices/climo"

# We have 153 days (Jun 1 to Oct 31). A 5x5 grid holds 25 plots.
# To match your timeline, we generate 25 consecutive 5-day blocks 
# starting on June 10th (Layer 10).
start_day <- 10 
blocks_5day <- lapply(0:24, function(i) (start_day + i * 5):(start_day + i * 5 + 4))

# Generate accurate Date strings for the facets
base_date <- as.Date("2021-06-01") # Non-leap year base
start_dates <- base_date + sapply(blocks_5day, min) - 1
end_dates   <- base_date + sapply(blocks_5day, max) - 1
layer_names <- sprintf("%s - %s", format(start_dates, "%d %b"), format(end_dates, "%d %b"))

# --- GLOBAL COLOR SCALES ---
feature_colors <- c(
  "Primary Ridge (H)" = "darkred",
  "TCWV > 25 kg/m2"   = "orange",
  "Q700 > 6 g/kg"     = "darkblue",
  "Z500 Heights"      = "gray30"
)
obs_pal <- colorRampPalette(c("wheat", "lightgreen", "deepskyblue3", "darkblue"))(100)

# -------------------------------------------------------------------------
# 2. LOAD, AGGREGATE & SMOOTH DATA
# -------------------------------------------------------------------------
cat("Loading daily data and aggregating into 25 5-day blocks...\n")
tcwv_climo <- rast(file.path(climo_dir, "climatology_tcwv_raw.tif"))
z500_climo <- rast(file.path(climo_dir, "climatology_z500_raw.tif"))
q700_climo <- rast(file.path(climo_dir, "climatology_q700_raw.tif"))

# Calculate means on the fly for each 5-day block
tcwv_sub <- rast(lapply(blocks_5day, function(idx) mean(tcwv_climo[[idx]], na.rm = TRUE)))
z500_sub <- rast(lapply(blocks_5day, function(idx) mean(z500_climo[[idx]], na.rm = TRUE)))
q700_sub <- rast(lapply(blocks_5day, function(idx) mean(q700_climo[[idx]], na.rm = TRUE)))

names(tcwv_sub) <- layer_names
names(z500_sub) <- layer_names
names(q700_sub) <- layer_names

# Spatial Smoothing for Q700 & TCWV contouring
spatial_smooth_mat <- matrix(1/25, nrow = 5, ncol = 5)
q700_sub <- focal(q700_sub, w = spatial_smooth_mat, fun = mean, na.rm = TRUE)
tcwv_smooth_for_contour <- focal(tcwv_sub, w = spatial_smooth_mat, fun = mean, na.rm = TRUE)

# -------------------------------------------------------------------------
# 3. ADVANCED "H" TRACKING & PLOTTING
# -------------------------------------------------------------------------
h_lon <- numeric(nlyr(z500_sub)); h_lat <- numeric(nlyr(z500_sub))
tracker_smooth_mat <- matrix(1/225, nrow = 15, ncol = 15); visual_smooth_mat <- matrix(1/25, nrow = 5, ncol = 5) 
monsoon_box <- ext(-114, -95, 15, 45)

for (i in 1:nlyr(z500_sub)) {
  z_search <- crop(z500_sub[[i]], monsoon_box)
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

h_pts <- data.frame(lon = h_lon, lat = h_lat, lyr = layer_names)
path_list <- list()
for (i in 1:length(layer_names)) {
  start_idx <- max(1, i - 2) 
  if (i - start_idx > 0) {
    path_list[[i]] <- data.frame(lon = h_lon[start_idx:i], lat = h_lat[start_idx:i], lyr = layer_names[i])
  }
}
h_paths <- do.call(rbind, path_list)

h_pts$lyr <- factor(h_pts$lyr, levels = layer_names)
if(!is.null(h_paths)) h_paths$lyr <- factor(h_paths$lyr, levels = layer_names)

# --- GGPLOT RENDER ---
world_map <- map_data("world"); state_map <- map_data("state"); b <- ext(tcwv_sub) 

p_grid <- ggplot() +
  geom_spatraster(data = tcwv_sub) +
  scale_fill_gradientn(colors = obs_pal, name = "Mean TCWV\n(kg/m2)", na.value = "transparent", limits = c(5, 65), oob = scales::squish) +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.5) +
  geom_polygon(data = state_map, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.3) +
  
  geom_spatraster_contour(data = z500_sub, aes(color = "Z500 Heights"), breaks = seq(5500, 6000, by = 20), linewidth = 0.3) +
  geom_spatraster_contour(data = q700_sub, aes(color = "Q700 > 6 g/kg"), breaks = c(6), linewidth = 0.7, linetype = "dashed") +
  
  # New smoothed TCWV 25 contour
  geom_spatraster_contour(data = tcwv_smooth_for_contour, aes(color = "TCWV > 25 kg/m2"), breaks = c(25), linewidth = 0.8, linetype = "solid") +
  
  #geom_path(data = h_paths, aes(x = lon, y = lat, color = "Primary Ridge (H)"), linewidth = 0.8) +
  geom_text(data = h_pts, aes(x = lon, y = lat, label = "H"), color = "blue1", fontface = "bold", size = 5) +
  scale_color_manual(name = "Features", values = feature_colors) +
  coord_sf(xlim = c(b[1], b[2]), ylim = c(b[3], b[4]), expand = FALSE) + 
  facet_wrap(~ factor(lyr, levels = layer_names), ncol = 5) + 
  labs(title = "North American Monsoon: Climatological Evolution (1991-2025)", x = "", y = "") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5, margin = margin(b = 20)),
    strip.text = element_text(face = "bold", size = 10, background = element_rect(fill = "gray90", color = NA)),
    axis.text = element_blank(), axis.ticks = element_blank(), panel.spacing = unit(0.5, "lines"),
    legend.position = "right", legend.key.height = unit(1.5, "cm")
  )

ggsave("output/climatology_monsoon_evolution_5x5.png", plot = p_grid, width = 16, height = 12, dpi = 300, bg = "white")
cat("Successfully saved sorted 5x5 climatology grid to: output/climatology_monsoon_evolution_5x5.png\n")