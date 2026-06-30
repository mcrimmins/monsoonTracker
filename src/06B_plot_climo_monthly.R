# =========================================================================
# CLIMATOLOGY GRID: MEAN MONTHLY MONSOON EVOLUTION (153 Daily Layers)
# =========================================================================

library(terra)
library(ggplot2)
library(tidyterra)
library(maps)

# -------------------------------------------------------------------------
# 1. SETUP PARAMETERS & MONTHLY AGGREGATION
# -------------------------------------------------------------------------
climo_dir <- "data/processed/4slices/climo"

# 153 Layers: Jun 1 to Oct 31.
# Map the exact daily layer indices to their respective calendar months.
monthly_indices <- list(
  "June"      = 1:30,     # 30 days
  "July"      = 31:61,    # 31 days
  "August"    = 62:92,    # 31 days
  "September" = 93:122,   # 30 days
  "October"   = 123:153   # 31 days
)
layer_names <- names(monthly_indices)

# --- GLOBAL COLOR SCALES ---
feature_colors <- c(
  "Primary Ridge (H)"         = "darkred",
  "TCWV > 25 kg/m2"           = "orange",
  "Q700 > 6 g/kg"             = "darkblue",
  "Z500 Heights"              = "gray30"
)
obs_pal <- colorRampPalette(c("wheat", "lightgreen", "deepskyblue3", "darkblue"))(100)

# -------------------------------------------------------------------------
# 2. LOAD, AGGREGATE & SMOOTH DATA
# -------------------------------------------------------------------------
cat("Loading daily climo data and aggregating into monthly means...\n")

tcwv_climo <- rast(file.path(climo_dir, "climatology_tcwv_raw.tif"))
z500_climo <- rast(file.path(climo_dir, "climatology_z500_raw.tif"))
q700_climo <- rast(file.path(climo_dir, "climatology_q700_raw.tif"))

# Calculate the mean across the exact daily layers for each month
tcwv_monthly <- rast(lapply(monthly_indices, function(idx) mean(tcwv_climo[[idx]], na.rm = TRUE)))
z500_monthly <- rast(lapply(monthly_indices, function(idx) mean(z500_climo[[idx]], na.rm = TRUE)))
q700_monthly <- rast(lapply(monthly_indices, function(idx) mean(q700_climo[[idx]], na.rm = TRUE)))

names(tcwv_monthly) <- layer_names
names(z500_monthly) <- layer_names
names(q700_monthly) <- layer_names

# Spatial Smoothing for Q700 & TCWV contouring
spatial_smooth_mat <- matrix(1/25, nrow = 5, ncol = 5)
q700_monthly <- focal(q700_monthly, w = spatial_smooth_mat, fun = mean, na.rm = TRUE)
tcwv_smooth_for_contour <- focal(tcwv_monthly, w = spatial_smooth_mat, fun = mean, na.rm = TRUE)

# -------------------------------------------------------------------------
# 3. ADVANCED "H" TRACKING & PLOTTING
# -------------------------------------------------------------------------
h_lon <- numeric(nlyr(z500_monthly))
h_lat <- numeric(nlyr(z500_monthly))
tracker_smooth_mat <- matrix(1/225, nrow = 15, ncol = 15) 
visual_smooth_mat  <- matrix(1/25, nrow = 5, ncol = 5) 
monsoon_box <- ext(-114, -95, 15, 45)

for (i in 1:nlyr(z500_monthly)) {
  z_search <- crop(z500_monthly[[i]], monsoon_box)
  z_macro <- focal(z_search, w = tracker_smooth_mat, fun = mean, na.rm = TRUE)
  macro_max <- global(z_macro, "max", na.rm = TRUE)[1, 1]
  macro_core <- ifel(z_macro >= (macro_max - 10), 1, NA)
  
  blobs <- patches(macro_core)
  blob_sizes <- freq(blobs)
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
  start_idx <- max(1, i - 1) 
  if (i - start_idx > 0) {
    path_list[[i]] <- data.frame(lon = h_lon[start_idx:i], lat = h_lat[start_idx:i], lyr = layer_names[i])
  }
}
h_paths <- do.call(rbind, path_list)

h_pts$lyr <- factor(h_pts$lyr, levels = layer_names)
if(!is.null(h_paths)) h_paths$lyr <- factor(h_paths$lyr, levels = layer_names)

# --- GGPLOT RENDER ---
world_map <- map_data("world")
state_map <- map_data("state")
b <- ext(tcwv_monthly) 

p_grid <- ggplot() +
  geom_spatraster(data = tcwv_monthly) +
  scale_fill_gradientn(colors = obs_pal, name = "Mean TCWV\n(kg/m2)", na.value = "transparent", limits = c(5, 65), oob = scales::squish) +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.5) +
  geom_polygon(data = state_map, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.3) +
  
  geom_spatraster_contour(data = z500_monthly, aes(color = "Z500 Heights"), breaks = seq(5500, 6000, by = 20), linewidth = 0.3) +
  geom_spatraster_contour(data = q700_monthly, aes(color = "Q700 > 6 g/kg"), breaks = c(6), linewidth = 0.6, linetype = "dashed") +
  
  # New TCWV 25 contour (using a slightly smoothed version so the line isn't too jagged)
  geom_spatraster_contour(data = tcwv_smooth_for_contour, aes(color = "TCWV > 25 kg/m2"), breaks = c(25), linewidth = 0.8, linetype = "solid") +
  
  #geom_path(data = h_paths, aes(x = lon, y = lat, color = "Primary Ridge (H)"), linewidth = 0.8) +
  geom_text(data = h_pts, aes(x = lon, y = lat, label = "H"), color = "blue1", fontface = "bold", size = 6) +
  scale_color_manual(name = "Features", values = feature_colors) +
  coord_sf(xlim = c(b[1], b[2]), ylim = c(b[3], b[4]), expand = FALSE) + 
  facet_wrap(~ factor(lyr, levels = layer_names), ncol = 5) + 
  labs(title = "North American Monsoon: Monthly Mean Evolution (1991-2025)", x = "", y = "") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5, margin = margin(b = 20)),
    strip.text = element_text(face = "bold", size = 14, background = element_rect(fill = "gray90", color = NA)),
    axis.text = element_blank(), axis.ticks = element_blank(), panel.spacing = unit(0.5, "lines"),
    legend.position = "bottom", legend.key.width = unit(2, "cm")
  )

ggsave("output/climatology_monthly_1x5.png", plot = p_grid, width = 20, height = 6, dpi = 300, bg = "white")
cat("Successfully saved monthly climatology grid to: output/climatology_monthly_1x5.png\n")