# -------------------------------------------------------------------------
# 1. SETUP & LOAD DATA
# -------------------------------------------------------------------------
library(terra)
library(ggplot2)
library(tidyterra)
library(maps)

slice_dir  <- "data/processed/4slices"
test_year  <- 2025
test_month <- "07"
test_day   <- 14 

map_date <- as.Date(sprintf("%d-%s-%02d", test_year, test_month, test_day))
map_title <- sprintf("Synoptic Anomaly\n%s", format(map_date, "%d %b %Y"))

# Load Daily Data
nc_single <- list.files(slice_dir, pattern = sprintf("single_%d_%s\\.nc$", test_year, test_month), full.names = TRUE)[1]
r_single  <- rast(nc_single)
tcwv_day  <- r_single[[test_day]] 

nc_press <- list.files(slice_dir, pattern = sprintf("pressure_%d_%s\\.nc$", test_year, test_month), full.names = TRUE)[1]
r_press  <- rast(nc_press)
days_in_month <- nlyr(r_press) / 8

q850_day <- r_press[[test_day]] * 1000                     
z500_day <- r_press[[test_day + (6 * days_in_month)]] / 9.80665 

# -------------------------------------------------------------------------
# 2. CALCULATE ANOMALIES (Mock logic updated to actually create gradients)
# -------------------------------------------------------------------------
# NOTE: Replace this with your real climatology subtractions!
# Here, I am subtracting the spatial mean of the daily map from itself to 
# force a clean gradient of positive and negative anomalies for the test plot.
tcwv_anom <- tcwv_day - global(tcwv_day, "mean", na.rm = TRUE)[1, 1]
z500_anom <- z500_day - global(z500_day, "mean", na.rm = TRUE)[1, 1]

# -------------------------------------------------------------------------
# 3. SPATIAL SMOOTHING
# -------------------------------------------------------------------------
visual_smooth_mat <- matrix(1/25, nrow = 5, ncol = 5)    

q850_day_smooth <- focal(q850_day, w = visual_smooth_mat, fun = mean, na.rm = TRUE)
z500_anom_smooth <- focal(z500_anom, w = visual_smooth_mat, fun = mean, na.rm = TRUE)

# -------------------------------------------------------------------------
# 4. BUILD THE ANOMALY FRAME
# -------------------------------------------------------------------------
out_file <- sprintf("clean_anom_%d%s%02d.png", test_year, test_month, test_day)
b <- ext(tcwv_day)

world_map <- map_data("world")
state_map <- map_data("state")

p <- ggplot() +
  # 1. Base Raster: TCWV Anomaly (Diverging Scale)
  geom_spatraster(data = tcwv_anom) +
  scale_fill_gradient2(
    low = "saddlebrown", 
    mid = "white", 
    high = "forestgreen", 
    midpoint = 0, # Forces 0 anomaly to be white
    name = "TCWV Anom\n(kg/m2)",
    limits = c(-15, 15),  # Adjust based on your real data range
    oob = scales::squish,
    na.value = "transparent"
  ) +
  
  # 2. Map Borders
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = NA, color = "black", linewidth = 0.8) +
  geom_polygon(data = state_map, aes(x = long, y = lat, group = group), 
               fill = NA, color = "black", linewidth = 0.5) +
  
  # 3. Z500 Anomalies (Cleanly separated by sign)
  # Negative Heights (Troughs) = Blue Dashed
  geom_spatraster_contour(data = z500_anom_smooth, aes(color = "Z500 Anom (-)"), 
                          breaks = seq(-300, -20, by = 20), linewidth = 0.8, linetype = "dashed") +
  # Positive Heights (Ridges) = Red Solid
  geom_spatraster_contour(data = z500_anom_smooth, aes(color = "Z500 Anom (+)"), 
                          breaks = seq(20, 300, by = 20), linewidth = 0.8, linetype = "solid") +
  
  # 4. The Single Anchor: Daily Moisture Boundary
  geom_spatraster_contour(data = q850_day_smooth, aes(color = "Obs Q850 > 10"), 
                          breaks = 10, linewidth = 1.5, linetype = "solid") +
  
  # 5. Legend & Theme Setup
  scale_color_manual(
    name = "Synoptic Features", 
    values = c(
      "Obs Q850 > 10" = "black", 
      "Z500 Anom (+)" = "darkred", 
      "Z500 Anom (-)" = "darkblue"
    )
  ) +
  coord_sf(xlim = c(b[1], b[2]), ylim = c(b[3], b[4]), expand = FALSE) + 
  
  labs(title = map_title, x = "", y = "") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 10)),
    axis.text = element_text(color = "black", size = 10),
    legend.position.inside = c(0.02, 0.02),
    legend.justification.inside = c(0, 0),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.margin = margin(5, 5, 5, 5)
  ) +
  guides(
    fill = guide_colorbar(order = 1, barheight = 15, title.vjust = 1),
    color = guide_legend(order = 2, position = "inside") 
  )

# 6. Save Plot
ggsave(out_file, plot = p, width = 9, height = 7, dpi = 150, bg = "white")
cat(sprintf("Successfully saved clean prototype: %s\n", out_file))