# src/05_examine_processed_data.R
library(terra)

# -------------------------------------------------------------------------
# 1. LOAD THE PROCESSED DATA
# -------------------------------------------------------------------------
# Point to the freshly processed June 1991 single-level file (PWAT)
file_sl <- "data/processed/4slices/daily_mean_4slices_era5_single_1991_06.nc"

cat("Loading processed data...\n")
pwat_daily <- rast(file_sl)

# -------------------------------------------------------------------------
# 2. INSPECT THE STRUCTURE
# -------------------------------------------------------------------------
cat("\n--- PROCESSED FILE METADATA ---\n")
print(pwat_daily)

cat("\n--- FIRST 5 LAYER NAMES ---\n")
print(head(names(pwat_daily), 5))

# -------------------------------------------------------------------------
# 3. VISUALIZATION
# -------------------------------------------------------------------------
# Set up a plotting layout: 1 row, 2 columns
par(mfrow = c(1, 2), mar = c(3, 3, 2, 5))

# Define a nice moisture color palette (Dry/Brown to Wet/Blue/Green)
moisture_pal <- colorRampPalette(c("#e5d5c1", "#a3c4dc", "#3288bd", "#5e4fa2"))(100)

# --- Plot A: A Single Specific Day ---
# Let's look at mid-month: June 15, 1991
day_to_plot <- 15
plot(pwat_daily[[day_to_plot]], 
     main = sprintf("Daily Mean TCWV (4-slices)\nJune %s, 1991", day_to_plot),
     col = moisture_pal,
     type = "continuous",
     plg = list(title = "kg/m²"))

# --- Plot B: The Entire Monthly Climatological Mean ---
# We can average all our daily averages to see the baseline for the whole month
cat("\nCalculating monthly mean for visualization...\n")
pwat_monthly_mean <- mean(pwat_daily)

plot(pwat_monthly_mean, 
     main = "Monthly Climatological Mean TCWV\nJune 1991",
     col = moisture_pal,
     type = "continuous",
     plg = list(title = "kg/m²"))

# Reset plot parameters
par(mfrow = c(1, 1))



##################
# src/07_plot_with_boundaries.R
library(terra)
library(rnaturalearth)

# -------------------------------------------------------------------------
# 1. LOAD THE PRESSURE LEVEL DATA
# -------------------------------------------------------------------------
file_pl <- "data/processed/4slices/daily_mean_4slices_era5_pressure_2025_06.nc"

cat("Loading 3D pressure level data...\n")
pl_daily <- rast(file_pl)

# -------------------------------------------------------------------------
# 2. EXTRACT SPECIFIC VARIABLES & LEVELS (BY INDEX)
# -------------------------------------------------------------------------
# Extract Specific Humidity (q) at 850 hPa (Block 1: Layers 1 to 30)
q_850_stack <- pl_daily[[1:30]]

# Extract Geopotential (z) at 500 hPa (Block 7: Layers 181 to 210)
z_500_stack <- pl_daily[[181:210]]

# -------------------------------------------------------------------------
# 3. DOWNLOAD GEOGRAPHIC BOUNDARIES
# -------------------------------------------------------------------------
cat("Fetching national and state boundaries...\n")

# Pull national borders (North America)
countries_sf <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
countries_vect <- vect(countries_sf) # Convert to terra's native SpatVector

# Pull state/province borders (US and Mexico)
states_sf <- ne_states(country = c("United States of America", "Mexico"), returnclass = "sf")
states_vect <- vect(states_sf)

# -------------------------------------------------------------------------
# 4. VISUALIZATION WITH BORDERS
# -------------------------------------------------------------------------
# Set up a plotting layout: 1 row, 2 columns
par(mfrow = c(1, 2), mar = c(3, 3, 3, 5))
day_to_plot <- 15

# --- Plot A: Geopotential Height (z) at 500 hPa ---
z_pal <- colorRampPalette(c("#313695", "#74add1", "#ffffbf", "#f46d43", "#a50026"))(100)
z_500_day <- z_500_stack[[day_to_plot]] / 9.80665 # Convert to meters

plot(z_500_day, 
     main = sprintf("500 hPa Geopotential Height\nJune %s, 1991", day_to_plot),
     col = z_pal,
     type = "continuous",
     plg = list(title = "Meters"))

# Add boundaries (National borders thicker, state borders thinner/dashed)
plot(states_vect, add = TRUE, border = "gray20", lwd = 0.5, lty = "dotted")
plot(countries_vect, add = TRUE, border = "black", lwd = 1.2)

# --- Plot B: Specific Humidity (q) at 850 hPa ---
q_pal <- colorRampPalette(c("#fff5eb", "#fd8d3c", "#41b6c4", "#081d58"))(100)
q_850_day <- q_850_stack[[day_to_plot]] * 1000 # Convert to g/kg

plot(q_850_day, 
     main = sprintf("850 hPa Specific Humidity\nJune %s, 1991", day_to_plot),
     col = q_pal,
     type = "continuous",
     plg = list(title = "g/kg"))

# Add boundaries
plot(states_vect, add = TRUE, border = "gray40", lwd = 0.5, lty = "dotted")
plot(countries_vect, add = TRUE, border = "black", lwd = 1.2)

# Reset plot parameters
par(mfrow = c(1, 1))

