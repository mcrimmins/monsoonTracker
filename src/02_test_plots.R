# src/02_test_plots_500mb.R
library(terra)
library(maps) # Adds geographic boundaries

# 1. Load the July 1991 data
nc_file <- "data/raw/era5_pressure_1991_07.nc"
message("Loading NetCDF: ", nc_file)
era5_data <- rast(nc_file)

# 2. Extract our variables
z_layers <- era5_data[[grep("^z", names(era5_data))]]
q_layers <- era5_data[[grep("^q", names(era5_data))]]

# 3. Target the 500 hPa level
# Assuming the order is 850, 700, 500, 250, we grab layer 3
z_500 <- z_layers[[3]] / 9.80665  # Convert to meters
q_500 <- q_layers[[3]] * 1000     # Convert to g/kg

# 4. Plot the results with map boundaries
par(mfrow = c(1, 2))

# --- Plot 1: 500mb Geopotential Height ---
plot(z_500, 
     main = "500 hPa Geopotential Height (m)\nJuly 1, 1991", 
     col = terrain.colors(50),
     mar = c(3, 3, 3, 4))

# Add country and state boundaries
map("world", add = TRUE, col = "black", lwd = 1.5)
map("state", add = TRUE, col = "gray20", lwd = 1, lty = 3)

# --- Plot 2: 500mb Specific Humidity ---
plot(q_500, 
     main = "500 hPa Specific Humidity (g/kg)\nJuly 1, 1991", 
     col = topo.colors(50),
     mar = c(3, 3, 3, 4))

# Add country and state boundaries
map("world", add = TRUE, col = "black", lwd = 1.5)
map("state", add = TRUE, col = "gray20", lwd = 1, lty = 3)