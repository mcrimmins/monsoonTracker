#### utility function...

library(terra)
library(maps)

# 1. Load the July 1991 data
nc_file <- "data/raw/era5_pressure_1991_07.nc"
message("Loading NetCDF: ", nc_file)
era5_data <- rast(nc_file)

# 2. Extract exactly the 500 hPa level for the 1st time step
# Using the exact name terra generated
z_500 <- era5_data[["z_pressure_level=500_1"]] / 9.80665  # Convert to meters
q_500 <- era5_data[["q_pressure_level=500_1"]] * 1000     # Convert to g/kg

# 3. Plot the results with map boundaries
par(mfrow = c(1, 2))

# --- Plot 1: 500mb Geopotential Height ---
plot(z_500, 
     main = "500 hPa Geopotential Height (m)\nJuly 1, 1991", 
     col = terrain.colors(50),
     mar = c(3, 3, 3, 4))

map("world", add = TRUE, col = "black", lwd = 1.5)
map("state", add = TRUE, col = "gray20", lwd = 1, lty = 3)

# --- Plot 2: 500mb Specific Humidity ---
plot(q_500, 
     main = "500 hPa Specific Humidity (g/kg)\nJuly 1, 1991", 
     col = topo.colors(50),
     mar = c(3, 3, 3, 4))

map("world", add = TRUE, col = "black", lwd = 1.5)
map("state", add = TRUE, col = "gray20", lwd = 1, lty = 3)