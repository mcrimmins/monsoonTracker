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


###### test plots 
# src/02_test_plots_500mb.R
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


####### test daily averaging....
# src/03_daily_averages.R
library(terra)
library(maps)

# 1. Load our July 1991 test data
nc_file <- "data/raw/era5_pressure_1991_07.nc"
message("Loading NetCDF: ", nc_file)
era5_data <- rast(nc_file)

# 2. Extract ALL hourly layers for 500 hPa Geopotential
# We use grep to find all layers that match the 500mb z-level
message("Extracting 500 hPa hourly layers...")
z_500_hourly <- era5_data[[grep("z_pressure_level=500", names(era5_data))]]

# Convert to Geopotential Height in meters
z_500_hourly <- z_500_hourly / 9.80665

# 3. Create a daily grouping index
# time() extracts the hourly datetime object for each layer
layer_times <- time(z_500_hourly)

# as.Date() strips the hours/minutes, leaving just the YYYY-MM-DD
daily_index <- as.Date(layer_times)

# 4. Calculate Daily Averages
# tapp() applies the 'mean' function based on our daily_index
message("Calculating daily averages (this may take a few seconds)...")
z_500_daily <- tapp(z_500_hourly, index = daily_index, fun = mean, na.rm = TRUE)

# 5. Verify the reduction in data volume
message("\nOriginal hourly layers: ", nlyr(z_500_hourly))
message("New daily layers: ", nlyr(z_500_daily))

# 6. Plot the first day to ensure spatial integrity is maintained
plot(z_500_daily[[1]], 
     main = "Daily Average 500 hPa Geopotential Height (m)\nJuly 1, 1991", 
     col = terrain.colors(50),
     mar = c(3, 3, 3, 4))

map("world", add = TRUE, col = "black", lwd = 1.5)
map("state", add = TRUE, col = "gray20", lwd = 1, lty = 3)


###### precip water calc testing...
# src/04_test_pwat.R
library(terra)
library(maps)

# 1. Load the data
nc_file <- "data/raw/era5_pressure_1991_07.nc"
message("Loading NetCDF: ", nc_file)
era5_data <- rast(nc_file)

# 2. Extract Specific Humidity (q) for the first hour (_1)
# Keep raw kg/kg units for the math
message("Extracting pressure levels...")
q_850 <- era5_data[["q_pressure_level=850_1"]]
q_700 <- era5_data[["q_pressure_level=700_1"]]
q_500 <- era5_data[["q_pressure_level=500_1"]]
q_250 <- era5_data[["q_pressure_level=250_1"]]

# 3. Perform Trapezoidal Vertical Integration
g <- 9.80665 # Gravity

message("Calculating layer moisture...")
# Layer 1: 850 hPa to 700 hPa (dp = 150 hPa = 15000 Pa)
layer1_mm <- ((q_850 + q_700) / 2) * (15000) / g

# Layer 2: 700 hPa to 500 hPa (dp = 200 hPa = 20000 Pa)
layer2_mm <- ((q_700 + q_500) / 2) * (20000) / g

# Layer 3: 500 hPa to 250 hPa (dp = 250 hPa = 25000 Pa)
layer3_mm <- ((q_500 + q_250) / 2) * (25000) / g

# Total Partial Column PWAT
pwat_mm <- layer1_mm + layer2_mm + layer3_mm

# 4. Plot the resulting PWAT field
# Using a custom color ramp typical for moisture (white -> light blue -> blue -> magenta)
plot(pwat_mm, 
     main = "Partial Column Precipitable Water (mm)\n850-250 hPa | July 1, 1991", 
     col = colorRampPalette(c("white", "lightblue", "blue", "magenta"))(50),
     mar = c(3, 3, 3, 4))

map("world", add = TRUE, col = "black", lwd = 1.5)
map("state", add = TRUE, col = "gray20", lwd = 1, lty = 3)
