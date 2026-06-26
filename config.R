# config.R - Shared parameters optimized for batch downloading
# MAC 06/26/26 Gemini Thinking/Pro

library(ecmwfr)

# 1. Directory Structure
dir_raw       <- "data/raw"
dir_processed <- "data/processed"
dir_temp      <- "data/temp"

# Create directories if missing
invisible(lapply(c(dir_raw, dir_processed, dir_temp), fs::dir_create, recurse = TRUE))

# 2. Spatial Bounding Box (Southwest US + Deep Mexico expansion)
aoi_southwest <- c(45, -125, 15, -95)

# 3. Temporal Parameters (June through October)
clim_years    <- as.character(1991:2025)
target_months <- c("06", "07", "08", "09", "10") 
target_hours  <- sprintf("%02d:00", 0:23)       

# 4. Vertical Atmospheric Levels & Elements
era5_pressure_dataset <- "reanalysis-era5-pressure-levels"
target_levels         <- list("250", "500", "700", "850")
target_variables      <- list("geopotential", "specific_humidity")

# 5. Concurrency Controls
MAX_CONCURRENT_JOBS   <- 3    # How many jobs can run on CDS simultaneously
POLL_INTERVAL_SECONDS <- 60   # How often to check status of running jobs