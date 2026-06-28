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

##### PWAT plotting test...
# src/05_test_daily_pwat.R
library(terra)
library(maps)

# 1. Load the raw hourly PWAT file
nc_file <- "data/raw/era5_single_1991_07.nc"
message("Loading NetCDF: ", nc_file)
pwat_hourly <- rast(nc_file)

# 2. Create a grouping index
# July has 31 days, so the file has 744 layers (31 * 24).
# We create a vector that repeats '1' twenty-four times, '2' twenty-four times, etc.
n_days <- nlyr(pwat_hourly) / 24
daily_index <- rep(1:n_days, each = 24)

# 3. Calculate daily averages
message("Averaging 24-hour blocks into daily means...")
pwat_daily <- tapp(pwat_hourly, index = daily_index, fun = mean)

# 4. Plot the daily average for July 1 (Layer 1)
plot(pwat_daily[[1]], 
     main = "Daily Average PWAT (mm)\nJuly 1, 1991", 
     col = colorRampPalette(c("white", "lightblue", "blue", "magenta"))(50),
     mar = c(3, 3, 3, 4))

map("world", add = TRUE, col = "black", lwd = 1.5)
map("state", add = TRUE, col = "gray20", lwd = 1, lty = 3)


###### generic daily averaging script ----
# src/03_daily_averages_benchmark.R
library(terra)
library(stringr)

# Source project configurations
source("config.R")

# =========================================================================
# 1. INITIALIZATION & ENGINE TUNING
# =========================================================================
cat("\n==================================================\n")
cat(" INITIALIZING PERFORMANCE BENCHMARK ENGINE\n")
cat("==================================================\n")

total_start_time <- Sys.time()

# Optimize hardware resource allocation
n_cores <- max(1, parallel::detectCores() - 1)
terraOptions(cores = n_cores)
terraOptions(memfrac = 0.8)

cat(sprintf("[INIT] Detected %s Logical Processors.\n", parallel::detectCores()))
cat(sprintf("[INIT] Threading set to %s cores (Intel TBB C++ mode active).\n", n_cores))
cat("[INIT] Memory allocation ceiling pushed to 80% RAM.\n")

# Aggregation metric choice ("mean" or "max")
aggr_metric <- "mean"

# Build the queue and isolate to our single test year
job_queue <- expand.grid(month = target_months, year = clim_years, stringsAsFactors = FALSE)
test_year <- 1991
job_queue <- job_queue[job_queue$year == test_year, ]

cat(sprintf("[INIT] Benchmark isolated to Year: %s (%s months total)\n", test_year, nrow(job_queue)))
cat("==================================================\n\n")

# =========================================================================
# 2. RUN PIPELINE WITH GRANULAR TIMERS
# =========================================================================
for (i in 1:nrow(job_queue)) {
  yr  <- job_queue$year[i]
  mon <- job_queue$month[i]
  jid <- paste0(yr, "_", mon)
  
  raw_pl_file  <- file.path(dir_raw, paste0("era5_pressure_", jid, ".nc"))
  raw_sl_file  <- file.path(dir_raw, paste0("era5_single_", jid, ".nc"))
  
  out_pl_file  <- file.path(dir_processed, sprintf("daily_%s_era5_pressure_%s.nc", aggr_metric, jid))
  out_sl_file  <- file.path(dir_processed, sprintf("daily_%s_era5_single_%s.nc", aggr_metric, jid))
  
  cat(sprintf("▶ STARTING PROCESSING BLOCK FOR: %s-%s\n", yr, mon))
  month_start_time <- Sys.time()
  
  # -------------------------------------------------------------------------
  # STEP A: 3D PRESSURE LEVELS (q & z)
  # -------------------------------------------------------------------------
  cat("  [Step 1/2] Checking 3D Pressure Levels...\n")
  
  if (file.exists(raw_pl_file) && !file.exists(out_pl_file)) {
    pl_start_time <- Sys.time()
    
    cat("    - Loading metadata from NetCDF disk... ")
    raw_stack   <- rast(raw_pl_file)
    layer_names <- names(raw_stack)
    
    vars_present   <- unique(str_extract(layer_names, "^[^_]+"))
    levels_present <- unique(str_extract(layer_names, "(?<=pressure_level=)\\d+"))
    levels_present <- levels_present[!is.na(levels_present)]
    cat("Done.\n")
    
    pl_layers <- list()
    
    for (v in vars_present) {
      for (lvl in levels_present) {
        loop_start <- Sys.time()
        cat(sprintf("    - Processing Variable '%s' at %s hPa:\n", v, lvl))
        
        sub_pattern <- sprintf("^%s_pressure_level=%s_", v, lvl)
        matching_indices <- which(str_detect(layer_names, sub_pattern))
        if (length(matching_indices) == 0) {
          cat("        * Skipped (No layers found).\n")
          next
        }
        
        sub_stack <- raw_stack[[matching_indices]]
        
        # --- THE RAM-FORCE TWEAK (CHUNKED) ---
        cat("        * Reading from disk to RAM cache (~2 mins)... ")
        ram_start <- Sys.time()
        sub_layer_names <- names(sub_stack)
        sub_stack <- sub_stack + 0 
        names(sub_stack) <- sub_layer_names
        ram_end <- Sys.time()
        cat(sprintf("Done (%s sec).\n", round(difftime(ram_end, ram_start, units="secs"), 2)))
        
        # --- CALCULATE DAILY AVERAGES ---
        cat(sprintf("        * Executing C++ '%s' reduction matrix... ", aggr_metric))
        n_days      <- nlyr(sub_stack) / 24
        daily_index <- rep(1:n_days, each = 24)
        
        daily_sub_stack <- tapp(sub_stack, index = daily_index, fun = aggr_metric)
        names(daily_sub_stack) <- sprintf("%s_%shPa_day_%s", v, lvl, 1:n_days)
        cat("Done.\n")
        
        pl_layers[[paste0(v, "_", lvl)]] <- daily_sub_stack
        
        loop_end <- Sys.time()
        cat(sprintf("        * Total sub-step time: %s sec.\n", round(difftime(loop_end, loop_start, units="secs"), 2)))
      }
    }
    
    # Write full 3D stack back to disk
    cat("    - Compiling variables and compiling output NetCDF... ")
    write_start <- Sys.time()
    master_pl_stack <- rast(pl_layers)
    writeCDF(master_pl_stack, out_pl_file, overwrite = TRUE)
    write_end <- Sys.time()
    cat(sprintf("Done (%s sec).\n", round(difftime(write_end, write_start, units="secs"), 2)))
    
    pl_end_time <- Sys.time()
    cat(sprintf("  ✔ Completed 3D Pressure Levels in: %s seconds.\n", 
                round(difftime(pl_end_time, pl_start_time, units="secs"), 1)))
    
  } else if (!file.exists(raw_pl_file)) {
    cat("  ⚠️  [Skip 3D] Raw source file missing from data/raw/.\n")
  } else {
    cat("  ℹ️  [Skip 3D] Processed target file already exists in data/processed/.\n")
  }
  
  # -------------------------------------------------------------------------
  # STEP B: SINGLE LEVEL PWAT (tcwv)
  # -------------------------------------------------------------------------
  cat("  [Step 2/2] Checking Single Level PWAT...\n")
  
  if (file.exists(raw_sl_file) && !file.exists(out_sl_file)) {
    sl_start_time <- Sys.time()
    
    cat("    - Loading and RAM-forcing PWAT layer... ")
    pwat_hourly <- rast(raw_sl_file)
    
    # RAM-Force PWAT
    pwat_names <- names(pwat_hourly)
    pwat_hourly <- pwat_hourly + 0
    names(pwat_hourly) <- pwat_names
    cat("Done.\n")
    
    cat(sprintf("    - Executing C++ '%s' reduction matrix... ", aggr_metric))
    n_days      <- nlyr(pwat_hourly) / 24
    daily_index <- rep(1:n_days, each = 24)
    
    pwat_daily  <- tapp(pwat_hourly, index = daily_index, fun = aggr_metric)
    names(pwat_daily) <- sprintf("tcwv_day_%s", 1:n_days)
    cat("Done.\n")
    
    cat("    - Writing processed PWAT NetCDF to disk... ")
    write_sl_start <- Sys.time()
    writeCDF(pwat_daily, out_sl_file, overwrite = TRUE)
    write_sl_end <- Sys.time()
    cat(sprintf("Done (%s sec).\n", round(difftime(write_sl_end, write_sl_start, units="secs"), 2)))
    
    sl_end_time <- Sys.time()
    cat(sprintf("  ✔ Completed PWAT Analysis in: %s seconds.\n", 
                round(difftime(sl_end_time, sl_start_time, units="secs"), 1)))
    
  } else if (!file.exists(raw_sl_file)) {
    cat("  ⚠️  [Skip PWAT] Raw source file missing from data/raw/.\n")
  } else {
    cat("  ℹ️  [Skip PWAT] Processed target file already exists in data/processed/.\n")
  }
  
  month_end_time <- Sys.time()
  cat(sprintf("🏁 FINISHED BLOCK %s-%s in total time of: %s seconds.\n", 
              yr, mon, round(difftime(month_end_time, month_start_time, units="secs"), 1)))
  cat("--------------------------------------------------\n\n")
}

# =========================================================================
# 3. GLOBAL PERFORMANCE METRIC REPORT
# =========================================================================
total_end_time <- Sys.time()
total_duration <- total_end_time - total_start_time

cat("==================================================\n")
cat(sprintf(" FINAL BENCHMARK PERFORMANCE REPORT (YEAR: %s)\n", test_year))
cat("==================================================\n")
cat("Total Clock Run Duration:\n")
print(total_duration)
cat("==================================================\n")


##### parallel testing...

# src/03_parallel_test_worker.R
library(terra)
library(stringr)

source("config.R")

# Accept specific years from the launcher
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {
  start_yr <- as.numeric(args[1])
  end_yr   <- as.numeric(args[2])
} else {
  stop("Missing year range arguments. Use the parallel launcher script.")
}

# Ultrabook-friendly resource configuration
n_cores <- 2  
terraOptions(cores = n_cores)
terraOptions(memfrac = 0.20) 

aggr_metric <- "mean"

# Filter the test queue for this worker's assigned block
job_queue <- expand.grid(month = target_months, year = clim_years, stringsAsFactors = FALSE)
job_queue <- job_queue[job_queue$year >= start_yr & job_queue$year <= end_yr, ]

message(sprintf("\n[WORKER START] Range: %s to %s | Processing %s months.\n", start_yr, end_yr, nrow(job_queue)))
worker_start_time <- Sys.time()

for (i in 1:nrow(job_queue)) {
  yr  <- job_queue$year[i]
  mon <- job_queue$month[i]
  jid <- paste0(yr, "_", mon)
  
  raw_pl_file  <- file.path(dir_raw, paste0("era5_pressure_", jid, ".nc"))
  raw_sl_file  <- file.path(dir_raw, paste0("era5_single_", jid, ".nc"))
  out_pl_file  <- file.path(dir_processed, sprintf("daily_%s_era5_pressure_%s.nc", aggr_metric, jid))
  out_sl_file  <- file.path(dir_processed, sprintf("daily_%s_era5_single_%s.nc", aggr_metric, jid))
  
  # Step A: 3D Pressure Levels
  if (file.exists(raw_pl_file) && !file.exists(out_pl_file)) {
    loop_start  <- Sys.time()
    raw_stack   <- rast(raw_pl_file)
    layer_names <- names(raw_stack)
    vars_present   <- unique(str_extract(layer_names, "^[^_]+"))
    levels_present <- unique(str_extract(layer_names, "(?<=pressure_level=)\\d+"))
    levels_present <- levels_present[!is.na(levels_present)]
    
    pl_layers <- list()
    for (v in vars_present) {
      for (lvl in levels_present) {
        sub_pattern <- sprintf("^%s_pressure_level=%s_", v, lvl)
        matching_indices <- which(str_detect(layer_names, sub_pattern))
        if (length(matching_indices) == 0) next
        
        sub_stack <- raw_stack[[matching_indices]]
        
        # Safe chunked RAM-force
        sub_layer_names <- names(sub_stack)
        sub_stack <- sub_stack + 0 
        names(sub_stack) <- sub_layer_names
        
        n_days      <- nlyr(sub_stack) / 24
        daily_index <- rep(1:n_days, each = 24)
        daily_sub_stack <- tapp(sub_stack, index = daily_index, fun = aggr_metric)
        names(daily_sub_stack) <- sprintf("%s_%shPa_day_%s", v, lvl, 1:n_days)
        
        pl_layers[[paste0(v, "_", lvl)]] <- daily_sub_stack
      }
    }
    master_pl_stack <- rast(pl_layers)
    writeCDF(master_pl_stack, out_pl_file, overwrite = TRUE)
    
    loop_end <- Sys.time()
    message(sprintf("  [%s-%s] Finished 3D Levels in %s sec.", yr, mon, round(difftime(loop_end, loop_start, units="secs"), 1)))
  }
  
  # Step B: PWAT
  if (file.exists(raw_sl_file) && !file.exists(out_sl_file)) {
    pwat_hourly <- rast(raw_sl_file)
    pwat_names  <- names(pwat_hourly)
    pwat_hourly <- pwat_hourly + 0
    names(pwat_hourly) <- pwat_names
    
    n_days      <- nlyr(pwat_hourly) / 24
    daily_index <- rep(1:n_days, each = 24)
    pwat_daily  <- tapp(pwat_hourly, index = daily_index, fun = aggr_metric)
    names(pwat_daily) <- sprintf("tcwv_day_%s", 1:n_days)
    
    writeCDF(pwat_daily, out_sl_file, overwrite = TRUE)
  }
}

worker_end_time <- Sys.time()
message(sprintf("\n[WORKER FINISHED] Range %s-%s completed in %s minutes.\n", 
                start_yr, end_yr, round(difftime(worker_end_time, worker_start_time, units="mins"), 1)))


