# src/15_calculate_daily_anomalies.R
library(terra)
library(dplyr)

# -------------------------------------------------------------------------
# 1. SETUP DIRECTORIES & DISCOVER YEARS
# -------------------------------------------------------------------------
slice_dir   <- "data/processed/4slices"
climo_dir   <- "data/processed/4slices/climo"
anom_dir    <- "data/processed/4slices/anomalies"

# Ensure the output folder exists
dir.create(anom_dir, showWarnings = FALSE)

months_to_process <- c("06", "07", "08", "09", "10")
month_names       <- c("June", "July", "August", "September", "October")

# Discover all available years in the slice directory dynamically
all_nc_files <- list.files(slice_dir, pattern = "\\.nc$")

# Use regex to extract the 4-digit year from the filenames (e.g., ..._2020_06.nc -> 2020)
available_years <- unique(as.numeric(sub(".*_([0-9]{4})_[0-9]{2}\\.nc$", "\\1", all_nc_files)))
available_years <- sort(available_years[!is.na(available_years)])

if (length(available_years) == 0) {
  stop("No valid NetCDF files with years found in the slice directory.")
}

cat(sprintf("Found %d years to process: %s\n", length(available_years), paste(available_years, collapse = ", ")))

# -------------------------------------------------------------------------
# 2. DEFINE SYSTEM ARCHITECTURE MATRIX (All 9 Variables)
# -------------------------------------------------------------------------
var_metadata <- tibble::tribble(
  ~var_id,  ~var_type,  ~block_pos, ~unit_scale,
  "tcwv",   "single",   1,          1,            # No scale adjustment needed
  "q850",   "q",        1,          1000,         # Convert to g/kg
  "q700",   "q",        2,          1000,
  "q500",   "q",        3,          1000,
  "q250",   "q",        4,          1000,
  "z850",   "z",        5,          1/9.80665,    # Convert to geopotential meters
  "z700",   "z",        6,          1/9.80665,
  "z500",   "z",        7,          1/9.80665,
  "z250",   "z",        8,          1/9.80665
)

# -------------------------------------------------------------------------
# 3. MASTER BATCH PROCESSING LOOP
# -------------------------------------------------------------------------
# Outer Loop: Variables
for (v in 1:nrow(var_metadata)) {
  target_var   <- var_metadata$var_id[v]
  v_type       <- var_metadata$var_type[v]
  b_pos        <- var_metadata$block_pos[v]
  scale_factor <- var_metadata$unit_scale[v]
  
  cat(sprintf("\n==================================================\n"))
  cat(sprintf("PROCESSING VARIABLE: %s\n", toupper(target_var)))
  cat(sprintf("==================================================\n"))
  
  # 1. Load the corresponding 153-day smoothed baseline climatology
  climo_path <- file.path(climo_dir, sprintf("climatology_%s_5day_smoothed.tif", target_var))
  if (!file.exists(climo_path)) {
    warning(sprintf("Skipping %s: Climatology baseline file not found at %s", target_var, climo_path))
    next
  }
  climo_raster <- rast(climo_path)
  
  # Inner Loop: Loop through all discovered years
  for (target_year in available_years) {
    cat(sprintf("\n -> Target Year: %d\n", target_year))
    
    actuals_month_list <- list()
    skip_year <- FALSE
    
    # 2. Extract and scale target year actuals month-by-month
    for (m in 1:length(months_to_process)) {
      m_num <- months_to_process[m]
      
      file_prefix <- if(v_type == "single") "single" else "pressure"
      pattern_str <- sprintf("%s_%d_%s\\.nc$", file_prefix, target_year, m_num)
      
      m_file <- list.files(slice_dir, pattern = pattern_str, full.names = TRUE)
      
      if (length(m_file) == 0) {
        warning(sprintf("    Missing file for %d-%s. Skipping year %d.", target_year, m_num, target_year))
        skip_year <- TRUE
        break # Abort month loop if data is missing for this year
      }
      
      r_month <- rast(m_file[1])
      
      divisor  <- if(v_type == "single") 1 else 8
      num_days <- nlyr(r_month) / divisor
      
      start_idx     <- ((b_pos - 1) * num_days) + 1
      end_idx       <- b_pos * num_days
      target_layers <- r_month[[start_idx:end_idx]]
      
      actuals_month_list[[m_num]] <- target_layers * scale_factor
    }
    
    # If any month was missing, safely skip to the next year
    if (skip_year) next
    
    actual_stack <- rast(actuals_month_list)
    names(actual_stack) <- names(climo_raster)
    
    # 3. Execute Matrix Subtraction (Actual - Baseline)
    anomaly_stack <- actual_stack - climo_raster
    
    # 4. Save Final Anomaly Stack
    out_path <- file.path(anom_dir, sprintf("anomaly_%s_%d.tif", target_var, target_year))
    writeRaster(anomaly_stack, out_path, overwrite = TRUE)
    cat(sprintf("    SUCCESS: Saved %s\n", basename(out_path)))
  }
}

cat("\n==================================================\n")
cat("BATCH PROCESSING COMPLETE: All Variables & Years\n")
cat("==================================================\n")