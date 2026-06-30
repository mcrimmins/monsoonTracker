# src/02_daily_stat_processing.R
library(terra)
library(stringr)

# Load global configuration (directories, months, years)
source("config.R")

# Optimize terra engine for ThinkPad X1 Carbon
terraOptions(cores = 2)
terraOptions(memfrac = 0.25)

#' ERA5 Flexible Daily Aggregation Pipeline
#'
#' @param target_years Vector of years to process (e.g., 1991:2020). If NULL, uses config.
#' @param target_months Vector of months to process (e.g., c("06", "07")). If NULL, uses config.
#' @param hours Vector of UTC hours to sample (0 to 23). 
#'              e.g., 0:23 for all, c(0,6,12,18) for 4-times, c(0) for 00Z slice.
#' @param stat Aggregation metric (e.g., "mean", "max").
#' @param overwrite Logical. If TRUE, overwrites existing processed files.
process_era5_pipeline <- function(target_years = NULL, 
                                  target_months = NULL, 
                                  hours = c(0, 6, 12, 18), 
                                  stat = "mean", 
                                  overwrite = FALSE) {
  
  # -------------------------------------------------------------------------
  # 1. VALIDATION & SETUP
  # -------------------------------------------------------------------------
  if (is.null(target_years)) target_years <- clim_years
  if (is.null(target_months)) target_months <- target_months
  if (any(hours < 0 | hours > 23)) stop("Hours must be between 0 and 23.")
  
  # Determine sub-directory and naming nomenclature based on selected hours
  n_hours <- length(hours)
  if (n_hours == 24) {
    method_dir <- "all_hours"
    file_prefix <- sprintf("daily_%s_all_hours", stat)
  } else if (n_hours == 1) {
    method_dir <- sprintf("slice_%02dZ", hours[1])
    file_prefix <- sprintf("daily_%s", method_dir) # Drop the stat name for a single slice
  } else {
    method_dir <- sprintf("%dslices", n_hours)
    file_prefix <- sprintf("daily_%s_%s", stat, method_dir)
  }
  
  # Create nested output directory
  out_dir <- file.path(dir_processed, method_dir)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  job_queue <- expand.grid(month = target_months, year = target_years, stringsAsFactors = FALSE)
  
  cat("\n==================================================\n")
  cat(" ERA5 REGIONAL CLIMATE PIPELINE\n")
  cat("==================================================\n")
  cat(sprintf(" Mode:       %s\n", method_dir))
  cat(sprintf(" Statistic:  %s\n", stat))
  cat(sprintf(" Output Dir: %s\n", out_dir))
  cat(sprintf(" Total Jobs: %s monthly blocks\n", nrow(job_queue)))
  cat("==================================================\n\n")
  
  total_start <- Sys.time()
  
  # Convert UTC hours to 1-based R indices for the daily 24-hour block
  base_indices <- hours + 1 
  
  # -------------------------------------------------------------------------
  # 2. MAIN PROCESSING LOOP
  # -------------------------------------------------------------------------
  for (i in 1:nrow(job_queue)) {
    yr  <- job_queue$year[i]
    mon <- job_queue$month[i]
    jid <- paste0(yr, "_", mon)
    
    raw_pl_file <- file.path(dir_raw, paste0("era5_pressure_", jid, ".nc"))
    raw_sl_file <- file.path(dir_raw, paste0("era5_single_", jid, ".nc"))
    out_pl_file <- file.path(out_dir, sprintf("%s_era5_pressure_%s.nc", file_prefix, jid))
    out_sl_file <- file.path(out_dir, sprintf("%s_era5_single_%s.nc", file_prefix, jid))
    
    cat(sprintf("▶ [%s/%s] PROCESSING %s-%s\n", i, nrow(job_queue), yr, mon))
    
    # === STEP A: 3D PRESSURE LEVELS ===
    if (file.exists(raw_pl_file)) {
      if (!file.exists(out_pl_file) || overwrite) {
        
        raw_stack <- rast(raw_pl_file)
        layer_names <- names(raw_stack)
        
        vars_present <- unique(str_extract(layer_names, "^[^_]+"))
        levels_present <- unique(str_extract(layer_names, "(?<=pressure_level=)\\d+"))
        levels_present <- levels_present[!is.na(levels_present)]
        
        # Calculate exactly which layers to extract for the entire month
        n_days <- nlyr(raw_stack) / length(vars_present) / length(levels_present) / 24
        extract_idx <- as.vector(sapply(0:(n_days - 1), function(d) (d * 24) + base_indices))
        
        pl_layers <- list()
        
        for (v in vars_present) {
          for (lvl in levels_present) {
            sub_pattern <- sprintf("^%s_pressure_level=%s_", v, lvl)
            matching_indices <- which(str_detect(layer_names, sub_pattern))
            if (length(matching_indices) == 0) next
            
            sub_stack <- raw_stack[[matching_indices]]
            sub_stack_sliced <- sub_stack[[extract_idx]]
            
            # RAM Cache
            sub_layer_names <- names(sub_stack_sliced)
            sub_stack_sliced <- sub_stack_sliced + 0
            names(sub_stack_sliced) <- sub_layer_names
            
            # --- THE METADATA & BYPASS FIX ---
            if (n_hours == 1) {
              daily_sub_stack <- sub_stack_sliced # Bypass tapp entirely for single slices
            } else {
              daily_index <- rep(1:n_days, each = n_hours)
              daily_sub_stack <- tapp(sub_stack_sliced, index = daily_index, fun = stat)
            }
            
            names(daily_sub_stack) <- sprintf("%s_%shPa_day_%s", v, lvl, 1:n_days)
            
            # Clear fragmented NetCDF tags to prevent writeCDF compilation errors
            time(daily_sub_stack) <- NULL
            depth(daily_sub_stack) <- NULL
            
            pl_layers[[paste0(v, "_", lvl)]] <- daily_sub_stack
          }
        }
        writeCDF(rast(pl_layers), out_pl_file, overwrite = TRUE)
        cat(sprintf("  ✔ 3D Pressure Levels saved to /%s/.\n", method_dir))
        
      } else {
        cat(sprintf("  ℹ️ 3D Levels already exist in /%s/. Skipping...\n", method_dir))
      }
    }
    
    # === STEP B: SINGLE LEVEL PWAT ===
    if (file.exists(raw_sl_file)) {
      if (!file.exists(out_sl_file) || overwrite) {
        
        pwat_hourly <- rast(raw_sl_file)
        n_days_sl <- nlyr(pwat_hourly) / 24
        extract_idx_sl <- as.vector(sapply(0:(n_days_sl - 1), function(d) (d * 24) + base_indices))
        
        pwat_sliced <- pwat_hourly[[extract_idx_sl]]
        
        pwat_names <- names(pwat_sliced)
        pwat_sliced <- pwat_sliced + 0
        names(pwat_sliced) <- pwat_names
        
        if (n_hours == 1) {
          pwat_daily <- pwat_sliced
        } else {
          daily_index_sl <- rep(1:n_days_sl, each = n_hours)
          pwat_daily <- tapp(pwat_sliced, index = daily_index_sl, fun = stat)
        }
        
        names(pwat_daily) <- sprintf("tcwv_day_%s", 1:n_days_sl)
        
        # Clear tags for consistency
        time(pwat_daily) <- NULL
        depth(pwat_daily) <- NULL
        
        writeCDF(pwat_daily, out_sl_file, overwrite = TRUE)
        cat(sprintf("  ✔ Single Level PWAT saved to /%s/.\n", method_dir))
        
      } else {
        cat(sprintf("  ℹ️ PWAT already exists in /%s/. Skipping...\n", method_dir))
      }
    }
  }
  
  total_end <- Sys.time()
  cat("\n==================================================\n")
  cat(sprintf("🏁 PIPELINE COMPLETE in %s minutes.\n", 
              round(difftime(total_end, total_start, units="mins"), 1)))
  cat("==================================================\n")
}

# =========================================================================
# USAGE EXAMPLES (Uncomment to run)
# =========================================================================

# --- 1. Test Run: Single Month, 00Z exact slice ---
# tictoc::tic()
# process_era5_pipeline(
#   target_years = 1991,
#   target_months = "06",
#   hours = c(0,6,12,18),            # Grabs exactly hour 00 (00:00 UTC)
#   stat = "mean",           # Tapp processes the single slice without altering data to preserve metadata
#   overwrite = TRUE
# )
# tictoc::toc()

# --- 2. Full Run: 4 Times Daily for the entire historical stack ---
process_era5_pipeline(
  target_years = clim_years,       # Explicitly pass the years from config.R
  target_months = target_months,   # Explicitly pass the months from config.R
  hours = c(0, 6, 12, 18),         # 4 slices per day
  stat = "mean", 
  overwrite = FALSE                # Will safely skip anything you've already processed
)


# --- 3. Benchmark Run: All 24 hours for a single year ---
# process_era5_pipeline(
#   target_years = 1992, 
#   target_months = NULL, 
#   hours = 0:23,            # Full 24-hour daily averaging
#   stat = "mean", 
#   overwrite = FALSE
# )