# src/01_sequential_download_era5.R
library(ecmwfr)
library(fs)
library(archive)
library(tictoc)

source("config.R")

message("--- Initializing Sequential Download Pipeline ---")

#job_queue <- expand.grid(year = clim_years, month = target_months, stringsAsFactors = FALSE)
job_queue <- expand.grid(month = target_months, year = clim_years, stringsAsFactors = FALSE)
job_queue$id <- paste0(job_queue$year, "_", job_queue$month)

message("Total months to check: ", nrow(job_queue))

for (i in 1:nrow(job_queue)) {
  yr  <- job_queue$year[i]
  mon <- job_queue$month[i]
  jid <- job_queue$id[i]
  
  final_nc_path_pl <- file.path(dir_raw, paste0("era5_pressure_", jid, ".nc"))
  final_nc_path_sl <- file.path(dir_raw, paste0("era5_single_", jid, ".nc"))
  
  # =========================================================================
  # 1. PRESSURE LEVELS (Geopotential & Specific Humidity)
  # =========================================================================
  if (!file.exists(final_nc_path_pl)) {
    zip_target_pl <- paste0("temp_era5_pl_", jid, ".zip")
    
    request_pl <- list(
      dataset_short_name = era5_pressure_dataset,
      product_type       = list("reanalysis"),
      variable           = target_variables,
      pressure_level     = target_levels,
      year               = list(yr),
      month              = list(mon),
      day                = as.list(sprintf("%02d", 1:31)),
      time               = as.list(target_hours),
      data_format        = "netcdf",
      download_format    = "zip",
      area               = aoi_southwest,
      target             = zip_target_pl
    )
    
    message("\n=======================================================")
    message("Processing Pressure Levels: ", jid)
    message("=======================================================")
    
    tictoc::tic()
    tryCatch({
      wf_request(request = request_pl, transfer = TRUE, path = dir_temp)
      
      zip_file <- file.path(dir_temp, zip_target_pl)
      if (file.exists(zip_file)) {
        archive::archive_extract(zip_file, dir = dir_temp)
        extracted_nc <- list.files(dir_temp, pattern = "\\.nc$", full.names = TRUE)
        
        if (length(extracted_nc) > 0) {
          fs::file_move(extracted_nc[1], final_nc_path_pl)
          message("  -> Stored permanently: ", final_nc_path_pl)
        }
        unlink(file.path(dir_temp, "*"), recursive = TRUE, force = TRUE)
      }
    }, error = function(e) {
      message("  !!! API Error on ", jid, " (PL): ", conditionMessage(e))
    })
    tictoc::toc()
  } else {
    message(sprintf("\nSkipping Pressure Levels %s - already downloaded.", jid))
  }
  
  # =========================================================================
  # 2. SINGLE LEVELS (Total Column Water Vapour / PWAT)
  # =========================================================================
  if (!file.exists(final_nc_path_sl)) {
    zip_target_sl <- paste0("temp_era5_sl_", jid, ".zip")
    
    request_sl <- list(
      dataset_short_name = dataset_single,   # from config
      product_type       = list("reanalysis"),
      variable           = list(var_single), # from config
      year               = list(yr),
      month              = list(mon),
      day                = as.list(sprintf("%02d", 1:31)),
      time               = as.list(target_hours),
      data_format        = "netcdf",
      download_format    = "zip",
      area               = aoi_southwest,
      target             = zip_target_sl
    )
    
    message("-------------------------------------------------------")
    message("Processing Single Levels (PWAT): ", jid)
    message("-------------------------------------------------------")
    
    tictoc::tic()
    tryCatch({
      wf_request(request = request_sl, transfer = TRUE, path = dir_temp)
      
      zip_file <- file.path(dir_temp, zip_target_sl)
      if (file.exists(zip_file)) {
        archive::archive_extract(zip_file, dir = dir_temp)
        extracted_nc <- list.files(dir_temp, pattern = "\\.nc$", full.names = TRUE)
        
        if (length(extracted_nc) > 0) {
          fs::file_move(extracted_nc[1], final_nc_path_sl)
          message("  -> Stored permanently: ", final_nc_path_sl)
        }
        unlink(file.path(dir_temp, "*"), recursive = TRUE, force = TRUE)
      }
    }, error = function(e) {
      message("  !!! API Error on ", jid, " (SL): ", conditionMessage(e))
    })
    tictoc::toc()
  } else {
    message(sprintf("Skipping Single Levels %s - already downloaded.", jid))
  }
}

message("\n--- Download Pipeline Complete! ---")