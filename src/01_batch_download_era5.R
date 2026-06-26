# src/01_batch_download_era5.R
library(ecmwfr)
library(fs)
library(archive)
library(tictoc)

source("config.R")

message("--- Initializing Monsoon ERA5 Climatology Batch Pipeline ---")

# 1. Build master queue
job_queue <- expand.grid(year = clim_years, month = target_months, stringsAsFactors = FALSE)
job_queue$id <- paste0(job_queue$year, "_", job_queue$month)

# Check data/raw/ and filter out chunks that are already finished
existing_files <- list.files(dir_raw, pattern = "\\.nc$")
job_queue <- job_queue[!paste0("era5_pressure_", job_queue$id, ".nc") %in% existing_files, ]

message("Total remaining chunks to process: ", nrow(job_queue))
if (nrow(job_queue) == 0) stop("All files are already downloaded!")

# Group remaining jobs into batches of 3 (Max concurrent jobs allowed by CDS)
batch_size <- 3
job_queue$batch <- rep(1:ceiling(nrow(job_queue)/batch_size), each = batch_size)[1:nrow(job_queue)]

# 2. Master Processing Loop
for (b in unique(job_queue$batch)) {
  current_batch <- job_queue[job_queue$batch == b, ]
  
  message("\n=======================================================")
  message("=== Starting Batch ", b, " of ", max(job_queue$batch), " ===")
  message("=======================================================")
  
  active_receipts <- list()
  
  # PHASE A: SUBMIT THE ENTIRE BATCH CONCURRENTLY
  for (i in 1:nrow(current_batch)) {
    yr  <- current_batch$year[i]
    mon <- current_batch$month[i]
    jid <- current_batch$id[i]
    
    zip_target <- paste0("temp_era5_", jid, ".zip")
    
    request <- list(
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
      target             = zip_target
    )
    
    message("-> Submitting Asynchronous Request for: ", jid)
    
    # transfer = FALSE returns a "receipt" immediately without waiting
    rcpt <- tryCatch({
      wf_request(request = request, transfer = FALSE, path = dir_temp)
    }, error = function(e) {
      message("  !!! API Error on ", jid, ": ", conditionMessage(e))
      NULL
    })
    
    if (!is.null(rcpt)) {
      active_receipts[[jid]] <- rcpt
    }
    
    # Pause 3 seconds so we don't trip the API rate-limiter while submitting
    Sys.sleep(3)
  }
  
  # PHASE B: HARVEST THE BATCH
  message("\n--- Batch submitted. Moving to Harvest Phase ---")
  
  for (jid in names(active_receipts)) {
    rcpt <- active_receipts[[jid]]
    zip_file <- file.path(dir_temp, paste0("temp_era5_", jid, ".zip"))
    
    message("\n[Waiting for CDS & Downloading: ", jid, "]")
    
    tictoc::tic()
    tryCatch({
      
      # 1. Manually update the status of this specific receipt
      rcpt$update_status()
      
      # 2. Build a bulletproof waiting loop
      while(rcpt$is_pending() || rcpt$is_running()) {
        message("  ... cloud job still processing, checking again in 30 seconds ...")
        Sys.sleep(30)
        rcpt$update_status() # Refresh status from server
      }
      
      # 3. Trigger download ONLY when we know it is definitively ready
      if (rcpt$is_success()) {
        message("  -> Cloud processing complete! Downloading...")
        rcpt$download() 
      } else if (rcpt$is_failed()) {
        message("  !!! Job failed on the Copernicus server.")
      }
      
    }, error = function(e) {
      message("  !!! Download interrupted/failed for ", jid, ": ", conditionMessage(e))
    })
    tictoc::toc()
    
    # Unzip, extract NetCDF, and clean up
    if (file.exists(zip_file)) {
      final_nc_path <- file.path(dir_raw, paste0("era5_pressure_", jid, ".nc"))
      
      archive::archive_extract(zip_file, dir = dir_temp)
      extracted_nc <- list.files(dir_temp, pattern = "\\.nc$", full.names = TRUE)
      
      if (length(extracted_nc) > 0) {
        fs::file_move(extracted_nc[1], final_nc_path)
        message("  -> Stored permanently: ", final_nc_path)
      } else {
        message("  !!! Unzip succeeded, but no NetCDF found inside.")
      }
      
      # Wipe the temp folder completely clean before processing the next month
      unlink(file.path(dir_temp, "*"), recursive = TRUE, force = TRUE)
    }
  }
}

message("\n--- Batch Download Pipeline Complete! ---")