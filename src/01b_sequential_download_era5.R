# src/01b_sequential_download_era5.R
library(ecmwfr)
library(fs)
library(archive)
library(tictoc)

source("config.R")

message("--- Initializing Sequential Download Pipeline ---")

job_queue <- expand.grid(year = clim_years, month = target_months, stringsAsFactors = FALSE)
job_queue$id <- paste0(job_queue$year, "_", job_queue$month)

# Filter out already downloaded files
existing_files <- list.files(dir_raw, pattern = "\\.nc$")
job_queue <- job_queue[!paste0("era5_pressure_", job_queue$id, ".nc") %in% existing_files, ]

message("Total remaining chunks to process: ", nrow(job_queue))
if (nrow(job_queue) == 0) stop("All files are already downloaded!")

for (i in 1:nrow(job_queue)) {
  yr  <- job_queue$year[i]
  mon <- job_queue$month[i]
  jid <- job_queue$id[i]
  
  zip_target <- paste0("temp_era5_", jid, ".zip")
  final_nc_path <- file.path(dir_raw, paste0("era5_pressure_", jid, ".nc"))
  
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
  
  message("\n=======================================================")
  message("Processing: ", jid)
  message("=======================================================")
  
  tictoc::tic()
  tryCatch({
    # transfer = TRUE makes R block until the file is 100% downloaded
    wf_request(request = request, transfer = TRUE, path = dir_temp)
    
    zip_file <- file.path(dir_temp, zip_target)
    if (file.exists(zip_file)) {
      archive::archive_extract(zip_file, dir = dir_temp)
      extracted_nc <- list.files(dir_temp, pattern = "\\.nc$", full.names = TRUE)
      
      if (length(extracted_nc) > 0) {
        fs::file_move(extracted_nc[1], final_nc_path)
        message("  -> Stored permanently: ", final_nc_path)
      }
      
      # Wipe temp directory clean for the next loop
      unlink(file.path(dir_temp, "*"), recursive = TRUE, force = TRUE)
    }
  }, error = function(e) {
    message("  !!! API Error on ", jid, ": ", conditionMessage(e))
  })
  tictoc::toc()
}

message("\n--- Download Pipeline Complete! ---")