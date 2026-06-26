# src/00_test_download.R
library(terra)
library(fs)
library(archive)
library(tictoc)
library(ecmwfr)

# 1. Source global configurations 
source("config.R")
source("./src/001_utilities.R")
# NOTE: Source the file where your custom `safe_wf_request` function lives here:
# source("path/to/your/utilities_script.R") 

message("--- Initiating Customized ERA5 Pressure Level Test ---")

zip_target <- "era5_pressure_test_snapshot.zip"

# Formulate the multi-layer pressure request using your syntax style
test_request <- list(
  dataset_short_name = era5_pressure_dataset,
  product_type       = list("reanalysis"),
  variable           = list("specific_humidity"),
  pressure_level     = list("700"),
  year               = list("2025"),
  month              = list("07"),
  day                = list("01"),
  time               = list("12:00"),
  data_format        = "netcdf",
  download_format    = "zip",
  area               = aoi_southwest,
  target             = zip_target
)

# Clear out any old debris in the temp directory
fs::dir_create(dir_temp, recurse = TRUE)
unlink(file.path(dir_temp, "*"), recursive = TRUE, force = TRUE)

message("Submitting test request via your safe_wf_request wrapper...")

# Run the download using your reliable custom wrapper
tictoc::tic()
tryCatch({
  zip_file <- safe_wf_request(
    request      = test_request,
    path         = dir_temp,
    max_attempts = 5,
    base_wait    = 180,
    jitter       = 60
  )
  tictoc::toc()
  
  # Extract the download package
  archive::archive_extract(zip_file, dir = dir_temp)
  
  # Identify the unzipped NetCDF file
  nc_files <- list.files(dir_temp, pattern = "\\.nc$", full.names = TRUE)
  
  if (length(nc_files) == 0) {
    stop("Extraction completed, but no NetCDF (.nc) file was found.")
  }
  
  # Load into terra to inspect the structural attributes
  r <- terra::rast(nc_files[1])
  
  message("\n--- TEST SUCCESSFUL ---")
  print(r) # Displays coordinate space, dimensions, and variable naming structure
  
  # === ADD THIS LINE RIGHT HERE ===
  terra::plot(r, main = "ERA5 700mb Specific Humidity Test")
  # ================================
  
  # Clean up the scratch workspace
  unlink(file.path(dir_temp, "*"), recursive = TRUE, force = TRUE)
  
}, error = function(e) {
  message("\n--- TEST FAILED ---")
  print(e$message)
})