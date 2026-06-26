# Southwest Monsoon Tracker

## Project Overview
This project is designed to track, analyze, and visualize the North American Monsoon (NAM) over the desert Southwest and Northern Mexico. By processing long-term atmospheric data, the pipeline constructs climatologies and tracks key monsoon metrics, including the 500 hPa monsoon ridge (Geopotential Height) and low-level moisture surges (Specific Humidity / Precipitable Water).

## Dataset: Copernicus ERA5 Reanalysis
Data is sourced from the Copernicus Climate Data Store (CDS) using the [ERA5 hourly data on pressure levels](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels) dataset.

* **Temporal Resolution:** Hourly data collapsed into daily averages.
* **Timeframe:** 1991 - 2025 (Monsoon Season: June through October).
* **Spatial Extent:** Bounding box covering the US Southwest and Mexico (Longitude: -125 to -95, Latitude: 15 to 45).
* **Target Variables:**
    * Geopotential (`z`)
    * Specific Humidity (`q`)
* **Pressure Levels:** 250, 500, 700, and 850 hPa.

## Directory Structure
    monsoonTracker/
    ├── data/
    │   ├── raw/      # Raw NetCDF files downloaded from CDS
    │   └── temp/     # Temporary processing directory for unzipping
    ├── src/          # R scripts for downloading and processing
    ├── config.R      # Master configuration variables
    └── README.md     # Project documentation

## Scripts & Pipeline

### 1. Configuration
* **`config.R`**: Holds the master parameters for the project, including the CDS API credentials setup, target years, months, bounding box, and variables. Keeps the main scripts clean and easily adjustable.

### 2. Data Acquisition
* **`src/01_sequential_download_era5.R`**: A robust, sequential download script utilizing the new `{ecmwfr}` (v2.0+) package and R6 object methods. It dynamically handles CDS Beta API rate-limiting, polls server status, downloads the `.zip` archives, extracts the `.nc` (NetCDF) files to `data/raw/`, and automatically skips previously downloaded months to prevent duplicate requests.

### 3. Data Inspection & Validation
* **`src/02_test_plots_500mb.R`** (and related test scripts): Utilizes the `{terra}` and `{maps}` packages to inspect the 3D NetCDF structures. These scripts verify that the data layers are correctly mapped to spatial coordinates, accurately converted to standard meteorological units (e.g., meters for geopotential height, g/kg for specific humidity), and perfectly aligned with the correct pressure levels based on parsed layer names.

## Requirements
To run this pipeline, the following R packages are required:
* `ecmwfr` (v2.0+ required for the CDS Beta API R6 objects)
* `terra` (for spatial raster processing and NetCDF handling)
* `fs` (for file system operations)
* `archive` (for unzipping CDS payloads)
* `tictoc` (for script timing)
* `maps` (for geographic boundary plotting)

## Current Status & Next Steps
* **Phase 1 (Complete):** Data pipeline established. Safe, automated downloads of hourly ERA5 NetCDF files for the 1991-2025 monsoon seasons are operational. Data validation confirms spatial extent and atmospheric levels are correct.
* **Phase 2 (In Development):** Developing `{terra}` processing scripts to collapse hourly NetCDF layers into daily averages, and calculating column-integrated Precipitable Water across the vertical pressure levels.