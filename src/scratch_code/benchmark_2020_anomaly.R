# src/14_benchmark_2020_anomaly.R
library(terra)
library(rnaturalearth)

# -------------------------------------------------------------------------
# 1. SETUP & LOAD BASELINE CLIMATOLOGY
# -------------------------------------------------------------------------
target_year <- 2020
climo_dir   <- "data/processed"          # Where your climo .tif files live
slice_dir   <- "data/processed/4slices"  # Where your raw monthly .nc files live
months_to_process <- c("06", "07", "08", "09", "10")

cat(sprintf("Benchmarking Anomaly Engine for %d...\n", target_year))

# Load our verified 153-day TCWV smoothed baseline (from the main processed folder)
climo_tcwv <- rast(file.path(climo_dir, "climatology_tcwv_5day_smoothed.tif"))

# -------------------------------------------------------------------------
# 2. EXTRACT TARGET YEAR ACTUALS (2020)
# -------------------------------------------------------------------------
actuals_list <- list()

for (m_num in months_to_process) {
  # Build our targeted filename regex
  pattern_str <- paste0("single_.*", target_year, "_", m_num, "\\.nc$")
  m_files <- list.files(slice_dir, pattern = pattern_str, full.names = TRUE)
  
  if (length(m_files) == 0) {
    stop(sprintf("Could not find a 'single' layer NetCDF for year %d, month %s in %s", 
                 target_year, m_num, slice_dir))
  }
  
  cat(sprintf(" -> Loading data from: %s (%d days)\n", basename(m_files[1]), nlyr(rast(m_files[1]))))
  
  # FIX: Since this single-level file is purely TCWV data, take all layers!
  actuals_list[[m_num]] <- rast(m_files[1])
}

# Combine into a single continuous 153-layer stack for the 2020 monsoon season
tcwv_2020_stack <- rast(actuals_list)
names(tcwv_2020_stack) <- names(climo_tcwv) # Force names to match climo ("June_01", etc.)



# -------------------------------------------------------------------------
# 3. CALCULATE THE ANOMALY GRID (Matrix Subtraction)
# -------------------------------------------------------------------------
cat("Calculating daily anomalies (Actual - Climatology)...\n")

# The magic happens here: cell-by-cell subtraction across all 153 days instantly
tcwv_anomaly_2020 <- tcwv_2020_stack - climo_tcwv

# -------------------------------------------------------------------------
# 4. VISUALIZE THE 2020 "NON-SOON" (August 15th Benchmark)
# -------------------------------------------------------------------------
cat("Rendering benchmark anomaly map...\n")

# Day 76 is August 15th - classically the peak of the monsoon
target_idx <- 76
target_name <- names(tcwv_anomaly_2020)[target_idx]

# Create a divergent palette: Brown (Dry/Negative) -> White (Normal) -> Teal (Wet/Positive)
divergent_pal <- colorRampPalette(c("#8c510a", "#d8b365", "#f6e8c3", "#f5f5f5", "#c7eae5", "#5ab4ac", "#01665e"))(100)

states_vect <- vect(ne_states(country = c("United States of America", "Mexico"), returnclass = "sf"))

par(mar = c(3, 3, 4, 6))

plot(tcwv_anomaly_2020[[target_idx]], 
     main = sprintf("TCWV Anomaly Benchmark: %s, %d\nThe 2020 'Non-soon'", target_name, target_year),
     col = divergent_pal,
     type = "continuous",
     range = c(-20, 20), # Lock the scale from -15 kg/m² to +15 kg/m² for visual symmetry
     plg = list(title = "Anomaly\n(kg/m²)", cex.title = 0.8))

plot(states_vect, add = TRUE, border = "gray20", lwd = 0.5)

par(mar = c(5, 4, 4, 2) + 0.1)