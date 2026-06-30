# src/16_verify_anomalies.R
library(terra)

# -------------------------------------------------------------------------
# 1. SETUP DIAGNOSTIC PARAMETERS
# -------------------------------------------------------------------------
anom_dir <- "data/processed/4slices/anomalies"

test_year  <- 2025      
test_layer <- "August_14" 

variables <- c("tcwv", "q850", "q700", "q500", "q250", "z850", "z700", "z500", "z250")

# Smooth diverging palette (Brown = Dry/Low, White = Normal, Teal = Wet/High)
div_pal <- colorRampPalette(c("blue", "white", "red"))(255)

cat(sprintf("Generating clean diagnostic dashboard for %s, %d...\n", test_layer, test_year))

# -------------------------------------------------------------------------
# 2. BUILD THE CLEAN 3x3 DASHBOARD
# -------------------------------------------------------------------------
# Set up a 3x3 plotting window with comfortable margins for clean legends
par(mfrow = c(3, 3), mar = c(1.5, 1.5, 3, 4.5)) 

for (var in variables) {
  file_path <- file.path(anom_dir, sprintf("anomaly_%s_%d.tif", var, test_year))
  
  if (!file.exists(file_path)) {
    plot.new()
    title(main = paste(toupper(var), "- MISSING"), col.main = "red")
    next
  }
  
  r_day <- rast(file_path)[[test_layer]]
  
  # Find the absolute maximum deviation to keep 0 locked right in the center
  max_val <- max(abs(minmax(r_day)[]), na.rm = TRUE)
  if (is.infinite(max_val) || max_val == 0) max_val <- 1 
  
  # Plot using standard continuous range mapping instead of discrete breaks
  plot(r_day, 
       col = div_pal, 
       type = "continuous",
       range = c(-max_val, max_val), # Symmetrical scale forces white at exactly 0
       main = sprintf("%s Anomaly\n(%s)", toupper(var), test_layer), 
       mar = c(1.5, 1.5, 3, 4.5),
       axes = FALSE, 
       box = TRUE,
       plg = list(shrink = 0.9, cex = 0.8)) # Cleans up and scales down legend text
}

# Reset plotting window
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
cat("Clean dashboard rendered successfully!\n")