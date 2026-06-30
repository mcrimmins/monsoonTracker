# src/13_plot_complete_climo_map.R
library(terra)
library(rnaturalearth)

# -------------------------------------------------------------------------
# 1. LOAD ALL SMOOTHED CLIMATOLOGY FILES
# -------------------------------------------------------------------------
tcwv_path <- "data/processed/climatology_tcwv_5day_smoothed.tif"
z500_path <- "data/processed/climatology_z500_5day_smoothed.tif"
q850_path <- "data/processed/climatology_q850_5day_smoothed.tif"
q700_path <- "data/processed/climatology_q700_5day_smoothed.tif"

cat("Loading smoothed climatology layers...\n")
tcwv <- rast(tcwv_path)
z500 <- rast(z500_path)
q850 <- rast(q850_path)
q700 <- rast(q700_path)

# Extract July 15th (Day 45)
target_idx  <- 90
target_name <- names(tcwv)[target_idx]

tcwv_day <- tcwv[[target_idx]]
z500_day <- z500[[target_idx]]
q850_day <- q850[[target_idx]]
q700_day <- q700[[target_idx]]

# Calculate Ridge High Center Coordinates
max_cell_idx <- which.max(values(z500_day))
high_coords  <- xyFromCell(z500_day, max_cell_idx)

# -------------------------------------------------------------------------
# 2. FETCH GEOGRAPHIC BOUNDARIES
# -------------------------------------------------------------------------
countries_vect <- vect(ne_countries(scale = "medium", continent = "North America", returnclass = "sf"))
states_vect    <- vect(ne_states(country = c("United States of America", "Mexico"), returnclass = "sf"))

# -------------------------------------------------------------------------
# 3. PLOT MASTER SYNOPTIC MAP
# -------------------------------------------------------------------------
cat("Generating complete synoptic map...\n")

# Use a vibrant, meteorologically intuitive moisture palette (Light to Deep Teal/Blue)
#moisture_pal <- colorRampPalette(c("#f7fbff", "#ebf3fb", "#bdd7e7", "#6baed6", "#3182bd", "#08519c"))(100)

moisture_pal <- colorRampPalette(c(
  "#7a5115",  # 1. Bone-Dry Desert Air (10 - 18 kg/m²)
  "#d8b365",  # 2. Modifying Continental Air / Pre-onset (18 - 25 kg/m²)
  "#f6e8c3",  # 3. The Monsoon Pivot/Threshold Zone (~25 - 28 kg/m²)
  "#c7eae5",  # 4. Shallow/Incoming Tropical Moisture (28 - 34 kg/m²)
  "#5ab4ac",  # 5. Active Convective Environment (34 - 42 kg/m²)
  "#01665e",  # 6. Deep Atmospheric Saturation (42 - 48 kg/m²)
  "#40004b"   # 7. Extreme Core Core / Tropical Surge Apex (> 48 kg/m²)
))(100)

# Expand margins slightly to prevent legend clipping in Base R
par(mar = c(4, 4, 4, 6))

# 1. Base Grid: Total Column Water Vapor
plot(tcwv_day, 
     main = sprintf("Monsoon Climatology Architecture: %s", target_name),
     col = moisture_pal,
     type = "continuous",
     range = c(5, 60),
     plg = list(title = "TCWV\n(kg/m²)", cex.title = 0.8))

# 2. Vector Layers
plot(states_vect, add = TRUE, border = "gray40", lwd = 0.4, lty = "dotted")
plot(countries_vect, add = TRUE, border = "gray10", lwd = 1.0)

# 3. 500 mb Upper-Level Steering Contours (Thin charcoal lines)
contour(z500_day, levels = seq(5700, 6000, by = 15), col = "gray20", lwd = 0.6, add = TRUE)

# 4. Moisture Boundary Isolines
#    - Red Solid Line = 850 mb specific humidity at 10 g/kg (The boundary of low-level surges)
#    - Blue Dashed Line = 700 mb specific humidity at 6 g/kg (The boundary of mountain convection)
contour(q850_day, levels = 10, col = "red", lwd = 3, add = TRUE, drawlabels = FALSE)
contour(q700_day, levels = 6, col = "blue", lwd = 3, lty = 2, add = TRUE, drawlabels = FALSE)

# 5. The Subtropical Ridge Core Engine
text(high_coords[1], high_coords[2], labels = "H", col = "darkblue", font = 2, cex = 2.0)

# 6. Comprehensive Legend
legend("bottomleft", 
       legend = c("850mb Moisture Boundary (10 g/kg)", 
                  "700mb Moisture Boundary (6 g/kg)", 
                  "500mb Height Contours",
                  "Subtropical Ridge Center"),
       col = c("red", "blue", "gray20", "darkblue"), 
       lwd = c(3, 3, 0.6, NA), 
       lty = c(1, 2, 1, NA),
       pch = c(NA, NA, NA, "H"),
       pt.cex = 1.4,
       bg = "white",
       cex = 0.75,
       box.col = "gray70")

# Reset plotting window
par(mar = c(5, 4, 4, 2) + 0.1)