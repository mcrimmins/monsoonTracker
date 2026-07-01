# =========================================================================
# GENERATE TABBED SEASON PAGES (Daily Player / Hovmoller / 5-Day Evolution)
# =========================================================================
# Rebuilds NAM_tracker/years/<year>.qmd for every year in `years`. Each page
# is a panel-tabset with three views:
#   - Daily Synoptic Player (the original frame slider, unchanged)
#   - Season Hovmoller (08_hovmoller.R output)
#   - 5-Day Evolution: absolute + anomaly (10_pentad_year.R output)
#
# Each tab checks for its underlying image(s) at RENDER time (inside an R
# chunk), not at generation time - so as you backfill more years with
# 08/10, the existing pages will pick up the new images on next render
# without needing to rerun this script. Missing data shows a callout-note
# instead of a broken image.
#
# Run from the project root (monsoonTracker/), consistent with 01-06.
# =========================================================================

years <- 1991:2025

if (!dir.exists("NAM_tracker/years")) {
  dir.create("NAM_tracker/years", recursive = TRUE)
}

for (y in years) {

  qmd_text <- c(
    "---",
    sprintf("title: \"%d Monsoon Season\"", y),
    "---",
    "",
    sprintf("Explore the %d monsoon season through three complementary views: a day-by-day synoptic player, a season-long Hovmöller diagram, and 5-day evolution snapshots. See [Methods](methods.qmd) for how each is built.", y),
    "",
    "::: {.panel-tabset}",
    "",
    "## Daily Synoptic Player",
    "",
    "Scrub day-by-day through the season's observed TCWV and 500 hPa ridge position (left) alongside that day's anomaly vs. the 1991–2025 baseline (right). Click any frame to enlarge.",
    "",
    "```{r, echo=FALSE, results='asis'}",
    sprintf("frame_dir <- \"../output/frames/%d\"", y),
    "has_frames <- dir.exists(frame_dir) && length(list.files(frame_dir)) > 0",
    "",
    "if (has_frames) {",
    "  cat(r\"--(",
    "<div class=\"player-container\" style=\"max-width: 1000px; margin: 0 auto; text-align: center; border: 1px solid #ddd; padding: 15px; border-radius: 8px; background: #fafafa;\">",
    "  <div style=\"width: 100%; min-height: 450px; background: #eee; display: flex; align-items: center; justify-content: center; margin-bottom: 15px;\">",
    sprintf("    <img id=\"monsoon-frame\" src=\"../output/frames/%d/dashboard_%d0601.png\" onclick=\"window.open(this.src, '_blank')\" title=\"Click to view full size\" style=\"width: 100%%; height: auto; border-radius: 4px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); cursor: pointer;\" />", y, y),
    "  </div>",
    "  <div style=\"display: flex; align-items: center; gap: 15px; justify-content: center;\">",
    "    <button onclick=\"prevFrame()\" class=\"btn btn-outline-primary btn-sm\">◀ Prev</button>",
    "    <input type=\"range\" id=\"frame-slider\" min=\"0\" max=\"152\" value=\"0\" style=\"width: 60%;\" oninput=\"updateFrame(this.value)\" />",
    "    <button onclick=\"nextFrame()\" class=\"btn btn-outline-primary btn-sm\">Next ▶</button>",
    "    <button onclick=\"togglePlay()\" id=\"play-btn\" class=\"btn btn-primary btn-sm\">Play</button>",
    "  </div>",
    "</div>",
    "",
    "<script>",
    "function pad(num) { return num.toString().padStart(2, '0'); }",
    "",
    "let currentFrame = 0;",
    "let playing = false;",
    "let playInterval;",
    "",
    "const frameImg = document.getElementById('monsoon-frame');",
    "const slider = document.getElementById('frame-slider');",
    "const playBtn = document.getElementById('play-btn');",
    "",
    "function updateFrame(index) {",
    "    currentFrame = parseInt(index);",
    "    slider.value = currentFrame;",
    "    ",
    sprintf("    let targetDate = new Date(%d, 5, 1);", y),
    "    targetDate.setDate(targetDate.getDate() + currentFrame);",
    "    ",
    "    let yy = targetDate.getFullYear();",
    "    let mm = pad(targetDate.getMonth() + 1);",
    "    let dd = pad(targetDate.getDate());",
    "    ",
    "    frameImg.src = `../output/frames/${yy}/dashboard_${yy}${mm}${dd}.png`;",
    "}",
    "",
    "function nextFrame() {",
    "    currentFrame = (currentFrame + 1) % 153;",
    "    updateFrame(currentFrame);",
    "}",
    "",
    "function prevFrame() {",
    "    currentFrame = (currentFrame - 1 + 153) % 153;",
    "    updateFrame(currentFrame);",
    "}",
    "",
    "function togglePlay() {",
    "    if(playing) {",
    "        clearInterval(playInterval);",
    "        playBtn.innerText = \"Play\";",
    "    } else {",
    "        playInterval = setInterval(nextFrame, 250);",
    "        playBtn.innerText = \"Pause\";",
    "    }",
    "    playing = !playing;",
    "}",
    "</script>",
    ")--\")",
    "} else {",
    "  cat(\"\\n::: {.callout-note}\\nDaily synoptic frames have not been generated yet for this year.\\n:::\\n\")",
    "}",
    "```",
    "",
    "## Season Hovmöller",
    "",
    "Zonal-mean TCWV across the Gulf of California corridor (115–105°W), June 1–Oct 31, with the 1991–2025 climatological mean overlaid for comparison. See [Methods](methods.qmd#sec-hovmoller) for how the corridor and threshold lines are defined.",
    "",
    "```{r, echo=FALSE, results='asis'}",
    sprintf("hov_path <- \"../output/hovmoller/hovmoller_%d.png\"", y),
    "if (file.exists(hov_path)) {",
    "  cat(sprintf('![](%s){fig-align=\"center\"}', hov_path))",
    "} else {",
    "  cat(\"\\n::: {.callout-note}\\nThe Hovmöller diagram has not been generated yet for this year.\\n:::\\n\")",
    "}",
    "```",
    "",
    "## 5-Day Evolution",
    "",
    "The season compressed into 25 five-day blocks (June 10–Oct 12): the absolute mean state, and the true pixel-level anomaly vs. 1991–2025. See [Methods](methods.qmd#sec-pentad) for how the blocks are defined.",
    "",
    "##### Absolute",
    "",
    "```{r, echo=FALSE, results='asis'}",
    sprintf("pentad_abs_path <- \"../output/pentad/pentad_absolute_%d.png\"", y),
    "if (file.exists(pentad_abs_path)) {",
    "  cat(sprintf('![](%s){fig-align=\"center\"}', pentad_abs_path))",
    "} else {",
    "  cat(\"\\n::: {.callout-note}\\nThe absolute 5-day evolution grid has not been generated yet for this year.\\n:::\\n\")",
    "}",
    "```",
    "",
    "##### Anomaly",
    "",
    "```{r, echo=FALSE, results='asis'}",
    sprintf("pentad_anom_path <- \"../output/pentad/pentad_anomaly_%d.png\"", y),
    "if (file.exists(pentad_anom_path)) {",
    "  cat(sprintf('![](%s){fig-align=\"center\"}', pentad_anom_path))",
    "} else {",
    "  cat(\"\\n::: {.callout-note}\\nThe anomaly 5-day evolution grid has not been generated yet for this year.\\n:::\\n\")",
    "}",
    "```",
    "",
    ":::"
  )

  # Define the exact file path
  file_path <- sprintf("NAM_tracker/years/%d.qmd", y)

  # Write the file directly to your working directory
  writeLines(qmd_text, file_path)

  # Print status to the console so you know it's working
  cat(sprintf("Generated: %s\n", file_path))
}

cat("Success! Tabbed season pages generated for 1991-2025.\n")
