# ===============================================================
# 05_ssom_plot_selected_years.R
# Purpose: Plot representative SSOM occupancy maps in a grid
# ===============================================================

# -----------------------------
# 1. Load packages
# -----------------------------
library(dplyr)
library(readr)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

# -----------------------------
# 2. Define years to plot
# -----------------------------
years_to_plot <- c(2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023)

# -----------------------------
# 3. File paths
# -----------------------------
input_dir  <- "outputData/SSOM/tables"
output_fig <- "outputData/SSOM/figures"

dir.create(output_fig, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 4. Load and combine predictions
# -----------------------------
occ_all <- lapply(years_to_plot, function(yr) {
  
  read_csv(
    file.path(input_dir,
              paste0("SSOM_Occupancy_Predictions_", yr, ".csv")),
    show_col_types = FALSE
  )
  
}) %>% bind_rows()

# Ensure Year is treated as a factor (important for faceting order)
occ_all$Year <- factor(occ_all$Year, levels = years_to_plot)

# -----------------------------
# 5. Determine common colour scale
# -----------------------------
# Using global limits ensures fair comparison
fill_limits <- range(occ_all$occupancy_prob, na.rm = TRUE)

# -----------------------------
# 6. Plot multi-panel grid
# -----------------------------
p <- ggplot(
  occ_all,
  aes(x = Longitude, y = Latitude, fill = occupancy_prob)
) +
  geom_tile() +
  coord_equal() +
  scale_fill_gradientn(
    colours = brewer.pal(5, "Spectral"),
    limits  = fill_limits,
    name    = "Occupancy Probability"
  ) +
  facet_wrap(~ Year, ncol = 4) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  guides(
    fill = guide_colourbar(
      direction = "horizontal",
      barheight = unit(1.2, "cm"),
      barwidth  = unit(8, "cm"),
      title.position = "top",
      title.hjust = 0.5,
      label.position = "bottom"
    )
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.spacing = unit(0.8, "lines"),
    strip.text    = element_text(size = 12, face = "bold"),
    plot.title    = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box.margin = margin(15, 10, 15, 10),
    legend.title = element_text(
      size = 13,
      face = "bold",
      margin = margin(b = 8)
    ),
    legend.text = element_text(
      size = 11,
      margin = margin(t = 6)
    )
  )

# -----------------------------
# 7. Save figure
# -----------------------------
ggsave(
  filename = file.path(
    output_fig,
    "SSOM_Occupancy_Representative_Years_PORTRAIT.png"
  ),
  plot   = p,
  width  = 8.3,    # A4 portrait width
  height = 11.7,   # A4 portrait height
  dpi    = 300
)

message("Representative SSOM occupancy maps saved successfully.")




