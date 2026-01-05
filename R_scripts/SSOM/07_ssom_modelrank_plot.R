# ===============================================================
# Plot: Top-3 SSOM model composition per year (stacked bars)
# ===============================================================

library(dplyr)
library(readr)
library(ggplot2)
library(purrr)

# -----------------------------
# 1. File paths
# -----------------------------
aic_dir <- "outputData/SSOM/tables"
output_fig <- "outputData/SSOM/figures"

if (!dir.exists(output_fig)) dir.create(output_fig, recursive = TRUE)

# -----------------------------
# 2. Read all AIC tables
# -----------------------------
aic_files <- list.files(
  aic_dir,
  pattern = "SSOM_AIC_\\d{4}\\.csv",
  full.names = TRUE
)

aic_all <- map_dfr(aic_files, read_csv, show_col_types = FALSE)

# -----------------------------
# 3. Extract top 3 models per year
# -----------------------------
library(dplyr)
library(tidyr)
library(dplyr)

top3_ranked <- aic_all %>%
  group_by(Year) %>%
  arrange(AIC, .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  mutate(
    Rank = case_when(
      row_number() == 1 ~ 3L,
      row_number() == 2 ~ 2L,
      row_number() == 3 ~ 1L
    )
  ) %>%
  ungroup()


table(top3_ranked$Rank)
# should show roughly equal counts of 1, 2, 3
head(top3_ranked)
library(dplyr)

top3_ranked <- top3_ranked %>%
  mutate(model = recode(model,
                        "wetland_natural" = "rangeland_waterbody"))

# -----------------------------
# 5. Stacked bar plot
# -----------------------------
p_ranked <- ggplot(
  top3_ranked,
  aes(
    x    = factor(Year),
    y    = Rank,
    fill = model
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.6,
    colour = "black"
  ) +
  scale_y_continuous(
    breaks = c(1, 2, 3),
    labels = c("3rd", "2nd", "Best"),
    limits = c(0, 3.3),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Year",
    y = "Model Rank",
    fill = "Model Type"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

p_ranked

# -----------------------------
# 6. Save figure (portrait)
# -----------------------------
ggsave(
  filename = file.path(
    output_fig,
    "SSOM_Top3_Model_Composition_By_Year.png"
  ),
  plot   = p_ranked,
  width  = 8,
  height = 10,
  dpi    = 300
)

message("Top-3 SSOM model composition plot saved successfully.")


##################################
#######################################
##########################################


library(dplyr)
library(readr)
library(purrr)

# Read all AIC tables
aic_files <- list.files(
  "outputData/SSOM/tables",
  pattern = "SSOM_AIC_\\d{4}\\.csv",
  full.names = TRUE
)

aic_all <- map_dfr(aic_files, read_csv, show_col_types = FALSE)

# Top 2 models per year
library(dplyr)

top2_models <- aic_all %>%
  group_by(Year) %>%
  arrange(AIC, .by_group = TRUE) %>%
  slice_head(n = 2) %>%
  mutate(
    Rank = row_number()
  ) %>%
  ungroup() %>%
  select(
    Year,
    Rank,
    Model = model,
    nPars,
    AIC,
    DeltaAIC = delta,
    AICwt
  )

# Save table
write.csv(
  top2_models,
  "outputData/SSOM/tables/SSOM_Top2_Models_2008_2023.csv",
  row.names = FALSE
)

