# -----------------------------
# Paths
# -----------------------------
effects_dir   <- "outputData/SSOM/marginal_Effect"
pred_dir      <- "outputData/SSOM/tables"
output_dir    <- "outputData/SSOM/diagnostics"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# -----------------------------
# Effect diagnostics from marginal-effect table
# -----------------------------
summarise_effect <- function(df, xvar) {
  
  df <- df %>% drop_na(fit, lo, hi)
  
  tibble(
    effect_size = max(df$fit) - min(df$fit),
    mean_ci_width = mean(df$hi - df$lo),
    slope_sign = sign(cor(df[[xvar]], df$fit))
  )
}


library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(tidyverse)

effect_files <- list.files(
  effects_dir,
  pattern = "SSOM_effect_.*\\.csv",
  full.names = TRUE
)

effect_summary <- map_dfr(effect_files, function(f) {
  
  df <- read_csv(f, show_col_types = FALSE)
  
  # var  <- str_extract(basename(f), "(NDVI|Tmax|Rf|AG|UB)")
  var <- str_extract(
    basename(f),
    regex("(ndvi|tmax|rf|ag|ub)", ignore_case = TRUE)
  )
  
  var <- toupper(var)
  
  year <- as.integer(str_extract(basename(f), "\\d{4}"))
  
  xvar <- names(df)[1]  # first column is x
  
  stats <- summarise_effect(df, xvar)
  
  tibble(
    Year = year,
    Variable = var
  ) %>% bind_cols(stats)
})

write_csv(effect_summary,
          file.path(output_dir, "SSOM_marginal_effect_diagnostics.csv"))


effect_summary <- effect_summary %>%
  mutate(
    effect_class = case_when(
      effect_size < 0.05 ~ "Weak",
      effect_size < 0.15 ~ "Moderate",
      TRUE               ~ "Strong"
    ),
    certainty = case_when(
      mean_ci_width < 0.10 ~ "High",
      mean_ci_width < 0.20 ~ "Medium",
      TRUE                 ~ "Low"
    )
  )

write_csv(effect_summary,
          file.path(output_dir, "SSOM_marginal_effect_diagnostics_classified.csv"))


############################
############################


pred_files <- list.files(
  pred_dir,
  pattern = "SSOM_Occupancy_Predictions_.*\\.csv",
  full.names = TRUE
)

prediction_summary <- map_dfr(pred_files, function(f) {
  
  df <- read_csv(f, show_col_types = FALSE)
  year <- as.integer(str_extract(basename(f), "\\d{4}"))
  
  tibble(
    Year = year,
    mean_psi = mean(df$occupancy_prob, na.rm = TRUE),
    sd_psi   = sd(df$occupancy_prob, na.rm = TRUE),
    prop_high_psi = mean(df$occupancy_prob > 0.7, na.rm = TRUE)
  )
})

write_csv(prediction_summary,
          file.path(output_dir, "SSOM_prediction_diagnostics.csv"))



final_diagnostics <- effect_summary %>%
  left_join(prediction_summary, by = "Year") %>%
  arrange(Year)

final_diagnostics <- final_diagnostics %>%
  filter(!is.na(Variable))


write_csv(final_diagnostics,
          file.path(output_dir, "SSOM_FINAL_DIAGNOSTICS.csv"))
