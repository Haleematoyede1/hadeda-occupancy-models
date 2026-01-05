# ===============================================================
# Script: 03b_ssom_gof_hpc.R
# Purpose: Goodness-of-fit diagnostics for best SSOM
# Method: MacKenzie–Bailey parametric bootstrap
# Compatible with: Local machine + UCT HPC
# ===============================================================

suppressPackageStartupMessages({
  library(unmarked)
  library(dplyr)
  library(tibble)
  library(readr)
})

# -----------------------------
# 1. Parse command-line args
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

YEAR <- if (length(args) >= 1) as.integer(args[1]) else 2008
NSIM <- if (length(args) >= 2) as.integer(args[2]) else 1000
SEED <- if (length(args) >= 3) as.integer(args[3]) else 123

set.seed(SEED)

cat("Running SSOM GOF\n")
cat("Year:", YEAR, " | nsim:", NSIM, " | seed:", SEED, "\n")

# -----------------------------
# 2. Resolve project root
# -----------------------------
project_root <- normalizePath(getwd())

cat("Project root:\n", project_root, "\n")

# -----------------------------
# 3. Define paths (portable)
# -----------------------------
input_best_models <- file.path(
  project_root, "outputData", "SSOM", "best_models"
)

input_ssom_models <- file.path(
  project_root, "outputData", "SSOM", "models"
)

output_gof <- file.path(
  project_root, "outputData", "SSOM_GOF"
)

dir.create(output_gof, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 4. Fail-fast sanity checks
# -----------------------------
stopifnot(
  dir.exists(input_best_models),
  dir.exists(input_ssom_models),
  file.exists(file.path(input_best_models, paste0("SSOM_best_", YEAR, ".rds"))),
  file.exists(file.path(input_ssom_models, paste0("ssom_inputs_", YEAR, ".rds")))
)

# -----------------------------
# 5. Load best model + inputs
# -----------------------------
best_model <- readRDS(
  file.path(input_best_models, paste0("SSOM_best_", YEAR, ".rds"))
)

ssom_inputs <- readRDS(
  file.path(input_ssom_models, paste0("ssom_inputs_", YEAR, ".rds"))
)

# (UMF not strictly required, but retained for consistency)
umf <- ssom_inputs$umf

# -----------------------------
# 6. Run MacKenzie–Bailey GOF
#    (capture warnings cleanly)
# -----------------------------
warning_messages <- character(0)

mb_gof <- withCallingHandlers(
  mb.gof.test(
    mod       = best_model,
    nsim      = NSIM,
    plot.hist = FALSE
  ),
  warning = function(w) {
    warning_messages <<- c(warning_messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

# -----------------------------
# 7. Save warnings log
# -----------------------------
warnings_file <- file.path(
  output_gof,
  paste0("SSOM_GOF_warnings_", YEAR, ".txt")
)

writeLines(unique(warning_messages), warnings_file)

# -----------------------------
# 8. Create clean GOF summary
# -----------------------------
gof_summary <- tibble(
  Year           = YEAR,
  Model_Type     = mb_gof$model.type,
  Test           = "MacKenzie–Bailey GOF",
  ChiSq_Observed = mb_gof$chi.square,
  ChiSq_Mean_Sim = mean(mb_gof$t.star),
  ChiSq_SD_Sim   = sd(mb_gof$t.star),
  p_value        = mb_gof$p.value,
  c_hat          = mb_gof$c.hat.est,
  nsim           = mb_gof$nsim,
  n_warnings     = length(unique(warning_messages))
) %>%
  mutate(
    GOF_Informative = ifelse(
      p_value > 0 & is.finite(c_hat) & c_hat < 10,
      "Yes",
      "No (sparse detection histories / extreme probabilities)"
    )
  )

# -----------------------------
# 9. Save outputs
# -----------------------------
write_csv(
  gof_summary,
  file.path(output_gof, paste0("SSOM_GOF_summary_", YEAR, ".csv"))
)

saveRDS(
  mb_gof,
  file.path(output_gof, paste0("SSOM_MBGof_raw_", YEAR, ".rds"))
)

# -----------------------------
# 10. Final message
# -----------------------------
cat("GOF completed successfully for YEAR =", YEAR, "\n")
cat("Summary saved to:\n",
    file.path(output_gof, paste0("SSOM_GOF_summary_", YEAR, ".csv")), "\n")
cat("Warnings log:\n", warnings_file, "\n")
