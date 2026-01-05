# ===============================================================
# Script: 03b_ssom_gof_hpc.R
# Purpose: MacKenzie–Bailey GOF for best SSOM (per year)
# Runs on: UCT HPC (batch/array)
# ===============================================================

# ---- 0) User library first (HPC) ----
.libPaths(c(Sys.getenv("R_LIBS_USER", "~/R/library"), .libPaths()))

suppressPackageStartupMessages({
  library(unmarked)
  library(AICcmodavg)   # provides mb.gof.test()
})

# ---- 1) Args ----
args <- commandArgs(trailingOnly = TRUE)
YEAR <- if (length(args) >= 1) as.integer(args[1]) else 2008
NSIM <- if (length(args) >= 2) as.integer(args[2]) else 1000
SEED <- if (length(args) >= 3) as.integer(args[3]) else 123
set.seed(SEED)

cat("=========================================\n")
cat("SSOM GOF | Year:", YEAR, "| nsim:", NSIM, "| seed:", SEED, "\n")
cat("=========================================\n")

# ---- 2) Project root (use current working directory) ----
project_root <- normalizePath(getwd())
cat("Working directory:\n", project_root, "\n")

# ---- 3) Paths ----
input_best_models <- file.path(project_root, "outputData", "SSOM", "best_models")
input_ssom_models <- file.path(project_root, "outputData", "SSOM", "models")
output_gof        <- file.path(project_root, "outputData", "SSOM_GOF")

if (!dir.exists(output_gof)) dir.create(output_gof, recursive = TRUE, showWarnings = FALSE)

best_file <- file.path(input_best_models, paste0("SSOM_best_", YEAR, ".rds"))
umf_file  <- file.path(input_ssom_models, paste0("ssom_inputs_", YEAR, ".rds"))

if (!dir.exists(input_best_models)) stop("Missing dir: ", input_best_models)
if (!dir.exists(input_ssom_models)) stop("Missing dir: ", input_ssom_models)
if (!file.exists(best_file)) stop("Missing file: ", best_file)
if (!file.exists(umf_file)) stop("Missing file: ", umf_file)

# ---- 4) Load ----
best_model  <- readRDS(best_file)
ssom_inputs <- readRDS(umf_file)

# ---- 5) Run GOF (capture warnings) ----
warning_messages <- character(0)

mb_gof <- withCallingHandlers(
  mb.gof.test(mod = best_model, nsim = NSIM, plot.hist = FALSE),
  warning = function(w) {
    warning_messages <<- c(warning_messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

# ---- 6) Save warnings ----
warn_file <- file.path(output_gof, paste0("SSOM_GOF_warnings_", YEAR, ".txt"))
writeLines(unique(warning_messages), warn_file)

# ---- 7) Save summary ----
stat_obs  <- mb_gof$chi.square
stat_mean <- mean(mb_gof$t.star)
p_val     <- mb_gof$p.value
chat      <- mb_gof$c.hat.est

out_df <- data.frame(
  Year      = YEAR,
  Stat_obs  = stat_obs,
  Stat_mean = stat_mean,
  p_value   = p_val,
  chat      = chat,
  nsim      = mb_gof$nsim,
  n_warn    = length(unique(warning_messages)),
  stringsAsFactors = FALSE
)

out_csv <- file.path(output_gof, paste0("SSOM_GOF_", YEAR, ".csv"))
write.csv(out_df, out_csv, row.names = FALSE)

out_rds <- file.path(output_gof, paste0("SSOM_GOF_object_", YEAR, ".rds"))
saveRDS(mb_gof, out_rds)

cat("GOF completed for YEAR =", YEAR, "\n")
cat("Saved:\n", out_csv, "\n", out_rds, "\n")
