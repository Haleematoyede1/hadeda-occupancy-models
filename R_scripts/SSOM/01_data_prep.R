# ===============================================================
# 01_data_prep.R
# Prepare SSOM inputs for Each year/season 
# ===============================================================

# -----------------------------
# 1. Load packages
# -----------------------------
library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(unmarked)
library(tibble)

# -----------------------------
# 2. User-defined year 2008-2023
# -----------------------------
YEAR <- 2008   # <<< CHANGE THIS ONLY (REPEAT FOR THE DURATION 2008-2023)

# -----------------------------
# 3. File paths
# -----------------------------
input_dir  <- "inputData/processedData"
output_dir <- "outputData/SSOM/models"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -----------------------------
# 4. Load and clean SABAP2 data
# -----------------------------
sabap <- read_csv(file.path(input_dir,
                            "Hadeda_SABAP2_DetectionData_2008_2023.csv"),
                  show_col_types = FALSE) %>%
  select(StartDate, EndDate, StartTime, CardNo, ObserverNo,
         Pentad, Spp, TotalSpp, TotalHours) %>%
  mutate(
    Spp = case_when(
      Spp == "84" ~ 1L,
      Spp == "-"  ~ 0L,
      TRUE        ~ NA_integer_
    ),
    StartDate  = as.Date(StartDate),
    EndDate    = as.Date(EndDate),
    SurveyDays = as.numeric(difftime(EndDate, StartDate, units = "days")),
    Year       = year(StartDate)
  ) %>%
  filter(
    Year == YEAR,
    TotalHours >= 2,
    SurveyDays <= 5,
    !is.na(TotalHours)
  ) %>%
  group_by(Pentad, StartDate, ObserverNo) %>%
  arrange(desc(TotalHours), desc(TotalSpp), desc(StartTime),
          .by_group = TRUE) %>%
  slice(1L) %>%
  ungroup() %>%
  drop_na()

# -----------------------------
# 5. Join site-level covariates
# -----------------------------
climate <- read_csv(file.path(input_dir,
                              "Hadeda_ClimateCovariates_Data_Annual.csv"),
                    show_col_types = FALSE) %>%
  filter(Year == YEAR)

road <- read_csv(file.path(input_dir,
                           "Hadeda_Accessibility_Data.csv"),
                 show_col_types = FALSE) %>%
  rename(Road = RoadDistance)

analysis_data <- sabap %>%
  left_join(climate, by = c("Pentad", "Year")) %>%
  left_join(road, by = "Pentad") %>%
  drop_na() %>%
  select(-StartTime, -EndDate, -CardNo, -ObserverNo)

# -----------------------------
# 6. Transform covariates
# -----------------------------
log_vars    <- c("TotalHours", "TotalSpp", "Rf", "Soil", "Road")
spline_vars <- c("Tmax", "Tmin", "Rf", "Soil", "ndvi")

numeric_vars <- analysis_data %>%
  select(TotalSpp, TotalHours,Road, Tmax, Tmin, Rf, Soil, ndvi, NV, WB, AG, UB, OTH) %>%
  names()

transform_info <- list()

for (v in numeric_vars) {
  x <- analysis_data[[v]]
  if (v %in% log_vars) {
    x <- log(x + 1)
    log_applied <- TRUE
  } else {
    log_applied <- FALSE
  }
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  analysis_data[[paste0(v, "_z")]] <- (x - m) / s
  transform_info[[v]] <- list(mean = m, sd = s, log = log_applied)
}

# ---- spline basis 
sortBase <- function(vec, n.knots = 2) {
  knots <- quantile(unique(vec),
                    probs = seq(0, 1, length.out = n.knots + 2)[-c(1, n.knots + 2)],
                    na.rm = TRUE)
  ZK <- abs(outer(vec, knots, "-"))^3
  OM <- abs(outer(knots, knots, "-"))^3
  sv <- svd(OM)
  sqrtOM <- t(sv$v %*% (t(sv$u) * sqrt(sv$d)))
  t(solve(sqrtOM, t(ZK)))
}

for (v in spline_vars) {
  z <- analysis_data[[paste0(v, "_z")]]
  base <- scale(sortBase(z))
  analysis_data[[paste0(v, "_spline1")]] <- base[, 1]
  analysis_data[[paste0(v, "_spline2")]] <- base[, 2]
}

# -----------------------------
# 7. Build detection & covariate matrices
# -----------------------------
y <- analysis_data %>%
  select(Pentad, StartDate, Spp) %>%
  arrange(Pentad, StartDate) %>%
  group_by(Pentad) %>%
  mutate(visit = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    id_cols   = Pentad,          # <<<<<< THIS IS THE FIX
    id_cols   = Pentad,          
    names_from = visit,
    values_from = Spp,
    values_fill = list(Spp = 0)
  ) %>%
  column_to_rownames("Pentad") %>%
  as.matrix()


make_obs_matrix <- function(var) {
  analysis_data %>%
    select(Pentad, StartDate, value = all_of(var)) %>%
    arrange(Pentad, StartDate) %>%
    group_by(Pentad) %>%
    mutate(visit = row_number()) %>%
    ungroup() %>%
    pivot_wider(
      id_cols   = Pentad,        # 
      id_cols   = Pentad,        
      names_from = visit,
      values_from = value
    ) %>%
    column_to_rownames("Pentad") %>%
    as.matrix()
}


obs_covs <- list(
  TotalHours_z = make_obs_matrix("TotalHours_z"),
  TotalSpp_z   = make_obs_matrix("TotalSpp_z")
)

site_covs <- analysis_data %>%
  select(Pentad,
         Tmax_spline1, Tmax_spline2, Tmin_spline1, Tmin_spline2,
         Rf_spline1, Rf_spline2, Soil_spline1, Soil_spline2,
         ndvi_spline1, ndvi_spline2,
         NV_z, AG_z, UB_z, WB_z, OTH_z, Road_z) %>%
  distinct(Pentad, .keep_all = TRUE) %>%
  column_to_rownames("Pentad")

# -----------------------------
# 8. Create unmarked frame
# -----------------------------
umf <- unmarkedFrameOccu(
  y        = y,
  siteCovs = site_covs,
  obsCovs  = obs_covs
)

summary(umf)

# -----------------------------
# 9. Save and stop
# -----------------------------
saveRDS(
  list(
    umf            = umf,
    y              = y,
    site_covs      = site_covs,
    obs_covs       = obs_covs,
    transform_info = transform_info,
    year           = YEAR
  ),
  file = file.path(output_dir,
                   paste0("ssom_inputs_", YEAR, ".rds"))
)

message("SSOM inputs prepared for year ", YEAR)
