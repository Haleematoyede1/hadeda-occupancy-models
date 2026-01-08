# ===================================================================
#   DYNAMIC OCCUPANCY + NONPARAMETRIC BOOTSTRAP SCRIPT
#   Hadeda Ibis (SABAP2) – UCT HPC
# ===================================================================
input_dir  <- "inputData/processedData"
output_dir <- "outputData/DOCM/models"
# ---------------------------
# 0. Load packages
# ---------------------------
library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(unmarked)
library(purrr)
library(tibble)
library(rlang)

set.seed(123)

# ===================================================================
# 1. Load SABAP2 data from CSV
# ===================================================================
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
    #Year == YEAR,
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


Ha_Ibis <- sabap

# Remove incomplete rows and problematic Pentad
Ha_Ibis <- Ha_Ibis %>% drop_na() %>% filter(Pentad != "2855_3155")

# ===================================================================
# 7. Subsample to max 20 visits per Pentad-Year
# ===================================================================

Ha_Ibis <- Ha_Ibis %>%
  group_by(Pentad, Year) %>%
  group_modify(~ slice_sample(.x, n = min(20, nrow(.x)))) %>%
  ungroup()

# ===================================================================
# 8. Add road covariate
# ===================================================================
roadCov <- read_csv(file.path(input_dir,"Hadeda_Accessibility_Data.csv"),
                    show_col_types = FALSE)
Ha_Ibis <- Ha_Ibis %>% left_join(roadCov, by = "Pentad") %>% drop_na()

# ===================================================================
# 9. Create detection covariate dataset BEFORE coding matrices
# ===================================================================
Ha_Ibis_Final <- Ha_Ibis %>%
  select(Pentad, Year,StartDate, Spp,
         TotalHours, TotalSpp, RoadDistance)

# ===================================================================
# 10. Stabilise detection covariates
# ===================================================================
Ha_Ibis_Final <- Ha_Ibis_Final %>%
  mutate(
    TotalSpp_log    = log(TotalSpp + 1),
    TotalSpp_log_Z  = as.numeric(scale(TotalSpp_log)),
    TotalHours_Z    = as.numeric(scale(TotalHours)),
    RoadDistance_Z  = as.numeric(scale(RoadDistance))
  )

# ===================================================================
# 11. Helper functions for matrices
# ===================================================================
expand_visits <- function(data, max_visits) {
  data %>%
    group_by(Pentad, Year) %>%
    mutate(Visit = row_number()) %>%
    summarise(MaxVisit = max(Visit), .groups = "drop") %>%
    rowwise() %>%
    do(tibble(
      Pentad = .$Pentad,
      Year = .$Year,
      Visit = 1:max_visits,
      ShouldBeNA = 1:max_visits > .$MaxVisit
    )) %>%
    ungroup()
}

prepare_matrix <- function(data, var, max_visits = 20) {
  var <- enquo(var)
  
  data <- data %>%
    arrange(Pentad, Year, StartDate) %>%
    group_by(Pentad, Year) %>%
    mutate(Visit = row_number()) %>%
    ungroup()
  
  all_sites <- sort(unique(data$Pentad))
  all_years <- sort(unique(data$Year))
  
  complete_grid <- expand.grid(
    Pentad = all_sites,
    Year   = all_years,
    Visit  = 1:max_visits,
    stringsAsFactors = FALSE
  )
  
  visit_grid <- expand_visits(data, max_visits)
  
  filled <- visit_grid %>%
    left_join(
      data %>% select(Pentad, Year, Visit, !!var),
      by = c("Pentad", "Year", "Visit")
    ) %>%
    mutate(!!var := ifelse(ShouldBeNA, NA, !!var)) %>%
    select(Pentad, Year, Visit, !!var)
  
  matrix_long <- complete_grid %>%
    left_join(filled, by = c("Pentad", "Year", "Visit")) %>%
    mutate(
      YearIndex = match(Year, all_years),
      ColIndex  = (YearIndex - 1) * max_visits + Visit
    ) %>%
    arrange(Pentad, Year, Visit)
  
  matrix_wide <- matrix_long %>%
    select(Pentad, ColIndex, !!var) %>%
    pivot_wider(names_from = ColIndex, values_from = !!var) %>%
    arrange(Pentad)
  
  out_matrix <- as.matrix(matrix_wide[, -1])
  rownames(out_matrix) <- matrix_wide$Pentad
  
  list(
    matrix = out_matrix,
    years  = all_years
  )
}

# ===================================================================
# 12. Apply prepare_matrix()
# ===================================================================
result_y <- prepare_matrix(Ha_Ibis_Final, var = Spp)
y <- result_y$matrix
years_vec <- result_y$years

result_spp  <- prepare_matrix(Ha_Ibis_Final, var = TotalSpp_log_Z)
TotalSpp_matrix <- result_spp$matrix

result_hours <- prepare_matrix(Ha_Ibis_Final, var = TotalHours_Z)
TotalHours_matrix <- result_hours$matrix

result_road <- prepare_matrix(Ha_Ibis_Final, var = RoadDistance_Z)
RoadDist_matrix <- result_road$matrix

# ===================================================================
# 13. Enforce NA pattern across y and obs covs
# ===================================================================
fix_obs_na <- function(x, y) {
  x[!is.na(y) & is.na(x)] <- 0
  x[is.na(y)] <- NA
  x
}

TotalSpp_matrix   <- fix_obs_na(TotalSpp_matrix,   y)
TotalHours_matrix <- fix_obs_na(TotalHours_matrix, y)
RoadDist_matrix   <- fix_obs_na(RoadDist_matrix,   y)

stopifnot(all(is.na(y) == is.na(TotalSpp_matrix)))
stopifnot(all(is.na(y) == is.na(TotalHours_matrix)))
stopifnot(all(is.na(y) == is.na(TotalHours_matrix)))
stopifnot(all(is.na(y) == is.na(RoadDist_matrix)))

obsCovs_list <- list(
  TotalSpp   = TotalSpp_matrix,
  TotalHours = TotalHours_matrix,
  Road       = RoadDist_matrix
)

# ===================================================================
# 14. Build static & dynamic site covariates
# ===================================================================
AnnualClimate_Data <- read_csv(file.path(input_dir,"Hadeda_ClimateCovariates_Data_Annual.csv"),
                               show_col_types = FALSE) %>%
  select(-OTH)
ClimateBaseline_Data <- read_csv(file.path(input_dir,"Hadeda_ClimateBaseline_Data_Static.csv"),
                                 show_col_types = FALSE)

siteCov <- AnnualClimate_Data %>%
  left_join(ClimateBaseline_Data, by = "Pentad")

# ===================================================================

# -----------------------------
# 6. Transform covariates
# -----------------------------
log_vars    <- c( "Rf","rain_base")
spline_vars <- c("tmax_base","rain_base",  "ndvi_base")

numeric_vars <- siteCov %>%
  select( tmax_base,rain_base, ndvi_base, cult_base, Tmax, Rf, ndvi, NV, WB, AG, UB) %>%
  names()

transform_info <- list()

for (v in numeric_vars) {
  x <- siteCov[[v]]
  if (v %in% log_vars) {
    x <- log(x + 1)
    log_applied <- TRUE
  } else {
    log_applied <- FALSE
  }
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  siteCov[[paste0(v, "_z")]] <- (x - m) / s
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
  z <- siteCov[[paste0(v, "_z")]]
  base <- scale(sortBase(z))
  siteCov[[paste0(v, "_spline1")]] <- base[, 1]
  siteCov[[paste0(v, "_spline2")]] <- base[, 2]
}

# FINAL siteCov_Final
siteCov_Final <- siteCov %>%
  select(Pentad, Year, ends_with("_Z"), ends_with("_base_spline1"), ends_with("_base_spline2") ) 

# ===================================================================
# 15. Extract static covariates
# ===================================================================
unique_pentads <- Ha_Ibis_Final %>% distinct(Pentad)

static_covars <- siteCov_Final %>%
  filter(Pentad %in% unique_pentads$Pentad) %>%
  distinct(Pentad, .keep_all = TRUE) %>%
  select(Pentad, ends_with("_base_spline1"), ends_with("_base_spline2"), ends_with("_base_z"))

static_covs_matrix <- static_covars %>%
  column_to_rownames("Pentad") %>%
  as.data.frame()

# ===================================================================
# 16. Extract dynamic yearly covariates
# ===================================================================
yearly_covars <- siteCov_Final %>%
  select(Pentad, Year,ends_with("_Z"), -ends_with("_base_z")) %>%
  filter(Pentad %in% unique_pentads$Pentad) %>%
  arrange(Pentad, Year)

dynamic_vars <- setdiff(names(yearly_covars), c("Pentad", "Year"))

yearly_siteCov_list <- lapply(dynamic_vars, function(v) {
  tmp <- yearly_covars %>%
    select(Pentad, Year, all_of(v)) %>%
    pivot_wider(
      names_from  = Year,
      values_from = all_of(v)
    ) %>%
    arrange(Pentad)
  
  M <- as.matrix(tmp[, -1])
  rownames(M) <- tmp$Pentad
  return(M)
})

names(yearly_siteCov_list) <- dynamic_vars

# ===================================================================
# 17. Build UMF
# ===================================================================
num_primary <- length(years_vec)
stopifnot(num_primary == 16)  # expected 2008–2023

umf <- unmarkedMultFrame(
  y              = y,
  siteCovs       = static_covs_matrix,
  yearlySiteCovs = yearly_siteCov_list,
  obsCovs        = obsCovs_list,
  numPrimary     = num_primary
)

# Save UMF Data
saveRDS(umf, "~/PROJECT_HADEDA/Outputs/umfDOCM.rds")
# ===================================================================
# 18. Fit FINAL MODEL (ndviCultMod)
# ===================================================================
# -------------------------
# Shared formulas
# -------------------------
psi_base <- ~ tmax_base_spline1 + tmax_base_spline2 +
  rain_base_spline1 + rain_base_spline2 +
  ndvi_base_spline1 + ndvi_base_spline2 +
  cult_base_z

psi_base2 <- ~ tmax_base_z +  rain_base_z +
  ndvi_base_z + cult_base_z 

p_effort <- ~ TotalSpp + TotalHours + Road
# ===================================================================

M0_null <- colext(
  psiformula = psi_base,
  gammaformula = ~ 1,
  epsilonformula = ~ 1,
  pformula = p_effort,
  data = umf
)


M1 <- colext(
  psiformula = psi_base,
  gammaformula   = ~ Tmax_z + Rf_z,
  epsilonformula = ~ Tmax_z + Rf_z,
  pformula = p_effort,
  data = umf
)#climate driven model


M2 <- colext(
  psiformula = psi_base,
  gammaformula   = ~ Rf_z,              # colonisation tracks rainfall opportunity
  epsilonformula = ~ Tmax_z + Rf_z,      # extinction tracks heat + drought
  pformula = p_effort,
  data = umf
)#climate driven model 2


M3 <- colext(
  psiformula = psi_base,
  gammaformula   = ~ ndvi_z,
  epsilonformula = ~ ndvi_z,
  pformula = p_effort,
  data = umf
)#productivity driven model


M4 <- colext(
  psiformula = psi_base,
  gammaformula   = ~ AG_z + UB_z,
  epsilonformula = ~ AG_z + UB_z,
  pformula = p_effort,
  data = umf
)#landuse driven model


M5 <- colext(
  psiformula = psi_base,
  gammaformula   = ~ ndvi_z + AG_z + UB_z,  # colonisation: resources + land-use
  epsilonformula = ~ Tmax_z + Rf_z,                 # extinction: climate stress
  pformula = p_effort,
  data = umf
)#human influenced ext_clim

M5b <- colext(
  psiformula = psi_base2,
  gammaformula   = ~ ndvi_z + AG_z + UB_z,  # colonisation: resources + land-use
  epsilonformula = ~ Tmax_z + Rf_z,                 # extinction: climate stress
  pformula = p_effort,
  data = umf
)#human influenced ext_clim

M6 <- colext(
  psiformula = psi_base,
  gammaformula   = ~ Tmax_z + Rf_z + ndvi_z,
  epsilonformula = ~ Tmax_z + Rf_z + ndvi_z,
  pformula = p_effort,
  data = umf
)#clim_ndvi_sym


Mod_list <- fitList(M0_null  = M0_null , M1 = M1, M2 = M2,  M3 = M3, M4 =M4, M5=M5, M5b=M5b, M6=M6 )
modSel(Mod_list)


summary(M5b)

# ====================================================================================
# ====================================================================================
# ====================================================================================

#model parameter estimates summary table 
library(unmarked)

# 1. Extract raw coefficients and Standard Errors
# ------------------------------------------------
# Get coefficients (Estimates)
coefs <- coef(M5b)
# Get Standard Errors
ses <- SE(M5b)

# 2. Create a Data Frame
# ------------------------------------------------
results_table <- data.frame(
  Parameter = names(coefs),
  Estimate = as.numeric(coefs),
  SE = as.numeric(ses)
)

# 3. Calculate Z-scores and P-values
# ------------------------------------------------
results_table$z_score <- results_table$Estimate / results_table$SE
results_table$P_value <- 2 * pnorm(-abs(results_table$z_score))

# 4. Add the "Transformed" Column (Probabilities & Odds Ratios)
# ------------------------------------------------
# Logic: If it's an Intercept, convert to Probability. If it's a Slope, convert to Odds Ratio.
results_table$Transformed <- NA # Initialize column

for(i in 1:nrow(results_table)) {
  if(grepl("Intercept", results_table$Parameter[i])) {
    # Inverse Logit for Intercepts (Probability)
    results_table$Transformed[i] <- plogis(results_table$Estimate[i])
  } else {
    # Exponentiate for Slopes (Odds Ratio)
    results_table$Transformed[i] <- exp(results_table$Estimate[i])
  }
}

# 5. Clean Up and Format for Publication
# ------------------------------------------------
# Round numerical columns for readability
results_table$Estimate <- round(results_table$Estimate, 3)
results_table$SE <- round(results_table$SE, 3)
results_table$z_score <- round(results_table$z_score, 2)
results_table$Transformed <- round(results_table$Transformed, 3)

# Format P-values: specific formatting for very small numbers
results_table$P_value_Clean <- ifelse(results_table$P_value < 0.001, "<0.001", round(results_table$P_value, 3))

# Add Significance Stars (Optional but helpful)
results_table$Significance <- ifelse(results_table$P_value < 0.001, "***",
                                     ifelse(results_table$P_value < 0.01, "**",
                                            ifelse(results_table$P_value < 0.05, "*", "ns")))

# 6. View and Export
# ------------------------------------------------
print(results_table)

# Save to CSV so you can copy-paste into Word/Excel
write.csv(results_table, "Final_Model_Parameters_M5b.csv", row.names = FALSE)
# ====================================================================================
# ====================================================================================
# ====================================================================================

# fitted relationship plot for model M5b
#occupancy plots

##marginal effect plots for occupancy ndvi
ndvi_seq <- seq(min(siteCov$ndvi_base, na.rm = TRUE), max(siteCov$ndvi_base, na.rm = TRUE),length.out = 100)
#backtransform 
scale_obj <- scale(siteCov$ndvi_base)
center_vals <- attr(scale_obj, "scaled:center")
scale_vals  <- attr(scale_obj, "scaled:scale")
ndvi_cov <-  (ndvi_seq - center_vals) / scale_vals

newdat_ndvi <- data.frame(
  tmax_base_z = mean(siteCov$tmax_base_z, na.rm = TRUE),
  rain_base_z  = mean(siteCov$rain_base_z, na.rm = TRUE),
  cult_base_z = mean(siteCov$cult_base_z, na.rm = TRUE),
  ndvi_base_z = ndvi_cov
)

pred_ndvi <- predict(M5b, type = "psi", newdata = newdat_ndvi, appendData = TRUE)

ndvi_plot <- data.frame(
  ndvi = ndvi_seq,
  fit  = pred_ndvi$Predicted,
  lo   = pmax(0, pred_ndvi$lower),
  hi   = pmin(1, pred_ndvi$upper)
)

# 6) Plot
p1 <- ggplot(ndvi_plot, aes(x = ndvi, y = fit)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3) +
  geom_line(linewidth = 1) +
  labs(
    x = "NDVI ",
    y = expression(paste("Occupancy (", psi, ")")),
    title = "NDVI"
  ) +
  theme_classic()+
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 15, face = "bold")
  )



##marginal effect plots for occupancy tmax
tmax_seq <- seq(min(siteCov$tmax_base, na.rm = TRUE), max(siteCov$tmax_base, na.rm = TRUE),length.out = 100)
#backtransform 
scale_obj <- scale(siteCov$tmax_base)
center_vals <- attr(scale_obj, "scaled:center")
scale_vals  <- attr(scale_obj, "scaled:scale")
tmax_cov <-  (tmax_seq - center_vals) / scale_vals

newdat_tmax <- data.frame(
  ndvi_base_z = mean(siteCov$ndvi_base_z, na.rm = TRUE),
  rain_base_z  = mean(siteCov$rain_base_z, na.rm = TRUE),
  cult_base_z = mean(siteCov$cult_base_z, na.rm = TRUE),
  tmax_base_z = tmax_cov
)

pred_tmax <- predict(M5b, type = "psi", newdata = newdat_tmax, appendData = TRUE)

tmax_plot <- data.frame(
  tmax = tmax_seq,
  fit  = pred_tmax$Predicted,
  lo   = pmax(0, pred_tmax$lower),
  hi   = pmin(1, pred_tmax$upper)
)

# 6) Plot
p2 <- ggplot(tmax_plot, aes(x = tmax, y = fit)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3) +
  geom_line(linewidth = 1) +
  labs(
    x = "Tmax ",
    y = expression(paste("Occupancy (", psi, ")")),
    title = "Temperature"
  ) +
  theme_classic()+
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 15, face = "bold")
  )

##marginal effect plots for occupancy cult
cult_seq <- seq(min(siteCov$cult_base, na.rm = TRUE), max(siteCov$cult_base, na.rm = TRUE),length.out = 100)
#backtransform 
scale_obj <- scale(siteCov$cult_base)
center_vals <- attr(scale_obj, "scaled:center")
scale_vals  <- attr(scale_obj, "scaled:scale")
cult_cov <-  (cult_seq - center_vals) / scale_vals

newdat_cult <- data.frame(
  ndvi_base_z = mean(siteCov$ndvi_base_z, na.rm = TRUE),
  rain_base_z  = mean(siteCov$rain_base_z, na.rm = TRUE),
  tmax_base_z = mean(siteCov$tmax_base_z, na.rm = TRUE),
  cult_base_z = cult_cov
)

pred_cult <- predict(M5b, type = "psi", newdata = newdat_cult, appendData = TRUE)

cult_plot <- data.frame(
  cult = cult_seq,
  fit  = pred_cult$Predicted,
  lo   = pmax(0, pred_cult$lower),
  hi   = pmin(1, pred_cult$upper)
)

# 6) Plot
p3 <- ggplot(cult_plot, aes(x = cult, y = fit)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3) +
  geom_line(linewidth = 1) +
  labs(
    x = "AG ",
    y = expression(paste("Occupancy (", psi, ")")),
    title = "Cultivated land use"
  ) +
  theme_classic()+
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 15, face = "bold")
  )


##marginal effect plots for occupancy Rainfall
# Raw rain sequence (for plotting later)
rain_seq <- seq(min(siteCov$rain_base, na.rm = TRUE), max(siteCov$rain_base, na.rm = TRUE),length.out = 100)

# Apply SAME transformation as model
rain_log <- log(rain_seq + 1)

# Extract scaling parameters
scale_obj  <- scale(log(siteCov$rain_base + 1))
center_val <- attr(scale_obj, "scaled:center")
scale_val  <- attr(scale_obj, "scaled:scale")

# Correct scaled covariate for prediction
rain_cov <- (rain_log - center_val) / scale_val

newdat_rain <- data.frame(
  ndvi_base_z = mean(siteCov$ndvi_base_z, na.rm = TRUE),
  cult_base_z  = mean(siteCov$cult_base_z, na.rm = TRUE),
  tmax_base_z = mean(siteCov$tmax_base_z, na.rm = TRUE),
  rain_base_z = rain_cov
)

pred_rain <- predict(M5b, type = "psi", newdata = newdat_rain, appendData = TRUE)

rain_plot <- data.frame(
  rain = rain_seq,
  fit = pred_rain$Predicted,
  lo   = pmax(0, pred_rain$lower),
  hi   = pmin(1, pred_rain$upper)
)

# 6) Plot
p4 <- ggplot(rain_plot, aes(x = rain, y = fit)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3) +
  geom_line(linewidth = 1) +
  labs(
    x = "Rf ",
    y = expression(paste("Occupancy (", psi, ")")),
    title = "Rainfall"
  ) +
  theme_classic()+
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 15, face = "bold")
  )


#colonisation plots
##marginal effect plots for colonisation AG
AG_seq <- seq(min(siteCov$AG_z, na.rm = TRUE), max(siteCov$AG_z, na.rm = TRUE),length.out = 100)
#backtransform 
scale_obj <- scale(siteCov$AG_z)
center_vals <- attr(scale_obj, "scaled:center")
scale_vals  <- attr(scale_obj, "scaled:scale")
AG_cov <-  (AG_seq - center_vals) / scale_vals

newdat_AG <- data.frame(
  ndvi_z = mean(siteCov$ndvi_z, na.rm = TRUE),
  UB_z  = mean(siteCov$UB_z, na.rm = TRUE),
  AG_z = AG_cov
)

pred_AG <- predict(M5b, type = "col", newdata = newdat_AG, appendData = TRUE)

AG_plot <- data.frame(
  AG = AG_seq,
  fit  = pred_AG$Predicted,
  lo   = pmax(0, pred_AG$lower),
  hi   = pmin(1, pred_AG$upper)
)

# 6) Plot
p5 <- ggplot(AG_plot, aes(x = AG, y = fit)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3) +
  geom_line(linewidth = 1) +
  labs(
    x = "AG ",
    y = expression(paste("Colonisation (", gamma, ")")),
    title = "Cultivated area"
  ) +
  theme_classic()+
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 15, face = "bold")
  )

##marginal effect plots for colonisation UB
UB_seq <- seq(min(siteCov$UB_z, na.rm = TRUE), max(siteCov$UB_z, na.rm = TRUE),length.out = 100)
#backtransform 
scale_obj <- scale(siteCov$UB_z)
center_vals <- attr(scale_obj, "scaled:center")
scale_vals  <- attr(scale_obj, "scaled:scale")
UB_cov <-  (UB_seq - center_vals) / scale_vals

newdat_UB <- data.frame(
  ndvi_z = mean(siteCov$ndvi_z, na.rm = TRUE),
  AG_z  = mean(siteCov$AG_z, na.rm = TRUE),
  UB_z = UB_cov
)

pred_UB <- predict(M5b, type = "col", newdata = newdat_UB, appendData = TRUE)

UB_plot <- data.frame(
  UB = UB_seq,
  fit  = pred_UB$Predicted,
  lo   = pmax(0, pred_UB$lower),
  hi   = pmin(1, pred_UB$upper)
)

# 6) Plot
p6 <- ggplot(UB_plot, aes(x = UB, y = fit)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3) +
  geom_line(linewidth = 1) +
  labs(
    x = "UB ",
    y = expression(paste("Colonisation (", gamma, ")")),
    title = "Urban area "
  ) +
  theme_classic()+
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 15, face = "bold")
  )



##marginal effect plots for colonisation ndvi
ndvi_seq <- seq(min(siteCov$ndvi_z, na.rm = TRUE), max(siteCov$ndvi_z, na.rm = TRUE),length.out = 100)
#backtransform 
scale_obj <- scale(siteCov$ndvi_z)
center_vals <- attr(scale_obj, "scaled:center")
scale_vals  <- attr(scale_obj, "scaled:scale")
ndvi_cov <-  (ndvi_seq - center_vals) / scale_vals

newdat_ndvi <- data.frame(
  UB_z = mean(siteCov$UB_z, na.rm = TRUE),
  AG_z  = mean(siteCov$AG_z, na.rm = TRUE),
  ndvi_z = ndvi_cov
)

pred_ndvi <- predict(M5b, type = "col", newdata = newdat_ndvi, appendData = TRUE)

ndvi_plot <- data.frame(
  ndvi = ndvi_seq,
  fit  = pred_ndvi$Predicted,
  lo   = pmax(0, pred_ndvi$lower),
  hi   = pmin(1, pred_ndvi$upper)
)

# 6) Plot
p7 <- ggplot(ndvi_plot, aes(x = ndvi, y = fit)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3) +
  geom_line(linewidth = 1) +
  labs(
    x = "ndvi ",
    y = expression(paste("Colonisation (", gamma, ")")),
    title = "NDVI"
  ) +
  theme_classic()+
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 15, face = "bold")
  )


##marginal effect plots for extinction Tmax & Rf
Tmax_seq <- seq(min(siteCov$Tmax_z, na.rm = TRUE), max(siteCov$Tmax_z, na.rm = TRUE),length.out = 100)
#backtransform 
scale_obj <- scale(siteCov$Tmax_z)
center_vals <- attr(scale_obj, "scaled:center")
scale_vals  <- attr(scale_obj, "scaled:scale")
Tmax_cov <-  (Tmax_seq - center_vals) / scale_vals

newdat_Tmax <- data.frame(
  Rf_z = mean(siteCov$Rf_z, na.rm = TRUE),
  Tmax_z = Tmax_cov
)

pred_Tmax <- predict(M5b, type = "ext", newdata = newdat_Tmax, appendData = TRUE)

Tmax_plot <- data.frame(
  Tmax = Tmax_seq,
  fit  = pred_Tmax$Predicted,
  lo   = pmax(0, pred_Tmax$lower),
  hi   = pmin(1, pred_Tmax$upper)
)

# 6) Plot
p8 <- ggplot(Tmax_plot, aes(x = Tmax, y = fit)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3) +
  geom_line(linewidth = 1) +
  labs(
    x = "Tmax ",
    y = expression(paste("Extinction (", gamma, ")")),
    title = "Temperature"
  ) +
  theme_classic()+
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 15, face = "bold")
  )


##marginal effect plots for extinction  Rf
Rf_seq <- seq(min(siteCov$Rf_z, na.rm = TRUE), max(siteCov$Rf_z, na.rm = TRUE),length.out = 100)
#backtransform 
scale_obj <- scale(siteCov$Rf_z)
center_vals <- attr(scale_obj, "scaled:center")
scale_vals  <- attr(scale_obj, "scaled:scale")
Rf_cov <-  (Rf_seq - center_vals) / scale_vals

newdat_Rf <- data.frame(
  Tmax_z = mean(siteCov$Tmax_z, na.rm = TRUE),
  Rf_z = Rf_cov
)

pred_Rf <- predict(M5b, type = "ext", newdata = newdat_Rf, appendData = TRUE)

Rf_plot <- data.frame(
  Rf = Rf_seq,
  fit  = pred_Rf$Predicted,
  lo   = pmax(0, pred_Rf$lower),
  hi   = pmin(1, pred_Rf$upper)
)

# 6) Plot
p9 <- ggplot(Rf_plot, aes(x = Rf, y = fit)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3) +
  geom_line(linewidth = 1) +
  labs(
    x = "Rf ",
    y = expression(paste("Extinction (", gamma, ")")),
    title = "Rainfall"
  ) +
  theme_classic()+
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 15, face = "bold")
  )


library(patchwork)

psi_panel <- (p1 | p2 | p3) /
  (p4 | p5 | p6) /
  (p7 | p8 | p9)


psi_panel

output_dir <- "outputData/DOCM/marginal_Effect_plot"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


ggsave(
  filename = file.path(output_dir, "DOCM_occupancy_effect_A4.png"),
  plot     = psi_panel,
  width    = 8.27,
  height   = 11.69,
  units    = "in",
  dpi      = 300
)

# #################################################################
# #################################################################
# #################################################################
# 
# 
# #save all data used in the analysis into their respective path
# models_dir <- "outputData/DOCM/models"
# dir.create(models_dir, showWarnings = FALSE, recursive = TRUE)
# 
# saveRDS(M0_null, file.path(models_dir, "M0_null.rds"))
# saveRDS(M1,      file.path(models_dir, "M1_climate.rds"))
# saveRDS(M2,      file.path(models_dir, "M2_climate_alt.rds"))
# saveRDS(M3,      file.path(models_dir, "M3_ndvi.rds"))
# saveRDS(M4,      file.path(models_dir, "M4_landuse.rds"))
# saveRDS(M5,      file.path(models_dir, "M5_full_spline.rds"))
# saveRDS(M5b,     file.path(models_dir, "M5b_linear_best.rds"))
# saveRDS(M6,      file.path(models_dir, "M6_clim_ndvi.rds"))
# 
# modsel_dir <- "outputData/DOCM/tables"
# dir.create(modsel_dir, showWarnings = FALSE, recursive = TRUE)
# 
# Mod_list <- fitList(M0_null  = M0_null , M1 = M1, M2 = M2,  M3 = M3, M4 =M4, M5=M5, M5b=M5b, M6=M6 )
# modsel_tab <- modSel(Mod_list)
# 
# 
# aic_table <- modsel_tab@Full %>%
#   select(model, nPars, AIC, delta, AICwt, cumltvWt)
# 
# write_csv(aic_table, file.path(modsel_dir, "DOCM_aicSummaryTable.csv"))
# saveRDS(Mod_list,   file.path(modsel_dir, "DOCM_model_list.rds"))
# saveRDS(modsel_tab, file.path(modsel_dir, "DOCM_model_selection.rds"))
# 
# write.csv(
#   as.data.frame(modsel_tab),
#   file.path(modsel_dir, "DOCM_model_selection.csv"),
#   row.names = FALSE
# )
# 
# best_dir <- "outputData/DOCM/best_models"
# dir.create(best_dir, showWarnings = FALSE, recursive = TRUE)
# 
# best_model_inputs <- list(
#   model          = M5b,
#   umf            = umf,
#   y              = y,
#   siteCovs       = siteCov,
#   obsCovs        = obsCovs_list,
#   yearlySiteCovs = yearly_siteCov_list)
# 
# saveRDS(best_model_inputs,
#         file.path(best_dir, "DOCM_best_model_with_inputs.rds"))
# 
# me_data_dir <- "outputData/DOCM/marginal_effect"
# dir.create(me_data_dir, showWarnings = FALSE, recursive = TRUE)
# 
# saveRDS(ndvi_plot, file.path(me_data_dir, "psi_ndvi_data.rds"))
# saveRDS(tmax_plot, file.path(me_data_dir, "psi_tmax_data.rds"))
# saveRDS(cult_plot, file.path(me_data_dir, "psi_cult_data.rds"))
# saveRDS(rain_plot, file.path(me_data_dir, "psi_rain_data.rds"))
# 
# saveRDS(AG_plot,   file.path(me_data_dir, "gamma_AG_data.rds"))
# saveRDS(UB_plot,   file.path(me_data_dir, "gamma_UB_data.rds"))
# saveRDS(ndvi_plot, file.path(me_data_dir, "gamma_ndvi_data.rds"))
# 
# saveRDS(Tmax_plot, file.path(me_data_dir, "epsilon_tmax_data.rds"))
# saveRDS(Rf_plot,   file.path(me_data_dir, "epsilon_rain_data.rds"))
# 
# 
# plot_dir <- "outputData/DOCM/marginal_Effect_plot"
# dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
# 
# ggsave(file.path(plot_dir, "psi_ndvi.png"), p1, width = 4, height = 4, dpi = 300)
# ggsave(file.path(plot_dir, "psi_tmax.png"), p2, width = 4, height = 4, dpi = 300)
# ggsave(file.path(plot_dir, "psi_cult.png"), p3, width = 4, height = 4, dpi = 300)
# ggsave(file.path(plot_dir, "psi_rain.png"), p4, width = 4, height = 4, dpi = 300)
# 
# ggsave(file.path(plot_dir, "gamma_AG.png"), p5, width = 4, height = 4, dpi = 300)
# ggsave(file.path(plot_dir, "gamma_UB.png"), p6, width = 4, height = 4, dpi = 300)
# ggsave(file.path(plot_dir, "gamma_ndvi.png"), p7, width = 4, height = 4, dpi = 300)
# 
# ggsave(file.path(plot_dir, "epsilon_tmax.png"), p8, width = 4, height = 4, dpi = 300)
# ggsave(file.path(plot_dir, "epsilon_rain.png"), p9, width = 4, height = 4, dpi = 300)
# 
# 
# 
# fig_dir <- "outputData/DOCM/figures"
# dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
# 
# ggsave(
#   filename = file.path(fig_dir, "DOCM_marginal_effects_A4.pdf"),
#   plot     = psi_panel,
#   width    = 8.27,
#   height   = 11.69,
#   units    = "in"
# )
# #################################################################
# #################################################################
# #################################################################

library(unmarked)
library(ggplot2)
library(dplyr)
library(patchwork) # For combining plots (install if needed: install.packages("patchwork"))

# ---------------------------------------------------------
# STEP 1: Extract the Projected Trajectory
# ---------------------------------------------------------
# The 'projected()' function in unmarked calculates the expected number 
# of occupied sites for each year based on your model parameters.
# ---------------------------------------------------------

# 1. Get the projected number of occupied sites (PAO: Projected Area Occupied)
# Note: 'mean' gives the expected value.
proj <- projected(M5b) 

# 2. Extract the values (Year 1 to Year 16)
# 'projected' returns a matrix; we want the second row (occupied) usually, 
# but let's use the explicit 'predict' method for the trajectory trend.
# A simpler, robust method for the *population level* trend:

# Get the raw estimates for each year
# (This sums the posterior probabilities of occupancy for all sites)
re <- ranef(M5b) 
# Sum up the posterior probabilities to get estimated N occupied per year
# This is often more accurate than using raw rates alone.
occupancy_per_year <- colSums(bup(re, stat="mean")) 

# 3. Convert to Proportion (Psi)
n_sites <- 14381 # From your data
psi_trend <- occupancy_per_year / n_sites

# ---------------------------------------------------------
# STEP 2: Calculate Growth Rate (Lambda)
# ---------------------------------------------------------
# Lambda = Psi(t+1) / Psi(t)
# Lambda > 1 = Expansion; Lambda < 1 = Contraction
# ---------------------------------------------------------

years <- 1:length(psi_trend)
lambda <- c(NA, psi_trend[2:length(psi_trend)] / psi_trend[1:(length(psi_trend)-1)])

# Create a clean data frame for plotting
plot_data <- data.frame(
  Year = years,
  Psi = psi_trend,
  Lambda = lambda
)

# ---------------------------------------------------------
# STEP 3: Create the Plots
# ---------------------------------------------------------
start_year <- 2008 

# Add a real year column to your existing plot_data
plot_data$RealYear <- start_year + (plot_data$Year - 1)


# PLOT A: Occupancy Trend (Psi over Time)
p1 <- ggplot(plot_data, aes(x = RealYear, y = Psi)) +
  geom_line(color = "steelblue", size = 1.2) + # Sea Green color
  geom_point(size = 3, color = "steelblue") +
  ylim(0, max(plot_data$Psi) * 1.2) + # Scale Y-axis nicely
  scale_x_continuous(breaks = seq(min(plot_data$RealYear), max(plot_data$RealYear), by = 2)) +
  labs(
    title = "A. Occupancy Trend",
    subtitle = "Proportion of sites occupied over 16 years",
    y = "Occupancy Probability (\u03A8)",
    x = "Year"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"))

# PLOT B: Growth Rate (Lambda over Time)
p2 <- ggplot(plot_data, aes(x = RealYear, y = Lambda)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) + # Reference line at 1
  geom_line(color = "steelblue", size = 1.2) +
  geom_point(size = 3, color = "steelblue") +
  scale_x_continuous(breaks = seq(min(plot_data$RealYear), max(plot_data$RealYear), by = 2)) +
  labs(
    title = "B. Population Growth Rate (\u03BB)",
    subtitle = "Values >1 indicate expansion; <1 indicate contraction",
    y = "Growth Rate (\u03BB)",
    x = "Year"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"))

# Combine plots side-by-side
combined_plot <- p1 / p2
print(combined_plot)


output_dir <- "outputData/DOCM/marginal_Effect_plot"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


ggsave(
  filename = file.path(output_dir, "DOCM_occupancytrendvsGrowthrate_A4.png"),
  plot     = combined_plot,
  width    = 10,
  height   = 11.69,
  units    = "in",
  dpi      = 300
)
# ---------------------------------------------------------
# STEP 4: Save the Data for your Report
# ---------------------------------------------------------
# This prints the numbers so you can put them in a table if needed
print(plot_data)

# ===============================================================

# ===================================================================
# 19. NONPARAMETRIC BOOTSTRAP (B = 100)
# ===================================================================
B <- 100

# function: refit model on bootstrap-resampled sites
fit_bootstrap_model <- function(umf, static_covs, yearly_covs, obs_covs) {
  sites <- rownames(static_covs)
  boot_sites <- sample(sites, length(sites), replace = TRUE)
  
  y_boot <- umf@y[boot_sites, , drop = FALSE]
  static_boot <- static_covs[boot_sites, , drop = FALSE]
  yearly_boot <- lapply(yearly_covs, function(x) x[boot_sites, , drop = FALSE])
  obs_boot <- lapply(obs_covs, function(x) x[boot_sites, , drop = FALSE])
  
  umf_boot <- unmarkedMultFrame(
    y              = y_boot,
    siteCovs       = static_boot,
    yearlySiteCovs = yearly_boot,
    obsCovs        = obs_boot,
    numPrimary     = umf@numPrimary
  )
  
  mod_boot <- try(
    colext(
      psiformula     = ~ tmax_base_1_Z + tmax_base_2_Z +
        rain_base_1_Z + rain_base_2_Z +
        cult_base_pct_Z + urban_base_pct_Z,
      gammaformula   = ~ Tmax_Z + ndvi_Z + AG_pct_Z,
      epsilonformula = ~ Tmax_Z + ndvi_Z + AG_pct_Z,
      pformula       = ~ TotalSpp + TotalHours + Road,
      data           = umf_boot
    ),
    silent = TRUE
  )
  
  mod_boot
}

# Matrix to store bootstrap projected mean occupancy per year
psi_boot <- matrix(NA_real_, nrow = B, ncol = num_primary)
colnames(psi_boot) <- paste0("Year", years_vec)

for (b in 1:B) {
  message("Bootstrap iteration: ", b, " of ", B)
  
  mb <- fit_bootstrap_model(
    umf         = umf,
    static_covs = static_covs_matrix,
    yearly_covs = yearly_siteCov_list,
    obs_covs    = obsCovs_list
  )
  
  if (inherits(mb, "try-error")) {
    message("  -> model failed to converge in this replicate; skipping.")
    next
  }
  
  # Use projected.mean slot (state 'occupied') for annual occupancy
  pm <- mb@projected.mean
  if (!("occupied" %in% rownames(pm))) {
    message("  -> 'occupied' row not found in projected.mean; skipping.")
    next
  }
  
  psi_b <- pm["occupied", ]
  
  # Basic sanity check on length
  if (length(psi_b) != num_primary) {
    message("  -> unexpected length of projected.mean; skipping.")
    next
  }
  
  psi_boot[b, ] <- as.numeric(psi_b)
}

# Save raw bootstrap matrix
saveRDS(psi_boot, "~/PROJECT_HADEDA/Outputs/psi_boot_ndviCultMod_B100.rds")

# ===================================================================
# 20. Summarise bootstrap (mean, 2.5%, 97.5%)
# ===================================================================
psi_mean <- apply(psi_boot, 2, mean,    na.rm = TRUE)
psi_lo   <- apply(psi_boot, 2, quantile, 0.025, na.rm = TRUE)
psi_hi   <- apply(psi_boot, 2, quantile, 0.975, na.rm = TRUE)

psi_boot_summary <- data.frame(
  Year    = years_vec,
  psi_hat = as.numeric(psi_hat),
  psi_mean_boot = as.numeric(psi_mean),
  psi_lo_boot   = as.numeric(psi_lo),
  psi_hi_boot   = as.numeric(psi_hi)
)

write_csv(psi_boot_summary,
          "~/PROJECT_HADEDA/Outputs/psi_boot_summary_ndviCultMod_B100.csv")

# ===================================================================
# END OF SCRIPT
# ===================================================================
