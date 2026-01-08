## assumes:
SSOM_best_2008 <- readRDS("C:/Users/OYDHAL001/Desktop/hadeda-occupancy-models/outputData/SSOM/best_models/SSOM_best_2008.rds")

 m   = SSOM_best_2008 
 umf = SSOM_best_2008@data
summary(umf)


## ---- Occupancy (ψ) on probability scale ----
psi_pred <- predict(m, type = "state")

psi_hat <- mean(psi_pred$Predicted)
psi_se  <- sd(psi_pred$Predicted)

## ---- Detection (p) on probability scale ----
p_pred <- predict(m, type = "det")

p_hat <- mean(p_pred$Predicted, na.rm = TRUE)
p_se  <- sd(p_pred$Predicted, na.rm = TRUE)
## ---- Output ----
psi_hat
psi_se
p_hat
p_se



###############################################

# Set working directory to folder with models
setwd("C:/Users/OYDHAL001/Desktop/hadeda-occupancy-models/outputData/SSOM/best_models")

# Vector of years
years <- 2008:2023

# Initialize results table
results <- data.frame(
  Year = years,
  psi_hat = NA,
  psi_se = NA,
  p_hat = NA,
  p_se = NA
)

# Loop through .rds files and extract predictions
for (i in seq_along(years)) {
  file <- paste0("SSOM_best_", years[i], ".rds")
  m <- readRDS(file)
  
  psi_pred <- predict(m, type = "state")
  p_pred   <- predict(m, type = "det")
  
  results$psi_hat[i] <- mean(psi_pred$Predicted, na.rm = TRUE)
  results$psi_se[i]  <- sd(psi_pred$Predicted, na.rm = TRUE)
  results$p_hat[i]   <- mean(p_pred$Predicted, na.rm = TRUE)
  results$p_se[i]    <- sd(p_pred$Predicted, na.rm = TRUE)
}

# Final output
results
