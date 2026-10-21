# ============================================================
# SIMPLE NPZ MODEL - TEXTBOOK VERSION
# ============================================================
#
# N = nutrients (DIN + DON)
# P = phytoplankton
# Z = zooplankton
#
# All concentrations in mmol N m^-3
# Time in days
#
# Includes:
# - Constant nutrient supply
# - Light-dependent phytoplankton growth
# - Nutrient limitation
# - Zooplankton grazing
# - Mortality
# - Remineralization
# - Export


# ============================================================
# 1. LOAD PACKAGES
# ============================================================

# Run this once if deSolve is not installed:

install.packages("deSolve")

library(deSolve)


# ============================================================
# 2. MODEL PARAMETERS
# ============================================================

parameters <- c(
  
  # Nutrient supply
  F0 = 0.05,          # mmol N m^-3 day^-1
  
  # Phytoplankton growth
  Vmax = 1.0,         # maximum growth rate (day^-1)
  K_N = 0.5,          # nutrient half-saturation
  
  # Light limitation (0 to 1)
  light_lim = 0.8,
  
  # Phytoplankton mortality
  lambda_P = 0.1,     # day^-1
  
  # Fraction of phytoplankton mortality remineralized
  mu_P = 0.7,
  
  # Zooplankton grazing
  g = 0.4,
  K_P = 1.0,
  
  # Zooplankton assimilation efficiency
  gamma_Z = 0.7,
  
  # Zooplankton mortality
  lambda_Z = 0.05,    # day^-1
  
  # Fraction of zooplankton losses remineralized
  mu_Z = 0.7
)


# ============================================================
# 3. NPZ MODEL EQUATIONS
# ============================================================

npz_model <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    # ----------------------------------------
    # Phytoplankton growth
    # ----------------------------------------
    
    nutrient_lim <- N / (K_N + N)
    
    growth <- Vmax * light_lim * nutrient_lim * P
    
    
    # ----------------------------------------
    # Zooplankton grazing
    # ----------------------------------------
    
    grazing <- g * Z * P / K_P
    
    
    # ----------------------------------------
    # Mortality
    # ----------------------------------------
    
    phyto_mortality <- lambda_P * P
    zoo_mortality <- lambda_Z * Z
    
    
    # ----------------------------------------
    # Remineralization
    # ----------------------------------------
    
    phyto_remin <- mu_P * phyto_mortality
    
    zoo_remin <- mu_Z * (
      (1 - gamma_Z) * grazing +
        zoo_mortality
    )
    
    
    # ----------------------------------------
    # Export
    # ----------------------------------------
    
    phyto_export <- (1 - mu_P) * phyto_mortality
    
    zoo_export <- (1 - mu_Z) * (
      (1 - gamma_Z) * grazing +
        zoo_mortality
    )
    
    
    # ----------------------------------------
    # Differential equations
    # ----------------------------------------
    
    dNdt <- F0 -
      growth +
      phyto_remin +
      zoo_remin
    
    dPdt <- growth -
      grazing -
      phyto_mortality
    
    dZdt <- gamma_Z * grazing -
      zoo_mortality
    
    
    return(list(c(dNdt, dPdt, dZdt)))
    
  })
}


# ============================================================
# 4. INITIAL CONDITIONS
# ============================================================

initial_conditions <- c(
  N = 10.0,
  P = 0.5,
  Z = 0.2
)


# ============================================================
# 5. TIME
# ============================================================

times <- seq(0, 150, by = 0.1)


# ============================================================
# 6. RUN MODEL
# ============================================================

solution <- ode(
  y = initial_conditions,
  times = times,
  func = npz_model,
  parms = parameters
)

# Convert to data frame
results <- as.data.frame(solution)


# ============================================================
# 7. PLOT NPZ MODEL
# ============================================================

plot(
  results$time,
  results$N,
  type = "l",
  col = "blue",
  lwd = 2,
  xlab = "Time (days)",
  ylab = "Concentration (mmol N m^-3)",
  ylim = range(results[, c("N", "P", "Z")])
)

lines(
  results$time,
  results$P,
  col = "green3",
  lwd = 2
)

lines(
  results$time,
  results$Z,
  col = "red",
  lwd = 2
)

legend(
  "topright",
  legend = c("Nutrients", "Phytoplankton", "Zooplankton"),
  col = c("blue", "green3", "red"),
  lwd = 2,
  bty = "n"
)


# ============================================================
# 8. CALCULATE EXPORT AND NITROGEN BUDGET
# ============================================================

results$grazing <- with(
  results,
  parameters["g"] * Z * P / parameters["K_P"]
)

results$phyto_export <- (
  (1 - parameters["mu_P"]) *
    parameters["lambda_P"] *
    results$P
)

results$zoo_export <- (
  (1 - parameters["mu_Z"]) *
    (
      (1 - parameters["gamma_Z"]) * results$grazing +
        parameters["lambda_Z"] * results$Z
    )
)

results$total_export <- (
  results$phyto_export +
    results$zoo_export
)

results$total_N <- (
  results$N +
    results$P +
    results$Z
)


# ============================================================
# 9. PLOT NITROGEN BUDGET AND EXPORT
# ============================================================

par(mfrow = c(2, 1))

# Total nitrogen
plot(
  results$time,
  results$total_N,
  type = "l",
  col = "black",
  lwd = 2,
  xlab = "Time (days)",
  ylab = "Total N (mmol N m^-3)",
  main = "Total nitrogen"
)


# Export fluxes
plot(
  results$time,
  results$total_export,
  type = "l",
  col = "black",
  lwd = 2,
  xlab = "Time (days)",
  ylab = "Export (mmol N m^-3 day^-1)",
  main = "Nitrogen export",
  ylim = range(
    c(
      results$phyto_export,
      results$zoo_export,
      results$total_export,
      parameters["F0"]
    )
  )
)

lines(
  results$time,
  results$phyto_export,
  col = "green3",
  lwd = 2
)

lines(
  results$time,
  results$zoo_export,
  col = "red",
  lwd = 2
)

abline(
  h = parameters["F0"],
  col = "blue",
  lty = 2,
  lwd = 2
)

legend(
  "topright",
  legend = c(
    "Total export",
    "Phytoplankton export",
    "Zooplankton export",
    "Nutrient supply"
  ),
  col = c("black", "green3", "red", "blue"),
  lty = c(1, 1, 1, 2),
  lwd = 2,
  bty = "n"
)

# Reset plotting layout
par(mfrow = c(1, 1))

