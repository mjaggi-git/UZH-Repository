# Author: Michael Jaggi
# Course: BCH 306 - Part "Binding equilibria and ITC" (Jelezarov)

# Load required libraries
library(ggplot2)
library(minpack.lm)

# Load data
Rel_Degree_Saturation <- c(0.025, 0.075, 0.155, 0.190, 0.200, 0.360, 0.425, 0.445, 0.494, 0.550, 0.600, 0.650, 0.600, 0.600, 0.650, 0.650, 0.700, 0.850, 0.900, 0.850, 0.895)

Free_Phd <- c(1.20E-09, 4.50E-09, 1.20E-08, 1.50E-08, 1.90E-08, 5.00E-08, 1.00E-07, 1.40E-07, 2.10E-07, 3.30E-07, 3.00E-07, 4.50E-07, 5.00E-07, 6.70E-07, 1.50E-06, 2.00E-06, 2.50E-06, 3.00E-06, 4.90E-06, 6.50E-06, 7.50E-06)

# Create data frame
data <- data.frame(Free_Phd, Rel_Degree_Saturation)

# Scatchard Plot
scatchard_plot <- ggplot(data, aes(x = Rel_Degree_Saturation, y = Rel_Degree_Saturation / Free_Phd)) +
  geom_point(size = 3, color = "blue") +
  labs(title = "Scatchard Plot: Phd Binding to Tandem Sites",
       x = "Relative Degree of Saturation (Y)",
       y = "Y / [Phd]") +
  theme_minimal()

print(scatchard_plot)

# Hill Plot
hill_data <- data.frame(logL = log10(Free_Phd), logY = log10(Rel_Degree_Saturation / (1 - Rel_Degree_Saturation)))
hill_plot <- ggplot(hill_data, aes(x = logL, y = logY)) +
  geom_point(size = 3, color = "red") +
  labs(title = "Hill Plot: Phd Binding to Tandem Sites",
       x = "log([Phd] µM)",
       y = "log(Y / (1-Y))") +
  theme_minimal()

print(hill_plot)

# Nonlinear Fitting using Cooperative Binding Model
binding_model <- function(params, Free_Phd, Rel_Degree_Saturation) {
  Kd = params[1]
  nH = params[2]
  predicted = (Free_Phd^nH) / (Kd^nH + Free_Phd^nH)
  return(predicted - Rel_Degree_Saturation)
}

# Initial guesses
params_start <- c(1e-8, 1.5)
fit <- nls.lm(par = params_start, fn = binding_model, Free_Phd = Free_Phd, Rel_Degree_Saturation = Rel_Degree_Saturation)
Kd_fit <- fit$par[1]
Hill_coeff <- fit$par[2]

cat("\nFitted Parameters:\n")
cat("Kd (Tandem Binding Sites):", Kd_fit, "M\n")
cat("Hill Coefficient (nH):", Hill_coeff, "\n")

# Comparison to Single Site KD
Kd_single <- 45e-6 # Given in problem statement
cat("\nComparison to Single-Site KD:\n")
cat("Single-Site KD:", Kd_single, "M\n")
cat("Ratio (Two-Site/Single-Site):", Kd_fit / Kd_single, "\n")
