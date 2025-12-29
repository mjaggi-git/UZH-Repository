# Author: Michael Jaggi
# Course: BCH 306 - Part "Binding equilibria and ITC" (Jelezarov)


# Load required libraries
library(minpack.lm)
library(ggplot2)

# Data from equilibrium dialysis
L_tot <- c(1.00E-07, 1.00E-06, 3.00E-06, 5.00E-06)
L_free <- c(5.00E-08, 6.20E-07, 2.30E-06, 4.20E-06)
Y <- c(0.05, 0.38, 0.69, 0.81)

# Data from fluorescence experiment
L_tot_spectra <- c(0, 1.00E-07, 5.00E-07, 1.00E-06, 2.50E-06, 4.00E-06, 8.00E-06, 1.20E-05, 1.50E-05, 2.00E-05, 3.00E-05)
S_spectra <- c(122.29, 177.29, 364.31, 540.26, 837.14, 957.97, 1089.76, 1133.74, 1144.73, 1166.72, 1188.71)


# Define variables
L_bound <- L_tot - L_free  # Bound ligand concentration
L_ratio <- L_bound / L_free  # New y-variable

# Create dataframe
data_linear <- data.frame(L_bound, L_ratio)

# Fit linear regression model: Y = (L_free / L_bound) ~ L_bound
linear_fit <- lm(L_bound ~ L_ratio, data = data_linear)  

# Extract x-intercept as P_tot estimate
P_tot <- (coef(linear_fit)[2] * (-1))



cat("Estimated Protein Concentration (P_tot):", P_tot, "\n")

# Plot linear transformation
p1 <- ggplot(data_linear, aes(x = L_bound, y = L_ratio)) +
      geom_point(color = "blue", size = 3) +
      geom_smooth(method = "lm", formula = y ~ x, color = "red", se = FALSE) +
      xlab(expression(L[bound])) +
      ylab(expression(L[bound]/L[free])) +
      ggtitle("Scatchard-Rosenthal Linear Transformation for P_tot Estimation") +
      theme_bw()

# print(p1)

# Define Binding Model (from script page 63)
binding_model <- function(KD, S_max, L_tot_spectra, P_tot) {
  term <- P_tot + L_tot_spectra + KD - sqrt((P_tot + L_tot_spectra + KD)^2 - 4 * P_tot * L_tot_spectra)
  Y_fit <- (S_max / (2 * P_tot)) * term
  return(Y_fit)
}

# Fit Non-Linear Regression
nls_fit <- nlsLM(S_spectra ~ binding_model(KD, S_max, L_tot_spectra, P_tot),
                 start = list(KD = 1e-6, S_max = max(S_spectra)), # KD value taken form xlsx file, approx. half of the integral value
                 control = nls.lm.control(maxiter = 500))

# Extract estimated parameters
KD_estimated <- coef(nls_fit)["KD"]
S_max_estimated <- coef(nls_fit)["S_max"]

cat("Estimated Dissociation Constant (K_D):", KD_estimated, "\n")
cat("Estimated Max Signal Change (S_max):", S_max_estimated, "\n")

# Step 4: Plot Data and Fit
L_smooth <- seq(min(L_tot_spectra), max(L_tot_spectra), length.out = 100)
S_fit <- binding_model(KD_estimated, S_max_estimated, L_smooth, P_tot_estimated)

p2 <- ggplot(data.frame(L_tot_spectra, S_spectra), aes(x = L_tot_spectra, y = S_spectra)) +
      geom_point(color = "blue", size = 3) +
      geom_line(data = data.frame(L_smooth, S_fit), aes(x = L_smooth, y = S_fit), 
                color = "red", linewidth = 1) +
      xlab(expression("Total Ligand Concentration " * L[tot] * " (M)")) +
      ylab("ΔS") +
      ggtitle("Non-linear regression fit") +
      theme_bw()

print(p2)
