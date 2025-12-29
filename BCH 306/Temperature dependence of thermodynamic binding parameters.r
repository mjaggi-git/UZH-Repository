# Author: Michael Jaggi
# Course: BCH 306 - Part "Binding Equilibria and ITC" (Jelezarov)

# Load necessary libraries
library(ggplot2)

# Step 1: Data Preparation
# Temperature in Kelvin
T <- c(278.15, 283.15, 293.15, 303.15, 313.15, 323.15)
# Binding constants from experiments (M^-1)
KA_Exp1 <- c(4.9e9, 3.0e9, 2.2e9, 8.5e8, 6.0e8, 2.5e8)
KA_Exp2 <- c(4.2e9, 3.7e9, 1.7e9, 1.2e9, 5.3e8, 3.3e8)
KA_Exp3 <- c(5.1e9, 2.9e9, 2.3e9, 9.0e8, 5.1e8, 3.5e8)
# Enthalpy from ITC in kJ/mol
DeltaH_ITC <- c(NA, -35.8, NA, -47.7, -54.1, NA)

# Calculate mean KA and ln(KA)
KA_mean <- rowMeans(cbind(KA_Exp1, KA_Exp2, KA_Exp3), na.rm = TRUE)
ln_KA <- log(KA_mean)

# Step 2: Van't Hoff Analysis
R <- 8.314 # J/mol*K
inv_T <- 1 / T

# Fit linear model ln(KA) = -DeltaH/R * (1/T) + DeltaS/R
fit <- lm(ln_KA ~ inv_T)
DeltaH <- -coef(fit)[2] * R / 1000 # in kJ/mol
DeltaS <- coef(fit)[1] * R / 1000 # in kJ/mol*K

# Step 3: Calculate Free Energy at 298.15 K and 356.15 K
DeltaG_25 <- DeltaH - 298.15 * DeltaS
DeltaG_83 <- DeltaH - 356.15 * DeltaS

# Step 4: Plot Van't Hoff Plot
plot_data <- data.frame(inv_T, ln_KA)
p1 <-    ggplot(plot_data, aes(x = inv_T, y = ln_KA)) +
          geom_point(color = "blue", size = 3) +
          geom_smooth(method = "lm", color = "red", se = FALSE) +
          labs(title = "Van't Hoff Plot",
              x = "1/T (1/K)",
              y = "ln(KA)") +
          theme_minimal()

print(p1)

# Step 5: Extrapolation and Results
T_optimal <- 356.15
KA_optimal <- exp(- (DeltaH * 1000 - DeltaS * 1000 * T_optimal) / (R * T_optimal))

cat("Results:\n")
cat("DeltaH (kJ/mol):", DeltaH, "\n")
cat("DeltaS (kJ/mol*K):", DeltaS, "\n")
cat("DeltaG at 25°C (kJ/mol):", DeltaG_25, "\n")
cat("DeltaG at 83°C (kJ/mol):", DeltaG_83, "\n")
cat("Binding constant at 83°C (M^-1):", KA_optimal, "\n")
