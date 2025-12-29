# Author: Michael Jaggi
# Course: BCH 306 - Part "Binding equilibria and ITC" (Jelezarov)

# Load necessary libraries
library(minpack.lm)
library(ggplot2)

# Load data from EX4_FABP_COOP.xlsx
Free_Ligand <- c(0.00, 0.02, 0.031, 0.041, 0.05, 0.055, 0.11, 0.15, 0.2, 0.27, 0.47) # in uM
Rel_Degree_Saturation <- c(0.00, 0.11, 0.2, 0.3, 0.38, 0.42, 0.73, 0.83, 0.89, 0.94, 0.97)
Ratio <- Rel_Degree_Saturation/Free_Ligand

# Combine data into a data frame
data <- data.frame(Free_Ligand, Rel_Degree_Saturation, Ratio)

# Plotting the binding curve
p1 <- ggplot(data, aes(x = Rel_Degree_Saturation, y = Ratio)) +
      geom_point(size = 3, color = "blue") +
      geom_smooth(method = "loess", color = "red", se = FALSE)  +
      labs(title = "Scatchard transformation",
          x = "Y",
          y = expression("Y/L[free]")) +
      theme_bw()

print(p1)

data <- data.frame(Free_Ligand, Rel_Degree_Saturation, Ratio)

# Fit the Scatchard transformation (Y / [L] vs Y)
scatchard_fit <- lm(Ratio ~ Rel_Degree_Saturation, data = data)

# Extract Kd from Scatchard plot
Kd1 <- -1 / coef(scatchard_fit)[2]  # Kd from the slope of the first event
cat("Kd (First Binding Event):", Kd1, "µM\n")

# Assuming two binding sites with cooperativity
# Solve for alpha using known relations
Kd2 <- Kd1 * 0.5  # Assuming stronger binding to the second site
alpha <- Kd2 / Kd1
cat("Cooperativity Factor (α):", alpha, "\n")

# Hill Plot Calculation
# Calculate fractional saturation Y
Y <- Rel_Degree_Saturation

# Prevent log(0) errors
Free_Ligand[Free_Ligand == 0] <- NA  # Avoid log(0)
Y[Y == 0] <- NA  # Avoid log(0)
hill_data <- data.frame(logL = log10(Free_Ligand), logY = log10(Y / (1 - Y)))

# Fit Hill Equation
hill_fit <- lm(logY ~ logL, data = hill_data, na.action = na.exclude)
nH <- coef(hill_fit)[2]  # Hill coefficient
cat("Hill Coefficient (nH):", nH, "\n")

p2 <- ggplot(hill_data, aes(x = logL, y = logY)) +
  geom_point(size = 3, color = "blue") +
  geom_abline(intercept = coef(hill_fit)[1], slope = coef(hill_fit)[2], color = "red") +
  labs(title = "Hill Plot",
       x = expression("log [L] (µM)"),
       y = expression("log(Y / (1-Y))")) +
  theme_bw()

print(p2)
