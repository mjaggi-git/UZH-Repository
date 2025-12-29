# Author: Michael Jaggi
# Course: BCH 306 - Part "Binding equilibria and ITC" (Jelezarov)

# Load required libraries
library(ggplot2)
library(minpack.lm)

# Load data from EX3_ADH.xlsx
L_tot <- c(158, 125, 105, 68.5, 34)    # Total NADH (µM)
L_free <- c(63, 41.2, 28.3, 13.5, 4.9) # Free NADH (µM)
L_bound <- L_tot - L_free              # Bound NADH (µM)
L_ratio <- L_bound / L_free            # Scatchard linear transformation

# Create dataframe
data1 <- data.frame(L_bound, L_ratio)

# Define total ADH concentration
ADH_tot <- 30  # µM

# Perform linear regression
linear_fit <- lm(L_bound ~ L_ratio, data = data1)

# Extract KD from slope (slope = -1/KD, so KD = -1/slope)
KD <- -1 / coef(linear_fit)[2]  # Negative inverse of slope
n <- (coef(linear_fit)[1])/ADH_tot # binding sites


# Print estimated KD
cat("Estimated KD:", KD, "µM\n")
cat("Estimated number of binding sites n:", n)

# Plot linear transformation
p1 <- ggplot(data1, aes(x = L_bound, y = L_ratio)) +
      geom_point(color = "blue", size = 3) +
      geom_smooth(method = "lm", color = "red", se = FALSE) +
      xlab(expression(L[bound] ~ "(µM)")) +
      ylab(expression(L[bound] / L[free])) +
      ggtitle(expression("Scatchard-Rosenthal Linear Transformation for" ~ K[D] ~ "Estimation")) +
      theme_bw()

print(p1)

# Method a: Directly from Data
Y_direct_rel <- (L_bound)/(ADH_tot * n)

#  Method b: Using KD
Y_KD_rel <- (KD*L_free)/(1 + KD*L_free)

# Compare and Discuss Discrepancies
discrepancy <- abs(Y_direct_rel - Y_KD_rel)
cat("\nDegree of Saturation (Method a):\n")
print(Y_direct_rel)
cat("\nDegree of Saturation (Method b):\n")
print(Y_KD_rel)

cat("\nDiscrepancies between methods:\n")
print(discrepancy)
