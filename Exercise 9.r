# Author: Michael Jaggi
# Course: BCH 306 - Part "Binding equilibria and ITC" (Jelezarov)

# Load required libraries
library(ggplot2)

# Provided experimental data
FGF_bound <- c(0.14, 0.45, 0.80, 1.06, 1.38, 1.64, 1.95, 2.46, 3.07, 3.08) # in nM
FGF_free <- c(0.05, 0.21, 0.49, 0.91, 1.98, 4.02, 8.86, 31, 300, 425) # in nM

# Non-linear regression (Langmuir Isotherm)
data <- data.frame(FGF_free, FGF_bound)
model <- nls(FGF_bound ~ Bmax * FGF_free / (KD + FGF_free), data = data,
             start = list(Bmax = max(FGF_bound), KD = 1))

# Summary of the model
summary(model)

# Extract KD and Bmax
KD <- coef(model)["KD"]
Bmax <- coef(model)["Bmax"]
cat("KD (nM):", KD, "\n")
cat("Bmax (nM):", Bmax, "\n")

# Estimating receptors per cell
Avogadro_Nr <- 6.022E+23  # Avogadro's number
receptors_per_cell <- (Bmax * Avogadro_Nr) / 1e6
cat("Number of high-affinity FGF receptors per cell:", receptors_per_cell, "\n")

# Plotting the binding curve
p1 <- ggplot(data, aes(x = FGF_free, y = FGF_bound)) +
      geom_point(color = "blue") +
      stat_function(fun = function(x) Bmax * x / (KD + x), color = "red") +
      labs(title = "Binding Curve of FGF to Receptor",
            x = "Free FGF (nM)",
            y = "Bound FGF (nM)") +
      theme_bw()

print(p1)