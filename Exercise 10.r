# Author: Michael Jaggi
# Course: BCH 306 - Part "Binding equilibria and ITC" (Jelezarov)

# Given Data
Q <- 2e-7 # Total accumulated heat (kJ)
V_cell <- 1e-3 # Volume of calorimetric cell (L)
P_tot <- 10e-6 # Total receptor concentration after injection j (M)
n <- 2 # Number of binding sites
Y <- 0.5 # Degree of saturation

# Calculate molar enthalpy of binding (Delta H)
Delta_H <- Q / (V_cell * P_tot * n * Y)


cat("Molar Enthalpy of Binding (ΔH):", Delta_H, "kJ/mol\n")
