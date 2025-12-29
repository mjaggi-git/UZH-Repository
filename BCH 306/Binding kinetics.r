# Author: Michael Jaggi
# Course: BCH 306 - Part "Binding equilibria and ITC" (Jelezarov)

# Load necessary libraries
library(ggplot2)
library(minpack.lm)

# Define fixed protein concentration
P_conc <- 10E-06  # 10 µM

# Time points (in seconds)
time <- c(0, 0.005, 0.007, 0.008, 0.01, 0.011, 0.013, 0.02, 0.03, 0.04, 
          0.05, 0.06, 0.084, 0.08, 0.069, 0.133, 0.117, 0.23, 0.33, 0.43, 
          0.53, 0.63, 0.73, 0.83, 0.93, 1.13, 1.33, 1.53, 1.73, 1.93, 2.13, 
          2.33, 2.53, 2.73, 2.93, 3.13, 3.33, 3.53)

# Fluorescence intensity for different ligand concentrations
L1 <- c(10, 12.469, 13.198, 13.921, 14.639, 25.508, 16.949, 19.516, 23.929, 
        28.127, 32.120, 35.918, 38.958, 42.968, 46.261, 50.133, 56.960, 
        78.336, 90.795, 98.352, 102.935, 103.493, 108.282, 108.424, 108.588, 
        110.082, 108.554, 110.388, 109.097, 110.863, 107.806, 112.222, 
        109.999, 108.248, 111.339, 109.131, 108.689, 110.897)

L2 <- c(10, 12.955, 13.825, 14.687, 15.541, 16.387, 17.226, 21.308, 26.473, 
        31.337, 35.834, 40.232, 44.087, 47.891, 51.725, 54.854, 63.278, 
        84.166, 96.190, 103.731, 104.070, 107.718, 108.747, 109.313, 109.623, 
        109.437, 110.863, 108.214, 109.097, 110.897, 108.690, 110.897, 
        107.806, 112.663, 108.248, 110.000, 108.248, 111.339)

L3 <- c(10, 13.439, 14.448, 15.446, 16.434, 17.411, 18.378, 23.064, 28.942, 
        34.422, 39.468, 44.295, 47.246, 52.879, 56.349, 59.983, 69.748, 
        90.011, 100.074, 105.071, 107.552, 108.785, 109.396, 109.700, 
        110.286, 110.422, 109.980, 110.456, 108.690, 109.999, 109.131, 
        110.897, 110.456, 109.131, 110.897, 110.897, 108.690, 110.456)

L4 <- c(10, 13.921, 15.067, 16.200, 17.318, 18.424, 19.516, 24.786, 31.337, 
        37.385, 41.976, 48.122, 52.879, 57.271, 59.754, 65.067, 72.581, 
        94.118, 101.973, 106.794, 108.559, 109.353, 109.709, 109.869, 
        109.941, 109.988, 109.196, 111.220, 108.791, 110.000, 109.601, 
        110.815, 110.000, 109.196, 110.000, 110.006, 110.815, 109.196)

# Ligand concentrations in µM
L_conc <- c(50, 60, 70, 80)

# Fit function for exponential decay
fit_function <- function(time, A, k_obs, C) {
  return(A * exp(-k_obs * time) + C)
}

# Store k_obs values
k_obs_values <- numeric(length(L_conc))

# Loop through each ligand dataset and fit the model
for (i in 1:length(L_conc)) {
  data <- data.frame(Time = time, Fluorescence = get(paste0("L", i)))

  # Fit the kinetic curve using nonlinear least squares
  fit <- nlsLM(Fluorescence ~ fit_function(Time, A, k_obs, C),
               data = data,
               start = list(A = max(data$Fluorescence) - min(data$Fluorescence),
                            k_obs = 0.1, 
                            C = min(data$Fluorescence)))
  
  # Extract k_obs
  k_obs_values[i] <- coef(fit)["k_obs"]
}

# Fit a linear regression to determine ka (association rate constant)
ka_fit <- lm(k_obs_values ~ L_conc)

# Extract ka (slope of the linear regression)
ka <- coef(ka_fit)[2]

# Print results
print(paste("Association rate constant (ka):", ka, "1/(M·s)"))

# Plot k_obs vs Ligand concentration
 p1 <-   ggplot(data.frame(L_conc, k_obs_values), aes(x = L_conc, y = k_obs_values)) +
         geom_point(color = "blue", size = 3) +
         geom_smooth(method = "lm", se = FALSE, color = "red") +
         labs(title = expression("Determination of Association Rate Constant " ~ (k[a])),
             x = "Ligand Concentration (µM)",
             y = expression("Observed Rate Constant" ~ (k[obs]))) +
         theme_bw()

print(p1)
