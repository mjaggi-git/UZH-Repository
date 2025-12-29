# Author: Michael Jaggi
# Course: BCH 306 - Part "Binding equilibria and ITC" (Jelezarov)

# Load required libraries
library(ggplot2)
library(minpack.lm)

# Load data
time <- c(0, 0.005, 0.05, 0.1, 0.15, 0.2, 0.23274, 0.3, 0.35851, 0.4, 0.45, 0.47881, 
          0.55, 0.57861, 0.64559, 0.7, 0.75, 0.811, 0.85, 0.87116, 0.95, 1.01743, 
          1.05, 1.07758, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 
          1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.05, 2.1, 2.15, 2.2, 2.25, 
          2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 
          2.95, 3, 3.05, 3.1, 3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 
          3.6, 3.65, 3.7, 3.75, 3.8, 3.85, 3.9, 3.95, 4, 4.05, 4.1, 4.15, 4.2, 
          4.25, 4.3, 4.35, 4.4, 4.45, 4.5, 4.55, 4.6, 4.65, 4.7, 4.75, 4.8, 4.85, 
          4.9, 4.95, 5) 

em_intensity <- c(10, 10.5, 14.7, 18.9, 22.6, 25.9, 29.4, 31.6, 33.4, 36.1, 38.1, 
                  39.8, 41.5, 42.9, 44.7, 45.6, 46.8, 47.2, 48.7, 49.7, 50.4, 50.6, 
                  51.8, 52.5, 53, 53.3, 54.5, 53.8, 55.8, 55.3, 54.7, 55.7, 56.5, 
                  56.1, 57.2, 56.6, 57.7, 56.8, 58.1, 57.8, 57.5, 59, 57.6, 58.8, 
                  58.3, 58.8, 58.5, 59.2, 58.4, 59, 59.9, 59.1, 58.7, 59.2, 59.1, 
                  59.6, 59.6, 59.2, 60.6, 59.1, 60.2, 59.4, 60.1, 59.4, 60.1, 59.9, 
                  59.3, 60.6, 59.3, 60.4, 59.1, 60.7, 60.3, 59.4, 60.5, 59.4, 60.3, 
                  59.2, 60.6, 58.8, 61, 59.2, 60.6, 59.5, 60.4, 59.5, 61, 58.6, 60.2, 
                  59.9, 60.4, 59.5, 60.8, 60, 59.3, 60.4, 59.3, 61.1, 59.5, 60.2, 
                  59.3, 60.9)

# Initial concentrations
P_init <- 10E-06 # in M
L_init <- 10E-06 # in M

# Store data in a dataframe
data1 <- data.frame(time, em_intensity)

# Define functions for equilibrium calculations
b_function <- function(kd, ka, L_init) {
    return(kd / (ka * L_init))
}

s_function <- function(b) {
    return(sqrt(4 * b + b^2))
}

z_function <- function(ka, L_init, b, s, t) {
    return(((2 + b - s) / (2 + b + s)) * exp(-ka * L_init * t))
}

get_P_and_L <- function(P_init, z, b, s) {
    return(((z * b + z * s - b + s) / (1 - z)) * (P_init / 2))
}

# Define fit function for fluorescence intensity
fit_function_S <- function(time, S_eq, S_max, ka, kd) {
    b <- b_function(kd, ka, L_init)
    s <- s_function(b)
    z <- z_function(ka, L_init, b, s, time)
    P_and_L <- get_P_and_L(P_init, z, b, s)
    return(S_eq + S_max * (1 - (P_and_L / P_init)))
}

# Fit the data using non-linear regression
nls_fit <- nlsLM(em_intensity ~ fit_function_S(time, S_eq, S_max, ka, kd),
                 start = list(S_eq = 60, S_max = 50.9, ka = 1e6, kd = 1e-3),
                 data = data1, control = nls.lm.control(maxiter = 500))

# Print the fitted parameters
print(summary(nls_fit))

# Generate fitted values for plotting
data1$fitted <- predict(nls_fit)

# Plot the data and fitted curve
p1 <- ggplot(data1, aes(x = time, y = em_intensity)) +
    geom_point(color = "blue") +
    geom_line(aes(y = fitted), color = "red") +
    labs(title = "Fluorescence Intensity vs. Time", x = "Time (s)", y = "Fluorescence Intensity") +
    theme_bw()

print(p1)
