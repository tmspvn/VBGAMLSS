#
#
#
#
#
#
#
#
#
#


# Define the plotting function
plot_gamlss_shash <- function(mu = 0, sigma = 1, nu = 1, tau = 1, x_range = c(-5, 5)) {
  # Create a sequence of x values
  x_vals <- seq(x_range[1], x_range[2], length.out = 500)

  # Create a grid of all parameter combinations passed to the function
  params <- expand.grid(mu = mu, sigma = sigma, nu = nu, tau = tau)

  # Calculate the density for each combination of parameters using gamlss.dist
  df_list <- lapply(1:nrow(params), function(i) {
    m <- params$mu[i]
    s <- params$sigma[i]
    n <- params$nu[i]
    t <- params$tau[i]

    # Using dSHASH from the gamlss.dist package
    y_vals <- dSHASH(x = x_vals, mu = m, sigma = s, nu = n, tau = t)

    # Create a descriptive label for the legend
    label <- sprintf("μ=%g, σ=%g, ν=%g, τ=%g", m, s, n, t)

    data.frame(x = x_vals, y = y_vals, label = as.factor(label))
  })

  # Combine the lists into a single data frame for ggplot
  df <- do.call(rbind, df_list)

  # Generate the overlaid plots
  ggplot(df, aes(x = x, y = y, color = label)) +
    geom_line(linewidth = 1) +
    labs(
      title = "SHASH Distribution",
      subtitle = "Comparing density curves across varying parameters",
      x = "x",
      y = "Density",
      color = "Parameters"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom", legend.direction = "vertical")
}

# --- Examples of how to use the function ---

# 1. Varying the Scale (sigma)
# plot_gamlss_shash(sigma = c(0.5, 1, 2))

# 2. Varying the Skewness (nu)
# plot_gamlss_shash(nu = c(-1, 0, 1, 2))

# 3. Varying the Tailweight/Kurtosis (tau)
# plot_gamlss_shash(tau = c(0.5, 0.8, 1, 1.5), x_range = c(-8, 8))

# 4. Combining multiple variations at once
# plot_gamlss_shash(nu = c(-1, 1), tau = c(0.5, 1.5), x_range = c(-6, 6))
