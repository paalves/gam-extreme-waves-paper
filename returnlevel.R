library(data.table)
library(POT)
library(evd)
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)

# Load data
data <- fread("C:/these_docs/data2024.csv")

# Remove missing values
hs <- na.omit(data$hs)

# ============================================================================
# PARAMETERS
# ============================================================================
# Chosen threshold (96th percentile)
chosen_threshold <- quantile(hs, 0.96)

# Run length from declustering analysis (hours)
run_length <- 72

# Return periods to estimate (years)
return_periods <- c(1, 2, 5, 10, 20, 50, 100, 200)

cat("========================================\n")
cat("GPD RETURN VALUE ANALYSIS\n")
cat("========================================\n")
cat(sprintf("Threshold: %.2f m (Q96)\n", chosen_threshold))
cat(sprintf("Run length: %d hours\n", run_length))
cat(sprintf("Total observations: %d\n", length(hs)))
cat(sprintf("Exceedances: %d (%.2f%%)\n\n", 
            sum(hs > chosen_threshold), 
            100 * mean(hs > chosen_threshold)))

# ============================================================================
# MANUAL DECLUSTERING
# ============================================================================
decluster_runs <- function(x, threshold, run_length) {
  exceed <- x > threshold
  cluster_max <- rep(FALSE, length(x))
  
  if (sum(exceed) == 0) return(list(n_clusters = 0, cluster_indices = integer(0)))
  
  exceed_indices <- which(exceed)
  if (length(exceed_indices) == 0) {
    return(list(n_clusters = 0, cluster_indices = integer(0)))
  }
  
  gaps <- diff(exceed_indices)
  new_cluster <- c(TRUE, gaps > run_length)
  cluster_id <- cumsum(new_cluster)
  cluster_indices <- integer(0)
  
  for (i in unique(cluster_id)) {
    cluster_members <- exceed_indices[cluster_id == i]
    max_idx <- cluster_members[which.max(x[cluster_members])]
    cluster_indices <- c(cluster_indices, max_idx)
  }
  
  return(list(
    n_clusters = length(cluster_indices),
    cluster_indices = cluster_indices
  ))
}

# Apply declustering
clusters <- decluster_runs(hs, chosen_threshold, run_length)
cluster_maxima <- hs[clusters$cluster_indices]

cat(sprintf("Independent events after declustering: %d\n", clusters$n_clusters))
cat(sprintf("Extremal index: %.4f\n\n", clusters$n_clusters / sum(hs > chosen_threshold)))

# ============================================================================
# FIT GPD MODEL - METHOD 1: DECLUSTERED DATA
# ============================================================================
fit_declust <- fpot(cluster_maxima, threshold = chosen_threshold, model = "gpd")

cat("========================================\n")
cat("METHOD 1: DECLUSTERED DATA\n")
cat("========================================\n")
cat("Fit GPD to cluster maxima only\n\n")
cat(sprintf("Shape (ξ):     %.4f (SE: %.4f)\n", 
            fit_declust$param["shape"], fit_declust$std.err["shape"]))
cat(sprintf("Scale (σ):     %.4f (SE: %.4f)\n", 
            fit_declust$param["scale"], fit_declust$std.err["scale"]))
cat(sprintf("95%% CI Shape: [%.4f, %.4f]\n",
            fit_declust$param["shape"] - 1.96 * fit_declust$std.err["shape"],
            fit_declust$param["shape"] + 1.96 * fit_declust$std.err["shape"]))
cat(sprintf("95%% CI Scale: [%.4f, %.4f]\n\n",
            fit_declust$param["scale"] - 1.96 * fit_declust$std.err["scale"],
            fit_declust$param["scale"] + 1.96 * fit_declust$std.err["scale"]))

# ============================================================================
# FIT GPD MODEL - METHOD 2: ALL EXCEEDANCES WITH EXTREMAL INDEX
# ============================================================================
all_exceedances <- hs[hs > chosen_threshold]
fit_all <- fpot(all_exceedances, threshold = chosen_threshold, model = "gpd")

# Calculate extremal index
extremal_index <- clusters$n_clusters / length(all_exceedances)

cat("========================================\n")
cat("METHOD 2: ALL EXCEEDANCES + EXTREMAL INDEX\n")
cat("========================================\n")
cat("Fit GPD to all exceedances, correct with extremal index\n\n")
cat(sprintf("Shape (ξ):     %.4f (SE: %.4f)\n", 
            fit_all$param["shape"], fit_all$std.err["shape"]))
cat(sprintf("Scale (σ):     %.4f (SE: %.4f)\n", 
            fit_all$param["scale"], fit_all$std.err["scale"]))
cat(sprintf("95%% CI Shape: [%.4f, %.4f]\n",
            fit_all$param["shape"] - 1.96 * fit_all$std.err["shape"],
            fit_all$param["shape"] + 1.96 * fit_all$std.err["shape"]))
cat(sprintf("95%% CI Scale: [%.4f, %.4f]\n",
            fit_all$param["scale"] - 1.96 * fit_all$std.err["scale"],
            fit_all$param["scale"] + 1.96 * fit_all$std.err["scale"]))
cat(sprintf("Extremal index (θ): %.4f\n\n", extremal_index))

cat("========================================\n")
cat("METHOD 3: BLOCK MAXIMA (GEV)\n")
cat("========================================\n")

# Monthly block maxima (hourly data assumed)
obs_per_year  <- 365.25 * 24
obs_per_month <- obs_per_year / 12
block_size    <- round(obs_per_month)  # ~730 hours per month

# Extract monthly block maxima
n_blocks <- floor(length(hs) / block_size)
block_maxima <- numeric(n_blocks)
for (i in 1:n_blocks) {
  start_idx <- (i - 1) * block_size + 1
  end_idx   <- min(i * block_size, length(hs))
  block_maxima[i] <- max(hs[start_idx:end_idx], na.rm = TRUE)
}

# Convenience
blocks_per_year <- obs_per_year / block_size  # ~12

cat(sprintf("Number of blocks: %d\n", n_blocks))
cat(sprintf("Block size: %d hours (%.2f months)\n\n", block_size, block_size / obs_per_month))

# Fit GEV distribution to monthly maxima
fit_gev <- fgev(block_maxima, std.err = TRUE)

cat(sprintf("Location (μ):  %.4f (SE: %.4f)\n",  fit_gev$estimate["loc"],   fit_gev$std.err["loc"]))
cat(sprintf("Scale (σ):     %.4f (SE: %.4f)\n",  fit_gev$estimate["scale"], fit_gev$std.err["scale"]))
cat(sprintf("Shape (ξ):     %.4f (SE: %.4f)\n",  fit_gev$estimate["shape"], fit_gev$std.err["shape"]))
cat(sprintf("95%% CI Shape: [%.4f, %.4f]\n\n",
            fit_gev$estimate["shape"] - 1.96 * fit_gev$std.err["shape"],
            fit_gev$estimate["shape"] + 1.96 * fit_gev$std.err["shape"]))

# ============================================================================
# RETURN VALUE ESTIMATION - ALL THREE METHODS
# ============================================================================
# Number of observations per year (hourly data)
obs_per_year <- 365.25 * 24

# METHOD 1: Rate based on declustered events
lambda_declust <- clusters$n_clusters / (length(hs) / obs_per_year)

# METHOD 2: Rate based on all exceedances times extremal index
lambda_all <- (length(all_exceedances) / (length(hs) / obs_per_year)) * extremal_index

cat("========================================\n")
cat("RETURN VALUE ESTIMATES COMPARISON\n")
cat("========================================\n")
cat(sprintf("Method 1 (GPD Declust) - Events per year: %.2f\n", lambda_declust))
cat(sprintf("Method 2 (GPD + Ext Idx) - Events per year: %.2f\n", lambda_all))
cat(sprintf("Method 3 (Block Maxima) - Blocks per year: %.2f\n\n", blocks_per_year))

# Calculate return values for all three methods
return_values_declust <- data.frame(
  period = return_periods,
  level = numeric(length(return_periods)),
  lower = numeric(length(return_periods)),
  upper = numeric(length(return_periods)),
  method = "GPD Declustered"
)

return_values_all <- data.frame(
  period = return_periods,
  level = numeric(length(return_periods)),
  lower = numeric(length(return_periods)),
  upper = numeric(length(return_periods)),
  method = "GPD Extremal Index"
)

return_periods_gev <- return_periods[return_periods >= 2]

return_values_gev <- data.frame(
  period = return_periods_gev,
  level  = numeric(length(return_periods_gev)),
  lower  = numeric(length(return_periods_gev)),
  upper  = numeric(length(return_periods_gev)),
  method = "Block Maxima (GEV)"
)

# METHOD 1: Declustered GPD
for (i in seq_along(return_periods)) {
  T_years <- return_periods[i]
  m <- lambda_declust * T_years
  
  xi <- fit_declust$param["shape"]
  sigma <- fit_declust$param["scale"]
  
  if (abs(xi) < 1e-6) {
    z_m <- chosen_threshold + sigma * log(m)
  } else {
    z_m <- chosen_threshold + (sigma / xi) * (m^xi - 1)
  }
  
  vcov <- fit_declust$var.cov
  
  if (abs(xi) < 1e-6) {
    grad <- c(log(m), 0)
  } else {
    grad <- c(
      (m^xi - 1) / xi,
      -sigma * (m^xi - 1) / xi^2 + sigma * m^xi * log(m) / xi
    )
  }
  
  se_z <- sqrt(t(grad) %*% vcov %*% grad)
  
  return_values_declust$level[i] <- z_m
  return_values_declust$lower[i] <- z_m - 1.96 * se_z
  return_values_declust$upper[i] <- z_m + 1.96 * se_z
}

# METHOD 2: All exceedances with extremal index
for (i in seq_along(return_periods)) {
  T_years <- return_periods[i]
  m <- lambda_all * T_years
  
  xi <- fit_all$param["shape"]
  sigma <- fit_all$param["scale"]
  
  if (abs(xi) < 1e-6) {
    z_m <- chosen_threshold + sigma * log(m)
  } else {
    z_m <- chosen_threshold + (sigma / xi) * (m^xi - 1)
  }
  
  vcov <- fit_all$var.cov
  
  if (abs(xi) < 1e-6) {
    grad <- c(log(m), 0)
  } else {
    grad <- c(
      (m^xi - 1) / xi,
      -sigma * (m^xi - 1) / xi^2 + sigma * m^xi * log(m) / xi
    )
  }
  
  se_z <- sqrt(t(grad) %*% vcov %*% grad)
  
  return_values_all$level[i] <- z_m
  return_values_all$lower[i] <- z_m - 1.96 * se_z
  return_values_all$upper[i] <- z_m + 1.96 * se_z
}

# METHOD 3: Block Maxima (GEV)
# METHOD 3: Block Maxima (GEV) with monthly blocks
for (i in seq_along(return_periods_gev)) {
  T_years <- return_periods_gev[i]

  mu        <- fit_gev$estimate["loc"]
  sigma_gev <- fit_gev$estimate["scale"]
  xi_gev    <- fit_gev$estimate["shape"]

  # Exceedance probability per MONTHLY block for a T-year return
  p <- 1 / (T_years * blocks_per_year)

  if (abs(xi_gev) < 1e-6) {
    z_p <- mu - sigma_gev * log(-log(1 - p))
  } else {
    z_p <- mu + (sigma_gev / xi_gev) * ((-log(1 - p))^(-xi_gev) - 1)
  }

  # Delta method SE
  if (abs(xi_gev) < 1e-6) {
    grad_gev <- c(1, -log(-log(1 - p)), 0)
  } else {
    y_p <- -log(1 - p)
    grad_gev <- c(
      1,
      (y_p^(-xi_gev) - 1) / xi_gev,
      -sigma_gev * (y_p^(-xi_gev) - 1) / xi_gev^2 + sigma_gev * y_p^(-xi_gev) * log(y_p) / xi_gev
    )
  }

  vcov_gev  <- fit_gev$var.cov
  se_z_gev  <- sqrt(t(grad_gev) %*% vcov_gev %*% grad_gev)

  return_values_gev$level[i] <- z_p
  return_values_gev$lower[i] <- z_p - 1.96 * se_z_gev
  return_values_gev$upper[i] <- z_p + 1.96 * se_z_gev
}

# Combine all methods
return_values <- rbind(return_values_declust, return_values_all, return_values_gev)

# Print comparison tables
cat("METHOD 1: GPD DECLUSTERED DATA\n")
cat("T (years) | Return Level (m) | 95% CI Lower | 95% CI Upper\n")
cat("----------|------------------|--------------|-------------\n")
for (i in 1:nrow(return_values_declust)) {
  cat(sprintf("%9.0f | %16.2f | %12.2f | %12.2f\n",
              return_values_declust$period[i],
              return_values_declust$level[i],
              return_values_declust$lower[i],
              return_values_declust$upper[i]))
}

cat("\nMETHOD 2: GPD ALL EXCEEDANCES + EXTREMAL INDEX\n")
cat("T (years) | Return Level (m) | 95% CI Lower | 95% CI Upper\n")
cat("----------|------------------|--------------|-------------\n")
for (i in 1:nrow(return_values_all)) {
  cat(sprintf("%9.0f | %16.2f | %12.2f | %12.2f\n",
              return_values_all$period[i],
              return_values_all$level[i],
              return_values_all$lower[i],
              return_values_all$upper[i]))
}

cat("\nMETHOD 3: BLOCK MAXIMA (GEV)\n")
cat("T (years) | Return Level (m) | 95% CI Lower | 95% CI Upper\n")
cat("----------|------------------|--------------|-------------\n")
for (i in 1:nrow(return_values_gev)) {
  cat(sprintf("%9.0f | %16.2f | %12.2f | %12.2f\n",
              return_values_gev$period[i],
              return_values_gev$level[i],
              return_values_gev$lower[i],
              return_values_gev$upper[i]))
}

# ============================================================================
# DIAGNOSTIC PLOTS
# ============================================================================

# Prepare excesses for both GPD fits
excesses_declust <- cluster_maxima - chosen_threshold
excesses_all <- all_exceedances - chosen_threshold

# 1. QQ Plot - Method 1 (Declustered)
n1 <- length(excesses_declust)
empirical_probs1 <- (1:n1 - 0.5) / n1
sorted_excesses1 <- sort(excesses_declust)

xi1 <- fit_declust$param["shape"]
sigma1 <- fit_declust$param["scale"]

if (abs(xi1) < 1e-6) {
  theoretical_quantiles1 <- -sigma1 * log(1 - empirical_probs1)
} else {
  theoretical_quantiles1 <- (sigma1 / xi1) * ((1 - empirical_probs1)^(-xi1) - 1)
}

qq_data1 <- data.frame(
  empirical = sorted_excesses1,
  theoretical = theoretical_quantiles1
)

p1 <- ggplot(qq_data1, aes(x = theoretical, y = empirical)) +
  geom_point(size = 2, alpha = 0.6, color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "QQ Plot - GPD Declustered",
    subtitle = sprintf("ξ=%.3f, σ=%.3f", xi1, sigma1),
    x = "Model Quantiles",
    y = "Empirical Quantiles"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.minor = element_blank()
  )

# 2. QQ Plot - Method 2 (All exceedances)
n2 <- length(excesses_all)
empirical_probs2 <- (1:n2 - 0.5) / n2
sorted_excesses2 <- sort(excesses_all)

xi2 <- fit_all$param["shape"]
sigma2 <- fit_all$param["scale"]

if (abs(xi2) < 1e-6) {
  theoretical_quantiles2 <- -sigma2 * log(1 - empirical_probs2)
} else {
  theoretical_quantiles2 <- (sigma2 / xi2) * ((1 - empirical_probs2)^(-xi2) - 1)
}

qq_data2 <- data.frame(
  empirical = sorted_excesses2,
  theoretical = theoretical_quantiles2
)

p2 <- ggplot(qq_data2, aes(x = theoretical, y = empirical)) +
  geom_point(size = 2, alpha = 0.6, color = "darkgreen") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "QQ Plot - GPD All Exceedances",
    subtitle = sprintf("ξ=%.3f, σ=%.3f", xi2, sigma2),
    x = "Model Quantiles",
    y = "Empirical Quantiles"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.minor = element_blank()
  )

# 3. Return Level Plot - All Methods
p3 <- ggplot(return_values, aes(x = period, y = level, color = method, fill = method)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = chosen_threshold, linetype = "dashed", 
             color = "black", alpha = 0.6) +
  annotate("text", x = max(return_periods) * 0.5, y = chosen_threshold, 
           label = sprintf("Threshold = %.2f m", chosen_threshold),
           vjust = -0.5, size = 3, color = "black") +
  scale_x_log10(breaks = return_periods) +
  scale_color_manual(values = c("GPD Declustered" = "blue", 
                                 "GPD Extremal Index" = "darkgreen",
                                 "Block Maxima (GEV)" = "red"),
                     name = "Method") +
  scale_fill_manual(values = c("GPD Declustered" = "blue", 
                                "GPD Extremal Index" = "darkgreen",
                                "Block Maxima (GEV)" = "red"),
                    name = "Method") +
  labs(
    title = "Return Level Plot - Three Methods Comparison",
    subtitle = "GPD (declustered), GPD (extremal index), and Block Maxima (GEV)",
    x = "Return Period (years)",
    y = "Return Level Hs (m)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# 4. Density Plot - Method 1 (Declustered)
hist_data1 <- data.frame(excesses = excesses_declust)
x_seq1 <- seq(0, max(excesses_declust), length.out = 200)

if (abs(xi1) < 1e-6) {
  fitted_density1 <- (1 / sigma1) * exp(-x_seq1 / sigma1)
} else {
  fitted_density1 <- (1 / sigma1) * (1 + xi1 * x_seq1 / sigma1)^(-(1/xi1 + 1))
}

density_data1 <- data.frame(x = x_seq1, density = fitted_density1)

p4 <- ggplot() +
  geom_histogram(data = hist_data1, aes(x = excesses, y = after_stat(density)),
                 bins = 30, fill = "lightblue", color = "black", alpha = 0.6) +
  geom_line(data = density_data1, aes(x = x, y = density),
            color = "blue", linewidth = 1.2) +
  labs(
    title = "Density Plot - GPD Declustered",
    subtitle = sprintf("ξ=%.3f, σ=%.3f", xi1, sigma1),
    x = "Excesses over Threshold (m)",
    y = "Density"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.minor = element_blank()
  )

# 5. Density Plot - Method 2 (All exceedances)
hist_data2 <- data.frame(excesses = excesses_all)
x_seq2 <- seq(0, max(excesses_all), length.out = 200)

if (abs(xi2) < 1e-6) {
  fitted_density2 <- (1 / sigma2) * exp(-x_seq2 / sigma2)
} else {
  fitted_density2 <- (1 / sigma2) * (1 + xi2 * x_seq2 / sigma2)^(-(1/xi2 + 1))
}

density_data2 <- data.frame(x = x_seq2, density = fitted_density2)

p5 <- ggplot() +
  geom_histogram(data = hist_data2, aes(x = excesses, y = after_stat(density)),
                 bins = 50, fill = "lightgreen", color = "black", alpha = 0.6) +
  geom_line(data = density_data2, aes(x = x, y = density),
            color = "darkgreen", linewidth = 1.2) +
  labs(
    title = "Density Plot - GPD All Exceedances",
    subtitle = sprintf("ξ=%.3f, σ=%.3f", xi2, sigma2),
    x = "Excesses over Threshold (m)",
    y = "Density"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.minor = element_blank()
  )

# 6. QQ Plot for GEV (Block Maxima)
n_gev <- length(block_maxima)
empirical_probs_gev <- (1:n_gev - 0.5) / n_gev
sorted_block_maxima <- sort(block_maxima)

mu_gev <- fit_gev$estimate["loc"]
sigma_gev <- fit_gev$estimate["scale"]
xi_gev <- fit_gev$estimate["shape"]

if (abs(xi_gev) < 1e-6) {
  theoretical_quantiles_gev <- mu_gev - sigma_gev * log(-log(empirical_probs_gev))
} else {
  theoretical_quantiles_gev <- mu_gev + (sigma_gev / xi_gev) * 
    ((-log(empirical_probs_gev))^(-xi_gev) - 1)
}

qq_data_gev <- data.frame(
  empirical = sorted_block_maxima,
  theoretical = theoretical_quantiles_gev
)

p6 <- ggplot(qq_data_gev, aes(x = theoretical, y = empirical)) +
  geom_point(size = 2, alpha = 0.6, color = "red") +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed", linewidth = 1) +
  labs(
    title = "QQ Plot - GEV Block Maxima",
    subtitle = sprintf("μ=%.3f, σ=%.3f, ξ=%.3f", mu_gev, sigma_gev, xi_gev),
    x = "Model Quantiles",
    y = "Empirical Quantiles"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.minor = element_blank()
  )

# 7. Density Plot - GEV (Block Maxima)
hist_data_gev <- data.frame(vals = block_maxima)
x_seq_gev <- seq(min(block_maxima), max(block_maxima), length.out = 200)

mu_gev2        <- fit_gev$estimate["loc"]
sigma_gev2     <- fit_gev$estimate["scale"]
xi_gev2        <- fit_gev$estimate["shape"]
fitted_density_gev <- evd::dgev(x_seq_gev, loc = mu_gev2, scale = sigma_gev2, shape = xi_gev2)

density_data_gev <- data.frame(x = x_seq_gev, density = fitted_density_gev)

p7 <- ggplot() +
  geom_histogram(data = hist_data_gev, aes(x = vals, y = after_stat(density)),
                 bins = 30, fill = "mistyrose", color = "black", alpha = 0.6) +
  geom_line(data = density_data_gev, aes(x = x, y = density),
            color = "red", linewidth = 1.2) +
  labs(
    title = "Density Plot - GEV (Monthly Maxima)",
    subtitle = sprintf("μ=%.3f, σ=%.3f, ξ=%.3f", mu_gev2, sigma_gev2, xi_gev2),
    x = "Monthly Block Maxima (m)",
    y = "Density"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.minor = element_blank()
  )

library(patchwork)
(p1+p2+p6) / (p4+p5+p7)

# ============================================================================
# COMBINE DIAGNOSTIC PLOTS
# ============================================================================
diagnostics_plot <- grid.arrange(
  p1, p2, p6,   # QQ plots: GPD-declust, GPD-all, GEV
  p4, p5, p7,   # Densities: GPD-declust, GPD-all, GEVs
  nrow = 2, ncol = 3,
  top = textGrob(
    "GPD and GEV Diagnostics (Monthly Block Maxima for GEV)",
    gp = gpar(fontsize = 14, fontface = "bold")
  )
)

p3 <- ggplot(return_values, aes(x = period, y = level, color = method, fill = method)) +
  #geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = chosen_threshold, linetype = "dashed", color = "black", alpha = 0.6) +
  annotate("text", x = max(return_values$period) * 0.5, y = chosen_threshold,
           label = sprintf("Threshold = %.2f m", chosen_threshold),
           vjust = -0.5, size = 3, color = "black") +
  scale_x_log10(breaks = sort(unique(return_values$period))) +
  scale_color_manual(values = c("GPD Declustered" = "blue",
                                "GPD Extremal Index" = "darkgreen",
                                "Block Maxima (GEV)" = "red"),
                     name = "Method") +
  scale_fill_manual(values = c("GPD Declustered" = "blue",
                               "GPD Extremal Index" = "darkgreen",
                               "Block Maxima (GEV)" = "red"),
                    name = "Method") +
  labs(
    title = "Return Level Plot - GPD vs. GEV",
    subtitle = "GEV return periods start at 2 years",
    x = "Return Period (years)",
    y = "Return Level Hs (m)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )


diagnostics_plot
p3
