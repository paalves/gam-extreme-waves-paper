# =============================================================================
# RETURN LEVEL TABLE FOR ULS COMBINATIONS
# =============================================================================
# This script calculates return values with 95% confidence intervals for:
# - 50 yr, waves (hs)
# - 5 yr and 50 yr, current speed (Cspd)
# - 50 yr, water level (wlv)
# - 5 yr and 50 yr, wind speed (Wspd)
#
# For wind speed and waves: declustering with threshold = 96th quantile, r = 72
# For all variables: GPD for return level calculations
# =============================================================================

library(data.table)
library(POT)
library(evd)
library(ggplot2)
library(scales)
library(tidyr)

# Load data
data <- fread("C:/these_docs/data2024.csv")

# Define declustering function (runs-based)
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

# Function to calculate return level with 95% CI using delta method
calculate_return_level <- function(fit, threshold, return_period, n_years, n_obs) {
  # Event rate (events per year)
  lambda <- n_obs / n_years

  # Number of events for return period
  m <- lambda * return_period

  xi <- fit$param["shape"]
  sigma <- fit$param["scale"]
  vcov <- fit$var.cov

  # Return level formula
  if (abs(xi) < 1e-6) {
    z_m <- threshold + sigma * log(m)
    grad <- c(log(m), 0)
  } else {
    z_m <- threshold + (sigma / xi) * (m^xi - 1)
    grad <- c(
      (m^xi - 1) / xi,
      -sigma * (m^xi - 1) / xi^2 + sigma * m^xi * log(m) / xi
    )
  }

  # Standard error using delta method
  se_z <- sqrt(as.numeric(t(grad) %*% vcov %*% grad))

  return(list(
    return_level = z_m,
    lower = z_m - 1.96 * se_z,
    upper = z_m + 1.96 * se_z
  ))
}

# =============================================================================
# 1. WAVES (hs) - 50 yr return level with declustering
# =============================================================================
hs <- na.omit(data$hs)
threshold_hs <- quantile(hs, 0.96)
run_length <- 72

# Decluster
clusters_hs <- decluster_runs(hs, threshold_hs, run_length)
cluster_maxima_hs <- hs[clusters_hs$cluster_indices]

# Fit GPD to declustered data
fit_hs <- fpot(cluster_maxima_hs, threshold = threshold_hs, model = "gpd")

# Number of observations per year (hourly data)
n_years_hs <- length(hs) / (365.25 * 24)

# Calculate 5-year and 50-year return levels
rl_hs_5 <- calculate_return_level(fit_hs, threshold_hs, 5, n_years_hs, clusters_hs$n_clusters)
rl_hs_50 <- calculate_return_level(fit_hs, threshold_hs, 50, n_years_hs, clusters_hs$n_clusters)

# =============================================================================
# 2. WIND SPEED (Wspd) - 5 yr and 50 yr return levels with declustering
# =============================================================================
Wspd <- na.omit(data$Wspd)
threshold_wspd <- quantile(Wspd, 0.96)

# Decluster
clusters_wspd <- decluster_runs(Wspd, threshold_wspd, run_length)
cluster_maxima_wspd <- Wspd[clusters_wspd$cluster_indices]

# Fit GPD to declustered data
fit_wspd <- fpot(cluster_maxima_wspd, threshold = threshold_wspd, model = "gpd")

# Number of observations per year
n_years_wspd <- length(Wspd) / (365.25 * 24)

# Calculate return levels
rl_wspd_5 <- calculate_return_level(fit_wspd, threshold_wspd, 5, n_years_wspd, clusters_wspd$n_clusters)
rl_wspd_50 <- calculate_return_level(fit_wspd, threshold_wspd, 50, n_years_wspd, clusters_wspd$n_clusters)

# =============================================================================
# 3. CURRENT SPEED (Cspd) - 5 yr and 50 yr return levels
# =============================================================================
Cspd <- na.omit(data$Cspd)
# Use 96th quantile threshold (user specified for waves/wind, applying consistently)
threshold_cspd <- quantile(Cspd, 0.96)

# For current speed, use all exceedances (no declustering specified)
excesses_cspd <- Cspd[Cspd > threshold_cspd]

# Fit GPD
fit_cspd <- fpot(excesses_cspd, threshold = threshold_cspd, model = "gpd")

# Number of observations per year
n_years_cspd <- length(Cspd) / (365.25 * 24)

# Calculate return levels
n_exc_cspd <- length(excesses_cspd)
rl_cspd_5 <- calculate_return_level(fit_cspd, threshold_cspd, 5, n_years_cspd, n_exc_cspd)
rl_cspd_50 <- calculate_return_level(fit_cspd, threshold_cspd, 50, n_years_cspd, n_exc_cspd)

# =============================================================================
# 4. WATER LEVEL (wlv) - 50 yr return level
# =============================================================================
wlv <- na.omit(data$wlv)
# Use 96th quantile threshold
threshold_wlv <- quantile(wlv, 0.96)

# For water level, use all exceedances (no declustering specified)
excesses_wlv <- wlv[wlv > threshold_wlv]

# Fit GPD
fit_wlv <- fpot(excesses_wlv, threshold = threshold_wlv, model = "gpd")

# Number of observations per year
n_years_wlv <- length(wlv) / (365.25 * 24)

# Calculate 50-year return level
n_exc_wlv <- length(excesses_wlv)
rl_wlv_50 <- calculate_return_level(fit_wlv, threshold_wlv, 50, n_years_wlv, n_exc_wlv)

# =============================================================================
# CREATE SUMMARY TABLE
# =============================================================================
results_table <- data.frame(
  Variable = c("Waves (hs)", "Wind Speed (Wspd)", "Wind Speed (Wspd)",
               "Current Speed (Cspd)", "Current Speed (Cspd)", "Water Level (wlv)"),
  Return_Period_yr = c("50", "5", "50", "5", "50", "50"),
  Return_Value = c(rl_hs_50$return_level, rl_wspd_5$return_level, rl_wspd_50$return_level,
                   rl_cspd_5$return_level, rl_cspd_50$return_level, rl_wlv_50$return_level),
  CI_95_Lower = c(rl_hs_50$lower, rl_wspd_5$lower, rl_wspd_50$lower,
                  rl_cspd_5$lower, rl_cspd_50$lower, rl_wlv_50$lower),
  CI_95_Upper = c(rl_hs_50$upper, rl_wspd_5$upper, rl_wspd_50$upper,
                  rl_cspd_5$upper, rl_cspd_50$upper, rl_wlv_50$upper),
  Threshold = c(threshold_hs, threshold_wspd, threshold_wspd,
                threshold_cspd, threshold_cspd, threshold_wlv),
  Declustered = c("Yes (r=72)", "Yes (r=72)", "Yes (r=72)", "No", "No", "No")
)

# Print table
cat("\n")
cat("================================================================================\n")
cat("RETURN LEVEL TABLE - 95% CONFIDENCE INTERVALS\n")
cat("================================================================================\n")
cat("Method: GPD fitted to declustered data (waves, wind) or all exceedances (current, water level)\n")
cat("Threshold: 96th quantile for all variables\n")
cat("Declustered with run length r = 72 hours for waves and wind speed\n")
cat("================================================================================\n\n")

print(results_table, row.names = FALSE)

cat("\n================================================================================\n")
cat("NOTES:\n")
cat("================================================================================\n")
cat("- Waves (hs) and Wind Speed (Wspd) were declustered using runs method (r=72 hours)\n")
cat("- Current Speed (Cspd) and Water Level (wlv) used all exceedances (no declustering)\n")
cat("- 95% CI calculated using delta method\n")
cat("- All thresholds set at 96th percentile\n")
cat("================================================================================\n\n")

# =============================================================================
# ULS COMBINATIONS PLOT
# =============================================================================
# Define the 4 ULS cases:
# Case 1: waves 50 yr, current 5yr, water level 50 yr, wind 5 yr
# Case 2: waves 5 yr, current 50 yr, water level 50 yr, wind 5 yr
# Case 3: waves 50 yr, current 5 yr, water level 50 yr, wind 50 yr
# Case 4: current 50 yr, water level 50 yr, wind 5 yr

# Create ULS combinations data frame
uls_cases <- data.frame(
  Case = factor(c(
    "Environmental combination 1", "Environmental combination 1", "Environmental combination 1", "Environmental combination 1",
    "Environmental combination 2", "Environmental combination 2", "Environmental combination 2", "Environmental combination 2",
    "Environmental combination 4", "Environmental combination 4", "Environmental combination 4", "Environmental combination 4",
    "Environmental combination 5", "Environmental combination 5", "Environmental combination 5"
  ), levels = c("Environmental combination 1", "Environmental combination 2", "Environmental combination 4", "Environmental combination 5")),
  Variable = factor(c(
    "Waves", "Current", "Water Level", "Wind Speed",
    "Waves", "Current", "Water Level", "Wind Speed",
    "Waves", "Current", "Water Level", "Wind Speed",
    "Current", "Water Level", "Wind Speed"
  ), levels = c("Waves", "Current", "Water Level", "Wind Speed")),
  Return_Period = c("50 yr", "5 yr", "50 yr", "5 yr",
                    "5 yr", "50 yr", "50 yr", "5 yr",
                    "50 yr", "5 yr", "50 yr", "50 yr",
                    "50 yr", "50 yr", "5 yr"),
  Return_Level = c(
    rl_hs_50$return_level, rl_cspd_5$return_level, rl_wlv_50$return_level, rl_wspd_5$return_level,
    rl_hs_5$return_level, rl_cspd_50$return_level, rl_wlv_50$return_level, rl_wspd_5$return_level,
    rl_hs_50$return_level, rl_cspd_5$return_level, rl_wlv_50$return_level, rl_wspd_50$return_level,
    rl_cspd_50$return_level, rl_wlv_50$return_level, rl_wspd_5$return_level
  ),
  CI_Lower = c(
    rl_hs_50$lower, rl_cspd_5$lower, rl_wlv_50$lower, rl_wspd_5$lower,
    rl_hs_5$lower, rl_cspd_50$lower, rl_wlv_50$lower, rl_wspd_5$lower,
    rl_hs_50$lower, rl_cspd_5$lower, rl_wlv_50$lower, rl_wspd_50$lower,
    rl_cspd_50$lower, rl_wlv_50$lower, rl_wspd_5$lower
  ),
  CI_Upper = c(
    rl_hs_50$upper, rl_cspd_5$upper, rl_wlv_50$upper, rl_wspd_5$upper,
    rl_hs_5$upper, rl_cspd_50$upper, rl_wlv_50$upper, rl_wspd_5$upper,
    rl_hs_50$upper, rl_cspd_5$upper, rl_wlv_50$upper, rl_wspd_50$upper,
    rl_cspd_50$upper, rl_wlv_50$upper, rl_wspd_5$upper
  )
)

# Color palette for variables
variable_colors <- c(
  "Waves" = "#1f77b4",       # blue
  "Current" = "#ff7f0e",     # orange
  "Water Level" = "#2ca02c", # green
  "Wind Speed" = "#d62728"   # red
)

# Create the plot
p_uls <- ggplot(uls_cases, aes(x = Variable, y = Return_Level, fill = Variable)) +
  geom_col(position = position_dodge(width = 0.9), alpha = 0.7) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper),
                width = 0.25, linewidth = 0.8) +
  geom_text(aes(label = paste0(round(Return_Level, 2), "\n(", round(CI_Lower, 2), "-", round(CI_Upper, 2), ")")),
            size = 2.5, vjust = -0.5) +
  facet_wrap(~ Case, ncol = 2) +
  scale_fill_manual(values = variable_colors) +
  labs(
    title = "ULS combinations of uncorrelated extreme events",
    subtitle = "GPD with 96th quantile threshold. 72-hour runs declustering applied for waves and wind.\nReturn levels with associated 95% confidence intervals.",
    x = "Variable",
    y = "Return level"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold")
  ) +
  coord_cartesian(ylim = c(0, max(uls_cases$CI_Upper) * 1.3))

# Print the plot
print(p_uls)

# =============================================================================
# CASE DESCRIPTIONS
# =============================================================================
case_labels <- data.frame(
  Case = c("Environmental combination 1", "Environmental combination 2", "Environmental combination 4", "Environmental combination 5"),
  Description = c(
    "Waves 50yr, Current 5yr, WL 50yr, Wind 5yr",
    "Waves 5yr, Current 50yr, WL 50yr, Wind 5yr",
    "Waves 50yr, Current 5yr, WL 50yr, Wind 50yr",
    "Current 50yr, WL 50yr, Wind 5yr (no Waves)"
  )
)

cat("\n")
cat("================================================================================\n")
cat("ULS CASES DESCRIPTION\n")
cat("================================================================================\n")
print(case_labels, row.names = FALSE)
cat("================================================================================\n\n")

# Return the results invisibly
invisible(list(
  results_table = results_table,
  uls_cases = uls_cases,
  plot = p_uls
))

