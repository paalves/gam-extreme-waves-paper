library(data.table)
library(POT)
library(evd)
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)

# ============================================================================
# SETUP
# ============================================================================
#data <- fread("C:/these_docs/data2024.csv")
data <- fread("C:/Users/minim/Documents/teletravail/data2024.csv")

hs   <- na.omit(data$hs)

chosen_threshold <- quantile(hs, 0.96)
thresholds       <- quantile(hs, probs = seq(0.90, 0.995, by = 0.005))

# Shared theme for all plots
theme_diag <- theme_bw(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 12, color = "grey30"),
    axis.title    = element_text(size = 13),
    axis.text     = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

# Shared vertical line + label for all plots
vline_threshold <- function(y_pos, label_text) {
  list(
    geom_vline(xintercept = chosen_threshold, linetype = "dashed",
               color = "blue", linewidth = 1),
    annotate("text", x = chosen_threshold, y = y_pos,
             label = label_text,
             hjust = -0.08, vjust = 1, size = 5, color = "blue",
             fontface = "bold")
  )
}

# ============================================================================
# 1. MEAN RESIDUAL LIFE PLOT
# ============================================================================
mrl_data <- do.call(rbind, lapply(thresholds, function(u) {
  excesses <- hs[hs > u] - u
  n  <- length(excesses)
  me <- mean(excesses)
  se <- sd(excesses) / sqrt(n)
  data.frame(threshold = u, mean_excess = me,
             lower = me - 1.96 * se, upper = me + 1.96 * se)
}))

p1 <- ggplot(mrl_data, aes(x = threshold, y = mean_excess)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  vline_threshold(max(mrl_data$upper),
                  sprintf("Chosen threshold\n(Q96 = %.2f m)", chosen_threshold)) +
  labs(
    title    = "Mean residual life plot",
    subtitle = "Approximate linearity indicates suitable threshold range",
    x        = expression(paste("Threshold ", italic(u), " (m)")),
    y        = "Mean excess (m)"
  ) +
  theme_diag

# ============================================================================
# 2. GPD PARAMETER STABILITY
# ============================================================================
param_data <- do.call(rbind, lapply(thresholds, function(u) {
  excesses <- hs[hs > u] - u
  if (length(excesses) <= 50) return(NULL)

  fit <- tryCatch(fpot(hs, threshold = u, model = "gpd"), error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  data.frame(
    threshold     = u,
    scale         = fit$param["scale"],
    shape         = fit$param["shape"],
    scale_lower   = fit$param["scale"] - 1.96 * fit$std.err["scale"],
    scale_upper   = fit$param["scale"] + 1.96 * fit$std.err["scale"],
    shape_lower   = fit$param["shape"] - 1.96 * fit$std.err["shape"],
    shape_upper   = fit$param["shape"] + 1.96 * fit$std.err["shape"],
    n_exceedances = length(excesses)
  )
}))

# Shape parameter
p2 <- ggplot(param_data, aes(x = threshold, y = shape)) +
  geom_ribbon(aes(ymin = shape_lower, ymax = shape_upper), fill = "grey80", alpha = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.6) +
  vline_threshold(max(param_data$shape_upper),
                  sprintf("Q96 = %.2f m", chosen_threshold)) +
  labs(
    title    = "Shape parameter stability",
    subtitle = "Stable values indicate appropriate threshold",
    x        = expression(paste("Threshold ", italic(u), " (m)")),
    y        = expression(paste("Shape parameter ", xi))
  ) +
  theme_diag

# Modified scale parameter (sigma* = sigma - xi * u)
p3 <- ggplot(param_data, aes(x = threshold, y = scale)) +
  geom_ribbon(aes(ymin = scale_lower, ymax = scale_upper), fill = "grey80", alpha = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  vline_threshold(max(param_data$scale_upper),
                  sprintf("Q96 = %.2f m", chosen_threshold)) +
  labs(
    title    = "Modified scale parameter stability",
    subtitle = expression(paste("Stability in ", sigma[u] - xi %.% u, " indicates good threshold")),
    x        = expression(paste("Threshold ", italic(u), " (m)")),
    y        = expression(paste("Modified scale ", sigma^"*"))
  ) +
  theme_diag

# ============================================================================
# 3. NUMBER OF EXCEEDANCES
# ============================================================================
p4 <- ggplot(param_data, aes(x = threshold, y = n_exceedances)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "darkred", alpha = 0.7) +
  annotate("text", x = min(param_data$threshold), y = 115,
           label = "Minimum recommended", hjust = 0, size = 4.5, color = "darkred") +
  vline_threshold(max(param_data$n_exceedances),
                  sprintf("Q96 = %.2f m", chosen_threshold)) +
  scale_y_log10(labels = comma) +
  labs(
    title    = "Number of threshold exceedances",
    subtitle = "Balance between bias (low threshold) and variance (high threshold)",
    x        = expression(paste("Threshold ", italic(u), " (m)")),
    y        = "Number of exceedances"
  ) +
  theme_diag

# ============================================================================
# COMBINE AND SAVE
# ============================================================================
combined_plot <- grid.arrange(
  p1, p2, p3,
  nrow = 1,
  top  = textGrob(
    "GPD threshold selection diagnostics",
    gp = gpar(fontsize = 17, fontface = "bold")
  )
)

# ============================================================================
# CONSOLE SUMMARIES
# ============================================================================
fit_chosen <- fpot(hs, threshold = chosen_threshold, model = "gpd")

cat("\n========================================\n")
cat("CHOSEN THRESHOLD ANALYSIS (Q96)\n")
cat("========================================\n")
cat(sprintf("Chosen threshold:     %.2f m\n",   chosen_threshold))
cat(sprintf("Number of exceedances: %d\n",       sum(hs > chosen_threshold)))
cat(sprintf("Percentage of data:   %.2f%%\n",   100 * mean(hs > chosen_threshold)))
cat(sprintf("Shape parameter:      %.3f (±%.3f)\n",
            fit_chosen$param["shape"], 1.96 * fit_chosen$std.err["shape"]))
cat(sprintf("Scale parameter:      %.3f (±%.3f)\n",
            fit_chosen$param["scale"], 1.96 * fit_chosen$std.err["scale"]))

# Suggest most stable threshold: smallest shape CI width with >= 100 exceedances
param_data$shape_ci_width <- param_data$shape_upper - param_data$shape_lower
stable <- subset(param_data, n_exceedances >= 100)

if (nrow(stable) > 0) {
  best <- stable[which.min(stable$shape_ci_width), ]
  cat("\n========================================\n")
  cat("ALTERNATIVE SUGGESTED THRESHOLD\n")
  cat("========================================\n")
  cat(sprintf("Suggested threshold:   %.2f m\n",  best$threshold))
  cat(sprintf("Number of exceedances: %d\n",       best$n_exceedances))
  cat(sprintf("Percentage of data:    %.2f%%\n",  100 * mean(hs > best$threshold)))
  cat(sprintf("Shape parameter:       %.3f (±%.3f)\n",
              best$shape, (best$shape_upper - best$shape) ))
  cat(sprintf("Scale parameter:       %.3f (±%.3f)\n",
              best$scale, (best$scale_upper - best$scale) ))
}

cat("\n========================================\n")
cat("DATA SUMMARY\n")
cat("========================================\n")
cat(sprintf("Total observations:  %d\n",   length(hs)))
cat(sprintf("Mean Hs:             %.2f m\n", mean(hs)))
cat(sprintf("Max Hs:              %.2f m\n", max(hs)))
cat(sprintf("95th percentile:     %.2f m\n", quantile(hs, 0.95)))
cat(sprintf("96th percentile:     %.2f m\n", quantile(hs, 0.96)))
cat(sprintf("99th percentile:     %.2f m\n", quantile(hs, 0.99)))