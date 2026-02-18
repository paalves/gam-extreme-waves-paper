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

# Define chosen threshold (96th percentile)
chosen_threshold <- quantile(hs, 0.96)

# ============================================================================
# 1. MEAN RESIDUAL LIFE PLOT
# ============================================================================
# Calculate mean excess function
thresholds <- quantile(hs, probs = seq(0.90, 0.995, by = 0.005))
mrl_data <- data.frame(
  threshold = numeric(),
  mean_excess = numeric(),
  lower = numeric(),
  upper = numeric()
)

for (u in thresholds) {
  excesses <- hs[hs > u] - u
  n <- length(excesses)
  me <- mean(excesses)
  se <- sd(excesses) / sqrt(n)
  
  mrl_data <- rbind(mrl_data, data.frame(
    threshold = u,
    mean_excess = me,
    lower = me - 1.96 * se,
    upper = me + 1.96 * se
  ))
}

p1 <- ggplot(mrl_data, aes(x = threshold, y = mean_excess)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.5) +
  geom_line(color = "black", linewidth = 0.8) +
  geom_point(size = 1.5, color = "black") +
  geom_vline(xintercept = chosen_threshold, linetype = "dashed", 
             color = "blue", linewidth = 1) +
  annotate("text", x = chosen_threshold, y = max(mrl_data$upper), 
           label = sprintf("Chosen threshold\n(Q96 = %.2f m)", chosen_threshold),
           hjust = -0.1, vjust = 1, size = 3, color = "blue") +
  labs(
    title = "Mean residual life plot",
    subtitle = "Approximate linearity indicates suitable threshold range",
    x = expression(paste("Threshold ", italic(u), " (m)")),
    y = "Mean excess (m)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.minor = element_blank()
  )

# ============================================================================
# 2. PARAMETER STABILITY PLOTS
# ============================================================================
# Fit GPD for range of thresholds and extract parameters
param_data <- data.frame(
  threshold = numeric(),
  scale = numeric(),
  shape = numeric(),
  scale_lower = numeric(),
  scale_upper = numeric(),
  shape_lower = numeric(),
  shape_upper = numeric(),
  n_exceedances = numeric()
)

for (u in thresholds) {
  excesses <- hs[hs > u] - u
  if (length(excesses) > 50) {  # Ensure sufficient exceedances
    fit <- tryCatch({
      fpot(hs, threshold = u, model = "gpd")
    }, error = function(e) NULL)
    
    if (!is.null(fit)) {
      param_data <- rbind(param_data, data.frame(
        threshold = u,
        scale = fit$param["scale"],
        shape = fit$param["shape"],
        scale_lower = fit$param["scale"] - 1.96 * fit$std.err["scale"],
        scale_upper = fit$param["scale"] + 1.96 * fit$std.err["scale"],
        shape_lower = fit$param["shape"] - 1.96 * fit$std.err["shape"],
        shape_upper = fit$param["shape"] + 1.96 * fit$std.err["shape"],
        n_exceedances = length(excesses)
      ))
    }
  }
}

# Shape parameter plot
p2 <- ggplot(param_data, aes(x = threshold, y = shape)) +
  geom_ribbon(aes(ymin = shape_lower, ymax = shape_upper), 
              fill = "grey80", alpha = 0.5) +
  geom_line(color = "black", linewidth = 0.8) +
  geom_point(size = 1.5, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.6) +
  geom_vline(xintercept = chosen_threshold, linetype = "dashed", 
             color = "blue", linewidth = 1) +
  annotate("text", x = chosen_threshold, y = max(param_data$shape_upper), 
           label = sprintf("Q96 = %.2f m", chosen_threshold),
           hjust = -0.1, vjust = 1, size = 3, color = "blue") +
  labs(
    title = "Shape parameter stability",
    subtitle = "Stable values indicate appropriate threshold range",
    x = expression(paste("Threshold ", italic(u), " (m)")),
    y = expression(paste("Shape parameter ", xi))
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.minor = element_blank()
  )

# Modified scale parameter plot
p3 <- ggplot(param_data, aes(x = threshold, y = scale)) +
  geom_ribbon(aes(ymin = scale_lower, ymax = scale_upper), 
              fill = "grey80", alpha = 0.5) +
  geom_line(color = "black", linewidth = 0.8) +
  geom_point(size = 1.5, color = "black") +
  geom_vline(xintercept = chosen_threshold, linetype = "dashed", 
             color = "blue", linewidth = 1) +
  annotate("text", x = chosen_threshold, y = max(param_data$scale_upper), 
           label = sprintf("Q96 = %.2f m", chosen_threshold),
           hjust = -0.1, vjust = 1, size = 3, color = "blue") +
  labs(
    title = "Modified scale parameter stability",
    subtitle = expression(paste("Stability in ", sigma[u] - xi %.% u, " indicates good threshold")),
    x = expression(paste("Threshold ", italic(u), " (m)")),
    y = expression(paste("Modified scale parameter ", sigma^"*" ))
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.minor = element_blank()
  )

# ============================================================================
# 3. NUMBER OF EXCEEDANCES
# ============================================================================
p4 <- ggplot(param_data, aes(x = threshold, y = n_exceedances)) +
  geom_line(color = "black", linewidth = 0.8) +
  geom_point(size = 1.5, color = "black") +
  geom_hline(yintercept = 100, linetype = "dashed", 
             color = "darkred", alpha = 0.6) +
  geom_vline(xintercept = chosen_threshold, linetype = "dashed", 
             color = "blue", linewidth = 1) +
  annotate("text", x = min(param_data$threshold), y = 110, 
           label = "Minimum recommended", hjust = 0, size = 3, color = "darkred") +
  annotate("text", x = chosen_threshold, y = max(param_data$n_exceedances), 
           label = sprintf("Q96 = %.2f m", chosen_threshold),
           hjust = -0.1, vjust = 1, size = 3, color = "blue") +
  labs(
    title = "Number of Threshold Exceedances",
    subtitle = "Balance between bias (low threshold) and variance (high threshold)",
    x = expression(paste("Threshold ", italic(u), " (m)")),
    y = "Number of Exceedances"
  ) +
  scale_y_log10(labels = comma) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.minor = element_blank()
  )

# ============================================================================
# COMBINE ALL PLOTS
# ============================================================================
# Create combined plot
combined_plot <- grid.arrange(p1, p2, p3, nrow = 1,
                              top = textGrob(
                                "GPD threshold selection diagnostics",
                                gp = gpar(fontsize = 14, fontface = "bold")
                              ))

# Save high-resolution figure
ggsave("gpd_threshold_diagnostics.png", combined_plot,
       width = 12, height = 10, dpi = 300, bg = "white")

# ============================================================================
# ANALYSIS OF CHOSEN THRESHOLD
# ============================================================================
# Fit GPD at chosen threshold
fit_chosen <- fpot(hs, threshold = chosen_threshold, model = "gpd")

cat("\n========================================\n")
cat("CHOSEN THRESHOLD ANALYSIS (Q96)\n")
cat("========================================\n")
cat(sprintf("Chosen threshold: %.2f m\n", chosen_threshold))
cat(sprintf("Number of exceedances: %d\n", sum(hs > chosen_threshold)))
cat(sprintf("Percentage of data: %.2f%%\n", 100 * mean(hs > chosen_threshold)))
cat(sprintf("Shape parameter: %.3f (±%.3f)\n",
            fit_chosen$param["shape"],
            1.96 * fit_chosen$std.err["shape"]))
cat(sprintf("Scale parameter: %.3f (±%.3f)\n",
            fit_chosen$param["scale"],
            1.96 * fit_chosen$std.err["scale"]))

# ============================================================================
# ADDITIONAL ANALYSIS: Suggest optimal threshold
# ============================================================================
# Look for stability in shape parameter (low variance in estimates)
shape_cv <- abs((param_data$shape_upper - param_data$shape_lower) / 
                  (2 * 1.96 * param_data$shape))
stable_indices <- which(shape_cv < median(shape_cv) & 
                          param_data$n_exceedances >= 100)

if (length(stable_indices) > 0) {
  suggested_threshold <- param_data$threshold[stable_indices[1]]
  cat("\n========================================\n")
  cat("ALTERNATIVE SUGGESTED THRESHOLD\n")
  cat("========================================\n")
  cat(sprintf("Suggested threshold: %.2f m\n", suggested_threshold))
  cat(sprintf("Number of exceedances: %d\n", 
              param_data$n_exceedances[stable_indices[1]]))
  cat(sprintf("Percentage of data: %.2f%%\n", 
              100 * mean(hs > suggested_threshold)))
  cat(sprintf("Shape parameter: %.3f (±%.3f)\n",
              param_data$shape[stable_indices[1]],
              1.96 * (param_data$shape_upper[stable_indices[1]] - 
                        param_data$shape[stable_indices[1]]) / 1.96))
  cat(sprintf("Scale parameter: %.3f (±%.3f)\n",
              param_data$scale[stable_indices[1]],
              1.96 * (param_data$scale_upper[stable_indices[1]] - 
                        param_data$scale[stable_indices[1]]) / 1.96))
}

# Print summary statistics
cat("\n========================================\n")
cat("DATA SUMMARY\n")
cat("========================================\n")
cat(sprintf("Total observations: %d\n", length(hs)))
cat(sprintf("Mean Hs: %.2f m\n", mean(hs)))
cat(sprintf("Max Hs: %.2f m\n", max(hs)))
cat(sprintf("95th percentile: %.2f m\n", quantile(hs, 0.95)))
cat(sprintf("96th percentile: %.2f m\n", quantile(hs, 0.96)))
cat(sprintf("99th percentile: %.2f m\n", quantile(hs, 0.99)))
cat("\nDiagnostic plots saved as 'gpd_threshold_diagnostics.png'\n")