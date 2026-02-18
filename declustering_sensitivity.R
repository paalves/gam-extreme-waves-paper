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
# MANUAL DECLUSTERING FUNCTION
# ============================================================================
# Simple runs-based declustering
decluster_runs <- function(x, threshold, run_length) {
  # Identify exceedances
  exceed <- x > threshold
  
  # Initialize cluster indicator
  cluster_max <- rep(FALSE, length(x))
  
  if (sum(exceed) == 0) return(list(n_clusters = 0, cluster_indices = integer(0)))
  
  # Find the start and end of each exceedance cluster
  exceed_indices <- which(exceed)
  
  if (length(exceed_indices) == 0) {
    return(list(n_clusters = 0, cluster_indices = integer(0)))
  }
  
  # Identify gaps larger than run_length
  gaps <- diff(exceed_indices)
  new_cluster <- c(TRUE, gaps > run_length)
  
  # Mark cluster maxima
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

# ============================================================================
# DECLUSTERING SENSITIVITY ANALYSIS
# ============================================================================

# Range of thresholds to test
thresholds <- quantile(hs, probs = seq(0.90, 0.995, by = 0.01))

# Range of run lengths (time separation between clusters in hours)
run_lengths <- c(12, 24, 48, 72, 96, 120)

# Storage for results
cluster_results <- data.frame()

cat("Analyzing declustering sensitivity...\n")
cat("Using manual runs-based declustering method\n\n")

pb <- txtProgressBar(min = 0, max = length(thresholds), style = 3)

for (i in seq_along(thresholds)) {
  u <- thresholds[i]
  n_exceedances <- sum(hs > u)
  
  for (r in run_lengths) {
    # Decluster the data
    clusters <- decluster_runs(hs, threshold = u, run_length = r)
    n_events <- clusters$n_clusters
    
    if (n_events >= 50 && n_exceedances > 0) {
      # Extract cluster maxima
      cluster_maxima <- hs[clusters$cluster_indices]
      excesses <- cluster_maxima - u
      
      # Fit GPD to declustered excesses
      fit <- tryCatch({
        # Use MLE fitting
        fit_result <- fpot(cluster_maxima, threshold = u, model = "gpd")
        fit_result
      }, error = function(e) {
        NULL
      })
      
      if (!is.null(fit) && !any(is.na(fit$param)) && !any(is.na(fit$std.err))) {
        cluster_results <- rbind(cluster_results, data.frame(
          threshold = u,
          run_length = r,
          n_exceedances = n_exceedances,
          n_clusters = n_events,
          cluster_ratio = n_events / n_exceedances,
          shape = fit$param["shape"],
          scale = fit$param["scale"],
          shape_se = fit$std.err["shape"],
          scale_se = fit$std.err["scale"]
        ))
      }
    }
  }
  setTxtProgressBar(pb, i)
}
close(pb)

cat("\n\nDeclustering analysis complete.\n")
cat(sprintf("Generated %d valid fits.\n\n", nrow(cluster_results)))

# Check if we have any results
if (nrow(cluster_results) == 0) {
  cat("ERROR: No valid fits generated.\n")
  cat("This may indicate issues with the data or threshold selection.\n")
  stop("No results to plot.")
}

# ============================================================================
# PLOT 1: Number of Independent Events vs Threshold
# ============================================================================
p1 <- ggplot(cluster_results, aes(x = threshold, y = n_clusters, 
                                   color = as.factor(run_length), 
                                   group = run_length)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  geom_vline(xintercept = chosen_threshold, linetype = "dashed", 
             color = "blue", linewidth = 1) +
  annotate("text", x = chosen_threshold, y = max(cluster_results$n_clusters), 
           label = sprintf("Q96 = %.2f m", chosen_threshold),
           hjust = -0.1, vjust = 1, size = 3.5) +
  scale_color_viridis_d(name = "Run length\n(hours)") +
  labs(
    title = "Number of independent storm events",
    subtitle = "Effect of threshold and run length on identified clusters",
    x = expression(paste("Threshold ", italic(u), " (m)")),
    y = "Number of independent events"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# ============================================================================
# PLOT 2: Cluster Ratio vs Threshold
# ============================================================================
p2 <- ggplot(cluster_results, aes(x = threshold, y = cluster_ratio, 
                                   color = as.factor(run_length), 
                                   group = run_length)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  geom_vline(xintercept = chosen_threshold, linetype = "dashed", 
             color = "blue", linewidth = 1) +
  annotate("text", x = chosen_threshold, y = max(cluster_results$cluster_ratio), 
           label = sprintf("Q96 = %.2f m", chosen_threshold),
           hjust = -0.1, vjust = 1, size = 3.5) +
  scale_color_viridis_d(name = "Run length\n(hours)") +
  labs(
    title = "Extremal index",
    subtitle = "Ratio of independent events to total exceedances",
    x = expression(paste("Threshold ", italic(u), " (m)")),
    y = "Extremal index (θ)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# ============================================================================
# PLOT 3: Shape Parameter Stability with Declustering
# ============================================================================
p3 <- ggplot(cluster_results, aes(x = threshold, y = shape, 
                                   color = as.factor(run_length), 
                                   group = run_length)) +
  geom_ribbon(aes(ymin = shape - 1.96 * shape_se, 
                  ymax = shape + 1.96 * shape_se,
                  fill = as.factor(run_length)), 
              alpha = 0.1, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_vline(xintercept = chosen_threshold, linetype = "dashed", 
             color = "blue", linewidth = 1) +
  annotate("text", x = chosen_threshold, 
           y = max(cluster_results$shape + 1.96 * cluster_results$shape_se, na.rm = TRUE), 
           label = sprintf("Q96 = %.2f m", chosen_threshold),
           hjust = -0.1, vjust = 1, size = 3.5) +
  scale_color_viridis_d(name = "Run length\n(hours)") +
  scale_fill_viridis_d(name = "Run length\n(hours)") +
  labs(
    title = "Shape parameter stabiility",
    subtitle = "Effect of run length on GPD shape parameter estimates",
    x = expression(paste("Threshold ", italic(u), " (m)")),
    y = expression(paste("Shape parameter ", xi))
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# ============================================================================
# PLOT 4: Scale Parameter Stability with Declustering
# ============================================================================
p4 <- ggplot(cluster_results, aes(x = threshold, y = scale, 
                                   color = as.factor(run_length), 
                                   group = run_length)) +
  geom_ribbon(aes(ymin = scale - 1.96 * scale_se, 
                  ymax = scale + 1.96 * scale_se,
                  fill = as.factor(run_length)), 
              alpha = 0.1, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  geom_vline(xintercept = chosen_threshold, linetype = "dashed", 
             color = "blue", linewidth = 1) +
  annotate("text", x = chosen_threshold, 
           y = max(cluster_results$scale + 1.96 * cluster_results$scale_se, na.rm = TRUE), 
           label = sprintf("Q96 = %.2f m", chosen_threshold),
           hjust = -0.1, vjust = 1, size = 3.5) +
  scale_color_viridis_d(name = "Run length\n(hours)") +
  scale_fill_viridis_d(name = "Run length\n(hours)") +
  labs(
    title = "Scale parameter with declustering",
    subtitle = "Effect of run length on GPD scale parameter estimates",
    x = expression(paste("Threshold ", italic(u), " (m)")),
    y = expression(paste("Scale parameter ", sigma))
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# ============================================================================
# COMBINE PLOTS
# ============================================================================
combined_plot <- grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2,
                              top = textGrob(
                                "Declustering sensitivity analysis",
                                gp = gpar(fontsize = 14, fontface = "bold")
                              ))

# Save high-resolution figure
ggsave("declustering_sensitivity.png", combined_plot,
       width = 14, height = 10, dpi = 300, bg = "white")

# ============================================================================
# ANALYSIS AT CHOSEN THRESHOLD
# ============================================================================
cat("\n========================================\n")
cat("DECLUSTERING ANALYSIS AT Q96 THRESHOLD\n")
cat("========================================\n")
cat(sprintf("Chosen threshold: %.2f m\n\n", chosen_threshold))

# Extract results for chosen threshold (find closest match)
threshold_tolerance <- 0.05
chosen_results <- cluster_results[abs(cluster_results$threshold - chosen_threshold) < threshold_tolerance, ]

if (nrow(chosen_results) > 0) {
  # Use the closest threshold
  closest_idx <- which.min(abs(unique(chosen_results$threshold) - chosen_threshold))
  actual_threshold <- unique(chosen_results$threshold)[closest_idx]
  chosen_results <- chosen_results[chosen_results$threshold == actual_threshold, ]
  
  cat(sprintf("Using threshold: %.2f m (closest to Q96)\n\n", actual_threshold))
  
  cat("Run Length | Events | Exceedances | θ      | ξ (shape) | σ (scale)\n")
  cat("-----------|--------|-------------|--------|-----------|----------\n")
  
  for (i in 1:nrow(chosen_results)) {
    cat(sprintf("%6d hrs | %6d | %11d | %.4f | %9.4f | %9.4f\n",
                chosen_results$run_length[i],
                chosen_results$n_clusters[i],
                chosen_results$n_exceedances[i],
                chosen_results$cluster_ratio[i],
                chosen_results$shape[i],
                chosen_results$scale[i]))
  }
  
  # Recommend run length
  cat("\n========================================\n")
  cat("RECOMMENDATIONS\n")
  cat("========================================\n")
  
  # Look for stability in extremal index
  theta_cv <- sd(chosen_results$cluster_ratio) / mean(chosen_results$cluster_ratio)
  
  # Typical storm duration for wave data is 48-72 hours
  recommended_runs <- chosen_results[chosen_results$run_length >= 48 & 
                                      chosen_results$run_length <= 72, ]
  
  if (nrow(recommended_runs) > 0) {
    cat(sprintf("Recommended run length: %d hours\n", 
                recommended_runs$run_length[1]))
    cat(sprintf("  - Independent events: %d\n", 
                recommended_runs$n_clusters[1]))
    cat(sprintf("  - Extremal index: %.4f\n", 
                recommended_runs$cluster_ratio[1]))
    cat(sprintf("  - Shape: %.4f (±%.4f)\n",
                recommended_runs$shape[1],
                1.96 * recommended_runs$shape_se[1]))
    cat(sprintf("  - Scale: %.4f (±%.4f)\n",
                recommended_runs$scale[1],
                1.96 * recommended_runs$scale_se[1]))
  }
  
  cat(sprintf("\nExtremal index variability (CV): %.2f%%\n", 100 * theta_cv))
  cat("\nNote: Extremal index θ represents the proportion of exceedances\n")
  cat("that are cluster maxima (independent events). Lower values indicate\n")
  cat("stronger temporal dependence in the data.\n")
} else {
  cat("No results available near the chosen threshold.\n")
  cat(sprintf("Available thresholds range from %.2f to %.2f m\n",
              min(cluster_results$threshold), max(cluster_results$threshold)))
}

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================
cat("\n========================================\n")
cat("SUMMARY STATISTICS\n")
cat("========================================\n")
cat(sprintf("Total observations: %d\n", length(hs)))
cat(sprintf("Thresholds tested: %d\n", length(thresholds)))
cat(sprintf("Run lengths tested: %d\n", length(run_lengths)))
cat(sprintf("Range of run lengths: %d - %d hours\n", 
            min(run_lengths), max(run_lengths)))
cat(sprintf("Successful fits: %d out of %d combinations\n",
            nrow(cluster_results), length(thresholds) * length(run_lengths)))

cat("\nPlots saved as 'declustering_sensitivity.png'\n")