library(data.table)
library(POT)
library(evd)
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)
library(viridis)

# Load data
#data <- fread("C:/these_docs/data2024.csv")
data <- fread("C:/Users/minim/Documents/teletravail/data2024.csv")

hs <- na.omit(data$hs)
chosen_threshold <- quantile(hs, 0.96)

# ============================================================================
# DECLUSTERING FUNCTION
# ============================================================================
decluster_runs <- function(x, threshold, run_length) {
  exceed <- x > threshold
  if (sum(exceed) == 0) return(list(n_clusters = 0, cluster_indices = integer(0)))
  
  exceed_indices <- which(exceed)
  gaps <- diff(exceed_indices)
  new_cluster <- c(TRUE, gaps > run_length)
  cluster_id <- cumsum(new_cluster)
  
  cluster_indices <- sapply(unique(cluster_id), function(i) {
    members <- exceed_indices[cluster_id == i]
    members[which.max(x[members])]
  })
  
  list(n_clusters = length(cluster_indices), cluster_indices = cluster_indices)
}

# ============================================================================
# SENSITIVITY ANALYSIS
# ============================================================================
thresholds <- quantile(hs, probs = seq(0.90, 0.995, by = 0.01))
run_lengths <- c(12, 24, 48, 72, 96, 120)
cluster_results <- data.frame()

for (u in thresholds) {
  n_exceedances <- sum(hs > u)
  
  for (r in run_lengths) {
    clusters <- decluster_runs(hs, threshold = u, run_length = r)
    n_events <- clusters$n_clusters
    
    if (n_events >= 50 && n_exceedances > 0) {
      cluster_maxima <- hs[clusters$cluster_indices]
      
      fit <- tryCatch({
        fpot(cluster_maxima, threshold = u, model = "gpd")
      }, error = function(e) NULL)
      
      if (!is.null(fit) && !any(is.na(fit$param)) && !any(is.na(fit$std.err))) {
        cluster_results <- rbind(cluster_results, data.frame(
          threshold = u, run_length = r, n_exceedances = n_exceedances,
          n_clusters = n_events, cluster_ratio = n_events / n_exceedances,
          shape = fit$param["shape"], scale = fit$param["scale"],
          shape_se = fit$std.err["shape"], scale_se = fit$std.err["scale"]
        ))
      }
    }
  }
}

# ============================================================================
# PLOTS
# ============================================================================
# Define theme with large blue text
plot_theme <- theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

# Plot 1: Number of Independent Events
p1 <- ggplot(cluster_results, aes(x = threshold, y = n_clusters, 
                                   color = as.factor(run_length), group = run_length)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_vline(xintercept = chosen_threshold, linetype = "dashed", color = "blue", linewidth = 1) +
  annotate("text", x = chosen_threshold, y = max(cluster_results$n_clusters), 
           label = sprintf("Q96 = %.2f m", chosen_threshold), hjust = -0.1, size = 6, fontface="bold", color = "blue") +
  scale_color_viridis_d(name = "Run length\n(hours)") +
  labs(title = "Number of independent storm events",
       subtitle = "Effect of threshold and run length on identified clusters",
       x = expression(paste("Threshold ", italic(u), " (m)")),
       y = "Number of independent events") +
  plot_theme

# Plot 2: Cluster Ratio
p2 <- ggplot(cluster_results, aes(x = threshold, y = cluster_ratio, 
                                   color = as.factor(run_length), group = run_length)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_vline(xintercept = chosen_threshold, linetype = "dashed", color = "blue", linewidth = 1) +
  annotate("text", x = chosen_threshold, y = max(cluster_results$cluster_ratio), 
           label = sprintf("Q96 = %.2f m", chosen_threshold), hjust = -0.1, size = 6, fontface="bold", color = "blue") +
  scale_color_viridis_d(name = "Run length\n(hours)") +
  labs(title = "Extremal index change",
       subtitle = "Ratio of independent events to total exceedances",
       x = expression(paste("Threshold ", italic(u), " (m)")),
       y = "Extremal index (θ)") +
  plot_theme

# Plot 3: Shape Parameter
p3 <- ggplot(cluster_results, aes(x = threshold, y = shape, 
                                   color = as.factor(run_length), group = run_length)) +
  geom_ribbon(aes(ymin = shape - 1.96 * shape_se, ymax = shape + 1.96 * shape_se, 
                  fill = as.factor(run_length)), alpha = 0.1, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_vline(xintercept = chosen_threshold, linetype = "dashed", color = "blue", linewidth = 1) +
  annotate("text", x = chosen_threshold, 
           y = max(cluster_results$shape + 1.96 * cluster_results$shape_se, na.rm = TRUE), 
           label = sprintf("Q96 = %.2f m", chosen_threshold), hjust = -0.1, size = 6, fontface="bold", color = "blue") +
  scale_color_viridis_d(name = "Run length\n(hours)") +
  scale_fill_viridis_d(name = "Run length\n(hours)") +
  labs(title = "Shape parameter stability",
       subtitle = "Effect of run length on GPD shape parameter estimates",
       x = expression(paste("Threshold ", italic(u), " (m)")),
       y = expression(paste("Shape parameter ", xi))) +
  plot_theme

# Plot 4: Scale Parameter
p4 <- ggplot(cluster_results, aes(x = threshold, y = scale, 
                                   color = as.factor(run_length), group = run_length)) +
  geom_ribbon(aes(ymin = scale - 1.96 * scale_se, ymax = scale + 1.96 * scale_se, 
                  fill = as.factor(run_length)), alpha = 0.1, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_vline(xintercept = chosen_threshold, linetype = "dashed", color = "blue", linewidth = 1) +
  annotate("text", x = chosen_threshold, 
           y = max(cluster_results$scale + 1.96 * cluster_results$scale_se, na.rm = TRUE), 
           label = sprintf("Q96 = %.2f m", chosen_threshold), hjust = -0.1, size = 6, fontface="bold", color = "blue") +
  scale_color_viridis_d(name = "Run length\n(hours)") +
  scale_fill_viridis_d(name = "Run length\n(hours)") +
  labs(title = "Scale parameter stability",
       subtitle = "Effect of run length on GPD scale parameter estimates",
       x = expression(paste("Threshold ", italic(u), " (m)")),
       y = expression(paste("Scale parameter ", sigma))) +
  plot_theme

# ============================================================================
# COMBINE PLOTS
# ============================================================================
grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2,
             top = textGrob("Declustering sensitivity analysis", 
                            gp = gpar(fontsize = 20, fontface = "bold")))
