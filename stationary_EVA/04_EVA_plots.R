library(ggplot2)
library(patchwork)
library(evd)

# Assumes 01_analysis.R has been sourced first

# ============================================================================
# THEME & PALETTE
# ============================================================================
col_gpd <- "#2166AC"
col_gev <- "#D6604D"

theme_rv <- function() {
  theme_bw(base_size = 14) +
    theme(
      plot.title       = element_text(face = "bold", size = 16, margin = margin(b = 4)),
      plot.subtitle    = element_text(size = 12, color = "grey40", margin = margin(b = 8)),
      axis.title       = element_text(size = 13),
      axis.text        = element_text(size = 12),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      panel.border     = element_rect(color = "grey70"),
      legend.position  = "none",
      plot.margin      = margin(10, 20, 8, 8)
    )
}

theme_diag <- function() {
  theme_bw(base_size = 13) +
    theme(
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(size = 11, color = "grey40"),
      axis.title       = element_text(size = 12),
      axis.text        = element_text(size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92")
    )
}

# ============================================================================
# OBSERVED DATA PLOTTING POSITIONS
# ============================================================================
n_c      <- length(cluster_maxima)
clust_df <- data.frame(
  period = sort((n_c + 1) / (n_c + 1 - rank(cluster_maxima))) / lambda_declust,
  level  = sort(cluster_maxima)
)

n_b   <- length(block_maxima)
bm_df <- data.frame(
  period = 1 / (1 - ((rank(block_maxima) - 0.44) / (n_b + 0.12))) / blocks_per_year,
  level  = sort(block_maxima)
)

y_min       <- floor(min(clust_df$level, bm_df$level) * 10) / 10
y_max       <- ceiling(max(max(rv_declust$upper),
                            quantile(rv_gev$upper, 0.85)) * 10) / 10
clust_plot  <- clust_df[clust_df$period >= 0.5, ]

# ============================================================================
# RETURN LEVEL PLOT
# ============================================================================
p_clean <- ggplot() +

  # CI ribbons
  geom_ribbon(data = rv_gev,
              aes(x = period, ymin = pmax(lower, y_min), ymax = pmin(upper, y_max + 1)),
              fill = col_gev, alpha = 0.13) +
  geom_ribbon(data = rv_declust,
              aes(x = period, ymin = lower, ymax = upper),
              fill = col_gpd, alpha = 0.18) +

  # Reference Hs line
  geom_hline(yintercept = hs_ref, linetype = "dotted",
             color = "grey20", linewidth = 0.7) +

  # Vertical drop lines
  geom_segment(aes(x = T_gev_ref,     xend = T_gev_ref,
                   y = y_min,         yend = hs_ref),
               linetype = "dotted", color = col_gev, linewidth = 0.7) +
  geom_segment(aes(x = T_declust_ref, xend = T_declust_ref,
                   y = y_min,         yend = hs_ref),
               linetype = "dotted", color = col_gpd, linewidth = 0.7) +

  # Fitted curves
  geom_line(data  = rv_gev,     aes(period, level), color = col_gev, linewidth = 1.2) +
  geom_line(data  = rv_declust, aes(period, level), color = col_gpd, linewidth = 1.2) +
  geom_point(data = rv_gev,     aes(period, level),
             color = col_gev, fill = "white", shape = 21, size = 3.2, stroke = 1.5) +
  geom_point(data = rv_declust, aes(period, level),
             color = col_gpd, fill = "white", shape = 21, size = 3.2, stroke = 1.5) +

  # T annotation labels — placed above the reference line, horizontally offset to avoid collision
  annotate("label", x = T_gev_ref * 0.70, y = hs_ref,
           label = sprintf("T = %.0f y", T_gev_ref),
           color = col_gev, fill = "white", label.size = 0.3,
           size = 6, fontface = "bold", vjust = -0.55) +
  annotate("label", x = T_declust_ref * 0.95, y = hs_ref,
           label = sprintf("T = %.0f y", T_declust_ref),
           color = col_gpd, fill = "white", label.size = 0.3,
           size = 6, fontface = "bold", vjust = 2) +

  # Hs reference label — top-left of the dotted line
  annotate("text", x = 1.05, y = hs_ref,
           label = sprintf("Hs = %.2f m", hs_ref),
           hjust = 0, vjust = -0.5, size = 6, color = "grey30") +

  # Direct curve labels — stacked clearly to the right of last points
  annotate("text", x = max(rv_gev$period) * 1.04,
           y = tail(rv_gev$level, 1) + 0.30,
           label = "GEV", color = col_gev, fontface = "bold", size = 4.2, hjust = 0) +
  annotate("text", x = max(rv_declust$period) * 1.04,
           y = tail(rv_declust$level, 1) - 0.30,
           label = "GPD", color = col_gpd, fontface = "bold", size = 4.2, hjust = 0) +

  scale_x_log10(
    breaks = c(1, 2, 5, 10, 20, 50, 100, 200),
    labels = c("1", "2", "5", "10", "20", "50", "100", "200"),
    expand = expansion(mult = c(0.02, 0.10))
  ) +
  coord_cartesian(ylim = c(y_min, y_max)) +
  labs(
    title    = "Return level plot",
    subtitle = sprintf("Shading = 95%% CI (delta method)  \u2022  Reference Hs = %.2f m", hs_ref),
    x = "Return period (years)", y = "Hs [m]"
  ) +
  theme_rv()

p_clean

# ============================================================================
# DIAGNOSTIC PLOTS — with CI in subtitles
# ============================================================================
gpd_quantiles <- function(probs, xi, sigma) {
  if (abs(xi) < 1e-6) -sigma * log(1 - probs)
  else                 (sigma / xi) * ((1 - probs)^(-xi) - 1)
}
gev_quantiles <- function(probs, mu, sigma, xi) {
  if (abs(xi) < 1e-6) mu - sigma * log(-log(probs))
  else                 mu + (sigma / xi) * ((-log(probs))^(-xi) - 1)
}

# Parameter estimates and CIs
xi1    <- fit_declust$param["shape"];  sigma1 <- fit_declust$param["scale"]
xi1_lo <- xi1    - 1.96 * fit_declust$std.err["shape"]
xi1_hi <- xi1    + 1.96 * fit_declust$std.err["shape"]
s1_lo  <- sigma1 - 1.96 * fit_declust$std.err["scale"]
s1_hi  <- sigma1 + 1.96 * fit_declust$std.err["scale"]

mu_g   <- fit_gev$estimate["loc"];    sg  <- fit_gev$estimate["scale"]
xi_g   <- fit_gev$estimate["shape"]
xi_g_lo <- xi_g - 1.96 * fit_gev$std.err["shape"]
xi_g_hi <- xi_g + 1.96 * fit_gev$std.err["shape"]
sg_lo   <- sg   - 1.96 * fit_gev$std.err["scale"]
sg_hi   <- sg   + 1.96 * fit_gev$std.err["scale"]
mu_lo   <- mu_g - 1.96 * fit_gev$std.err["loc"]
mu_hi   <- mu_g + 1.96 * fit_gev$std.err["loc"]

excesses_declust <- cluster_maxima - chosen_threshold

# QQ — GPD
n1  <- length(excesses_declust)
qq1 <- ggplot(data.frame(emp = sort(excesses_declust),
                          th  = gpd_quantiles((1:n1 - 0.5)/n1, xi1, sigma1)),
              aes(th, emp)) +
  geom_abline(slope = 1, intercept = 0, color = "grey60", linetype = "dashed") +
  geom_point(color = col_gpd, size = 2.5, alpha = 0.7) +
  labs(
    title    = "QQ - GPD (declustered)",
    subtitle = sprintf(
      "\u03be = %.3f [%.3f, %.3f]   \u03c3 = %.3f [%.3f, %.3f]",
      xi1, xi1_lo, xi1_hi, sigma1, s1_lo, s1_hi),
    x = "Theoretical quantiles", y = "Sample quantiles"
  ) +
  theme_diag()

# QQ — GEV
n_g <- length(block_maxima)
qq3 <- ggplot(data.frame(emp = sort(block_maxima),
                          th  = gev_quantiles((1:n_g - 0.5)/n_g, mu_g, sg, xi_g)),
              aes(th, emp)) +
  geom_abline(slope = 1, intercept = 0, color = "grey60", linetype = "dashed") +
  geom_point(color = col_gev, size = 2.5, alpha = 0.7) +
  labs(
    title    = "QQ - GEV (monthly maxima)",
    subtitle = sprintf(
      "\u03bc = %.3f [%.3f, %.3f]   \u03c3 = %.3f [%.3f, %.3f]   \u03be = %.3f [%.3f, %.3f]",
      mu_g, mu_lo, mu_hi, sg, sg_lo, sg_hi, xi_g, xi_g_lo, xi_g_hi),
    x = "Theoretical quantiles", y = "Sample quantiles"
  ) +
  theme_diag()

# Density — GPD
x1    <- seq(0, max(excesses_declust), length.out = 300)
dens1 <- (1/sigma1) * (1 + xi1 * x1 / sigma1)^(-(1/xi1 + 1))
d1 <- ggplot() +
  geom_histogram(data = data.frame(x = excesses_declust),
                 aes(x, after_stat(density)),
                 breaks = seq(min(excesses_declust), max(excesses_declust), length.out = 26),
                 fill = col_gpd, alpha = 0.25, color = "white") +
  geom_line(data = data.frame(x = x1, d = dens1), aes(x, d),
            color = col_gpd, linewidth = 1.2) +
  labs(
    title    = "Density - GPD (declustered)",
    subtitle = sprintf(
      "\u03be = %.3f [%.3f, %.3f]   \u03c3 = %.3f [%.3f, %.3f]",
      xi1, xi1_lo, xi1_hi, sigma1, s1_lo, s1_hi),
    x = "Excess over threshold (m)", y = "Density"
  ) +
  theme_diag()

# Density — GEV
x_g <- seq(min(block_maxima), max(block_maxima), length.out = 300)
d3  <- ggplot() +
  geom_histogram(data = data.frame(x = block_maxima),
                 aes(x, after_stat(density)),
                 breaks = seq(min(block_maxima), max(block_maxima), length.out = 26),
                 fill = col_gev, alpha = 0.25, color = "white") +
  geom_line(data = data.frame(x = x_g,
                               d = dgev(x_g, loc = mu_g, scale = sg, shape = xi_g)),
            aes(x, d), color = col_gev, linewidth = 1.2) +
  labs(
    title    = "Density - GEV (monthly maxima)",
    subtitle = sprintf(
      "\u03bc = %.3f [%.3f, %.3f]   \u03c3 = %.3f [%.3f, %.3f]   \u03be = %.3f [%.3f, %.3f]",
      mu_g, mu_lo, mu_hi, sg, sg_lo, sg_hi, xi_g, xi_g_lo, xi_g_hi),
    x = "Monthly block maximum (m)", y = "Density"
  ) +
  theme_diag()

p_diag <- (qq1 + qq3) / (d1 + d3) +
  plot_annotation(
    title = "Model diagnostics",
    theme = theme(plot.title = element_text(face = "bold", size = 15))
  )

# ============================================================================
# DISPLAY
# ============================================================================
p_diag
p_clean