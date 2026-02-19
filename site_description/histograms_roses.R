library(data.table)
library(ggplot2)
library(gridExtra)
library(dplyr)

# ============================================================================
# DATA
# ============================================================================
data <- fread("C:/Users/minim/Documents/teletravail/data2024.csv")

data[, TideCondition := fifelse(fje == "flot",   "Flood",
                         fifelse(fje == "jusant", "Ebb", "Slack"))]
data$TideCondition <- factor(data$TideCondition, levels = c("Slack", "Flood", "Ebb"))

data <- data %>%
  mutate(dir_bin = as.numeric(as.character(cut(
    Cdirdeg %% 360,
    breaks = seq(0, 360, by = 10),
    include.lowest = TRUE,
    labels = seq(5, 355, by = 10)
  ))))

# ============================================================================
# THEMES
# ============================================================================
theme_dist <- theme_bw(base_size = 18) +
  theme(
    plot.title       = element_text(face = "bold", size = 20),
    axis.title       = element_text(size = 16),
    axis.text        = element_text(size = 14),
    panel.grid.minor = element_blank()
  )

theme_rose <- theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_text(face = "bold", size = 16, color = "gray20", hjust = 0.5),
    axis.title       = element_blank(),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.text.x      = element_text(color = "gray50", size = 10),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.position  = "none"
  )

tide_colors <- c("Slack" = "#66BB6A", "Flood" = "#5C6BC0", "Ebb" = "#EF5350")

# ============================================================================
# DISTRIBUTION HISTOGRAMS
# ============================================================================
make_hist <- function(var, fill_col, title, xlab, adjust = 1.5) {
  brk <- seq(0, max(data[[var]], na.rm = TRUE), length.out = 50)
  ggplot(data, aes(x = .data[[var]])) +
    geom_histogram(aes(y = after_stat(density)), breaks = brk,
                   fill = fill_col, color = "black", linewidth = 0.2) +
    geom_density(linewidth = 0.75, color = "red", adjust = adjust) +
    labs(title = title, x = xlab, y = "Density") +
    theme_dist
}

p1 <- make_hist("hs",               "lightblue",  "Significant wave height distribution", expression(H[s]~"[m]"))
p2 <- make_hist("tp",               "lightgreen", "Wave peak period distribution",        expression(T[p]~"[s]"))
p3 <- make_hist("Cspd",             "lavender",   "Current speed distribution",           "Current speed [m/s]", adjust = 2)
p4 <- make_hist("OrbSpeed_current", "orange",     "Orbital speed distribution",           "Orbital speed [m/s]")

# ============================================================================
# DIRECTION ROSES
# ============================================================================
rose_scale_x <- scale_x_continuous(
  limits = c(0, 360),
  breaks = seq(0, 315, 45),
  labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
)

p_rose_cur <- ggplot(data, aes(x = dir_bin)) +
  geom_bar(aes(fill = after_stat(count)), width = 10) +
  coord_polar(start = 0) +
  rose_scale_x +
  scale_fill_gradient(low = "#E0F7FA", high = "purple") +
  labs(title = "Current direction rose") +
  theme_rose

p_rose_wave <- ggplot(data, aes(x = dp)) +
  geom_histogram(bins = 16, center = 11.25, fill = "darkblue",
                 color = "white", linewidth = 0.1) +
  coord_polar(start = 0) +
  rose_scale_x +
  labs(title = "Wave direction rose") +
  theme_rose

# ============================================================================
# ANGULAR DIFFERENCE HISTOGRAM
# ============================================================================
p_hist_diff <- ggplot(data, aes(x = angleDiff, fill = TideCondition, color = TideCondition)) +
  geom_histogram(breaks = seq(0, 360, by = 5), linewidth = 0.1, alpha = 0.45, aes(y=..density..)) +
  geom_density(linewidth = 0.8, fill = NA, adjust = 1.5, key_glyph = "path") +
  scale_fill_manual(values = tide_colors) +
  scale_color_manual(values = tide_colors) +
  scale_x_continuous(limits = c(0, 360), breaks = seq(0, 360, by = 45)) +
  labs(
    title = "Angular difference between wave and current directions",
    subtitle = "180° means waves are opposing currents. 0°/360° means waves are following currents.\n 90°/270° means waves are perpendicular to currents.",
    x     = "Angular difference [°]",
    y     = "Density",
    fill  = "Tide Condition",
    color = "Tide Condition"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title            = element_text(face = "bold", size = 16, color = "gray20"),
    panel.grid.minor      = element_blank(),
    panel.border          = element_rect(color = "gray80"),
    legend.position       = "top",
    legend.box.background = element_rect(color = "gray90", fill = "white")
  )


# ============================================================================
# LAYOUT
# ============================================================================
grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)

grid.arrange(
  arrangeGrob(p_rose_cur, p_rose_wave, ncol = 2),
  p_hist_diff,
  nrow = 2,
  heights = c(1, 0.8)
)
