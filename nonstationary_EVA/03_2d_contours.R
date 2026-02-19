library(ggplot2)
library(viridis)

# Requires objects from 01_evgam_fit.R: data

# -----------------------------------------------------------------------------
# 1. CORE FUNCTIONS
# -----------------------------------------------------------------------------
calculate_C_theta <- function(data, theta, Pe, var1 = "Cspd", var2 = "hs") {
  Y <- data[[var1]] * cos(theta) + data[[var2]] * sin(theta)
  n <- length(Y)
  k <- max(1, min(ceiling(n * (1 - Pe)), n))
  sort(Y, partial = k)[k]
}

check_supporting_hyperplane <- function(C_values, idx, delta, tolerance = 1e-6) {
  n         <- length(C_values)
  idx_minus <- ((idx - 2) %% n) + 1
  idx_plus  <- (idx %% n) + 1
  lhs       <- 0.5 * (C_values[idx_minus] + C_values[idx_plus])
  rhs       <- cos(delta) * C_values[idx]
  lhs > (rhs - tolerance)
}

smooth_C_values <- function(angles, C_values, span = 0.1) {
  angles_ext <- c(angles - 2*pi, angles, angles + 2*pi)
  C_ext      <- rep(C_values, 3)
  fit        <- loess(C_ext ~ angles_ext, span = span, degree = 1)
  predict(fit, newdata = angles)
}

# -----------------------------------------------------------------------------
# 2. HUSEBY METHOD 1
# -----------------------------------------------------------------------------
huseby_method1_improved <- function(data, Pe, n_angles = 180,
                                    var1 = "Cspd", var2 = "hs",
                                    remove_artifacts = TRUE,
                                    smooth_C = TRUE) {

  angles  <- seq(0, 2*pi, length.out = n_angles + 1)[1:n_angles]
  delta   <- angles[2] - angles[1]

  C_values <- sapply(angles, function(theta) calculate_C_theta(data, theta, Pe, var1, var2))

  if (smooth_C) C_values <- smooth_C_values(angles, C_values, span = 0.05)

  if (remove_artifacts) {
    is_supporting <- rep(TRUE, n_angles)

    for (test_delta in c(delta, 2*delta, 3*delta)) {
      n_test <- floor(2*pi / test_delta)
      if (n_test < 3) next
      test_angles <- seq(0, 2*pi, length.out = n_test + 1)[1:n_test]
      test_C      <- approx(angles, C_values, test_angles, rule = 2)$y

      for (i in 1:n_test) {
        if (!check_supporting_hyperplane(test_C, i, test_delta)) {
          nearest_idx <- which.min(abs(angles - test_angles[i]))
          is_supporting[nearest_idx] <- FALSE
        }
      }
    }

    for (i in 1:n_angles) {
      idx_prev <- ((i - 2) %% n_angles) + 1
      idx_next <- (i %% n_angles) + 1
      if (!is_supporting[i] && is_supporting[idx_prev] && is_supporting[idx_next]) {
        if (check_supporting_hyperplane(C_values, i, delta, tolerance = 1e-5)) {
          is_supporting[i] <- TRUE
        }
      }
    }

    valid_indices <- which(is_supporting)

    if (length(valid_indices) < 4) {
      warning("Too few supporting hyperplanes (", length(valid_indices), "). Using all points.")
      valid_indices <- 1:n_angles
    } else {
      cat(sprintf("Removed %d non-supporting hyperplanes (%.1f%%)\n",
                  n_angles - length(valid_indices),
                  100 * (n_angles - length(valid_indices)) / n_angles))
    }
  } else {
    valid_indices <- 1:n_angles
  }

  n_valid        <- length(valid_indices)
  contour_points <- matrix(0, nrow = n_valid, ncol = 2)

  for (i in 1:n_valid) {
    j     <- valid_indices[i]
    j_next <- valid_indices[(i %% n_valid) + 1]

    theta_j   <- angles[j];   theta_jp1 <- angles[j_next]
    C_j       <- C_values[j]; C_jp1     <- C_values[j_next]
    delta_theta <- theta_jp1 - theta_j
    if (delta_theta < 0) delta_theta <- delta_theta + 2*pi

    sin_delta <- sin(delta_theta)
    if (abs(sin_delta) > 1e-10) {
      x1 <- (sin(theta_jp1) * C_j - sin(theta_j) * C_jp1) / sin_delta
      x2 <- (-cos(theta_jp1) * C_j + cos(theta_j) * C_jp1) / sin_delta
    } else {
      x1 <- 0.5 * (C_j * cos(theta_j) + C_jp1 * cos(theta_jp1))
      x2 <- 0.5 * (C_j * sin(theta_j) + C_jp1 * sin(theta_jp1))
    }
    contour_points[i, ] <- c(x1, x2)
  }

  if (remove_artifacts && n_valid > 10) {
    centroid  <- colMeans(contour_points)
    distances <- sqrt(rowSums((contour_points - matrix(centroid, n_valid, 2, byrow = TRUE))^2))
    keep      <- distances < median(distances) * 3

    if (sum(!keep) > 0 && sum(keep) > 10) {
      cat(sprintf("Removed %d outlier points in post-processing\n", sum(!keep)))
      contour_points <- contour_points[keep, , drop = FALSE]
    }
  }

  contour_df <- data.frame(
    var1 = c(contour_points[, 1], contour_points[1, 1]),
    var2 = c(contour_points[, 2], contour_points[1, 2])
  )
  names(contour_df)[1:2] <- c(var1, var2)
  contour_df
}

# -----------------------------------------------------------------------------
# 3. CALCULATE CONTOURS FOR MULTIPLE RETURN PERIODS
# -----------------------------------------------------------------------------
calculate_environmental_contours <- function(data,
                                             return_periods = c(1, 2, 5, 10, 20, 50, 100),
                                             n_angles = 180,
                                             var1 = "Cspd", var2 = "hs",
                                             remove_artifacts = TRUE) {

  Pe_values <- 1 / (return_periods * 365.25 * 24)

  cat("=== Calculating Environmental Contours ===\n")
  cat(sprintf("Data points: %d | Angles: %d | Remove artifacts: %s\n\n",
              nrow(data), n_angles, remove_artifacts))

  contours_list <- list()

  for (i in seq_along(return_periods)) {
    rp <- return_periods[i]; Pe <- Pe_values[i]
    cat(sprintf("Calculating %3d-year contour (Pe = %.2e)... ", rp, Pe))

    contour            <- huseby_method1_improved(data, Pe, n_angles, var1, var2, remove_artifacts, smooth_C = TRUE)
    contour$return_period <- rp
    contours_list[[paste0("RP", rp)]] <- contour

    max_val <- max(contour[[var2]], na.rm = TRUE)
    max_idx <- which.max(contour[[var2]])
    cat(sprintf("max %s = %.2f at %s = %.2f\n", var2, max_val, var1, contour[[var1]][max_idx]))
  }

  list(
    contours_list  = contours_list,
    all_contours   = do.call(rbind, contours_list),
    data           = data,
    return_periods = return_periods
  )
}

# -----------------------------------------------------------------------------
# 4. PLOT CONTOURS
# -----------------------------------------------------------------------------
plot_environmental_contours <- function(contour_result,
                                        sample_size = nrow(contour_result$data),
                                        point_alpha = 1,
                                        var1_label = "Current Speed (m/s)",
                                        var2_label = "Significant Wave Height (m)") {

  data           <- contour_result$data
  all_contours   <- contour_result$all_contours
  return_periods <- contour_result$return_periods
  var1 <- names(all_contours)[1]; var2 <- names(all_contours)[2]

  if (!is.null(sample_size) && nrow(data) > sample_size) {
    set.seed(42)
    data <- data[sample(1:nrow(data), sample_size), ]
  }

  ggplot() +
    geom_point(data = data, aes(x = .data[[var1]], y = .data[[var2]]),
               color = "black", size = 0.5, alpha = point_alpha) +
    geom_path(data = all_contours,
              aes(x = .data[[var1]], y = .data[[var2]],
                  color = as.factor(return_period), group = return_period),
              linewidth = 1.2) +
    scale_color_viridis_d(option = "plasma", name = "Return Period (years)",
                          breaks = return_periods, labels = return_periods) +
    labs(title    = "Environmental Contours - Huseby Method 1",
         subtitle = paste("Return periods:", paste(return_periods, collapse = ", "), "years"),
         x = var1_label, y = var2_label) +
    theme_minimal(base_size = 12) +
    theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          legend.position  = "right",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "gray90")) +
    coord_cartesian(xlim = range(data[[var1]]), ylim = range(data[[var2]]))
}

print_contour_summary <- function(contour_result) {
  contours_list  <- contour_result$contours_list
  return_periods <- contour_result$return_periods
  var1 <- names(contours_list[[1]])[1]; var2 <- names(contours_list[[1]])[2]

  cat("\n=== Contour Summary ===\n")
  for (rp in return_periods) {
    contour <- contours_list[[paste0("RP", rp)]]
    max_idx <- which.max(contour[[var2]])
    cat(sprintf("\n%d-year:  %s in [%.2f, %.2f]  |  %s in [%.2f, %.2f]  |  max %s = %.2f at %s = %.2f\n",
                rp,
                var1, min(contour[[var1]]), max(contour[[var1]]),
                var2, min(contour[[var2]]), max(contour[[var2]]),
                var2, contour[[var2]][max_idx], var1, contour[[var1]][max_idx]))
  }
  cat("\n=== COMPLETE ===\n")
}

# -----------------------------------------------------------------------------
# 5. RUN
# -----------------------------------------------------------------------------
result <- calculate_environmental_contours(
  data,
  return_periods   = c(1, 2, 5, 10, 20, 50, 100),
  n_angles         = 180,
  var1             = "tp",
  var2             = "hs",
  remove_artifacts = TRUE
)

result$all_contours$hs <- pmax(0, result$all_contours$hs)

p <- plot_environmental_contours(result, var1_label = "Tp (s)", var2_label = "Hs (m)")
print(p)

print_contour_summary(result)