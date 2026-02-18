library(plotly)
library(geometry)
library(akima)
library(abind)

# Requires objects from 01_threshold_gpd_fit.R: data
# Requires objects from 02_simulation.R: simulate_structured_storms, best_fit,
#   fit_gpd, fit_ald, exceedance_rate, mean_cluster_size

# -----------------------------------------------------------------------------
# 1. VANEM (2019) 3D CONTOUR FUNCTION
# -----------------------------------------------------------------------------
calculate_vanem_3d_contour <- function(data, var_names, return_period, n_steps = 50) {

  X_raw <- data[[var_names[1]]]; Y_raw <- data[[var_names[2]]]; Z_raw <- data[[var_names[3]]]
  n_samples <- nrow(data)
  Pe <- 1 / (return_period * 365.25 * 24)
  if (n_samples * Pe < 1) Pe <- 1 / n_samples
  k <- floor(n_samples * (1 - Pe))

  mean_X <- mean(X_raw); sd_X <- sd(X_raw)
  mean_Y <- mean(Y_raw); sd_Y <- sd(Y_raw)
  mean_Z <- mean(Z_raw); sd_Z <- sd(Z_raw)
  X <- (X_raw - mean_X) / sd_X
  Y <- (Y_raw - mean_Y) / sd_Y
  Z <- (Z_raw - mean_Z) / sd_Z

  t1_seq <- seq(0, 2*pi, length.out = n_steps + 1)
  t2_seq <- seq(0, 2*pi, length.out = n_steps + 1)
  d1 <- t1_seq[2] - t1_seq[1]; d2 <- t2_seq[2] - t2_seq[1]

  C_mat <- matrix(0, nrow = length(t1_seq), ncol = length(t2_seq))

  cat(paste("Calculating projections for", return_period, "yr contour...\n"))
  for (j in seq_along(t2_seq)) {
    t2 <- t2_seq[j]; sin_t2 <- sin(t2); cos_t2 <- cos(t2)
    for (i in seq_along(t1_seq)) {
      t1 <- t1_seq[i]
      R  <- X*cos(t1)*cos_t2 + Y*sin(t1)*cos_t2 + Z*sin_t2
      C_mat[i, j] <- -sort(-R, partial = n_samples - k + 1)[n_samples - k + 1]
    }
  }

  # Smooth C-matrix
  C_pad <- rbind(C_mat[nrow(C_mat), ], C_mat, C_mat[1, ])
  C_pad <- cbind(C_pad[, ncol(C_pad)], C_pad, C_pad[, 1])
  C_smooth <- C_mat
  for (i in 1:nrow(C_mat)) for (j in 1:ncol(C_mat)) C_smooth[i, j] <- mean(C_pad[i:(i+2), j:(j+2)])
  C_mat <- C_smooth

  points_list <- list(); counter <- 1

  for (j in 1:n_steps) {
    t2 <- t2_seq[j]
    if (abs(cos(t2)) < 1e-4) next

    for (i in 1:n_steps) {
      t1       <- t1_seq[i]
      idx_i    <- i; idx_i_next <- if (i == n_steps) 1 else i + 1
      idx_j    <- j; idx_j_next <- if (j == n_steps) 1 else j + 1

      C1 <- C_mat[idx_i, idx_j]; C2 <- C_mat[idx_i_next, idx_j]; C3 <- C_mat[idx_i, idx_j_next]

      D_val <- (cos(t1)*sin(t1+d1) - sin(t1)*cos(t1+d1)) * cos(t2) * (cos(t2)*sin(t2+d2) - sin(t2)*cos(t2+d2))
      if (abs(D_val) < 1e-5) next

      Ax <- sin(t1+d1)*cos(t2)*sin(t2+d2) - sin(t1)*sin(t2)*cos(t2+d2)
      Bx <- sin(t1) * (sin(t2)*cos(t2+d2) - cos(t2)*sin(t2+d2))
      Cx <- -(sin(t1+d1) - sin(t1)) * cos(t2) * sin(t2)

      Ay <- cos(t1)*sin(t2)*cos(t2+d2) - cos(t1+d1)*cos(t2)*sin(t2+d2)
      By <- cos(t1) * (cos(t2)*sin(t2+d2) - sin(t2)*cos(t2+d2))
      Cy <- (cos(t1+d1) - cos(t1)) * cos(t2) * sin(t2)

      val_x <- (Ax*C1 + Bx*C2 + Cx*C3) / D_val
      val_y <- (Ay*C1 + By*C2 + Cy*C3) / D_val
      val_z <- (cos(t2+d2)*C1 - cos(t2)*C3) / (sin(t2)*cos(t2+d2) - cos(t2)*sin(t2+d2))

      points_list[[counter]] <- c(val_x, val_y, val_z)
      counter <- counter + 1
    }
  }

  pts_norm  <- do.call(rbind, points_list)
  dist_sq   <- rowSums(pts_norm^2)
  pts_norm  <- pts_norm[dist_sq < median(dist_sq) * 10, ]

  final_pts <- pts_norm
  final_pts[, 1] <- final_pts[, 1] * sd_X + mean_X
  final_pts[, 2] <- final_pts[, 2] * sd_Y + mean_Y
  final_pts[, 3] <- final_pts[, 3] * sd_Z + mean_Z

  ts <- convhulln(final_pts)

  list(
    x = final_pts[, 1], y = final_pts[, 2], z = final_pts[, 3],
    i = ts[, 1] - 1, j = ts[, 2] - 1, k = ts[, 3] - 1
  )
}

# -----------------------------------------------------------------------------
# 2. PLOT SINGLE 3D CONTOUR (ORIGINAL DATA)
# -----------------------------------------------------------------------------
vars           <- c("Cspd", "angularDifference", "hs")
contour_mesh   <- calculate_vanem_3d_contour(data, vars, return_period = 100, n_steps = 50)
contour_mesh$x <- pmax(0, contour_mesh$x)
contour_mesh$y <- pmax(0, pmin(360, contour_mesh$y))
contour_mesh$z <- pmax(0, contour_mesh$z)

max_idx   <- which.max(contour_mesh$z)
max_label <- paste0("Hs: ", round(contour_mesh$z[max_idx], 2), " m",
                    "<br>Current speed: ", round(contour_mesh$x[max_idx], 2), " m.s-1",
                    "<br>Angular difference: ", round(contour_mesh$y[max_idx], 1), "°")

fig <- plot_ly() %>%
  add_trace(x = data$Cspd, y = data$angularDifference, z = data$hs,
            type = "scatter3d", mode = "markers",
            marker = list(size = 2, opacity = 0.8, color = "black"),
            name = "Observations") %>%
  add_trace(type = "mesh3d",
            x = contour_mesh$x, y = contour_mesh$y, z = contour_mesh$z,
            i = contour_mesh$i, j = contour_mesh$j, k = contour_mesh$k,
            intensity = contour_mesh$z,
            colorscale = list(c(0, "green"), c(1, "red")),
            opacity = 0.4, showscale = TRUE,
            colorbar = list(title = "Hs (m)"),
            name = "100-yr Contour") %>%
  add_trace(x = c(contour_mesh$x[max_idx]), y = c(contour_mesh$y[max_idx]), z = c(contour_mesh$z[max_idx]),
            type = "scatter3d", mode = "markers+text",
            marker = list(size = 8, color = "red", symbol = "diamond"),
            text = c(max_label), textposition = "top center",
            name = "Max Hs", hoverinfo = "text") %>%
  layout(title  = list(text = "3D smoothed environmental contour (original data)", y = 0.9, yanchor = "bottom"),
         scene  = list(xaxis = list(title = "Current speed [m.s-1]"),
                       yaxis = list(title = "Angular difference [°]"),
                       zaxis = list(title = "Hs [m]")),
         margin = list(t = 50))
fig

# -----------------------------------------------------------------------------
# 3. PLOT SINGLE 3D CONTOUR (SIMULATED DATA)
# -----------------------------------------------------------------------------
contour_mesh_sim <- calculate_vanem_3d_contour(simulated_data, vars, return_period = 100, n_steps = 50)
contour_mesh_sim$x <- pmax(0, contour_mesh_sim$x)
contour_mesh_sim$y <- pmax(0, pmin(360, contour_mesh_sim$y))
contour_mesh_sim$z <- pmax(0, contour_mesh_sim$z)

max_idx_s   <- which.max(contour_mesh_sim$z)
max_label_s <- paste0("Hs: ", round(contour_mesh_sim$z[max_idx_s], 2), " m",
                      "<br>Current speed: ", round(contour_mesh_sim$x[max_idx_s], 2), " m/s",
                      "<br>Angular difference: ", round(contour_mesh_sim$y[max_idx_s], 1), "°")

fig_sim <- plot_ly() %>%
  add_trace(x = simulated_data$Cspd, y = simulated_data$angularDifference, z = simulated_data$hs,
            type = "scatter3d", mode = "markers",
            marker = list(size = 2, opacity = 0.8, color = "black"),
            name = "Simulated data") %>%
  add_trace(type = "mesh3d",
            x = contour_mesh_sim$x, y = contour_mesh_sim$y, z = contour_mesh_sim$z,
            i = contour_mesh_sim$i, j = contour_mesh_sim$j, k = contour_mesh_sim$k,
            intensity = contour_mesh_sim$z,
            colorscale = list(c(0, "blue"), c(1, "gold")),
            opacity = 0.4, showscale = TRUE,
            colorbar = list(title = "Simulated Hs (m)"),
            name = "100-yr Simulated Contour") %>%
  add_trace(x = c(contour_mesh_sim$x[max_idx_s]),
            y = c(contour_mesh_sim$y[max_idx_s]),
            z = c(contour_mesh_sim$z[max_idx_s]),
            type = "scatter3d", mode = "markers+text",
            marker = list(size = 8, color = "red", symbol = "diamond"),
            text = c(max_label_s), textposition = "top center",
            name = "Max Contour Point", hoverinfo = "text") %>%
  layout(title = "3D smoothed environmental contour (simulated data)",
         scene = list(xaxis = list(title = "Current speed [m.s-1]"),
                      yaxis = list(title = "Angular difference [°]"),
                      zaxis = list(title = "Hs [m]")))
fig_sim

# -----------------------------------------------------------------------------
# 4. UNCERTAINTY QUANTIFICATION (N SIMULATIONS)
# -----------------------------------------------------------------------------
N_sims   <- 50
grid_res <- 100
set.seed(2024)

max_cspd_grid <- max(data$Cspd) * 1.2
x_grid <- seq(0, max_cspd_grid, length.out = grid_res)
y_grid <- seq(0, 360, length.out = grid_res)

z_stack <- array(0, dim = c(grid_res, grid_res, N_sims))

n_years <- 100
n_sim   <- 8766 * n_years

cat(paste("Running", N_sims, "simulations...\n"))

for (i in 1:N_sims) {
  if (i %% 10 == 0) cat(paste0(i, "/", N_sims, "...\n"))

  sim_iter  <- simulate_structured_storms(n_sim, data, best_fit, fit_gpd, fit_ald, exceedance_rate, mean_cluster_size)
  mesh_iter <- calculate_vanem_3d_contour(sim_iter, vars, return_period = 100, n_steps = 50)

  df_surf   <- data.frame(x = pmax(0, mesh_iter$x), y = mesh_iter$y, z = pmax(0, mesh_iter$z))
  pad_low   <- df_surf[df_surf$y > 300, ]; pad_low$y  <- pad_low$y - 360
  pad_high  <- df_surf[df_surf$y < 60,  ]; pad_high$y <- pad_high$y + 360
  df_interp <- rbind(df_surf, pad_low, pad_high)

  interp_res <- interp(x = df_interp$x, y = df_interp$y, z = df_interp$z,
                       xo = x_grid, yo = y_grid,
                       duplicate = "mean", linear = TRUE, extrap = FALSE)

  z_vals          <- interp_res$z
  z_vals[is.na(z_vals)] <- 0
  z_stack[, , i]  <- z_vals
}

# -----------------------------------------------------------------------------
# 5. AGGREGATE STATISTICS
# -----------------------------------------------------------------------------
z_mean  <- apply(z_stack, c(1, 2), mean,     na.rm = TRUE)
z_lower <- apply(z_stack, c(1, 2), quantile, probs = 0.025, na.rm = TRUE)
z_upper <- apply(z_stack, c(1, 2), quantile, probs = 0.975, na.rm = TRUE)
z_med   <- apply(z_stack, c(1, 2), quantile, probs = 0.50,  na.rm = TRUE)

# -----------------------------------------------------------------------------
# 6. PLOTTING FUNCTION FOR AGGREGATED SURFACES
# -----------------------------------------------------------------------------
create_smooth_mesh_plot <- function(z_matrix, title_str, color_scale,
                                    xlab = "Current speed [m.s-1]",
                                    ylab = "Angular difference [°]",
                                    zlab = "Hs [m]") {

  grid_coords     <- expand.grid(x = x_grid, y = y_grid)
  grid_coords$z   <- as.vector(z_matrix)
  pts             <- grid_coords[grid_coords$z > 0.01, ]

  pts_norm <- as.matrix(cbind(
    (pts$x - mean(pts$x)) / sd(pts$x),
    (pts$y - mean(pts$y)) / sd(pts$y),
    (pts$z - mean(pts$z)) / sd(pts$z)
  ))
  ts <- convhulln(pts_norm)

  max_idx   <- which.max(pts$z)
  max_pt    <- pts[max_idx, ]
  max_label <- paste0("Max Hs: ", round(max_pt$z, 2), " m",
                      "<br>Current speed: ", round(max_pt$x, 2), " m/s",
                      "<br>Angular difference: ", round(max_pt$y, 1), "°")

  plot_ly() %>%
    add_trace(type = "mesh3d",
              x = pts$x, y = pts$y, z = pts$z,
              i = ts[, 1] - 1, j = ts[, 2] - 1, k = ts[, 3] - 1,
              intensity = pts$z, colorscale = color_scale,
              opacity = 0.5, flatshading = FALSE, showscale = TRUE,
              colorbar = list(title = "Hs (m)"),
              name = "Contour Surface") %>%
    add_trace(x = c(max_pt$x), y = c(max_pt$y), z = c(max_pt$z),
              type = "scatter3d", mode = "markers+text",
              marker = list(size = 6, color = "red", symbol = "diamond"),
              text = c(max_label), textposition = "top center",
              name = "Max Point", hoverinfo = "text") %>%
    layout(title  = list(text = title_str, y = 0.9, yanchor = "bottom"),
           scene  = list(xaxis = list(title = xlab),
                         yaxis = list(title = ylab),
                         zaxis = list(title = zlab),
                         aspectmode = "cube"),
           margin = list(t = 50))
}

# -----------------------------------------------------------------------------
# 7. GENERATE OUTPUTS
# -----------------------------------------------------------------------------
fig_mean  <- create_smooth_mesh_plot(z_mean,  paste("Mean 100-year 3D contour (", N_sims, "simulations)"),  list(c(0, "blue"),    c(1, "cyan")))
fig_lower <- create_smooth_mesh_plot(z_lower, paste("2.5% quantile of 100-year 3D contour (", N_sims, "simulations)"),  list(c(0, "darkred"), c(1, "orange")))
fig_upper <- create_smooth_mesh_plot(z_upper, paste("97.5% quantile of 100-year 3D contour (", N_sims, "simulations)"), list(c(0, "darkred"), c(1, "orange")))

fig_mean
fig_lower
fig_upper