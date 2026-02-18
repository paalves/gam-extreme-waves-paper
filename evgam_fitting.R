  rm(list=ls())
  library(data.table)
  library(evgam)
  library(extRemes)

  # 1. Load Data
  #data <- fread("C:/Users/minim/Documents/teletravail/data2024.csv")
  data <- fread("C:/these_docs/data2024.csv")

  # Ensure covariate exists in the main dataframe for the ALD fit
  data$angularDifference <- data$angleDiff

  # 2. Fit Dynamic Threshold (ALD)
  # We use the same smooth structure on the location parameter (quantile)
  # The scale parameter is kept constant (~ 1)
  fit_ald <- evgam(
    formula = list(
      hs ~ s(Cspd, bs = "tp", k = 5) + s(angularDifference, bs = "cc", k = 5, m = 1),
      ~ 1
    ),
    data = data,
    family = "ald",
    ald.args = list(tau = 0.96)
  )

  plot(fit_ald,scheme=2)
  predict(fit_ald, newdata = data.frame(Cspd=3.01, angularDifference=102), type = "response")


  # 3. Predict the dynamic threshold
  # predict(..., type="response") returns a matrix where column 1 is the location (quantile)
  thresh_pred <- predict(fit_ald, newdata = data, type = "response")
  data$threshold <- thresh_pred[, 1]

  # 4. Decluster using the vector threshold
  # extRemes::decluster accepts a vector for the 'threshold' argument
  decl <- decluster(data$hs, threshold = data$threshold, r = 72)

  # 5. Calculate Excesses
  data$excess <- pmax(decl - data$threshold, 0)
  excess_df <- subset(data, excess > 0)

  # 6. Fit GPD on Excesses
  # Note: 'angularDifference' is already defined in excess_df from the subset
  fit_gpd <- evgam(
    formula = list(
      excess ~ s(Cspd, bs = "tp", k = 5, m = 1) +
        s(angularDifference, bs = "cc", k = 5, m=2) +
        ti(Cspd, angularDifference, bs=c("tp","cc"), k = c(5,5), m=c(1,2)),
      ~ 1
    ),
    family = "gpd",
    data = excess_df
  )

  AIC(fit_gpd) ; BIC(fit_gpd)
  summary(fit_gpd)
  plot(fit_gpd, scheme = 2, shade = FALSE)

  par(mfrow=c(1,1))
  predict(fit_gpd, newdata = excess_df, type = "qqplot")
  #predict(fit_gpd, newdata = data.frame(Cspd=3.01, angularDifference=180), type = "response")

  # --------------------------------------------------------------------------------------

library(ggplot2)
  library(MASS)
  library(fitdistrplus)
  library(KernSmooth)
  library(viridis)

  body_df <- subset(data, excess == 0)

  kde_bandwidth <- c(MASS::bandwidth.nrd(data$Cspd), 
                    MASS::bandwidth.nrd(data$angularDifference))

  fit_weibull <- fitdist(body_df$hs, "weibull")
  fit_gamma <- fitdist(body_df$hs, "gamma")
  fit_lnorm <- fitdist(body_df$hs, "lnorm")
  fit_norm <- fitdist(body_df$hs, "norm")

  comparison <- data.frame(
    Distribution = c("Weibull", "Gamma", "Lognormal", "Normal"),
    AIC = c(fit_weibull$aic, fit_gamma$aic, fit_lnorm$aic, fit_norm$aic),
    BIC = c(fit_weibull$bic, fit_gamma$bic, fit_lnorm$bic, fit_norm$bic)
  )
  print(comparison)

  best_fit <- list(fit_weibull, fit_gamma, fit_lnorm, fit_norm)[[which.min(comparison$AIC)]]
  best_name <- comparison$Distribution[which.min(comparison$AIC)]

  excess_runs <- rle(data$excess > 0)
  mean_cluster_size <- sum(data$hs > data$threshold) / sum(data$excess > 0) 
  exceedance_rate <- sum(data$excess > 0) / nrow(data)

  n_years <- 31
  n_sim <- 8760 * n_years

simulate_structured_storms <- function(n_total, data_ref, body_fit, gpd_fit, ald_fit, exc_rate, mean_dur) {
  
  # Split Data for Specific Sampling
  storm_ref <- subset(data_ref, excess > 0)
  calm_ref  <- subset(data_ref, excess == 0)
  
  # Bandwidth calculation
  bw_c <- MASS::bandwidth.nrd(data_ref$Cspd) 
  bw_a <- MASS::bandwidth.nrd(data_ref$angularDifference)
  
  # For circular variable, we need special handling
  # Use von Mises or wrapped normal, but for simplicity with existing tools:
  # We'll use a bandwidth that accounts for wrapping
  bw_c <- bw_c * 0.5
  bw_a_circular <- bw_a * 0.05 # You may want to adjust this
  
  total_storm_steps <- round(n_total * exc_rate)
  num_storms <- round(total_storm_steps) 
  
  storm_data_list <- list()
  
  # ==============================================================================
  # PART A: SIMULATE EXTREME (STORM) EVENTS - CONTINUOUS KDE SAMPLING
  # ==============================================================================
  if(num_storms > 0) {
    # Sample directly from the storm data with Gaussian kernel perturbation
    sampled_indices <- sample(1:nrow(storm_ref), size = num_storms, replace = TRUE)
    
    # Add Gaussian noise for continuous KDE
    cspd_storm <- storm_ref$Cspd[sampled_indices] + rnorm(num_storms, 0, bw_c)
    ang_storm_raw <- storm_ref$angularDifference[sampled_indices] + rnorm(num_storms, 0, bw_a_circular)
    
    # Handle circular wrapping for angular difference
    ang_storm <- ang_storm_raw %% 360
    
    # Clip to data range for Cspd
    cspd_storm <- pmax(pmin(cspd_storm, max(data_ref$Cspd)), min(data_ref$Cspd))
    
    storm_covs <- data.frame(Cspd = cspd_storm, angularDifference = ang_storm)
    
    # Predict Parameters
    storm_thresh <- predict(ald_fit, newdata = storm_covs, type = "response")[, 1]
    gpd_params <- predict(gpd_fit, newdata = storm_covs, type = "response")
    scale_vec <- gpd_params[, 1]
    shape_vec <- gpd_params[, 2]
    
    for(i in 1:num_storms) {
      dur <- max(1, rpois(1, mean_dur))
      
      # Inverse Transform Sampling for GPD
      u <- runif(1)
      if(abs(shape_vec[i]) < 1e-6) {
        peak_excess <- -scale_vec[i] * log(1 - u)
      } else {
        peak_excess <- (scale_vec[i] / shape_vec[i]) * ((1 - u)^(-shape_vec[i]) - 1)
      }
      peak_excess <- max(peak_excess, 0)
      
      # Storm Profile
      if(dur == 1) {
        profile <- 1
      } else {
        t_peak <- ceiling(dur / 2)
        rise <- seq(0.2, 1, length.out = t_peak + 1)[-1] 
        seq_fall <- seq(1, 0.2, length.out = (dur - t_peak) + 2)
        fall <- seq_fall[-c(1, length(seq_fall))]        
        
        if (length(fall) == 0 && (dur - t_peak) > 0) {
          fall <- 0.1 
        }
        profile <- c(rise, fall) 
        
        if(length(profile) != dur) {
           profile <- approx(seq(0, 1, along.with = profile), profile, n = dur)$y
        }
      }
      
      hs_storm <- storm_thresh[i] + (peak_excess * profile)
      
      # =========================================================================
      # NEW: Generate varying Cspd and angularDifference during the storm
      # =========================================================================
      # Add small perturbations around the peak values
      # Scale the perturbation by the storm intensity profile to have more variation at peak
      
      # For Cspd: add noise proportional to bandwidth and profile intensity
      cspd_variation <- cspd_storm[i] + rnorm(dur, 0, bw_c * 0.3 * (0.5 + 0.5 * profile))
      cspd_variation <- pmax(pmin(cspd_variation, max(data_ref$Cspd)), min(data_ref$Cspd))
      
      # For angularDifference: add circular noise
      ang_variation_raw <- ang_storm[i] + rnorm(dur, 0, bw_a_circular * 2 * (0.5 + 0.5 * profile))
      ang_variation <- ang_variation_raw %% 360
      
      # Recalculate threshold and GPD parameters for each time step
      varying_covs <- data.frame(Cspd = cspd_variation, angularDifference = ang_variation)
      varying_thresh <- predict(ald_fit, newdata = varying_covs, type = "response")[, 1]
      varying_gpd_params <- predict(gpd_fit, newdata = varying_covs, type = "response")
      varying_scale <- varying_gpd_params[, 1]
      varying_shape <- varying_gpd_params[, 2]
      
      storm_data_list[[i]] <- data.frame(
        Cspd = cspd_variation,
        angularDifference = ang_variation,
        hs = hs_storm,
        extreme = TRUE,
        threshold = varying_thresh,
        gpd_scale = varying_scale,
        gpd_shape = varying_shape
      )
    }
  }
  
  # Calculate remaining steps needed for calm
  simulated_storm_steps <- sum(sapply(storm_data_list, nrow))
  n_calm_needed <- n_total - simulated_storm_steps
  
  # ==============================================================================
  # PART B: SIMULATE BODY (CALM) EVENTS - CONTINUOUS KDE SAMPLING
  # ==============================================================================
  if(n_calm_needed > 0) {
    # Sample directly from calm data with Gaussian kernel perturbation
    sampled_indices_calm <- sample(1:nrow(calm_ref), size = n_calm_needed, replace = TRUE)
    
    # Add Gaussian noise for continuous KDE
    cspd_calm <- calm_ref$Cspd[sampled_indices_calm] + rnorm(n_calm_needed, 0, bw_c)
    ang_calm_raw <- calm_ref$angularDifference[sampled_indices_calm] + rnorm(n_calm_needed, 0, bw_a_circular)
    
    # Handle circular wrapping
    ang_calm <- ang_calm_raw %% 360
    
    # Clip to data range for Cspd
    cspd_calm <- pmax(pmin(cspd_calm, max(data_ref$Cspd)), min(data_ref$Cspd))
    
    calm_covs <- data.frame(Cspd = cspd_calm, angularDifference = ang_calm)
    
    calm_thresh <- predict(ald_fit, newdata = calm_covs, type = "response")[, 1]
    
    # Generate Body Values
    dist_name <- body_fit$distname
    
    if(dist_name == "weibull") {
      raw_body <- rweibull(n_calm_needed, shape = body_fit$estimate["shape"], scale = body_fit$estimate["scale"])
    } else if(dist_name == "gamma") {
      raw_body <- rgamma(n_calm_needed, shape = body_fit$estimate["shape"], rate = body_fit$estimate["rate"])
    } else if(dist_name == "lnorm") {
      raw_body <- rlnorm(n_calm_needed, meanlog = body_fit$estimate["meanlog"], sdlog = body_fit$estimate["sdlog"])
    } else {
      raw_body <- rnorm(n_calm_needed, mean = body_fit$estimate["mean"], sd = body_fit$estimate["sd"])
    }
    
    # Rejection Sampling
    invalid_indices <- which(raw_body > calm_thresh)
    max_iter <- 100
    iter <- 0
    
    while(length(invalid_indices) > 0 && iter < max_iter) {
      n_fix <- length(invalid_indices)
      if(dist_name == "weibull") {
        new_vals <- rweibull(n_fix, shape = body_fit$estimate["shape"], scale = body_fit$estimate["scale"])
      } else if(dist_name == "gamma") {
        new_vals <- rgamma(n_fix, shape = body_fit$estimate["shape"], rate = body_fit$estimate["rate"])
      } else if(dist_name == "lnorm") {
        new_vals <- rlnorm(n_fix, meanlog = body_fit$estimate["meanlog"], sdlog = body_fit$estimate["sdlog"])
      } else {
        new_vals <- rnorm(n_fix, mean = body_fit$estimate["mean"], sd = body_fit$estimate["sd"])
      }
      
      raw_body[invalid_indices] <- new_vals
      invalid_indices <- invalid_indices[raw_body[invalid_indices] > calm_thresh[invalid_indices]]
      iter <- iter + 1
    }
    
    hs_calm <- raw_body
    
    calm_df <- data.frame(
      Cspd = cspd_calm,
      angularDifference = ang_calm,
      hs = hs_calm,
      extreme = FALSE,
      threshold = calm_thresh,
      gpd_scale = NA, 
      gpd_shape = NA
    )
    
    # Combine and Interleave
    combined_list <- list()
    calm_chunks <- split(calm_df, sort(sample(1:max(1, num_storms), nrow(calm_df), replace = TRUE)))
    
    max_len <- max(length(calm_chunks), length(storm_data_list))
    for(k in 1:max_len) {
      if(k <= length(calm_chunks)) combined_list[[length(combined_list)+1]] <- calm_chunks[[k]]
      if(k <= length(storm_data_list)) combined_list[[length(combined_list)+1]] <- storm_data_list[[k]]
    }
    
    final_df <- do.call(rbind, combined_list)
    return(final_df)
    
  } else {
    return(do.call(rbind, storm_data_list))
  }
}

simulated_data <- simulate_structured_storms(n_sim, data, best_fit, fit_gpd, fit_ald, exceedance_rate, mean_cluster_size)

  par(mfrow = c(2, 2))
  hist(data$hs, prob = TRUE, breaks = 50, main = "Original Hs", xlab = "Hs", col = "lightblue")
  hist(simulated_data$hs, prob = TRUE, breaks = 50, main = "Simulated Hs", xlab = "Hs", col = "lightcoral")

  qqplot(data$hs, simulated_data$hs, main = "Q-Q Plot", xlab = "Original", ylab = "Simulated")
  abline(0, 1, col = "red")

  plot(density(data$hs), main = "Density Comparison", xlab = "Hs", lwd = 2, col = "blue")
  lines(density(simulated_data$hs), lwd = 2, col = "red")
  legend("topright", legend = c("Original", "Simulated"), col = c("blue", "red"), lwd = 2)
  par(mfrow = c(1, 1))

  print(paste("Original data - Mean:", round(mean(data$hs), 3), "SD:", round(sd(data$hs), 3)))
  print(paste("Simulated data - Mean:", round(mean(simulated_data$hs), 3), "SD:", round(sd(simulated_data$hs), 3)))
  print(paste("Original exceedance rate:", round(exceedance_rate, 4)))
  print(paste("Simulated exceedance rate:", round(mean(simulated_data$extreme), 4)))
  print(paste("Mean Cluster Duration Used:", round(mean_cluster_size, 2)))

  ggplot(simulated_data, aes(x=1:nrow(simulated_data), y = hs, color = extreme)) +
    geom_line() +
    theme_minimal() +
    labs(title = "Subset of Simulated Time Series (First 500 hours)", x = "Time", y = "Hs")


  simulated_data[which.max(simulated_data$hs),]


######################################################################
######## ENVIRONMENTAL CONTOUR #######################################
######################################################################

# Improved Huseby Environmental Contours with Artifact Removal
# Based on Huseby et al. (2015) - Structural Safety - METHOD 1

library(ggplot2)
library(viridis)

#' Calculate C(theta) function for Huseby's method
#' 
#' @param data Data frame with environmental variables
#' @param theta Angle in radians
#' @param Pe Exceedance probability
#' @param var1 Name of first variable (default: "Cspd")
#' @param var2 Name of second variable (default: "hs")
#' @return Value of C(theta)
calculate_C_theta <- function(data, theta, Pe, var1 = "Cspd", var2 = "hs") {
  # Calculate projections Y = X1*cos(theta) + X2*sin(theta)
  Y <- data[[var1]] * cos(theta) + data[[var2]] * sin(theta)
  
  # Find k such that (n-k)/n ≈ Pe
  n <- length(Y)
  k <- ceiling(n * (1 - Pe))  # Use ceiling to be conservative
  k <- max(1, min(k, n))  # Ensure k is in valid range
  
  # C(theta) is the k-th order statistic
  C_theta <- sort(Y, partial = k)[k]  # More efficient partial sort
  
  return(C_theta)
}


#' Check if hyperplane supports the contour (Theorem 2.6)
#' 
#' Implements condition: 0.5*[C(theta-delta) + C(theta+delta)] > cos(delta)*C(theta)
#' 
#' @param C_values Vector of C(theta) values
#' @param idx Index to check
#' @param delta Angular spacing
#' @param tolerance Numerical tolerance for comparison
#' @return TRUE if hyperplane supports the contour
check_supporting_hyperplane <- function(C_values, idx, delta, tolerance = 1e-6) {
  n <- length(C_values)
  
  # Get indices with wraparound
  idx_minus <- ((idx - 2) %% n) + 1
  idx_plus <- (idx %% n) + 1
  
  # Check condition from Theorem 2.6 with tolerance
  lhs <- 0.5 * (C_values[idx_minus] + C_values[idx_plus])
  rhs <- cos(delta) * C_values[idx]
  
  # Add small tolerance to handle numerical precision
  return(lhs > (rhs - tolerance))
}


#' Smooth C values using local regression
#' 
#' @param angles Vector of angles
#' @param C_values Vector of C(theta) values
#' @param span Smoothing parameter (0-1)
#' @return Smoothed C values
smooth_C_values <- function(angles, C_values, span = 0.1) {
  # Use loess for smoothing but preserve overall shape
  # Duplicate data at boundaries for circular smoothing
  angles_ext <- c(angles - 2*pi, angles, angles + 2*pi)
  C_ext <- rep(C_values, 3)
  
  fit <- loess(C_ext ~ angles_ext, span = span, degree = 1)
  smoothed <- predict(fit, newdata = angles)
  
  return(smoothed)
}


#' Huseby Method 1: Improved with robust artifact removal
#' 
#' @param data Data frame with environmental variables
#' @param Pe Exceedance probability
#' @param n_angles Number of angles to use
#' @param var1 Name of first variable
#' @param var2 Name of second variable
#' @param remove_artifacts Whether to remove non-supporting hyperplanes
#' @param smooth_C Whether to smooth C function before artifact removal
#' @param adaptive_angles Use adaptive angle spacing in problematic regions
#' @return Data frame with contour points
huseby_method1_improved <- function(data, Pe, n_angles = 180, 
                                    var1 = "Cspd", var2 = "hs",
                                    remove_artifacts = TRUE,
                                    smooth_C = TRUE,
                                    adaptive_angles = FALSE) {
  
  # Generate evenly spaced angles
  angles <- seq(0, 2*pi, length.out = n_angles + 1)[1:n_angles]
  delta <- angles[2] - angles[1]
  
  # Calculate C(theta) for each angle
  C_values <- sapply(angles, function(theta) {
    calculate_C_theta(data, theta, Pe, var1, var2)
  })
  
  # Optional: Smooth C values before artifact detection
  if (smooth_C) {
    C_values_smooth <- smooth_C_values(angles, C_values, span = 0.05)
    # Keep original for comparison but use smoothed for artifact detection
    C_values <- C_values_smooth
  }
  
  # Identify supporting hyperplanes with improved logic
  if (remove_artifacts) {
    # Check each hyperplane multiple times with different deltas
    is_supporting <- rep(TRUE, n_angles)
    
    # Test with multiple delta values for robustness
    test_deltas <- c(delta, 2*delta, 3*delta)
    
    for (test_delta in test_deltas) {
      n_test <- floor(2*pi / test_delta)
      if (n_test < 3) next
      
      test_angles <- seq(0, 2*pi, length.out = n_test + 1)[1:n_test]
      test_C <- approx(angles, C_values, test_angles, rule = 2)$y
      
      for (i in 1:n_test) {
        passes <- check_supporting_hyperplane(test_C, i, test_delta)
        if (!passes) {
          # Find corresponding original angles and mark as non-supporting
          angle_val <- test_angles[i]
          nearest_idx <- which.min(abs(angles - angle_val))
          is_supporting[nearest_idx] <- FALSE
        }
      }
    }
    
    # Additional check: remove isolated non-supporting points
    # (likely numerical artifacts)
    for (i in 1:n_angles) {
      idx_prev <- ((i - 2) %% n_angles) + 1
      idx_next <- (i %% n_angles) + 1
      
      # If neighbors are supporting but this isn't, check more carefully
      if (!is_supporting[i] && is_supporting[idx_prev] && is_supporting[idx_next]) {
        # Recheck with tighter tolerance
        if (check_supporting_hyperplane(C_values, i, delta, tolerance = 1e-5)) {
          is_supporting[i] <- TRUE
        }
      }
    }
    
    # Keep only supporting hyperplanes
    valid_indices <- which(is_supporting)
    
    if (length(valid_indices) < 4) {
      warning("Too few supporting hyperplanes found (", length(valid_indices), 
              "). Using all points.")
      valid_indices <- 1:n_angles
    } else {
      cat(sprintf("Removed %d non-supporting hyperplanes (%.1f%%)\n", 
                  n_angles - length(valid_indices),
                  100 * (n_angles - length(valid_indices)) / n_angles))
    }
  } else {
    valid_indices <- 1:n_angles
  }
  
  # Calculate intersection points only for valid hyperplanes
  n_valid <- length(valid_indices)
  contour_points <- matrix(0, nrow = n_valid, ncol = 2)
  
  for (i in 1:n_valid) {
    j <- valid_indices[i]
    j_next <- valid_indices[(i %% n_valid) + 1]
    
    theta_j <- angles[j]
    theta_jp1 <- angles[j_next]
    C_j <- C_values[j]
    C_jp1 <- C_values[j_next]
    
    # Calculate delta (handle wraparound)
    delta_theta <- theta_jp1 - theta_j
    if (delta_theta < 0) {
      delta_theta <- delta_theta + 2*pi
    }
    
    # Solve for intersection point
    sin_delta <- sin(delta_theta)
    
    if (abs(sin_delta) > 1e-10) {
      x1 <- (sin(theta_jp1) * C_j - sin(theta_j) * C_jp1) / sin_delta
      x2 <- (-cos(theta_jp1) * C_j + cos(theta_j) * C_jp1) / sin_delta
    } else {
      # Handle parallel lines (use midpoint)
      x1 <- 0.5 * (C_j * cos(theta_j) + C_jp1 * cos(theta_jp1))
      x2 <- 0.5 * (C_j * sin(theta_j) + C_jp1 * sin(theta_jp1))
    }
    
    contour_points[i, 1] <- x1
    contour_points[i, 2] <- x2
  }
  
  # Post-process: remove any obvious outliers
  if (remove_artifacts && n_valid > 10) {
    # Calculate distances from centroid
    centroid <- colMeans(contour_points)
    distances <- sqrt(rowSums((contour_points - matrix(centroid, n_valid, 2, byrow = TRUE))^2))
    
    # Remove points that are much farther than median (likely artifacts)
    median_dist <- median(distances)
    outlier_threshold <- median_dist * 3
    keep <- distances < outlier_threshold
    
    if (sum(!keep) > 0 && sum(keep) > 10) {
      cat(sprintf("Removed %d outlier points in post-processing\n", sum(!keep)))
      contour_points <- contour_points[keep, , drop = FALSE]
    }
  }
  
  # Create data frame with closed contour
  contour_df <- data.frame(
    var1 = c(contour_points[, 1], contour_points[1, 1]),
    var2 = c(contour_points[, 2], contour_points[1, 2])
  )
  names(contour_df)[1:2] <- c(var1, var2)
  
  return(contour_df)
}


#' Calculate contours for multiple return periods
#' 
#' @param data Data frame with environmental variables
#' @param return_periods Vector of return periods in years
#' @param n_angles Number of angles for contour calculation
#' @param var1 Name of first variable
#' @param var2 Name of second variable
#' @param remove_artifacts Whether to remove artifacts
#' @return List with all contours and combined data frame
calculate_environmental_contours <- function(data, 
                                             return_periods = c(1, 2, 5, 10, 20, 50, 100),
                                             n_angles = 180,
                                             var1 = "Cspd",
                                             var2 = "hs",
                                             remove_artifacts = TRUE) {
  
  # For hourly data: Pe = 1/(N*365.25*24)
  Pe_values <- 1 / (return_periods * 365.25 * 24)
  
  cat("=== Calculating Environmental Contours ===\n")
  cat(sprintf("Data points: %d\n", nrow(data)))
  cat(sprintf("Return periods: %s years\n", paste(return_periods, collapse=", ")))
  cat(sprintf("Angles: %d\n", n_angles))
  cat(sprintf("Remove artifacts: %s\n\n", remove_artifacts))
  
  # Calculate contours
  contours_list <- list()
  
  for (i in 1:length(return_periods)) {
    rp <- return_periods[i]
    Pe <- Pe_values[i]
    
    cat(sprintf("Calculating %3d-year contour (Pe = %.6e)... ", rp, Pe))
    
    contour <- huseby_method1_improved(
      data, 
      Pe, 
      n_angles = n_angles,
      var1 = var1, 
      var2 = var2,
      remove_artifacts = remove_artifacts,
      smooth_C = TRUE
    )
    
    contour$return_period <- rp
    contours_list[[paste0("RP", rp)]] <- contour
    
    max_val <- max(contour[[var2]], na.rm = TRUE)
    max_idx <- which.max(contour[[var2]])
    cat(sprintf("max %s = %.2f at %s = %.2f\n", 
                var2, max_val, var1, contour[[var1]][max_idx]))
  }
  
  # Combine all contours
  all_contours <- do.call(rbind, contours_list)
  
  return(list(
    contours_list = contours_list,
    all_contours = all_contours,
    data = data,
    return_periods = return_periods
  ))
}


#' Plot environmental contours with data
#' 
#' @param contour_result Result from calculate_environmental_contours
#' @param sample_size Number of data points to plot (NULL for all)
#' @param point_alpha Transparency of data points
#' @param var1_label Label for x-axis
#' @param var2_label Label for y-axis
#' @return ggplot object
plot_environmental_contours <- function(contour_result, 
                                        sample_size = nrow(data),
                                        point_alpha = 1,
                                        var1_label = "Current Speed (m/s)",
                                        var2_label = "Significant Wave Height (m)") {
  
  data <- contour_result$data
  all_contours <- contour_result$all_contours
  return_periods <- contour_result$return_periods
  
  # Sample data points for plotting
  if (!is.null(sample_size) && nrow(data) > sample_size) {
    set.seed(42)
    plot_indices <- sample(1:nrow(data), sample_size)
    data_plot <- data[plot_indices, ]
  } else {
    data_plot <- data
  }
  
  # Get variable names from contours
  var1 <- names(all_contours)[1]
  var2 <- names(all_contours)[2]
  
  # Create the plot
  p <- ggplot() +
    geom_point(data = data_plot, 
               aes(x = .data[[var1]], y = .data[[var2]]), 
               color = "black", 
               size = 0.5, 
               alpha = point_alpha) +
    geom_path(data = all_contours, 
              aes(x = .data[[var1]], y = .data[[var2]], 
                  color = as.factor(return_period), 
                  group = return_period),
              linewidth = 1.2) +
    scale_color_viridis_d(option = "plasma",
                          name = "Return Period (years)",
                          breaks = return_periods,
                          labels = return_periods) +
    labs(title = "Environmental Contours - Huseby Method 1 (Improved)",
         subtitle = sprintf("Return periods: %s years | Artifacts removed", 
                           paste(return_periods, collapse=", ")),
         x = var1_label,
         y = var2_label) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = "gray80"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90")
    ) +
    coord_cartesian(xlim = range(data_plot[[var1]]), 
                    ylim = range(data_plot[[var2]]))
  
  return(p)
}


#' Print summary of contours
#' 
#' @param contour_result Result from calculate_environmental_contours
print_contour_summary <- function(contour_result) {
  contours_list <- contour_result$contours_list
  return_periods <- contour_result$return_periods
  var1 <- names(contours_list[[1]])[1]
  var2 <- names(contours_list[[1]])[2]
  
  cat("\n=== Contour Summary ===\n")
  for (rp in return_periods) {
    contour <- contours_list[[paste0("RP", rp)]]
    cat(sprintf("\n%d-year return period:\n", rp))
    cat(sprintf("  Points: %d\n", nrow(contour) - 1))
    cat(sprintf("  %s range: [%.2f, %.2f]\n", 
                var1, min(contour[[var1]]), max(contour[[var1]])))
    cat(sprintf("  %s range: [%.2f, %.2f]\n", 
                var2, min(contour[[var2]]), max(contour[[var2]])))
    max_idx <- which.max(contour[[var2]])
    cat(sprintf("  Maximum: %s = %.2f at %s = %.2f\n", 
                var2, contour[[var2]][max_idx], var1, contour[[var1]][max_idx]))
  }
  cat("\n=== COMPLETE ===\n")
}

result <- calculate_environmental_contours(
  data, 
  return_periods = c(1, 2, 5, 10, 20, 50, 100),
  n_angles = 180,
  var1 = "tp", # Set x-variable
  var2 = "hs",                # Set y-variable
  remove_artifacts = TRUE
)

result$all_contours$hs <- pmax(0, result$all_contours$hs)

p <- plot_environmental_contours(
  result, 
  var1_label = "Tp (s)", 
  var2_label = "Hs (m)"
)
print(p)

print_contour_summary(result)


######################################################################
######## 3D ENVIRONMENTAL CONTOURS (VANEM, 2019) #####################
######################################################################

library(plotly)
library(geometry) # Required for convhulln

# 1. Implementation of the Vanem (2019) 3D Direct Sampling Method
calculate_vanem_3d_contour <- function(data, var_names, return_period, n_steps = 50) {
  
  # --- [Setup and C-Calculation remains identical to previous version] ---
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
  d1 <- t1_seq[2] - t1_seq[1]
  d2 <- t2_seq[2] - t2_seq[1]
  
  C_mat <- matrix(0, nrow = length(t1_seq), ncol = length(t2_seq))
  
  cat(paste("Calculating projections for", return_period, "yr contour...\n"))
  for (j in 1:length(t2_seq)) {
    t2 <- t2_seq[j]; sin_t2 <- sin(t2); cos_t2 <- cos(t2)
    for (i in 1:length(t1_seq)) {
      t1 <- t1_seq[i]
      R <- X*cos(t1)*cos_t2 + Y*sin(t1)*cos_t2 + Z*sin_t2
      C_mat[i, j] <- -sort(-R, partial = n_samples - k + 1)[n_samples - k + 1]
    }
  }
  
  # --- Smoothing C-matrix (Essential for clean hull) ---
  C_pad <- rbind(C_mat[nrow(C_mat),], C_mat, C_mat[1,])
  C_pad <- cbind(C_pad[,ncol(C_pad)], C_pad, C_pad[,1])
  C_smooth <- C_mat
  for(i in 1:nrow(C_mat)) {
    for(j in 1:ncol(C_mat)) {
       C_smooth[i,j] <- mean(C_pad[i:(i+2), j:(j+2)])
    }
  }
  C_mat <- C_smooth
  
  # --- Solve Intersections ---
  # Store as list of points instead of matrix grid
  points_list <- list()
  counter <- 1
  
  for (j in 1:n_steps) {
    t2 <- t2_seq[j]
    if (abs(cos(t2)) < 1e-4) next 
    
    for (i in 1:n_steps) {
      t1 <- t1_seq[i]
      
      idx_i <- i; idx_i_next <- if(i == n_steps) 1 else i + 1
      idx_j <- j; idx_j_next <- if(j == n_steps) 1 else j + 1
      
      C1 <- C_mat[idx_i, idx_j]; C2 <- C_mat[idx_i_next, idx_j]; C3 <- C_mat[idx_i, idx_j_next]
      
      # Denominator [cite: 141]
      D_val <- (cos(t1)*sin(t1+d1) - sin(t1)*cos(t1+d1)) * cos(t2) * (cos(t2)*sin(t2+d2) - sin(t2)*cos(t2+d2))
      
      if (abs(D_val) < 1e-5) next # Skip singularities
      
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
  
  # Combine to matrix
  pts_norm <- do.call(rbind, points_list)
  
  # --- Filter Outliers (Fixes the "Artefacts") ---
  # Remove points that are extremely far from the center (numerical explosions)
  dist_sq <- rowSums(pts_norm^2)
  threshold <- median(dist_sq) * 10 # Robust outlier detection
  valid_idx <- which(dist_sq < threshold)
  pts_norm <- pts_norm[valid_idx, ]
  
  # Destandardize
  final_pts <- pts_norm
  final_pts[,1] <- final_pts[,1] * sd_X + mean_X
  final_pts[,2] <- final_pts[,2] * sd_Y + mean_Y
  final_pts[,3] <- final_pts[,3] * sd_Z + mean_Z
  
  # --- Convex Hull Construction ---
  # This creates a clean triangular mesh enclosing the points
  ts <- convhulln(final_pts)
  
  return(list(
    x = final_pts[,1], 
    y = final_pts[,2], 
    z = final_pts[,3],
    i = ts[,1]-1, # Plotly uses 0-based indexing
    j = ts[,2]-1,
    k = ts[,3]-1
  ))
}

##################################################################
########### PLOTTING FOR THE ORIGINAL DATA  ###################
##################################################################

# 1. Calculate
vars <- c("Cspd", "angularDifference", "hs")
# We can use fewer steps now because the hull interpolates well
contour_mesh <- calculate_vanem_3d_contour(data, vars, return_period = 100, n_steps = 50)

# Clip negative Z (physical constraint)
contour_mesh$z <- pmax(0, contour_mesh$z)

# 2. Plot
fig <- plot_ly()

# Add Data
sub_idx <- sample(1:nrow(data), nrow(data))
fig <- fig %>% add_trace(
  x = data$Cspd[sub_idx], 
  y = data$angularDifference[sub_idx], 
  z = data$hs[sub_idx],
  type = "scatter3d", mode = "markers",
  marker = list(size = 2, opacity = 0.8, color="black"),
  name = "Observations"
)

contour_mesh$x <- pmax(0, contour_mesh$x)
contour_mesh$y <- pmin(360, contour_mesh$y)
contour_mesh$y <- pmax(0, contour_mesh$y)
contour_mesh$z <- pmax(0, contour_mesh$z)

# Add Convex Hull Mesh
fig <- fig %>% add_trace(
  type = "mesh3d",
  x = contour_mesh$x,
  y = contour_mesh$y,
  z = contour_mesh$z,
  i = contour_mesh$i,
  j = contour_mesh$j,
  k = contour_mesh$k,
  intensity = contour_mesh$z,
  colorscale = list(c(0, "green"), c(1, "red")),
  opacity = 0.4,
  showscale = TRUE,
  colorbar = list(title = "Hs (m)"),
  name = "100-yr Contour (Convex Hull)"
)

# --- UPDATED LAYOUT SECTION ---
fig <- fig %>% layout(
  title = list(
    text = "3D smoothed environmental contour",
    y = 0.9,             # Moves title lower (0.9 = 90% height)
    yanchor = "bottom"   # Anchors text bottom to the Y position
  ),
  scene = list(
    xaxis = list(title = "Current speed [m.s-1]"),
    yaxis = list(title = "Angular difference [°]"),
    zaxis = list(title = "Hs [m]")
  ),
  margin = list(t = 50)  # Reduces the top padding margin
)

# 1. Identify the index of the maximum Hs value within the contour mesh
max_idx <- which.max(contour_mesh$z)

# 2. Extract the coordinates
max_hs   <- contour_mesh$z[max_idx]
max_cspd <- contour_mesh$x[max_idx]
max_ang  <- contour_mesh$y[max_idx]

# 3. Create a label string for the hover text/annotation
max_label <- paste0("Hs: ", round(max_hs, 2), " m",
                    "<br>Current speed: ", round(max_cspd, 2), " m.s-1",
                    "<br>Angular difference: ", round(max_ang, 1), "°")

# 4. Add the red point to your existing fig
fig <- fig %>% add_trace(
  x = c(max_cspd),
  y = c(max_ang),
  z = c(max_hs),
  type = "scatter3d",
  mode = "markers+text",
  marker = list(size = 8, color = "red", symbol = "diamond"),
  text = c(max_label),
  textposition = "top center",
  name = "Max Hs contour point",
  hoverinfo = "text"
)

# Display the updated figure
fig


##################################################################
########### PLOTTING FOR THE SIMULATED DATA  ###################
##################################################################

# 1. Define variables and calculate the 3D contour using Simulated Data
# Note: Ensure the variable names match the columns in your simulated_data object
vars_sim <- c("Cspd", "angularDifference", "hs")

# Calculate the mesh using the Vanem (2019) method on simulated points
contour_mesh_sim <- calculate_vanem_3d_contour(
  data = simulated_data, 
  var_names = vars_sim, 
  return_period = 100, 
  n_steps = 50
)

# Apply physical constraints (Clips negative values and circular bounds)
contour_mesh_sim$x <- pmax(0, contour_mesh_sim$x)
contour_mesh_sim$y <- pmax(0, pmin(360, contour_mesh_sim$y))
contour_mesh_sim$z <- pmax(0, contour_mesh_sim$z)

# 2. Build the Plotly Visualization
fig_sim <- plot_ly()

fig_sim <- fig_sim %>% add_trace(
  x = simulated_data$Cspd, 
  y = simulated_data$angularDifference, 
  z = simulated_data$hs,
  type = "scatter3d", 
  mode = "markers",
  marker = list(size = 2, opacity = 0.8, color = "black"),
  name = "Simulated data"
)

# Add the Convex Hull Mesh based on Simulated Extremes
fig_sim <- fig_sim %>% add_trace(
  type = "mesh3d",
  x = contour_mesh_sim$x,
  y = contour_mesh_sim$y,
  z = contour_mesh_sim$z,
  i = contour_mesh_sim$i,
  j = contour_mesh_sim$j,
  k = contour_mesh_sim$k,
  intensity = contour_mesh_sim$z, 
  colorscale = list(c(0, "blue"), c(1, "gold")), # Different color to distinguish from original
  opacity = 0.4,
  showscale = TRUE,
  colorbar = list(title = "Simulated Hs (m)"),
  name = "100-yr Simulated Contour"
)

# Layout and Labels
fig_sim <- fig_sim %>% layout(
  title = "3D smoothed environmental contour (simulated data)",
  scene = list(
    xaxis = list(title = "Current speed [m.s-1]"),
    yaxis = list(title = "Angular difference [°]"),
    zaxis = list(title = "Hs [m])")
  )
)

# 1. Identify the index of the maximum Hs value within the contour mesh
max_idx <- which.max(contour_mesh_sim$z)

# 2. Extract the coordinates
max_hs   <- contour_mesh_sim$z[max_idx]
max_cspd <- contour_mesh_sim$x[max_idx]
max_ang  <- contour_mesh_sim$y[max_idx]

# 3. Create a label string for the hover text/annotation
max_label <- paste0("Hs: ", round(max_hs, 2), " m",
                    "<br>Current speed: ", round(max_cspd, 2), " m/s",
                    "<br>Angular difference: ", round(max_ang, 1), "°")

# 4. Add the red point to your existing fig_sim
fig_sim <- fig_sim %>% add_trace(
  x = c(max_cspd),
  y = c(max_ang),
  z = c(max_hs),
  type = "scatter3d",
  mode = "markers+text",
  marker = list(size = 8, color = "red", symbol = "diamond"),
  text = c(max_label),
  textposition = "top center",
  name = "Max Contour Point",
  hoverinfo = "text"
)

# Display the updated figure
fig_sim

#################################################################################
#################################################################################


library(plotly)
library(akima)    # For gridding/interpolation
library(geometry) # For convex hull (convhulln)
library(abind)    # For array handling

# -----------------------------------------------------------------------------
# 1. SETUP SIMULATION PARAMETERS
# -----------------------------------------------------------------------------
N_sims <- 50          # Number of simulations (Increase to 100+ for smoother quantiles)
grid_res <- 75       # High resolution for smoother interpolation (100x100)

# Define a fixed grid to aggregate all simulations onto
# We add a small buffer to the max Cspd to ensure we capture the full tail
max_cspd_grid <- max(data$Cspd) * 1.2
x_grid <- seq(0, max_cspd_grid, length.out = grid_res)
y_grid <- seq(0, 360, length.out = grid_res)

# Array to store the "Height" (Hs) surfaces: [x, y, simulation_index]
z_stack <- array(0, dim = c(grid_res, grid_res, N_sims))

print(paste("Running", N_sims, "simulations..."))

# -----------------------------------------------------------------------------
# 2. SIMULATION LOOP
# -----------------------------------------------------------------------------
set.seed(2024) 
n_years <- 100
n_sim <- 8766 * n_years

for (i in 1:N_sims) {
  if(i %% 10 == 0) cat(paste0(i, "/", N_sims, "...\n"))
  
  # A. Generate Synthetic Data
  sim_iter <- simulate_structured_storms(
    n_sim, data, best_fit, fit_gpd, fit_ald, exceedance_rate, mean_cluster_size
  )
  
  # B. Calculate 3D Contour Points (Vertices only)
  vars_iter <- c("Cspd", "angularDifference", "hs")
  
  # We reuse your existing function logic but just need the points, not the plot indices
  mesh_iter <- calculate_vanem_3d_contour(
    data = sim_iter, 
    var_names = vars_iter, 
    return_period = 100, 
    n_steps = 50
  )
  
  # C. Prepare Data for Gridding
  df_surf <- data.frame(x = mesh_iter$x, y = mesh_iter$y, z = mesh_iter$z)
  
  # Physical constraints
  df_surf$x <- pmax(0, df_surf$x)
  df_surf$z <- pmax(0, df_surf$z)
  
  # Circular Padding: Duplicate points near boundaries to ensure smooth interpolation across 0/360
  pad_low  <- df_surf[df_surf$y > 300, ]; pad_low$y  <- pad_low$y - 360
  pad_high <- df_surf[df_surf$y < 60,  ]; pad_high$y <- pad_high$y + 360
  df_interp_input <- rbind(df_surf, pad_low, pad_high)
  
  # D. Interpolate onto the fixed common grid
  # interp() fills the gaps between the specific simulation points
  interp_res <- interp(
    x = df_interp_input$x,
    y = df_interp_input$y,
    z = df_interp_input$z,
    xo = x_grid,
    yo = y_grid,
    duplicate = "mean",
    linear = TRUE,
    extrap = FALSE
  )
  
  # Store. NAs (outside the hull) are treated as 0 height for aggregation
  z_vals <- interp_res$z
  z_vals[is.na(z_vals)] <- 0
  z_stack[, , i] <- z_vals
}

# -----------------------------------------------------------------------------
# 3. STATISTICAL AGGREGATION
# -----------------------------------------------------------------------------
# Calculate the surfaces element-wise
z_mean  <- apply(z_stack, c(1, 2), mean, na.rm = TRUE)
z_lower <- apply(z_stack, c(1, 2), quantile, probs = 0.025, na.rm = TRUE)
z_upper <- apply(z_stack, c(1, 2), quantile, probs = 0.975, na.rm = TRUE)

# -----------------------------------------------------------------------------
# 4. PLOTTING FUNCTION (Using Convex Hull Mesh)
# -----------------------------------------------------------------------------
# This function converts the grid back into a mesh3d object exactly like fig_sim

create_smooth_mesh_plot <- function(z_matrix, title_str, color_scale, xlab="Current speed [m.s-1]", ylab="Angular difference [°]", zlab="Hs [m]") {
  
  # 1. Convert Grid Matrix back to Point Cloud (XYZ)
  grid_coords <- expand.grid(x = x_grid, y = y_grid)
  grid_coords$z <- as.vector(z_matrix)
  
  pts <- grid_coords[grid_coords$z > 0.01, ]
  
  # 2. Re-calculate Convex Hull
  mx <- mean(pts$x); sx <- sd(pts$x)
  my <- mean(pts$y); sy <- sd(pts$y)
  mz <- mean(pts$z); sz <- sd(pts$z)
  
  pts_norm <- as.matrix(cbind(
    (pts$x - mx)/sx,
    (pts$y - my)/sy,
    (pts$z - mz)/sz
  ))
  
  ts <- convhulln(pts_norm)
  
  # 3. Identify Max Point
  max_idx <- which.max(pts$z)
  max_pt <- pts[max_idx, ]
  max_label <- paste0("Max Hs: ", round(max_pt$z, 2), " m",
                      "<br>Current speed: ", round(max_pt$x, 2), " m/s",
                      "<br>Angular difference: ", round(max_pt$y, 1), "°")
  
  # 4. Build Plot
  fig <- plot_ly()
  
  fig <- fig %>% add_trace(
    type = "mesh3d",
    x = pts$x,
    y = pts$y,
    z = pts$z,
    i = ts[,1]-1, 
    j = ts[,2]-1,
    k = ts[,3]-1,
    intensity = pts$z,
    colorscale = color_scale,
    opacity = 0.5, 
    flatshading = FALSE, 
    showscale = TRUE,
    colorbar = list(title = "Hs (m)"),
    name = "Contour Surface"
  )
  
  fig <- fig %>% add_trace(
    x = c(max_pt$x),
    y = c(max_pt$y),
    z = c(max_pt$z),
    type = "scatter3d",
    mode = "markers+text",
    marker = list(size = 6, color = "red", symbol = "diamond"),
    text = c(max_label),
    textposition = "top center",
    name = "Max Point",
    hoverinfo = "text"
  )
  
  # Formatting - UPDATED TITLE POSITION
  fig <- fig %>% layout(
    title = list(
      text = title_str,
      y = 0.9,            # Moves title down (1.0 is top, 0 is bottom)
      yanchor = "bottom"  # Anchors the bottom of the text to the y position
    ),
    scene = list(
      xaxis = list(title = xlab),
      yaxis = list(title = ylab),
      zaxis = list(title = zlab),
      aspectmode = "cube" 
    ),
    margin = list(t = 50) # Reduces top margin whitespace
  )
  
  return(fig)
}

# -----------------------------------------------------------------------------
# 5. GENERATE OUTPUTS
# -----------------------------------------------------------------------------

# Plot 1: Mean Contour
fig_mean_smooth <- create_smooth_mesh_plot(
  z_mean, 
  "Mean 100 year 3D contour (100 simulations)",
  list(c(0, "blue"), c(1, "cyan"))
)
fig_mean_smooth

# Plot 2: Lower Bound (2.5%)
fig_lower_smooth <- create_smooth_mesh_plot(
  z_lower, 
  "2.5% quantile of 100 year 3D contour (100 simulations)", 
  list(c(0, "darkred"), c(1, "orange"))
)
fig_lower_smooth

# Plot 3: Upper Bound (97.5%)
fig_upper_smooth <- create_smooth_mesh_plot(
  z_upper, 
  "97.5% quantile of 100 year 3D contour (100 simulations)", 
  list(c(0, "darkred"), c(1, "orange"))
)
fig_upper_smooth

