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


# ===============================================================================
# QUICK TEST VERSION - Return Value Analysis
# ===============================================================================
# 
# This is a reduced version for quick testing:
# - N = 10 simulations (instead of 50)
# - T = 20 years (instead of 100)
# - Coarser covariate grid
# 
# Runtime: ~2-5 minutes
# ===============================================================================

library(data.table)
library(ggplot2)
library(viridis)

# -------------------------------------------------------------------------------
# CONFIGURATION
# -------------------------------------------------------------------------------

N_simulations <- 15      # Reduced from 50
T_years <- 100            # Reduced from 100
n_sim_total <- 8760 * T_years

set.seed(456)  # Different seed than full analysis

# -------------------------------------------------------------------------------
# RUN SIMULATIONS
# -------------------------------------------------------------------------------

all_simulations_quick <- vector("list", N_simulations)

cat("QUICK TEST MODE\n")
cat("===============\n")
cat("Running", N_simulations, "simulations of", T_years, "years each...\n\n")

for (i in 1:N_simulations) {
  cat("  Simulation", i, "of", N_simulations, "\n")
  
  all_simulations_quick[[i]] <- simulate_structured_storms(
    n_total = n_sim_total,
    data_ref = data,
    body_fit = best_fit,
    gpd_fit = fit_gpd,
    ald_fit = fit_ald,
    exc_rate = exceedance_rate,
    mean_dur = mean_cluster_size
  )
  
  all_simulations_quick[[i]]$sim_id <- i
}

cat("\nSimulations complete!\n\n")

# -------------------------------------------------------------------------------
# COMPUTE RETURN VALUES
# -------------------------------------------------------------------------------

return_periods <- c(1, 5, 10, 20, 50, 100)  # Reduced set

# Extract annual maxima
annual_maxima_quick <- lapply(1:N_simulations, function(i) {
  sim_data <- all_simulations_quick[[i]]
  sim_data$year <- ceiling(1:nrow(sim_data) / 8760)
  annual_max <- aggregate(hs ~ year, data = sim_data, FUN = max)
  annual_max$sim_id <- i
  return(annual_max)
})

all_annual_maxima_quick <- do.call(rbind, annual_maxima_quick)

# Compute return values
return_values_matrix_quick <- matrix(NA, nrow = N_simulations, ncol = length(return_periods))
colnames(return_values_matrix_quick) <- paste0("RP_", return_periods)

for (i in 1:N_simulations) {
  am_sim <- all_annual_maxima_quick$hs[all_annual_maxima_quick$sim_id == i]
  
  for (j in 1:length(return_periods)) {
    rp <- return_periods[j]
    prob <- 1 - 1/rp
    return_values_matrix_quick[i, j] <- quantile(am_sim, probs = prob, type = 8)
  }
}

# Statistics
return_value_stats_quick <- data.frame(
  ReturnPeriod = return_periods,
  Median = apply(return_values_matrix_quick, 2, median),
  Mean = apply(return_values_matrix_quick, 2, mean),
  CI_lower = apply(return_values_matrix_quick, 2, quantile, probs = 0.025),
  CI_upper = apply(return_values_matrix_quick, 2, quantile, probs = 0.975),
  SD = apply(return_values_matrix_quick, 2, sd)
)

cat("\nQUICK TEST - RETURN VALUE ESTIMATES:\n")
cat("=" , rep("=", 70), "\n", sep = "")
print(return_value_stats_quick, row.names = FALSE, digits = 3)
cat("\n")
cat("NOTE: These are quick estimates with limited simulations.\n")
cat("For final results, run the full analysis with N=50, T=100.\n\n")

# -------------------------------------------------------------------------------
# PRETTY VISUALIZATION (GGPLOT2) - WITH RIBBON
# -------------------------------------------------------------------------------
library(ggplot2)
library(data.table)

# 1. Prepare Data
# ----------------------------

# Convert individual simulations matrix to long format
sim_dt <- as.data.table(return_values_matrix_quick)
sim_dt[, SimID := .I]
sim_long <- melt(sim_dt, id.vars = "SimID", variable.name = "RP_Label", value.name = "Hs")
sim_long[, ReturnPeriod := as.numeric(gsub("RP_", "", RP_Label))]

# Ensure summary stats are a data.frame
stats_df <- as.data.frame(return_value_stats_quick)

# 2. Create the Plot
# ------------------
p <- ggplot() +
  # A. Individual Simulations (Background "Spaghetti")
  geom_line(data = sim_long, 
            aes(x = ReturnPeriod, y = Hs, group = SimID), 
            color = "gray85", size = 0.4, alpha = 0.6) +
  
  # B. Confidence Interval Ribbon (The change you requested)
  geom_ribbon(data = stats_df, 
              aes(x = ReturnPeriod, ymin = CI_lower, ymax = CI_upper),
              fill = "#C0392B", alpha = 0.2) + # Light red transparent fill
  
  # C. Median Line
  geom_line(data = stats_df, 
            aes(x = ReturnPeriod, y = Median), 
            color = "#2980B9", size = 1.2) + 
  
  # D. Median Points
  geom_point(data = stats_df, 
             aes(x = ReturnPeriod, y = Median), 
             color = "#2980B9", fill = "white", size = 3, shape = 21, stroke = 1.5) +
  
  # E. Text Labels
  # 1. Median Labels
  geom_text(data = stats_df,
            aes(x = ReturnPeriod, y = Median, label = sprintf("%.2f", Median)),
            vjust = -1.2, color = "#2980B9", fontface = "bold", size = 3.5) +
  
  # 2. Upper CI Labels (Sitting on top of ribbon)
  geom_text(data = stats_df,
            aes(x = ReturnPeriod, y = CI_upper, label = sprintf("%.2f", CI_upper)),
            vjust = -0.5, color = "#922B21", size = 3, fontface = "italic") +
  
  # 3. Lower CI Labels (Sitting below ribbon)
  geom_text(data = stats_df,
            aes(x = ReturnPeriod, y = CI_lower, label = sprintf("%.2f", CI_lower)),
            vjust = 1.5, color = "#922B21", size = 3, fontface = "italic") +
  
  # F. Scales and Theme
  scale_x_log10(breaks = return_periods, labels = return_periods) +
  annotation_logticks(sides = "b") +
  labs(
    title = "Return Value Analysis: Significant Wave Height",
    subtitle = paste0("Estimates based on ", N_simulations, " simulations"),
    x = "Return Period (Years)",
    y = "Significant Wave Height Hs (m)",
    caption = "Shaded Area: 95% Confidence Interval | Blue Line: Median"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "gray50", size = 11),
    axis.text = element_text(color = "gray30")
  )

# Display the plot
print(p)

# -------------------------------------------------------------------------------
# SIMPLIFIED COVARIATE ANALYSIS
# -------------------------------------------------------------------------------

cat("Computing simplified covariate analysis...\n")

# Coarser grid for quick testing
Cspd_quick <- quantile(data$Cspd, c(0.25, 0.5, 0.75))
angle_quick <- c(0, 90, 180, 270)

covariate_grid_quick <- expand.grid(
  Cspd = Cspd_quick,
  angularDifference = angle_quick
)

# Function for quick conditional RV
compute_conditional_RV_quick <- function(target_Cspd, target_ang, all_sims) {
  rv_vals <- numeric(length(all_sims))
  
  for (i in 1:length(all_sims)) {
    sim_data <- all_sims[[i]]
    
    # Wider windows for quick test
    in_window <- (abs(sim_data$Cspd - target_Cspd) <= 0.8) &
                 (abs(((sim_data$angularDifference - target_ang + 180) %% 360) - 180) <= 60)
    
    filtered <- sim_data$hs[in_window]
    
    if (length(filtered) >= 50) {
      n_chunks <- floor(length(filtered) / 50)
      chunks <- split(filtered, ceiling(seq_along(filtered) / 50)[1:(n_chunks*50)])
      chunk_max <- sapply(chunks, max)
      rv_vals[i] <- quantile(chunk_max, 1 - 1/20)  # 20-year RV
    } else {
      rv_vals[i] <- NA
    }
  }
  
  rv_vals <- rv_vals[!is.na(rv_vals)]
  
  if (length(rv_vals) >= 5) {
    return(list(
      median = median(rv_vals),
      ci_lower = quantile(rv_vals, 0.025),
      ci_upper = quantile(rv_vals, 0.975)
    ))
  } else {
    return(list(median = NA, ci_lower = NA, ci_upper = NA))
  }
}

covariate_grid_quick$RV_20_median <- NA
covariate_grid_quick$RV_20_ci_lower <- NA
covariate_grid_quick$RV_20_ci_upper <- NA

for (i in 1:nrow(covariate_grid_quick)) {
  result <- compute_conditional_RV_quick(
    covariate_grid_quick$Cspd[i],
    covariate_grid_quick$angularDifference[i],
    all_simulations_quick
  )
  covariate_grid_quick$RV_20_median[i] <- result$median
  covariate_grid_quick$RV_20_ci_lower[i] <- result$ci_lower
  covariate_grid_quick$RV_20_ci_upper[i] <- result$ci_upper
}

cat("\nCOVARIATE-CONDITIONAL RETURN VALUES (20-year):\n")
cat("=" , rep("=", 70), "\n", sep = "")
print(covariate_grid_quick[!is.na(covariate_grid_quick$RV_20_median), ], 
      row.names = FALSE, digits = 3)
cat("\n")

# -------------------------------------------------------------------------------
# SAVE RESULTS
# -------------------------------------------------------------------------------

cat("\n===============================================================================\n")
cat("QUICK TEST COMPLETE!\n")
cat("===============================================================================\n")
cat("\nThis was a quick test with:\n")
cat("  - Only", N_simulations, "simulations (full analysis uses 50)\n")
cat("  - Only", T_years, "years per simulation (full analysis uses 100)\n")
cat("  - Coarser covariate grid\n")
cat("\nFor publication-quality results, run the full analysis:\n")
cat("  source('return_value_analysis.R')\n")
cat("\n===============================================================================\n")