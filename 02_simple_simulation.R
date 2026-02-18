library(MASS)
library(fitdistrplus)
library(ggplot2)

# Requires objects from 01_threshold_gpd_fit.R:
#   data, fit_ald, fit_gpd, excess_df

# -----------------------------------------------------------------------------
# 1. FIT BODY DISTRIBUTIONS
# -----------------------------------------------------------------------------
body_df <- subset(data, excess == 0)

fit_weibull <- fitdist(body_df$hs, "weibull")
fit_gamma   <- fitdist(body_df$hs, "gamma")
fit_lnorm   <- fitdist(body_df$hs, "lnorm")
fit_norm    <- fitdist(body_df$hs, "norm")

comparison <- data.frame(
  Distribution = c("Weibull", "Gamma", "Lognormal", "Normal"),
  AIC = c(fit_weibull$aic, fit_gamma$aic, fit_lnorm$aic, fit_norm$aic),
  BIC = c(fit_weibull$bic, fit_gamma$bic, fit_lnorm$bic, fit_norm$bic)
)
print(comparison)

best_fit  <- list(fit_weibull, fit_gamma, fit_lnorm, fit_norm)[[which.min(comparison$AIC)]]
best_name <- comparison$Distribution[which.min(comparison$AIC)]

# -----------------------------------------------------------------------------
# 2. COMPUTE SIMULATION PARAMETERS
# -----------------------------------------------------------------------------
exceedance_rate   <- sum(data$excess > 0) / nrow(data)
mean_cluster_size <- sum(data$hs > data$threshold) / sum(data$excess > 0)

# -----------------------------------------------------------------------------
# 3. SIMULATION FUNCTION
# -----------------------------------------------------------------------------
simulate_structured_storms <- function(n_total, data_ref, body_fit, gpd_fit, ald_fit, exc_rate, mean_dur) {

  storm_ref <- subset(data_ref, excess > 0)
  calm_ref  <- subset(data_ref, excess == 0)

  bw_c          <- MASS::bandwidth.nrd(data_ref$Cspd) * 0.5
  bw_a_circular <- MASS::bandwidth.nrd(data_ref$angularDifference) * 0.05

  total_storm_steps <- round(n_total * exc_rate)
  num_storms        <- round(total_storm_steps)

  storm_data_list <- list()

  # ---------------------------------------------------------------------------
  # PART A: STORM EVENTS
  # ---------------------------------------------------------------------------
  if (num_storms > 0) {
    sampled_indices <- sample(1:nrow(storm_ref), size = num_storms, replace = TRUE)

    cspd_storm    <- storm_ref$Cspd[sampled_indices] + rnorm(num_storms, 0, bw_c)
    ang_storm_raw <- storm_ref$angularDifference[sampled_indices] + rnorm(num_storms, 0, bw_a_circular)
    ang_storm     <- ang_storm_raw %% 360
    cspd_storm    <- pmax(pmin(cspd_storm, max(data_ref$Cspd)), min(data_ref$Cspd))

    storm_covs  <- data.frame(Cspd = cspd_storm, angularDifference = ang_storm)
    storm_thresh <- predict(ald_fit, newdata = storm_covs, type = "response")[, 1]
    gpd_params  <- predict(gpd_fit, newdata = storm_covs, type = "response")
    scale_vec   <- gpd_params[, 1]
    shape_vec   <- gpd_params[, 2]

    for (i in 1:num_storms) {
      dur <- max(1, rpois(1, mean_dur))
      u   <- runif(1)

      if (abs(shape_vec[i]) < 1e-6) {
        peak_excess <- -scale_vec[i] * log(1 - u)
      } else {
        peak_excess <- (scale_vec[i] / shape_vec[i]) * ((1 - u)^(-shape_vec[i]) - 1)
      }
      peak_excess <- max(peak_excess, 0)

      if (dur == 1) {
        profile <- 1
      } else {
        t_peak <- ceiling(dur / 2)
        rise <- seq(0.2, 1, length.out = t_peak + 1)[-1]
        seq_fall <- seq(1, 0.2, length.out = (dur - t_peak) + 2)
        fall <- seq_fall[-c(1, length(seq_fall))]
        if (length(fall) == 0 && (dur - t_peak) > 0) fall <- 0.1
        profile <- c(rise, fall)
        if (length(profile) != dur) {
          profile <- approx(seq(0, 1, along.with = profile), profile, n = dur)$y
        }
      }

      hs_storm <- storm_thresh[i] + (peak_excess * profile)

      cspd_variation    <- cspd_storm[i] + rnorm(dur, 0, bw_c * 0.3 * (0.5 + 0.5 * profile))
      cspd_variation    <- pmax(pmin(cspd_variation, max(data_ref$Cspd)), min(data_ref$Cspd))
      ang_variation_raw <- ang_storm[i] + rnorm(dur, 0, bw_a_circular * 2 * (0.5 + 0.5 * profile))
      ang_variation     <- ang_variation_raw %% 360

      varying_covs   <- data.frame(Cspd = cspd_variation, angularDifference = ang_variation)
      varying_thresh <- predict(ald_fit, newdata = varying_covs, type = "response")[, 1]
      varying_gpd    <- predict(gpd_fit, newdata = varying_covs, type = "response")

      storm_data_list[[i]] <- data.frame(
        Cspd             = cspd_variation,
        angularDifference = ang_variation,
        hs               = hs_storm,
        extreme          = TRUE,
        threshold        = varying_thresh,
        gpd_scale        = varying_gpd[, 1],
        gpd_shape        = varying_gpd[, 2]
      )
    }
  }

  simulated_storm_steps <- sum(sapply(storm_data_list, nrow))
  n_calm_needed         <- n_total - simulated_storm_steps

  # ---------------------------------------------------------------------------
  # PART B: CALM EVENTS
  # ---------------------------------------------------------------------------
  if (n_calm_needed > 0) {
    sampled_indices_calm <- sample(1:nrow(calm_ref), size = n_calm_needed, replace = TRUE)

    cspd_calm    <- calm_ref$Cspd[sampled_indices_calm] + rnorm(n_calm_needed, 0, bw_c)
    ang_calm_raw <- calm_ref$angularDifference[sampled_indices_calm] + rnorm(n_calm_needed, 0, bw_a_circular)
    ang_calm     <- ang_calm_raw %% 360
    cspd_calm    <- pmax(pmin(cspd_calm, max(data_ref$Cspd)), min(data_ref$Cspd))

    calm_covs  <- data.frame(Cspd = cspd_calm, angularDifference = ang_calm)
    calm_thresh <- predict(ald_fit, newdata = calm_covs, type = "response")[, 1]

    dist_name <- body_fit$distname
    sample_body <- function(n) {
      if (dist_name == "weibull") rweibull(n, shape = body_fit$estimate["shape"], scale = body_fit$estimate["scale"])
      else if (dist_name == "gamma") rgamma(n, shape = body_fit$estimate["shape"], rate = body_fit$estimate["rate"])
      else if (dist_name == "lnorm") rlnorm(n, meanlog = body_fit$estimate["meanlog"], sdlog = body_fit$estimate["sdlog"])
      else rnorm(n, mean = body_fit$estimate["mean"], sd = body_fit$estimate["sd"])
    }

    raw_body       <- sample_body(n_calm_needed)
    invalid_indices <- which(raw_body > calm_thresh)
    iter <- 0
    while (length(invalid_indices) > 0 && iter < 100) {
      raw_body[invalid_indices] <- sample_body(length(invalid_indices))
      invalid_indices <- invalid_indices[raw_body[invalid_indices] > calm_thresh[invalid_indices]]
      iter <- iter + 1
    }

    calm_df <- data.frame(
      Cspd              = cspd_calm,
      angularDifference = ang_calm,
      hs                = raw_body,
      extreme           = FALSE,
      threshold         = calm_thresh,
      gpd_scale         = NA,
      gpd_shape         = NA
    )

    combined_list <- list()
    calm_chunks   <- split(calm_df, sort(sample(1:max(1, num_storms), nrow(calm_df), replace = TRUE)))
    max_len       <- max(length(calm_chunks), length(storm_data_list))

    for (k in 1:max_len) {
      if (k <= length(calm_chunks))    combined_list[[length(combined_list) + 1]] <- calm_chunks[[k]]
      if (k <= length(storm_data_list)) combined_list[[length(combined_list) + 1]] <- storm_data_list[[k]]
    }

    return(do.call(rbind, combined_list))
  } else {
    return(do.call(rbind, storm_data_list))
  }
}

# -----------------------------------------------------------------------------
# 4. RUN SIMULATION & DIAGNOSTICS
# -----------------------------------------------------------------------------
n_years       <- 31
n_sim         <- 8760 * n_years
simulated_data <- simulate_structured_storms(n_sim, data, best_fit, fit_gpd, fit_ald, exceedance_rate, mean_cluster_size)

par(mfrow = c(2, 2))
hist(data$hs, prob = TRUE, breaks = 50, main = "Original Hs",   xlab = "Hs", col = "lightblue")
hist(simulated_data$hs, prob = TRUE, breaks = 50, main = "Simulated Hs", xlab = "Hs", col = "lightcoral")

qqplot(data$hs, simulated_data$hs, main = "Q-Q Plot", xlab = "Original", ylab = "Simulated")
abline(0, 1, col = "red")

plot(density(data$hs), main = "Density Comparison", xlab = "Hs", lwd = 2, col = "blue")
lines(density(simulated_data$hs), lwd = 2, col = "red")
legend("topright", legend = c("Original", "Simulated"), col = c("blue", "red"), lwd = 2)
par(mfrow = c(1, 1))

print(paste("Original  - Mean:", round(mean(data$hs), 3), "SD:", round(sd(data$hs), 3)))
print(paste("Simulated - Mean:", round(mean(simulated_data$hs), 3), "SD:", round(sd(simulated_data$hs), 3)))
print(paste("Original exceedance rate:", round(exceedance_rate, 4)))
print(paste("Simulated exceedance rate:", round(mean(simulated_data$extreme), 4)))
print(paste("Mean cluster duration:", round(mean_cluster_size, 2)))

ggplot(simulated_data, aes(x = 1:nrow(simulated_data), y = hs, color = extreme)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Simulated Time Series", x = "Time", y = "Hs")