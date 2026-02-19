library(data.table)
library(POT)
library(evd)
library(extRemes)

# ============================================================================
# SETUP
# ============================================================================
#data <- fread("C:/these_docs/data2024.csv")
data <- fread("C:/Users/minim/Documents/teletravail/data2024.csv")

hs   <- na.omit(data$hs)

chosen_threshold <- quantile(hs, 0.96)
run_length       <- 72
return_periods   <- c(1, 2, 5, 10, 20, 50, 100, 200)
obs_per_year     <- 8766
hs_ref           <- 7.59

cat(sprintf("Threshold: %.2f m (Q96) | Total obs: %d | Exceedances: %d (%.2f%%)\n",
            chosen_threshold, length(hs),
            sum(hs > chosen_threshold), 100 * mean(hs > chosen_threshold)))

# ============================================================================
# DECLUSTERING
# ============================================================================
clusters       <- decluster(data$hs, r=72, "runs", threshold=chosen_threshold)
excesses <- clusters - chosen_threshold

print(excesses)

# ============================================================================
# MODEL FITTING
# ============================================================================

# --- GPD on declustered cluster maxima ---
fit_declust <- fpot(excesses, threshold = 0, model = "gpd")

cat("\n--- GPD Declustered (delta method CIs) ---\n")
cat(sprintf("Scale: %.4f (SE: %.4f) | 95%% CI: [%.4f, %.4f]\n",
            fit_declust$param["scale"], fit_declust$std.err["scale"],
            fit_declust$param["scale"] - 1.96 * fit_declust$std.err["scale"],
            fit_declust$param["scale"] + 1.96 * fit_declust$std.err["scale"]))
cat(sprintf("Shape: %.4f (SE: %.4f) | 95%% CI: [%.4f, %.4f]\n",
            fit_declust$param["shape"], fit_declust$std.err["shape"],
            fit_declust$param["shape"] - 1.96 * fit_declust$std.err["shape"],
            fit_declust$param["shape"] + 1.96 * fit_declust$std.err["shape"]))

# --- GEV on monthly block maxima ---
obs_per_month   <- obs_per_year / 12
block_size      <- round(obs_per_month)
n_blocks        <- floor(length(hs) / block_size)
block_maxima    <- sapply(1:n_blocks, function(i)
  max(hs[((i-1)*block_size + 1):min(i*block_size, length(hs))], na.rm = TRUE))
blocks_per_year <- obs_per_year / block_size

fit_gev <- fgev(block_maxima, std.err = TRUE)

cat("\n--- GEV Block Maxima (delta method CIs) ---\n")
cat(sprintf("Loc:   %.4f (SE: %.4f) | 95%% CI: [%.4f, %.4f]\n",
            fit_gev$estimate["loc"],   fit_gev$std.err["loc"],
            fit_gev$estimate["loc"]   - 1.96 * fit_gev$std.err["loc"],
            fit_gev$estimate["loc"]   + 1.96 * fit_gev$std.err["loc"]))
cat(sprintf("Scale: %.4f (SE: %.4f) | 95%% CI: [%.4f, %.4f]\n",
            fit_gev$estimate["scale"], fit_gev$std.err["scale"],
            fit_gev$estimate["scale"] - 1.96 * fit_gev$std.err["scale"],
            fit_gev$estimate["scale"] + 1.96 * fit_gev$std.err["scale"]))
cat(sprintf("Shape: %.4f (SE: %.4f) | 95%% CI: [%.4f, %.4f]\n",
            fit_gev$estimate["shape"], fit_gev$std.err["shape"],
            fit_gev$estimate["shape"] - 1.96 * fit_gev$std.err["shape"],
            fit_gev$estimate["shape"] + 1.96 * fit_gev$std.err["shape"]))

# ============================================================================
# RETURN VALUE ESTIMATION — delta method with all three fixes applied
# ============================================================================
n_years        <- length(hs) / obs_per_year
n_clusters     <- length(excesses[excesses > 0])
lambda_declust <- n_clusters / n_years
# FIX 2: Poisson variance of lambda estimate
var_lambda <- lambda_declust / n_years

return_periods_gev <- return_periods[return_periods >= 2]

rv_declust <- data.frame(period = return_periods,
                         level  = numeric(length(return_periods)),
                         lower  = numeric(length(return_periods)),
                         upper  = numeric(length(return_periods)),
                         method = "GPD Declustered")

rv_gev     <- data.frame(period = return_periods_gev,
                         level  = numeric(length(return_periods_gev)),
                         lower  = numeric(length(return_periods_gev)),
                         upper  = numeric(length(return_periods_gev)),
                         method = "Block Maxima (GEV)")

# ---------------------------------------------------------------------------
# GPD return levels — delta method now includes lambda uncertainty (FIX 2)
# ---------------------------------------------------------------------------
xi    <- fit_declust$param["shape"]
sigma <- fit_declust$param["scale"]
vcov  <- fit_declust$var.cov   # 2x2 covariance matrix for (sigma, xi)

for (i in seq_along(return_periods)) {
  T_i <- return_periods[i]
  m   <- lambda_declust * T_i   # expected number of events in T years

  z_m <- if (abs(xi) < 1e-6) chosen_threshold + sigma * log(m)
         else                 chosen_threshold + (sigma / xi) * (m^xi - 1)

  # Gradient w.r.t. (sigma, xi) — unchanged from original (correct)
  grad_theta <- if (abs(xi) < 1e-6) {
    c(log(m), 0)
  } else {
    c((m^xi - 1) / xi,
      -sigma * (m^xi - 1) / xi^2 + sigma * m^xi * log(m) / xi)
  }

  # FIX 2: gradient w.r.t. lambda, then propagate lambda variance
  grad_lambda <- if (abs(xi) < 1e-6) {
    sigma * T_i / m                          # d/d_lambda [sigma * log(lambda*T)]
  } else {
    sigma * T_i * m^(xi - 1)                # d/d_lambda [(sigma/xi)*(m^xi - 1)]
  }

  var_theta  <- as.numeric(t(grad_theta) %*% vcov %*% grad_theta)
  var_lam    <- grad_lambda^2 * var_lambda
  se_z       <- sqrt(var_theta + var_lam)

  rv_declust$level[i] <- z_m
  rv_declust$lower[i] <- z_m - 1.96 * se_z
  rv_declust$upper[i] <- z_m + 1.96 * se_z
}

# ---------------------------------------------------------------------------
# GEV return levels — FIX 1: corrected sign in d(z_p)/d(xi)
# ---------------------------------------------------------------------------
mu       <- fit_gev$estimate["loc"]
sigma_g  <- fit_gev$estimate["scale"]
xi_g     <- fit_gev$estimate["shape"]
vcov_gev <- fit_gev$var.cov   # 3x3 covariance matrix for (loc, scale, shape)

for (i in seq_along(return_periods_gev)) {
  p   <- 1 / (return_periods_gev[i] * blocks_per_year)
  y_p <- -log(1 - p)   # reduced variate

  z_p <- if (abs(xi_g) < 1e-6) mu - sigma_g * log(y_p)
         else                   mu + (sigma_g / xi_g) * (y_p^(-xi_g) - 1)

  # FIX 1: corrected gradient — note the MINUS sign before the log(y_p) term
  grad <- if (abs(xi_g) < 1e-6) {
    c(1, -log(y_p), 0)
  } else {
    c(1,
      (y_p^(-xi_g) - 1) / xi_g,
      -sigma_g * (y_p^(-xi_g) - 1) / xi_g^2 - sigma_g * y_p^(-xi_g) * log(y_p) / xi_g)
      #                                       ^
      #                          FIXED: was '+', must be '-'
      #   because d/dxi [y^{-xi}] = -log(y) * y^{-xi}
  }

  se_z <- sqrt(as.numeric(t(grad) %*% vcov_gev %*% grad))
  rv_gev$level[i] <- z_p
  rv_gev$lower[i] <- z_p - 1.96 * se_z
  rv_gev$upper[i] <- z_p + 1.96 * se_z
}

return_values <- rbind(rv_declust, rv_gev)

# ============================================================================
# RETURN PERIOD AT REFERENCE Hs — log-linear inverse interpolation
# ============================================================================
interp_return_period <- function(rv_df, hs_target) {
  below <- rv_df[rv_df$level <= hs_target, ]
  above <- rv_df[rv_df$level >  hs_target, ]
  if (nrow(below) == 0 || nrow(above) == 0) return(NA)
  x1 <- below$period[nrow(below)]; y1 <- below$level[nrow(below)]
  x2 <- above$period[1];           y2 <- above$level[1]
  exp(log(x1) + (hs_target - y1) / (y2 - y1) * (log(x2) - log(x1)))
}

T_declust_ref <- interp_return_period(rv_declust, hs_ref)
T_gev_ref     <- interp_return_period(rv_gev,     hs_ref)

cat(sprintf("\nReturn period for Hs = %.2f m:\n", hs_ref))
cat(sprintf("  GPD Declustered:    %.1f years\n", T_declust_ref))
cat(sprintf("  Block Maxima (GEV): %.1f years\n", T_gev_ref))

# ============================================================================
# PRINT RETURN VALUE TABLES
# ============================================================================
print_rv_table <- function(df, title) {
  cat(sprintf("\n--- %s ---\n", title))
  cat("T (years) | Return Level (m) | 95% CI\n")
  for (i in 1:nrow(df))
    cat(sprintf("%9.0f | %16.2f | [%.2f, %.2f]\n",
                df$period[i], df$level[i], df$lower[i], df$upper[i]))
}

print_rv_table(rv_declust, "GPD Declustered")
print_rv_table(rv_gev,     "Block Maxima (GEV)")

# ============================================================================
# NOTE ON METHODOLOGY
# ============================================================================
# The delta method produces symmetric CIs which can be poorly calibrated for
# extreme value parameters, especially shape xi. If precision matters, consider
# replacing these with profile likelihood CIs via:
#   confint(fit_declust)   # POT
#   profile(fit_gev)       # evd
# ============================================================================