library(data.table)
library(fitdistrplus)
library(ggplot2)
library(mgcv) # For GAM

# ==============================================================================
# 1. FUNCTION: FIT STATISTICAL MODELS
# ==============================================================================
# Fits a Marginal Weibull (var1) and a Conditional Lognormal GAM (var2 | var1)
fit_iform_models <- function(data, var_marginal="hs", var_conditional="tp") {
  
  # 1. Prepare Data
  df_clean <- na.omit(data[, c(var_marginal, var_conditional), with=FALSE])
  names(df_clean) <- c("v1", "v2") # Internal renaming for consistency
  df_clean <- df_clean[v1 > 0 & v2 > 0]
  
  # 2. Fit Marginal (Weibull) to Variable 1
  fit_marg <- fitdist(df_clean$v1, "weibull")
  
  # 3. Fit Conditional (Lognormal via GAM) to Variable 2 given Variable 1
  # We model log(v2) because v2 is lognormal
  df_clean$log_v2 <- log(df_clean$v2)
  
  # A. Mean trend
  model_mu <- gam(log_v2 ~ s(v1, k=5), data = df_clean)
  
  # B. Variance trend (Residuals)
  df_clean$resid_sq <- (df_clean$log_v2 - predict(model_mu, df_clean))^2
  model_var <- gam(resid_sq ~ s(v1, k=5), family = Gamma(link="log"), data = df_clean)
  
  # Return list of models and the name mapping
  return(list(
    marg_dist = fit_marg,
    cond_mu = model_mu,
    cond_var = model_var,
    var_names = c(marginal=var_marginal, conditional=var_conditional)
  ))
}

# ==============================================================================
# 2. FUNCTION: CALCULATE CONTOURS
# ==============================================================================
calculate_contours <- function(models, 
                               return_periods = c(1, 10, 50, 100), 
                               sea_state_duration_hours = 1,
                               omission_factor_alpha2 = 0.15) {
  
  results_list <- list()
  
  # Extract model parameters
  shape_w <- models$marg_dist$estimate["shape"]
  scale_w <- models$marg_dist$estimate["scale"]
  
  # Helper: Transform U1 -> Marginal Variable (Weibull)
  get_v1 <- function(u1) {
    qweibull(pnorm(u1), shape = shape_w, scale = scale_w)
  }
  
  # Helper: Transform U2 -> Conditional Variable (Lognormal GAM)
  get_v2 <- function(v1_val, u2) {
    # Predict Mean and Sigma for the specific v1 value
    pred_data <- data.frame(v1 = v1_val)
    mu <- predict(models$cond_mu, newdata = pred_data)
    var <- predict(models$cond_var, newdata = pred_data, type="response")
    exp(mu + u2 * sqrt(var))
  }
  
  # LOOP over Return Periods
  for (rp in return_periods) {
    
    # Calculate Beta
    n_events <- (rp * 365 * 24) / sea_state_duration_hours
    beta_target <- qnorm(1 - 1/n_events)
    beta_inflated <- beta_target / sqrt(1 - omission_factor_alpha2)
    
    # Generate Circle
    theta <- seq(0, 2*pi, length.out = 500)
    u1 <- beta_inflated * cos(theta)
    u2 <- beta_inflated * sin(theta)
    
    # Transform
    v1_vec <- get_v1(u1)
    v2_vec <- get_v2(v1_vec, u2)
    
    # Store
    tmp_df <- data.frame(
      RP = as.factor(rp), # Factor for plotting order
      V1 = v1_vec,
      V2 = v2_vec
    )
    
    # Rename columns back to original variable names for clarity
    names(tmp_df)[2:3] <- c(models$var_names["marginal"], models$var_names["conditional"])
    
    results_list[[paste0(rp)]] <- tmp_df
  }
  
  return(do.call(rbind, results_list))
}

# ==============================================================================
# 3. FUNCTION: VISUALIZE
# ==============================================================================
plot_environmental_contours <- function(raw_data, 
                                        contour_data, 
                                        x_col, 
                                        y_col,
                                        title = "Environmental Contours",
                                        subtitle = NULL,
                                        xlab = NULL,
                                        ylab = NULL,
                                        xlims = NULL,
                                        ylims = NULL) {
  
  # Defaults if labels not provided
  if(is.null(xlab)) xlab <- x_col
  if(is.null(ylab)) ylab <- y_col
  
  # --- 1. Identifier les Points de Conception (Max Hs pour chaque RP) ---
  # On groupe par Période de Retour (RP) et on prend la ligne où x_col est maximal
  design_points <- contour_data %>%
    group_by(RP) %>%
    filter(.data[[x_col]] == max(.data[[x_col]])) %>%
    ungroup()
  
  # --- 2. Préparer l'étiquette spécifique pour 100 ans ---
  # On filtre pour RP = 100 (en gérant le fait que RP puisse être un facteur ou numérique)
  dp_100 <- design_points %>%
    filter(RP == 100 | RP == "100") %>%
    mutate(label_txt = paste0("Point de conception\n (100 ans) :\n",
                              "Hs = ", round(.data[[x_col]], 2), " m\n",
                              "Tp = ", round(.data[[y_col]], 2), " s"))
  
  # --- 3. Construction du graphique ---
  p <- ggplot() +
    # A. Données brutes
    geom_point(data = raw_data, aes_string(x = x_col, y = y_col), 
               color = "black", alpha = 0.5, size = 0.8) +
    
    # B. Contours (Lignes)
    geom_path(data = contour_data, aes_string(x = x_col, y = y_col, color = "factor(RP)"), 
              size = 1.2) +
    
    # C. Points de conception (Rouges) - Sur tous les contours
    geom_point(data = design_points, aes_string(x = x_col, y = y_col),
               color = "red", size = 3, shape = 19) +
    
    # D. Bulle d'info pour le 100 ans
    geom_label(data = dp_100, aes_string(x = x_col, y = y_col, label = "label_txt"),
               fill = "white", color = "black", fontface = "bold", size = 3.5,
               hjust = 1, # Décale le texte vers la gauche du point pour ne pas le cacher
               vjust = -1,
               alpha = 0.9) +
    
    # E. Styling
    scale_color_viridis_d(name = "Période de retour\n(années)", option = "plasma", end = 0.9) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    theme_bw() +
    theme(legend.position = "right")
  
  # 4. Apply Limits if provided
  if(!is.null(xlims)) p <- p + scale_x_continuous(limits = xlims)
  if(!is.null(ylims)) p <- p + scale_y_continuous(limits = ylims)
  
  return(p)
}


    
# --- Step A: Load Data ---
data <- fread("C:/these_docs/data2024.csv")
# Ensure columns are numeric
data$hs <- as.numeric(data$hs)
data$tp <- as.numeric(data$tp)

# --- Step B: Fit Models ---
# "hs" is the Marginal (V1), "tp" is the Conditional (V2)
my_models <- fit_iform_models(data, var_marginal = "hs", var_conditional = "tp")

# --- Step C: Compute Contours (Multiple Years) ---
# Calculate 1, 10, 50, and 100 year contours
my_contours <- calculate_contours(my_models, 
                                  return_periods = c(1, 2, 5, 10, 20, 50, 100), 
                                  sea_state_duration_hours = 1, # Change if your data is 3-hourly
                                  omission_factor_alpha2 = 0.15)

# --- Step D: Visualize ---
# You want Hs on X-axis, Tp on Y-axis.
# You want Y-axis limit 0 to 30.
plot_environmental_contours(
  raw_data = data,
  contour_data = my_contours,
  x_col = "hs", 
  y_col = "tp",
  title = "Contours IFORM (Hs, Tp), Raz Blanchard",
  subtitle = "Hs: Weibull | Tp|Hs: Lognormal-GAM | Facteur d'omission: 0.15",
  xlab = "Hs [m]",
  ylab = "Tp [s]",
  ylims = c(0, 30) # YOUR CUSTOM LIMIT
)

  

###################diags##########################
###################diags##########################
###################diags##########################

library(gridExtra)

check_model_fit <- function(models, data) {
  
  # Extract variable names and raw data
  var_marg <- models$var_names["marginal"]
  var_cond <- models$var_names["conditional"]
  
  # Filter data exactly as done in fitting (remove NAs and <=0)
  df <- na.omit(data[, c(var_marg, var_cond), with=FALSE])
  names(df) <- c("v1", "v2")
  df <- df[v1 > 0 & v2 > 0]
  
  # ============================================================
  # PART 1: MARGINAL DIAGNOSTICS (Weibull)
  # ============================================================
  
  # Parameters
  shape <- models$marg_dist$estimate["shape"]
  scale <- models$marg_dist$estimate["scale"]
  
  # A. QQ-Plot Data Preparation
  # Sort data and calculate empirical probabilities
  df <- df[order(df$v1), ]
  n <- nrow(df)
  df$prob_emp <- (1:n) / (n + 1)
  
  # Calculate Theoretical Quantiles based on probabilities
  df$theo_quant <- qweibull(df$prob_emp, shape, scale)
  
  # Plot 1: Marginal QQ-Plot
  p1 <- ggplot(df, aes(x = theo_quant, y = v1)) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    geom_point(alpha = 0.5, size = 0.8) +
    labs(title = paste("Marginal QQ-Plot:", var_marg),
         subtitle = "Points should follow the red line",
         x = "Theoretical Quantiles (Weibull)", 
         y = "Empirical Quantiles (Data)") +
    theme_bw()
  
  # Plot 2: Marginal CDF Comparison
  # Create a smooth theoretical line
  x_seq <- seq(min(df$v1), max(df$v1), length.out = 200)
  theo_cdf <- data.frame(x = x_seq, y = pweibull(x_seq, shape, scale))
  
  p2 <- ggplot() +
    stat_ecdf(data = df, aes(x = v1), geom = "step", color = "black", size = 0.8) +
    geom_line(data = theo_cdf, aes(x = x, y = y), color = "red", size = 1) +
    labs(title = paste("Marginal CDF:", var_marg),
         subtitle = "Black: Empirical Data | Red: Fitted Model",
         x = var_marg, y = "Cumulative Probability") +
    theme_bw()
  
  # ============================================================
  # PART 2: CONDITIONAL DIAGNOSTICS (Lognormal GAM)
  # ============================================================
  
  # To check the conditional fit, we calculate Z-scores (Standardized Residuals).
  # If the model is correct, Z should be Standard Normal N(0,1).
  
  # 1. Predict Mean and Sigma for every data point
  df$log_v2 <- log(df$v2)
  df$pred_mu <- predict(models$cond_mu, newdata = data.frame(v1 = df$v1))
  df$pred_var <- predict(models$cond_var, newdata = data.frame(v1 = df$v1), type = "response")
  df$pred_sigma <- sqrt(df$pred_var)
  
  # 2. Calculate Z-score: (Actual - Predicted_Mean) / Predicted_Sigma
  df$z_score <- (df$log_v2 - df$pred_mu) / df$pred_sigma
  
  # Plot 3: Residual QQ-Plot (Check for Normality)
  p3 <- ggplot(df, aes(sample = z_score)) +
    stat_qq(alpha = 0.5, size = 0.8) +
    stat_qq_line(color = "red", linetype = "dashed") +
    labs(title = paste("Conditional Residuals:", var_cond, "|", var_marg),
         subtitle = "Standardized Residuals should align with Normal(0,1)",
         x = "Theoretical Quantiles (Normal)", 
         y = "Standardized Residuals") +
    theme_bw()
  
  # Plot 4: Residuals vs Covariate (Check for Bias/Heteroscedasticity)
  p4 <- ggplot(df, aes(x = v1, y = z_score)) +
    geom_point(alpha = 0.1) +
    geom_hline(yintercept = c(-2, 0, 2), linetype = "dashed", color = "blue") +
    labs(title = "Residuals vs Covariate",
         subtitle = "Should be a flat cloud centered at 0 (no trend)",
         x = var_marg, 
         y = "Standardized Residuals (Z-score)") +
    theme_bw() +
    coord_cartesian(ylim = c(-4, 4))
  
  # ============================================================
  # OUTPUT
  # ============================================================
  
  # Arrange in a 2x2 Grid
  grid.arrange(p1, p3, nrow=1)
}


check_model_fit(my_models, data)
