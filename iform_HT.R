rm(list=ls())
library(texmex)
library(data.table)
library(ggplot2)
library(dplyr)

data <- fread("C:/these_docs/data2024.csv")

# ==============================================================================
# 1. FIT HT04 MODELS (BIDIRECTIONAL) - DIRECT LAPLACE
# ==============================================================================

fit_ht04_fast <- function(data, var1="hs", var2="tp", threshold_quantile=0.8, sep_hours = 72) {
  
  # 1. Clean Data
  df_clean <- na.omit(data[, c(var1, var2), with = FALSE])
  names(df_clean) <- c("v1", "v2")
  df_clean <- df_clean[v1 > 0 & v2 > 0]
  
  # 2. Every-x max thinning
  g <- ceiling(seq_len(nrow(df_clean)) / sep_hours)
  df_fit <- df_clean[, .(v1 = max(v1), v2 = max(v2)), by = g][, g := NULL]
  
  # 3. Empirical CDF for marginals (full data)
  # Store sorted values and their empirical probabilities
  ecdf_v1 <- ecdf(df_clean$v1)
  ecdf_v2 <- ecdf(df_clean$v2)
  
  # 4. HT04 models (on df_fit)
  m1 <- mex(data = df_fit, which = "v1",
            mqu = threshold_quantile, dqu = threshold_quantile)
  
  m2 <- mex(data = df_fit, which = "v2",
            mqu = threshold_quantile, dqu = threshold_quantile)
  
  return(list(
    df_fit = df_fit,
    model1 = m1,
    model2 = m2,
    ecdf_v1 = ecdf_v1,
    ecdf_v2 = ecdf_v2,
    sorted_v1 = sort(df_clean$v1),
    sorted_v2 = sort(df_clean$v2),
    names = c(v1 = var1, v2 = var2)
  ))
}

# ==============================================================================
# 2. PREDICTION HELPER (ROBUST & SAFE) - DIRECT LAPLACE
# ==============================================================================
predict_ht04_value <- function(driver_val, u_resid, mex_model, 
                               driver_ecdf, response_sorted) {
  
  # 1. Extract HT04 parameters
  dep <- mex_model$dependence
  alpha <- dep$coefficients[1]
  beta  <- dep$coefficients[2]
  residuals_z <- na.omit(dep$Z)
  
  # 2. Physical Driver -> Probability (using ECDF)
  p_driver <- driver_ecdf(driver_val)
  
  # Clamp probabilities
  epsilon <- 1e-6
  p_driver <- pmin(pmax(p_driver, epsilon), 1 - epsilon)
  
  # 3. Probability -> Laplace Scale (Driver)
  lap_driver <- ifelse(p_driver < 0.5, log(2 * p_driver), -log(2 * (1 - p_driver)))
  
  # 4. Get Residual (Z)
  prob_z <- pnorm(u_resid)
  
  # Robust Quantile Lookup
  n_res <- length(residuals_z)
  prob_z_clamped <- pmin(pmax(prob_z, 1/(n_res+1)), n_res/(n_res+1))
  z_star <- quantile(residuals_z, probs = prob_z_clamped, type = 8)
  
  # 5. Regression Equation (HT04)
  lap_response <- alpha * lap_driver + (abs(lap_driver)^beta) * z_star
  
  # 6. Laplace Scale (Response) -> Probability
  prob_response <- ifelse(lap_response < 0, 
                          0.5 * exp(lap_response), 
                          1 - 0.5 * exp(-lap_response))
  
  prob_response <- pmin(pmax(prob_response, epsilon), 1 - epsilon)
  
  # 7. Probability -> Physical Response (using empirical quantile)
  n <- length(response_sorted)
  idx <- prob_response * n
  response_val <- approx(x = seq_len(n), y = response_sorted, xout = idx, 
                         method = "linear", rule = 2)$y
  
  return(as.numeric(response_val))
}

# ==============================================================================
# 3. CALCULATE CONTOURS (SMOOTH BLENDING) - DIRECT LAPLACE
# ==============================================================================
calculate_contours_bidirectional <- function(models, 
                                       return_periods = c(1,2,5,10,20, 50, 100), 
                                       sea_state_duration_hours = 1) {
  
  results_list <- list()
  
  # Extract sorted values for fast access
  sorted_v1 <- models$sorted_v1
  sorted_v2 <- models$sorted_v2
  n_v1 <- length(sorted_v1)
  n_v2 <- length(sorted_v2)
  
  # Helper: Vectorized HT04 Prediction (Direct Laplace)
  predict_vec <- function(driver_vals, u_resid, mex_model, response_sorted) {
    
    # A. Get Model Parameters
    dep <- mex_model$dependence
    a <- dep$coefficients[1] # alpha
    b <- dep$coefficients[2] # beta
    Z <- na.omit(dep$Z)
    
    # B. Driver values are already probabilities (p_d)
    p_d <- pmin(pmax(driver_vals, 1e-7), 1 - 1e-7)
    lap_d <- ifelse(p_d < 0.5, log(2 * p_d), -log(2 * (1 - p_d)))
    
    # C. Get Z quantile (Vectorized)
    p_z <- pmin(pmax(pnorm(u_resid), 1/(length(Z)+1)), length(Z)/(length(Z)+1))
    z_star <- quantile(Z, probs = p_z, type = 8, names = FALSE)
    
    # D. HT04 Equation
    lap_resp <- a * lap_d + (abs(lap_d)^b) * z_star
    
    # E. Back transform to Probability
    p_resp <- ifelse(lap_resp < 0, 0.5*exp(lap_resp), 1 - 0.5*exp(-lap_resp))
    p_resp <- pmin(pmax(p_resp, 1e-7), 1 - 1e-7)
    
    # F. Probability to Physical (empirical quantile)
    n <- length(response_sorted)
    idx <- p_resp * n
    approx(x = seq_len(n), y = response_sorted, xout = idx, 
           method = "linear", rule = 2)$y
  }
  
  omission_factor_alpha2 <- 0
  
  for (rp in return_periods) {
    
    # 1. Define Circular Probability Space (U-space)
    n_events <- (rp * 365 * 24) / sea_state_duration_hours
    beta <- qnorm(1 - 1/n_events)
    beta <- beta / sqrt(1 - omission_factor_alpha2)

    theta <- seq(0, 2*pi, length.out = 2000)
    u1 <- beta * cos(theta)
    u2 <- beta * sin(theta)
    
    v1_out <- numeric(length(theta))
    v2_out <- numeric(length(theta))
    
    # 2. SECTOR SWITCHING LOGIC
    mask_s1 <- abs(u1) >= abs(u2)
    mask_s2 <- !mask_s1
    
    # --- Process Sector 1 (Driver = V1/Hs) ---
    if(any(mask_s1)) {
      # V1 is EXACTLY the marginal value (empirical quantile)
      p_v1 <- pnorm(u1[mask_s1])
      p_v1 <- pmin(pmax(p_v1, 1e-7), 1 - 1e-7)
      idx_v1 <- p_v1 * n_v1
      v1_vals <- approx(x = seq_len(n_v1), y = sorted_v1, xout = idx_v1, 
                        method = "linear", rule = 2)$y
      
      # V2 is Predicted using Model 1
      v2_vals <- predict_vec(p_v1, u2[mask_s1], models$model1, sorted_v2)
      
      v1_out[mask_s1] <- v1_vals
      v2_out[mask_s1] <- v2_vals
    }
    
    # --- Process Sector 2 (Driver = V2/Tp) ---
    if(any(mask_s2)) {
      # V2 is EXACTLY the marginal value (empirical quantile)
      p_v2 <- pnorm(u2[mask_s2])
      p_v2 <- pmin(pmax(p_v2, 1e-7), 1 - 1e-7)
      idx_v2 <- p_v2 * n_v2
      v2_vals <- approx(x = seq_len(n_v2), y = sorted_v2, xout = idx_v2, 
                        method = "linear", rule = 2)$y
      
      # V1 is Predicted using Model 2
      v1_vals <- predict_vec(p_v2, u1[mask_s2], models$model2, sorted_v1)
      
      v1_out[mask_s2] <- v1_vals
      v2_out[mask_s2] <- v2_vals
    }
    
    # 3. Smoothing
    w_size <- 20
    v1_smooth <- stats::filter(c(v1_out, v1_out[1:w_size]), rep(1/w_size, w_size), sides=1)[(w_size+1):(length(v1_out)+w_size)]
    v2_smooth <- stats::filter(c(v2_out, v2_out[1:w_size]), rep(1/w_size, w_size), sides=1)[(w_size+1):(length(v2_out)+w_size)]
    
    res_df <- data.frame(RP = as.factor(rp), V1 = v1_smooth, V2 = v2_smooth)
    names(res_df)[2:3] <- c(models$names["v1"], models$names["v2"])
    results_list[[paste0(rp)]] <- res_df
  }
  
  return(do.call(rbind, results_list))
}

# ==============================================================================
# 4. VISUALIZATION FUNCTIONS
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
  
  if(is.null(xlab)) xlab <- x_col
  if(is.null(ylab)) ylab <- y_col
  
  design_points <- contour_data %>%
    group_by(RP) %>%
    filter(.data[[x_col]] == max(.data[[x_col]])) %>%
    ungroup()
  
  dp_100 <- design_points %>%
    filter(RP == 100 | as.character(RP) == "100") %>%
    mutate(label_text = paste0("Point de conception\n (100 ans) :\n",
                               "Hs = ", round(.data[[x_col]], 2), "\n",
                               "Tp = ", round(.data[[y_col]], 2)))
  
  p <- ggplot() +
    geom_point(data = raw_data, aes_string(x = x_col, y = y_col), 
               color = "black", alpha = 0.5, size = 0.8) +
    geom_point(data = contour_data, aes_string(x = x_col, y = y_col, color = "factor(RP)"), 
              size = 1.2) +
    geom_point(data = design_points, aes_string(x = x_col, y = y_col),
               color = "red", size = 3, shape = 19) +
    geom_label(data = dp_100, aes_string(x = x_col, y = y_col, label = "label_text"),
               fill = "white", color = "black", fontface = "bold", size = 3.5,
               hjust = 1, vjust = -1, alpha = 0.9, show.legend = FALSE) +
    scale_color_viridis_d(name = "Période de retour\n(années)", option = "plasma", end = 0.9) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    theme_bw() +
    theme(legend.position = "right")
  
  if(!is.null(xlims)) p <- p + scale_x_continuous(limits = xlims)
  if(!is.null(ylims)) p <- p + scale_y_continuous(limits = ylims)
  
  return(p)
}

plot_environmental_contours_noetiquette <- function(raw_data, 
                                        contour_data, 
                                        x_col, 
                                        y_col,
                                        title = "Environmental Contours",
                                        subtitle = NULL,
                                        xlab = NULL,
                                        ylab = NULL,
                                        xlims = NULL,
                                        ylims = NULL) {
  
  if(is.null(xlab)) xlab <- x_col
  if(is.null(ylab)) ylab <- y_col
  
  x_coll <- "hs"
  design_points <- contour_data %>%
    group_by(RP) %>%
    filter(.data[[x_coll]] == max(.data[[x_coll]])) %>%
    ungroup()
  
  p <- ggplot() +
    geom_point(data = raw_data, aes_string(x = x_col, y = y_col), 
               color = "black", alpha = 0.5, size = 0.8) +
    geom_path(data = contour_data, aes_string(x = x_col, y = y_col, color = "factor(RP)"), 
              size = 0.9) +
    geom_point(data = design_points, aes_string(x = x_col, y = y_col),
               color = "red", size = 3, shape = 19) +
    scale_color_viridis_d(name = "Return period\n(years)", option = "plasma", end = 0.9) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    theme_bw() +
    theme(legend.position = "right")
  
  if(!is.null(xlims)) p <- p + scale_x_continuous(limits = xlims)
  if(!is.null(ylims)) p <- p + scale_y_continuous(limits = ylims)
  
  return(p)
}

# ==============================================================================
# 5. EXECUTION & PLOTTING
# ==============================================================================

# Fit Bidirectional Models
ht04_models <- fit_ht04_fast(data, var1="hs", var2="tp", threshold_quantile=0.95, sep_hours = 24)

# Calculate Contours
contours_ht04 <- calculate_contours_bidirectional(ht04_models, 
                                                  return_periods = c(1, 2, 5, 10, 20, 50, 100), 
                                                  sea_state_duration_hours = 1)

plot_environmental_contours_noetiquette(
  raw_data = data,
  contour_data = contours_ht04,
  x_col = "hs",
  y_col = "tp",
  title = "IFORM contours (Hs, Tp)",
  subtitle = "Joint distribution : Heffernan & Tawn (2004) conditional extremes model",
  xlab = "Hs [m]",
  ylab = "Tp [s]",
  ylims = c(0, 25)
)







# Fit Bidirectional Models
ht04_models <- fit_ht04_fast(data, var1="Cspd", var2="hs", threshold_quantile=0.95, sep_hours = 24)

# Calculate Contours
contours_ht04 <- calculate_contours_bidirectional(ht04_models, 
                                                  return_periods = c(1, 2, 5, 10, 20, 50, 100), 
                                                  sea_state_duration_hours = 1)

plot_environmental_contours_noetiquette(
  raw_data = data,
  contour_data = contours_ht04,
  x_col = "Cspd",
  y_col = "hs",
  title = "IFORM contours (Hs, Cs)",
  subtitle = "Joint distribution : Heffernan & Tawn (2004) conditional extremes model",
  xlab = "Cs [m.s-1]",
  ylab = "Hs [m]",
  ylims = c(0, 8)
)


# Fit Bidirectional Models
ht04_models <- fit_ht04_fast(data, var1="angleDiff", var2="hs", threshold_quantile=0.95, sep_hours = 24)

# Calculate Contours
contours_ht04 <- calculate_contours_bidirectional(ht04_models, 
                                                  return_periods = c(1, 2, 5, 10, 20, 50, 100), 
                                                  sea_state_duration_hours = 1)

plot_environmental_contours_noetiquette(
  raw_data = data,
  contour_data = contours_ht04,
  x_col = "angleDiff",
  y_col = "hs",
  title = "IFORM contours (Hs, angular difference)",
  subtitle = "Joint distribution : Heffernan & Tawn (2004) conditional extremes model",
  xlab = "Angular difference [°]",
  ylab = "Hs [m]",
  ylims = c(0, 8)
)
