library(shiny)
library(resourcecode)
library(resourcecodedata)
library(dplyr)
library(lubridate)
library(ggplot2)
library(plotly)
library(leaflet)
library(mgcv)
library(evgam)
library(quantreg)
library(splines)
library(plot3D)
library(ggnewscale)
library(tidyr)
library(rlang)
library(gratia)
library(extRemes)
library(DT)
library(car)
library(data.table) # Essentiel pour la rapidité de calcul sur les séries temporelles
library(fitdistrplus) # Pour les marginales paramétriques
library(texmex)       # Pour le modèle Heffernan & Tawn
library(viridis)


# --- FONCTIONS GLOBALES POUR CONTOURS ENVIRONNEMENTAUX ---
# --- NOUVELLES FONCTIONS : BOOTSTRAP & GÉOMÉTRIE ---

# Fonction pour interpoler un contour sur une grille angulaire fixe (0 à 360°)
# Nécessaire pour faire des stats "point par point" entre des contours de formes légèrement différentes
interpolate_contour_polar <- function(df_contour, center, n_angles = 360) {
  # df_contour doit avoir les colonnes v1, v2 (ou noms des variables)
  # center : c(mean_x, mean_y)
  
  # 1. Conversion en polaire
  dx <- df_contour[[1]] - center[1]
  dy <- df_contour[[2]] - center[2]
  
  r <- sqrt(dx^2 + dy^2)
  theta <- atan2(dy, dx) # -pi à pi
  
  # Gestion des doublons d'angles pour l'interpolation (garder le r max si conflit, conservateur)
  df_pol <- data.frame(theta = theta, r = r)
  df_pol <- df_pol[order(df_pol$theta), ]
  
  # On s'assure que la boucle est fermée pour l'approx
  df_pol <- rbind(df_pol, data.frame(theta = df_pol$theta[1] + 2*pi, r = df_pol$r[1]))
  
  # 2. Grille cible
  target_theta <- seq(-pi, pi, length.out = n_angles)
  
  # 3. Interpolation linéaire du rayon
  # approx peut échouer si duplicates exacts, on agrège avant si besoin (ici supposé propre via IFORM)
  r_interp <- approx(df_pol$theta, df_pol$r, xout = target_theta, rule = 2)$y
  
  # 4. Retour en Cartésien
  x_new <- center[1] + r_interp * cos(target_theta)
  y_new <- center[2] + r_interp * sin(target_theta)
  
  return(data.frame(x = x_new, y = y_new, theta = target_theta))
}

# Fonction wrapper pour exécuter une itération de contour
run_contour_iteration <- function(data, method, v_x, v_y, period, params) {
  res <- NULL
  
  if(method == "direct") {
    res <- run_direct_iform(
      data = data, vars = c(v_x, v_y), return_period = period,
      spacing = params$spacing, threshold_prob = params$prob, peak_sep_h = params$sea_state
    )
  } else if (method == "wl") {
    # Gestion d'erreur robuste pour le bootstrap (si le fit échoue sur un resample)
    try({
      mod <- fit_iform_wl(as.data.table(data), v_x, v_y)
      res <- calc_contours_wl(mod, period, params$sea_state)
    }, silent = TRUE)
  } else if (method == "ht04") {
    try({
      mod <- fit_iform_ht04(as.data.table(data), v_x, v_y, th_quant = params$ht_threshold)
      res <- calc_contours_ht04(mod, period, params$sea_state)
    }, silent = TRUE)
  }
  
  # Harmonisation des noms de colonnes pour la sortie
  if(!is.null(res)) {
    # On s'assure que les 2 premières colonnes sont les variables
    if(!all(c(v_x, v_y) %in% names(res))) names(res)[1:2] <- c(v_x, v_y)
    return(res[, c(v_x, v_y)]) # Retourne juste les coords
  }
  return(NULL)
}

# 1. Génération des vecteurs de direction (Hypersphère)
generate_directions <- function(ndim, spacing = 0.1) {
  get_simplex_points <- function(n, s) {
    if (n == 1) return(matrix(s, ncol = 1))
    res <- list()
    vals <- seq(0, s, by = spacing)
    for (v in vals) {
      rem <- round(s - v, 10)
      suffix <- get_simplex_points(n - 1, rem)
      prefix <- matrix(v, nrow = nrow(suffix), ncol = 1)
      res[[length(res) + 1]] <- cbind(prefix, suffix)
    }
    return(do.call(rbind, res))
  }
  
  points_pos <- get_simplex_points(ndim, 1.0)
  signs <- as.matrix(expand.grid(rep(list(c(1, -1)), ndim)))
  
  all_dirs <- list()
  for (i in 1:nrow(signs)) {
    orthant_points <- sweep(points_pos, 2, signs[i, ], "*")
    all_dirs[[i]] <- orthant_points
  }
  
  raw_vecs <- unique(do.call(rbind, all_dirs))
  norms <- sqrt(rowSums(raw_vecs^2))
  unit_vecs <- raw_vecs / norms
  
  if (ndim == 2) {
    angles <- atan2(unit_vecs[,2], unit_vecs[,1])
    unit_vecs <- unit_vecs[order(angles), ]
  }
  return(unit_vecs)
}

# 2. Ajustement GPD Contraint (évite les formes irréalistes)
fit_constrained_gpd <- function(x, threshold) {
  nll <- function(theta) {
    scale <- exp(theta[1])
    xi <- theta[2]
    # Contraintes physiques pour éviter des queues trop lourdes ou absurdes
    if (xi >= -1e-4 || xi < -1) return(1e9)
    limit <- threshold - scale/xi
    if (max(x + threshold) > limit) return(1e9)
    
    term <- 1 + xi * (x / scale)
    if (any(term <= 0)) return(1e9)
    
    val <- length(x) * log(scale) + (1/xi + 1) * sum(log(term))
    return(val)
  }
  
  x_bar <- mean(x); s2 <- var(x)
  xi_init <- 0.5 * (1 - (x_bar^2 / s2))
  xi_init <- max(min(xi_init, -0.1), -0.9)
  sigma_init <- x_bar * (1 - xi_init)
  
  tryCatch({
    opt <- optim(c(log(sigma_init), xi_init), nll, method = "Nelder-Mead")
    return(list(scale = exp(opt$par[1]), shape = opt$par[2]))
  }, error = function(e) {
    return(list(scale = sigma_init, shape = xi_init))
  })
}

# 3. Reconstruction du contour 2D à partir des niveaux de retour projetés
reconstruct_contour_2d <- function(U, R, center_mean, center_sd) {
  n <- nrow(U)
  vertices <- matrix(NA, nrow = n, ncol = 2)
  for (i in 1:n) {
    idx1 <- i; idx2 <- if (i == n) 1 else i + 1 
    if (is.na(R[idx1]) || is.na(R[idx2])) next
    u1 <- U[idx1, 1]; v1 <- U[idx1, 2]; r1 <- R[idx1]
    u2 <- U[idx2, 1]; v2 <- U[idx2, 2]; r2 <- R[idx2]
    det <- u1*v2 - v1*u2
    if (abs(det) > 1e-6) {
      x_int <- (r1*v2 - r2*v1) / det
      y_int <- (u1*r2 - u2*r1) / det
      vertices[i, ] <- c(x_int, y_int)
    }
  }
  vertices <- vertices[complete.cases(vertices), ]
  if(nrow(vertices) > 0) {
    vertices[,1] <- vertices[,1] * center_sd[1] + center_mean[1]
    vertices[,2] <- vertices[,2] * center_sd[2] + center_mean[2]
    vertices <- rbind(vertices, vertices[1,]) # Fermer la boucle
  }
  return(vertices)
}

# 4. Fonction Principale Direct IFORM (CORRIGÉE)
run_direct_iform <- function(data, vars, return_period, spacing=0.1, threshold_prob=0.1, peak_sep_h=72) {
  # Conversion et validation
  dt <- as.data.table(data)[, c("time", vars), with = FALSE]
  # Nettoyage initial des données
  if(any(is.na(dt))) dt <- na.omit(dt)
  
  if(nrow(dt) < 100) return(NULL) # Sécurité si trop peu de données
  
  # Normalisation
  stats_list <- lapply(vars, function(v) list(med = median(dt[[v]]), sd = sd(dt[[v]])))
  Y <- matrix(NA, nrow = nrow(dt), ncol = length(vars))
  for (i in seq_along(vars)) Y[, i] <- (dt[[vars[i]]] - stats_list[[i]]$med) / stats_list[[i]]$sd
  
  # Directions
  U <- generate_directions(length(vars), spacing)
  return_values <- numeric(nrow(U))
  
  # Temps total en années
  total_years <- as.numeric(difftime(max(dt$time), min(dt$time), units="days"))/365.25
  
  # Taille fenêtre (en pas de temps)
  dt_hours <- median(as.numeric(diff(dt$time), units="hours"), na.rm=TRUE)
  if(is.na(dt_hours) || dt_hours == 0) dt_hours <- 1 # Sécurité division par 0
  
  w_size <- round(peak_sep_h / dt_hours)
  if (w_size < 1) w_size <- 1
  if (w_size %% 2 == 0) w_size <- w_size + 1
  
  # Calcul par direction
  for (k in 1:nrow(U)) {
    r_proj <- Y %*% U[k,]
    
    # Fast Declustering
    roll_max <- frollapply(as.numeric(r_proj), n=w_size, FUN=max, align="center", fill=NA)
    
    # Identification des pics (Gestion robuste des NA)
    # Si roll_max est NA, l'égalité renvoie NA. On force ces cas à FALSE.
    is_peak <- (r_proj == roll_max)
    is_peak[is.na(is_peak)] <- FALSE 
    
    peak_vals <- r_proj[is_peak]
    
    # Nettoyage ultime avant calculs statistiques
    peak_vals <- peak_vals[!is.na(peak_vals)]
    
    if (length(peak_vals) < 20) { return_values[k] <- NA; next }
    
    # --- CORRECTION DE L'ERREUR ICI (na.rm = TRUE) ---
    u_thresh <- quantile(peak_vals, probs = 1 - threshold_prob, na.rm = TRUE)
    
    excesses <- peak_vals[peak_vals > u_thresh]
    
    if(length(excesses) < 5) { return_values[k] <- NA; next }
    
    fit <- fit_constrained_gpd(excesses - u_thresh, u_thresh)
    lambda <- length(excesses) / total_years
    
    # Calcul Niveau de retour dans l'espace réduit
    xi <- fit$shape; sigma <- fit$scale
    
    # Calcul robuste (éviter les NaN dans le log ou la puissance)
    if (abs(xi) < 1e-6) {
      rv <- u_thresh + sigma * log(lambda * return_period)
    } else {
      # Sécurité contre les nombres complexes ou NaN si le terme est négatif
      term <- lambda * return_period
      if(term <= 0) rv <- NA 
      else rv <- u_thresh + (sigma/xi) * ((term)^xi - 1)
    }
    
    return_values[k] <- rv
  }
  
  # Reconstruction
  if(length(vars) == 2) {
    # On filtre les directions où le calcul a échoué (NA)
    valid_idx <- !is.na(return_values)
    if(sum(valid_idx) < 3) return(NULL) # Pas assez de points pour un contour
    
    U_clean <- U[valid_idx, , drop=FALSE]
    R_clean <- return_values[valid_idx]
    
    coords <- reconstruct_contour_2d(U_clean, R_clean, c(stats_list[[1]]$med, stats_list[[2]]$med), c(stats_list[[1]]$sd, stats_list[[2]]$sd))
    
    if(nrow(coords) > 0) {
        hull_idx <- chull(coords)
        coords <- coords[c(hull_idx, hull_idx[1]), ]
        df_res <- as.data.frame(coords); names(df_res) <- vars
        df_res$Period <- as.factor(return_period)
        return(df_res)
    }
  }
  return(NULL)
}

# --- NOUVELLES FONCTIONS IFORM (Weibull-Lognormal & HT04) ---

# A. IFORM WEIBULL-LOGNORMAL (C-GAM)
fit_iform_wl <- function(data, var_marg, var_cond) {
  df_clean <- na.omit(data[, c(var_marg, var_cond), with=FALSE])
  names(df_clean) <- c("v1", "v2")
  df_clean <- df_clean[v1 > 0 & v2 > 0]
  
  # 1. Marginale V1 (Weibull)
  fit_marg <- fitdist(df_clean$v1, "weibull")
  
  # 2. Conditionnelle V2 | V1 (Lognormal via GAM)
  # On modélise log(v2)
  df_clean$log_v2 <- log(df_clean$v2)
  
  # Moyenne
  model_mu <- gam(log_v2 ~ s(v1, k=5), data = df_clean)
  
  # Variance (sur les résidus au carré)
  df_clean$resid_sq <- (df_clean$log_v2 - predict(model_mu, df_clean))^2
  # Gamma link log pour la variance positive
  model_var <- gam(resid_sq ~ s(v1, k=5), family = Gamma(link="log"), data = df_clean)
  
  return(list(
    type = "WL",
    marg_dist = fit_marg,
    cond_mu = model_mu,
    cond_var = model_var,
    vars = c(var_marg, var_cond)
  ))
}

calc_contours_wl <- function(models, periods, sea_state_h=1) {
  res <- data.frame()
  shape_w <- models$marg_dist$estimate["shape"]
  scale_w <- models$marg_dist$estimate["scale"]
  
  for(rp in periods) {
    # Cercle dans l'espace U
    n_events <- (rp * 365.25 * 24) / sea_state_h
    beta <- qnorm(1 - 1/n_events)
    # Inflation due à l'omission factor (ex: 0.15 standard)
    beta <- beta / sqrt(1 - 0.15) 
    
    theta <- seq(0, 2*pi, length.out=360)
    u1 <- beta * cos(theta)
    u2 <- beta * sin(theta)
    
    # Transformation Inverse
    # 1. U1 -> V1 (Weibull)
    v1_vec <- qweibull(pnorm(u1), shape=shape_w, scale=scale_w)
    
    # 2. U2 -> V2 (Lognormal cond)
    # Pour chaque v1, on prédit mu et sigma
    pred_data <- data.frame(v1 = v1_vec)
    mu_vec <- predict(models$cond_mu, newdata=pred_data)
    var_vec <- predict(models$cond_var, newdata=pred_data, type="response")
    sigma_vec <- sqrt(var_vec)
    
    # V2 = exp(mu + U2 * sigma)
    v2_vec <- exp(mu_vec + u2 * sigma_vec)
    
    tmp <- data.frame(v1=v1_vec, v2=v2_vec, Period=factor(rp))
    names(tmp)[1:2] <- models$vars
    res <- rbind(res, tmp)
  }
  return(res)
}

# B. IFORM HEFFERNAN & TAWN (HT04)
fit_iform_ht04 <- function(data, var1, var2, th_quant=0.8) {
  df_clean <- na.omit(data[, c(var1, var2), with=FALSE])
  names(df_clean) <- c("v1", "v2")
  df_clean <- df_clean[v1 > 0 & v2 > 0]
  
  # Marginals (Weibull for simplicity/robustness in this context)
  m_v1 <- fitdist(df_clean$v1, "weibull")
  m_v2 <- fitdist(df_clean$v2, "weibull")
  
  # HT04 Models (Bidirectional)
  # Model 1: V2 | V1 extreme
  mod_v2_v1 <- mex(data=df_clean, which="v1", mqu=th_quant, dqu=th_quant)
  # Model 2: V1 | V2 extreme
  mod_v1_v2 <- mex(data=df_clean, which="v2", mqu=th_quant, dqu=th_quant)
  
  return(list(
    type = "HT04",
    m1 = mod_v2_v1, m2 = mod_v1_v2,
    d1 = m_v1, d2 = m_v2,
    vars = c(var1, var2)
  ))
}

# Helper pour HT04 Prediction
pred_ht04_val <- function(driver_val, u_resid, mex_mod, d_driver, d_resp) {
  dep <- mex_mod$dependence
  a <- dep$coefficients[1]; b <- dep$coefficients[2]
  Z <- na.omit(dep$Z)
  
  # 1. Driver Val -> Prob -> Laplace
  p_d <- pweibull(driver_val, shape=d_driver$estimate["shape"], scale=d_driver$estimate["scale"])
  p_d <- pmin(pmax(p_d, 1e-6), 1-1e-6) # Clamp
  lap_d <- ifelse(p_d < 0.5, log(2*p_d), -log(2*(1-p_d)))
  
  # 2. U_resid -> Z residual
  p_z <- pnorm(u_resid)
  n_z <- length(Z)
  p_z <- pmin(pmax(p_z, 1/(n_z+1)), n_z/(n_z+1))
  z_star <- quantile(Z, probs=p_z, type=8)
  
  # 3. Model
  lap_r <- a * lap_d + (abs(lap_d)^b) * z_star
  
  # 4. Laplace -> Prob -> Response Val
  p_r <- ifelse(lap_r < 0, 0.5*exp(lap_r), 1 - 0.5*exp(-lap_r))
  p_r <- pmin(pmax(p_r, 1e-6), 1-1e-6)
  
  val_r <- qweibull(p_r, shape=d_resp$estimate["shape"], scale=d_resp$estimate["scale"])
  return(as.numeric(val_r))
}

calc_contours_ht04 <- function(models, periods, sea_state_h=1) {
  res <- data.frame()
  
  for(rp in periods) {
    n_events <- (rp * 365.25 * 24) / sea_state_h
    beta <- qnorm(1 - 1/n_events) / sqrt(1 - 0.15)
    
    theta <- seq(0, 2*pi, length.out=360)
    u1 <- beta * cos(theta)
    u2 <- beta * sin(theta)
    
    v1_out <- numeric(length(theta))
    v2_out <- numeric(length(theta))
    
    for(i in seq_along(theta)) {
      # Weighting (Sigmoid blending around 45 deg)
      diff_u <- abs(u1[i]) - abs(u2[i])
      w <- 1 / (1 + exp(-5 * diff_u))
      
      # M1: V1 Driver
      v1_m1 <- qweibull(pnorm(u1[i]), shape=models$d1$estimate["shape"], scale=models$d1$estimate["scale"])
      v2_m1 <- pred_ht04_val(v1_m1, u2[i], models$m1, models$d1, models$d2)
      
      # M2: V2 Driver
      v2_m2 <- qweibull(pnorm(u2[i]), shape=models$d2$estimate["shape"], scale=models$d2$estimate["scale"])
      v1_m2 <- pred_ht04_val(v2_m2, u1[i], models$m2, models$d2, models$d1)
      
      v1_out[i] <- w * v1_m1 + (1-w) * v1_m2
      v2_out[i] <- w * v2_m1 + (1-w) * v2_m2
    }
    
    tmp <- data.frame(v1=v1_out, v2=v2_out, Period=factor(rp))
    names(tmp)[1:2] <- models$vars
    res <- rbind(res, tmp)
  }
  return(res)
}


var_groups <- list(
  "Hauteur significative des vagues (m)" = "hs",
  "Période de pic (s)" = "tp",
  "Cambrure" = "steepness",
  "Vitesse et direction du courant" = c("ucur", "vcur"),
  "Vitesse du vent (m/s)" = c("uwnd", "vwnd"),
  "Direction du pic des vagues (°)" = "dp",
  "Hauteur d'eau (m)" = "wlv",
  "Vitesse orbitale au fond (RMS) (m/s)" = c("uubr", "vubr")
)

vars_dir <- c(
  "Direction du pic des vagues (°)" = "dp",
  "Direction du courant (°)" = "Cdirdeg",
  "Ecart angulaire vagues-courant (°)" = "angleDiff",
  "Direction du vent (°)" = "Wdirdeg"
)

vars_cat <- c(
  "Provenance cardinale des vagues" = "dpCardinale",
  "Profil vagues-courant" = "pvc",
  "Profil vagues-courant (avec étale)" = "pvce",
  "Interaction Vagues-Courant" = "dirvc",
  "Séparation flot-jusant" = "fj",
  "Séparation flot-jusant-étale" = "fje",
  "Mois" = "mois",
  "Année" = "annee"
)

vars_num_pure <- c(
  "Hauteur significative des vagues (m)" = "hs",
  "Période du pic (s)" = "tp",
  "Cambrure" = "steepness",
  "Vitesse du courant (m/s)" = "Cspd",
  "Vitesse orbitale au fond (m/s)" = "OrbSpeed",
  "Vitesse du vent (m/s)" = "Wspd",
  "Hauteur d'eau (m)" = "wlv",
  "Cycle annuel" = "cycle",
  "Composante U Orbitale (m/s)" = "uubr",
  "Composante V Orbitale (m/s)" = "vubr"
)

vars_all_num <- c(vars_num_pure, vars_dir)
var_labels <- c(vars_num_pure, vars_dir, vars_cat)

var_colors <- c(
  "hs" = "#0073C2", "tp" = "#6699CC", "steepness" = "#333333",
  "OrbSpeed" = "orange", "Cspd" = "#009E73", "angleDiff" = "#CC79A7",
  "dirvc" = "#882255", "Wspd" = "#94C973", "dp" = "#F0E442",
  "wlv" = "#56B4E9", "Cdirdeg" = "red", "Wdirdeg" = "#A0522D"
)

get_selected_variables <- function(selected_groups) {
  vars <- unique(unlist(var_groups[selected_groups]))
  if (any(c("Vitesse et direction du courant", "Direction du pic des vagues (°)") %in% selected_groups)) {
    vars <- unique(c(vars, "dp", "ucur", "vcur"))
  }
  return(vars)
}

create_location_map <- function() {
  leaflet(rscd_field) %>% 
    addTiles() %>% 
    setView(-4, 48.5, 7) %>%
    addCircleMarkers(lng=~longitude, lat=~latitude, layerId=~node, 
                     radius=3, color="blue", opacity = 0.8, fillOpacity = 0.5,
                     popup=~paste("Node:", node))
}

transform_dataset <- function(data, variables) {
  opposite_dir <- c("N"="S", "NE"="SO", "E"="O", "SE"="NO", "S"="N", "SO"="NE", "O"="E", "NO"="SE")
  
  df <- data %>%
    mutate(mois=month(time), annee=year(time), jour=day(time), heure=hour(time))
  
  if(all(c("uubr", "vubr") %in% names(df))) df <- df %>% mutate(OrbSpeed = sqrt(uubr^2 + vubr^2))
  if(all(c("hs", "tp") %in% names(df))) df <- df %>% mutate(steepness = hs / tp)
  
  df <- df %>%
    mutate(
      current_datetime = as.POSIXct(sprintf("%04d-%02d-%02d %02d:00:00", annee, mois, jour, heure), tz = "UTC"),
      year_start = as.POSIXct(sprintf("%04d-01-01 00:00:00", annee), tz = "UTC"),
      cycle = as.numeric(difftime(current_datetime, year_start, units = "hours"))
    ) %>%
    group_by(annee) %>%
    mutate(cycle = cycle / max(cycle, na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::select(-current_datetime, -year_start)
  
  if(all(c("ucur", "vcur") %in% names(df))) {
    df <- df %>%
      mutate(
        Cspd = sqrt(ucur^2 + vcur^2), Cdir = atan2(vcur, ucur), Cdirdeg = (90 - Cdir * (180/pi)) %% 360, 
        fje = factor(case_when(Cspd < 1 ~ "étale", Cdirdeg %/% 180 == 1 ~ "jusant", TRUE ~ "flot")),
        fj = factor(if_else((Cdirdeg %/% 180) == 1, "jusant", "flot"))
      )
  }
  
  if("wlv" %in% names(df)) df <- df %>% mutate(wlvSep = c(NA, ifelse(diff(wlv) > 0, "Marée montante", "Marée descendante")))
  if(all(c("uwnd", "vwnd") %in% names(df))) df <- df %>% mutate(Wspd = sqrt(uwnd^2 + vwnd^2), Wdir = atan2(vwnd, uwnd), Wdirdeg = (90 - Wdir * (180/pi)) %% 360)
  
  if(all(c("dp", "Cdirdeg") %in% names(df))) {
    df <- df %>%
      mutate(
        dpAngleDiff = case_when(dp>0 & dp<180 ~ dp+180, dp>180 & dp<360 ~ dp-180, dp==0 | dp==360 ~ 180, TRUE ~ 0),
        angleDiff = (dpAngleDiff - Cdirdeg) %% 360
      )
  }
  
  if("dp" %in% names(df)) {
    df <- df %>%
      mutate(
        dpCardinale = case_when(
          dp >= 337.5 | dp < 22.5 ~ "N", dp >= 22.5 & dp < 67.5 ~ "NE", dp >= 67.5 & dp < 112.5 ~ "E",
          dp >= 112.5 & dp < 157.5 ~ "SE", dp >= 157.5 & dp < 202.5 ~ "S", dp >= 202.5 & dp < 247.5 ~ "SO",
          dp >= 247.5 & dp < 292.5 ~ "O", dp >= 292.5 & dp < 337.5 ~ "NO", TRUE ~ NA_character_
        ),
        pvc = case_when(!is.na(dpCardinale) & fj == "flot" ~ paste0("vagues direction ", opposite_dir[dpCardinale], ", courant N"),
                        !is.na(dpCardinale) & fj != "flot" ~ paste0("vagues direction ", opposite_dir[dpCardinale], ", courant S"), TRUE ~ NA_character_),
        pvce = case_when(fje == "étale" ~ paste0("vagues direction ", opposite_dir[dpCardinale], ", courant à l'étale"),
                         fje == "flot" ~ paste0("vagues direction ", opposite_dir[dpCardinale], ", courant N"),
                         fje == "jusant" ~ paste0("vagues direction ", opposite_dir[dpCardinale], ", courant S"), TRUE ~ NA_character_)
      ) %>%
      mutate(dpCardinale = factor(dpCardinale, levels = c("N", "NE", "E", "SE", "S", "SO", "O", "NO")), pvc = as.factor(pvc), pvce = as.factor(pvce))
  }
  
  if("angleDiff" %in% names(df)) {
    df <- df %>%
      mutate(
        dirvc = case_when(angleDiff >= 140 & angleDiff <= 225 ~ "contre-courant",
                          (angleDiff >= 315 & angleDiff <= 360) | (angleDiff >= 0 & angleDiff <= 45) ~ "co-courant", TRUE ~ "perpendiculaire"),
        dirvc = factor(dirvc, levels = c("co-courant", "perpendiculaire", "contre-courant"))
      )
  }
  return(df)
}

create_bivariate_plot <- function(df, var_x, var_y, facet_var = "none", bins=70, show_q=FALSE, q_vals="") {
  name_x <- names(var_labels)[var_labels == var_x]
  name_y <- names(var_labels)[var_labels == var_y]
  
  df_clean <- df[!is.na(df[[var_x]]) & !is.na(df[[var_y]]), ]
  if(nrow(df_clean) == 0) return(ggplot() + annotate("text", x=0, y=0, label="Pas de données"))

  x_range <- range(df_clean[[var_x]])
  y_range <- range(df_clean[[var_y]])
  x_breaks <- seq(x_range[1], x_range[2], length.out = bins + 1)
  y_breaks <- seq(y_range[1], y_range[2], length.out = bins + 1)
  x_mids <- x_breaks[-length(x_breaks)] + diff(x_breaks)/2
  y_mids <- y_breaks[-length(y_breaks)] + diff(y_breaks)/2
  
  df_clean$x_bin_idx <- findInterval(df_clean[[var_x]], x_breaks, all.inside = TRUE)
  df_clean$y_bin_idx <- findInterval(df_clean[[var_y]], y_breaks, all.inside = TRUE)
  
  grp_cols <- c("x_bin_idx", "y_bin_idx")
  if(facet_var != "none") grp_cols <- c(grp_cols, facet_var)
  
  df_agg <- df_clean %>%
    group_by(across(all_of(grp_cols))) %>%
    summarise(Count = n(), .groups = "drop") %>%
    mutate(x_coord = x_mids[x_bin_idx], y_coord = y_mids[y_bin_idx])

  p <- ggplot() + 
    geom_tile(data = df_agg, aes(x = x_coord, y = y_coord, fill = Count)) + 
    scale_fill_viridis_c(name="Compte") + 
    theme_bw() +
    ggtitle(paste("Densité bivariée :", name_x, "vs", name_y))
  
  if(facet_var != "none") p <- p + facet_wrap(as.formula(paste("~", facet_var)))
  
  if(show_q && q_vals != "") {
    qs <- as.numeric(trimws(strsplit(q_vals, ",")[[1]]))
    if(!any(is.na(qs))) {
      p <- p + new_scale_color()
      p <- p + geom_quantile(
        data = df_clean,
        aes_string(x = var_x, y = var_y, color = "factor(after_stat(quantile))"),
        quantiles = qs, formula = y ~ x, size = 0.8
      ) + scale_color_viridis_d(option = "plasma", name = "Quantiles", end = 0.9) 
    }
  }
  p + labs(x = name_x, y = name_y)
}

# --- DANS LA SECTION DES FONCTIONS GLOBALES ---

create_histogram_plot <- function(df, var, stats_to_show = character(0), facet_var = "none", show_q=FALSE, q_vals="") {
  name_v <- names(var_labels)[var_labels == var]
  p <- ggplot(df, aes_string(x = var)) + 
    geom_histogram(fill = var_colors[var], alpha = 0.7, bins = 30) +
    theme_bw() +
    ggtitle(paste("Distribution de :", name_v))
  
  compute_stats <- function() {
    qs <- numeric(0)
    if(show_q && q_vals != "") qs <- as.numeric(trimws(strsplit(q_vals, ",")[[1]]))
    grp_cols <- if(facet_var != "none") facet_var else character(0)
    
    stats_df <- df %>%
      group_by(across(all_of(grp_cols))) %>%
      reframe(
        Type = c(
          if("basic" %in% stats_to_show) "Moyenne" else NULL,
          if("min" %in% stats_to_show) "Min" else NULL,   # AJOUT
          if("max" %in% stats_to_show) "Max" else NULL,   # AJOUT
          if(length(qs) > 0) paste0("Q", qs*100, "%") else NULL
        ),
        Value = c(
          if("basic" %in% stats_to_show) mean(!!sym(var), na.rm=TRUE) else NULL,
          if("min" %in% stats_to_show) min(!!sym(var), na.rm=TRUE) else NULL,   # AJOUT
          if("max" %in% stats_to_show) max(!!sym(var), na.rm=TRUE) else NULL,   # AJOUT
          if(length(qs) > 0) quantile(!!sym(var), probs=qs, na.rm=TRUE) else NULL
        )
      )
    return(stats_df)
  }
  
  # Modification de la condition pour inclure min ou max
  if(length(stats_to_show) > 0 || (show_q && q_vals != "")) {
    sdata <- compute_stats()
    p <- p + 
      geom_vline(data = sdata, aes(xintercept = Value, color = Type), linetype = "dashed", size = 0.6) +
      scale_color_viridis_d(option = "plasma", name = "Statistiques", end = 0.9) + 
      geom_text(data = sdata, aes(x = Value, y = Inf, label = round(Value, 2), color = Type), 
                angle = 90, vjust = -0.5, hjust = 1.1, size = 3, show.legend = FALSE)
  }
  if (facet_var != "none") p <- p + facet_wrap(as.formula(paste("~", facet_var)))
  p + labs(x = name_v, y = "Fréquence")
}

create_timeseries_plot <- function(df, var, ts_agg, quantile_values = NULL, facet_var = "none", show_q=FALSE, q_vals="") {
  name_v <- names(var_labels)[var_labels == var]
  
  if (ts_agg == "none") {
    df$date <- as.Date(with(df, paste(annee, mois, jour, sep = "-")))
    p <- ggplot(df, aes(x = date, y = !!sym(var))) + 
      geom_line(color = var_colors[var], alpha = 0.5, size=0.3) +
      theme_bw() +
      ggtitle(paste("Série temporelle :", name_v))
    
    if (facet_var != "none") p <- p + facet_wrap(as.formula(paste("~", facet_var)))
    
    if(show_q && q_vals != "") {
      qs <- as.numeric(trimws(strsplit(q_vals, ",")[[1]]))
      if(!any(is.na(qs))) {
        p <- p + geom_quantile(quantiles = qs, formula = y ~ bs(x, df=5), 
                               aes(color = factor(after_stat(quantile))), size=0.6) +
          scale_color_viridis_d(option = "plasma", name = "Quantiles", end = 0.9) 
      }
    }
    p <- p + labs(x = "Date", y = name_v)
    
  } else {
    group_var <- if (grepl("month", ts_agg)) "mois" else "annee"
    x_lab <- if(group_var=="mois") "Mois" else "Année"
    
    if(grepl("quantile", ts_agg)) {
      quants <- as.numeric(trimws(strsplit(quantile_values, ",")[[1]]))
      df_agg <- df %>% 
        group_by(!!sym(group_var), across(all_of(if(facet_var!="none") facet_var else character(0)))) %>%
        reframe(quantile_type = paste0("Q", round(quants * 100, 0), "%"),
                value = quantile(!!sym(var), probs = quants, na.rm = TRUE)) %>% ungroup()
      
      p <- ggplot(df_agg, aes(x = !!sym(group_var), y = value, color = quantile_type)) + 
        geom_line(size=0.8) + geom_point() +
        scale_color_viridis_d(option = "plasma", name = "Quantiles", end = 0.9) + 
        theme_bw() + ggtitle(paste("Evolution temporelle (quantiles) :", name_v))
      
    } else {
      type_agg <- if(grepl("mean", ts_agg)) "Moyenne" else "Maximum"
      df_agg <- df %>% 
        group_by(!!sym(group_var), across(all_of(if(facet_var!="none") facet_var else character(0)))) %>%
        summarise(value = if(grepl("mean", ts_agg)) mean(!!sym(var), na.rm=T) else max(!!sym(var), na.rm=T), .groups="drop")
      
      p <- ggplot(df_agg, aes(x = !!sym(group_var), y = value)) + 
        geom_line(color=var_colors[var], size=0.8) + geom_point(color=var_colors[var]) +
        theme_bw() + ggtitle(paste("Evolution temporelle (", type_agg, ") :", name_v))
    }
    
    if(facet_var != "none") p <- p + facet_wrap(as.formula(paste("~", facet_var)))
    p <- p + labs(x = x_lab, y = name_v)
  }
  return(p)
}

create_wind_rose <- function(data, direction_var, speed_var, n_cuts = 5, color_scale = "Blues", facet_var = "none") {
  if(!direction_var %in% names(data) || !speed_var %in% names(data)) return(NULL)
  
  breaks <- seq(0, max(data[[speed_var]], na.rm = TRUE), length.out = n_cuts + 1)
  labels <- sprintf("%.2f-%.2f", breaks[-length(breaks)], breaks[-1])
  
  base_data <- data %>% 
    mutate(
      speed_cuts = cut(!!sym(speed_var), breaks = breaks, labels = labels, include.lowest = TRUE),
      direction_bins = cut(!!sym(direction_var), breaks = seq(0, 360, by = 20), include.lowest = T, labels = seq(10, 350, by = 20))
    )
  
  create_single_rose <- function(d) {
    d %>% group_by(direction_bins, speed_cuts) %>% summarise(count = n(), .groups = "drop") %>%
      mutate(direction_bins = as.numeric(as.character(direction_bins))) %>%
      plot_ly(r = ~count, theta = ~direction_bins, color = ~speed_cuts, colors = color_scale, type = "barpolar") %>%
      layout(polar = list(angularaxis = list(direction = "clockwise", rotation = 90)), showlegend = FALSE)
  }
  
  if (facet_var != "none") {
    facet_data <- split(base_data, base_data[[facet_var]])
    plots <- lapply(names(facet_data), function(nm) create_single_rose(facet_data[[nm]]) %>% layout(title = nm))
    subplot(plots, nrows = ceiling(sqrt(length(plots))), shareX = FALSE, shareY = FALSE)
  } else { 
    create_single_rose(base_data) %>% layout(showlegend = TRUE, title = paste("Rose :", names(var_labels)[var_labels == direction_var]))
  }
}

create_quantile_plot <- function(data, x_var, y_var, quantiles, df = 5, conf_level = 0.95, separation = "none", bins = 40) {
  x_seq <- seq(min(data[[x_var]], na.rm=T), max(data[[x_var]], na.rm=T), length.out = 100)
  all_preds <- data.frame()
  fit_predict <- function(sub, grp) {
    res <- data.frame()
    for (q in quantiles) {
        fit <- rq(as.formula(paste(y_var, "~ bs(", x_var, ", df =", df, ")")), tau = q, data = sub)
        pred <- predict(fit, newdata = setNames(data.frame(x_seq), x_var), interval = "confidence", level = conf_level)
        res <- rbind(res, data.frame(x_val = x_seq, fit = pred[,1], lower=pred[,2], upper=pred[,3], quantile=factor(sprintf("Q%.0f%%", q*100)), separation=grp))
    }
    if(nrow(res) > 0 && separation != "none") names(res)[6] <- separation
    res
  }
  if(separation == "none") all_preds <- fit_predict(data, "Overall")
  else for(lvl in unique(data[[separation]])) all_preds <- rbind(all_preds, fit_predict(subset(data, data[[separation]]==lvl), lvl))
  
  p <- ggplot() + 
    geom_bin2d(data=data, aes_string(x=x_var, y=y_var), bins=bins) + 
    scale_fill_viridis_c(name="Densité") + 
    new_scale_fill() +
    ggtitle(paste("Régression quantile :", names(var_labels)[var_labels == y_var], "vs", names(var_labels)[var_labels == x_var]))
  
  if(nrow(all_preds) > 0) {
    p <- p + 
      geom_ribbon(data=all_preds, aes_string(x="x_val", ymin="lower", ymax="upper", fill="quantile"), alpha=0.5) +
      geom_line(data=all_preds, aes_string(x="x_val", y="fit", color="quantile"), size=0.8) +
      scale_color_viridis_d(option = "plasma", name = "Quantiles", end = 0.9) +
      scale_fill_viridis_d(option = "plasma", name = "IC 95%", end = 0.9)
  }
  p <- p + theme_bw()
  if(separation != "none") p <- p + facet_wrap(as.formula(paste("~", separation)))
  p + labs(x = names(var_labels)[var_labels == x_var], y = names(var_labels)[var_labels == y_var])
}

heatmap_with_progress <- function(data, xvar, yvar, statvar, threshold, nbins = 10) {
  x_breaks <- seq(min(data[[xvar]], na.rm=T), max(data[[xvar]], na.rm=T), length.out = nbins)
  y_breaks <- seq(min(data[[yvar]], na.rm=T), max(data[[yvar]], na.rm=T), length.out = nbins)
  data$x_bin <- cut(data[[xvar]], breaks = x_breaks, include.lowest = TRUE)
  data$y_bin <- cut(data[[yvar]], breaks = y_breaks, include.lowest = TRUE)
  df_freq <- as.data.frame(table(data$x_bin, data$y_bin)); names(df_freq) <- c("x_interval", "y_interval", "Frequency")
  df_stats <- data %>% group_by(x_bin, y_bin) %>% summarise(Mean = mean(!!sym(statvar), na.rm=T), Max = max(!!sym(statvar), na.rm=T), Prop = mean(!!sym(statvar) > threshold, na.rm=T)*100, .groups="drop") %>% rename(x_interval=x_bin, y_interval=y_bin)
  df_freq <- merge(df_freq, df_stats, by = c("x_interval", "y_interval"), all.x = TRUE)
  df_freq$Label <- with(df_freq, paste0("n=", Frequency, "\nAvg:", round(Mean,1), "\n>", threshold, ":", round(Prop,1), "%"))
  
  ggplot(df_freq, aes(x=x_interval, y=y_interval, fill=Frequency)) + 
    geom_tile(color="white") + 
    geom_text(aes(label=Label), size=3) +
    scale_fill_gradient(low="white", high="red") + 
    theme_minimal() + 
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle(paste("Heatmap Fréquentielle :", names(var_labels)[var_labels == statvar])) +
    labs(x = names(var_labels)[var_labels == xvar], y = names(var_labels)[var_labels == yvar])
}

fit_gev_model <- function(data, response, covariates, separation="none", spline_params=list()) {
  cov_terms <- character(0)
  if (!is.null(covariates) && length(covariates) > 0) {
    cov_terms <- sapply(covariates, function(cov) {
      bs_val <- if(!is.null(spline_params[[cov]]$bs)) spline_params[[cov]]$bs else "cs"
      k_val <- if(!is.null(spline_params[[cov]]$k)) spline_params[[cov]]$k else 5
      sprintf("s(%s, bs='%s', k=%d %s)", cov, bs_val, k_val, 
              if(separation != "none") paste0(", by=", separation) else "")
    })
  }

  rhs_loc <- if(length(cov_terms) > 0) paste(cov_terms, collapse="+") else "1"
  if(separation != "none" && length(cov_terms) > 0) rhs_loc <- paste(separation, "+", rhs_loc)
  else if(separation != "none") rhs_loc <- separation
  
  fmla <- list(
    location = as.formula(paste(response, "~", rhs_loc)),
    scale = if(separation != "none") as.formula(paste("~", separation)) else ~1,
    shape = ~1
  )

  evgam(fmla, data=data, family="gev")
}

fit_gpd_model <- function(data, covariates, threshold, separation="none", stationary="stationary", threshold_var="cycle", spline_params=list()) {
  cov_terms <- character(0)
  if (!is.null(covariates) && length(covariates) > 0) {
    cov_terms <- sapply(covariates, function(cov) {
      bs_val <- if(!is.null(spline_params[[cov]]$bs)) spline_params[[cov]]$bs else "cs"
      k_val <- if(!is.null(spline_params[[cov]]$k)) spline_params[[cov]]$k else 5
      sprintf("s(%s, bs='%s', k=%d %s)", cov, bs_val, k_val, 
              if(separation != "none") paste0(", by=", separation) else "")
    })
  }
  if(stationary!="stationary") {
      bs_val_thresh <- "cc"
      k_val_thresh <- 5
      if(!is.null(spline_params[[threshold_var]])) {
          bs_val_thresh <- spline_params[[threshold_var]]$bs
          k_val_thresh <- spline_params[[threshold_var]]$k
      }
      cov_terms <- c(cov_terms, sprintf("s(%s, bs='%s', k=%d %s)", threshold_var, bs_val_thresh, k_val_thresh, if(separation!="none") paste0(", by=", separation) else ""))
  }
  rhs_scale <- if(length(cov_terms) > 0) paste(cov_terms, collapse="+") else "1"
  evgam(list(scale=as.formula(paste("excess ~", rhs_scale)), shape=~1), data=data, family="gpd")
}

compute_nonstationary_threshold <- function(data, var, quantile, time_var) {
  
  fmla_string <- paste(var, "~ s(", time_var, ", bs = 'cc', k = 5)")
  
  fmla_list <- list(as.formula(fmla_string), ~ 1)
  
  fit <- evgam(
    formula = fmla_list, 
    data = data, 
    family = "ald", 
    ald.args = list(tau = quantile)
  )
  
  data$threshold <- predict(fit, newdata = data)$location
  
  return(data)
}

plot_gev_rl_stationary <- function(model, return_periods = c(2, 5, 10, 20, 50, 100)) {
  coefs <- coef(model)
  mu <- as.numeric(coefs[grep("location", names(coefs))])[1]
  sig <- exp(as.numeric(coefs[grep("scale", names(coefs))])[1])
  xi <- as.numeric(coefs[grep("shape", names(coefs))])[1]
  
  rls <- sapply(return_periods, function(T_val) { 
    p <- 1 - 1/T_val
    if(abs(xi)<1e-6) mu - sig*log(-log(p)) else mu - (sig/xi)*(1-(-log(p))^(-xi)) 
  })
  
  df_rl <- data.frame(Period = return_periods, Level = rls)
  
  seq_T <- seq(1.1, max(return_periods)+10, length.out=100)
  seq_rl <- sapply(seq_T, function(T_val) {
    p <- 1 - 1/T_val
    if(abs(xi)<1e-6) mu - sig*log(-log(p)) else mu - (sig/xi)*(1-(-log(p))^(-xi))
  })
  df_curve <- data.frame(Period = seq_T, Level = seq_rl)
  
  ggplot() +
    geom_line(data=df_curve, aes(x=Period, y=Level), color="blue") +
    geom_point(data=df_rl, aes(x=Period, y=Level), color="red", size=3) +
    scale_x_log10(breaks=return_periods, labels=return_periods) +
    theme_bw() +
    labs(x="Période de retour (années) [Echelle Log]", y="Niveau de retour", title="GEV - Niveaux de retour (Stationnaire)")
}

plot_gpd_rl_stationary <- function(model, threshold, n_exc, total_years, return_periods = c(2, 5, 10, 20, 50, 100)) {
  lambda <- n_exc / total_years
  coefs <- coef(model)
  scale_val <- exp(as.numeric(coefs[grep("scale", names(coefs))])[1])
  shape_val <- as.numeric(coefs[grep("shape", names(coefs))])[1]
  u <- threshold
  
  rls <- sapply(return_periods, function(T_val) { 
    if(abs(shape_val)<1e-6) u + scale_val * log(T_val * lambda) else u + (scale_val/shape_val) * ( (T_val * lambda)^shape_val - 1 ) 
  })
  
  df_rl <- data.frame(Period = return_periods, Level = rls)
  
  seq_T <- seq(1.1, max(return_periods)+10, length.out=100)
  seq_rl <- sapply(seq_T, function(T_val) {
    if(abs(shape_val)<1e-6) u + scale_val * log(T_val * lambda) else u + (scale_val/shape_val) * ( (T_val * lambda)^shape_val - 1 )
  })
  df_curve <- data.frame(Period = seq_T, Level = seq_rl)
  
  ggplot() +
    geom_line(data=df_curve, aes(x=Period, y=Level), color="darkgreen") +
    geom_point(data=df_rl, aes(x=Period, y=Level), color="orange", size=3) +
    scale_x_log10(breaks=return_periods, labels=return_periods) +
    theme_bw() +
    labs(x="Période de retour (années) [Echelle Log]", y="Niveau de retour", title="GPD - Niveaux de retour (Stationnaire)")
}

ui <- fluidPage(
  navbarPage("Analyse par point Resourcecode",
    
    tabPanel("Sélection des données",
      fluidRow(
        column(8, leafletOutput("map", height = "600px")),
        column(4,
          checkboxInput("show_grid_points", "Afficher tous les points de la grille", value = FALSE),
          numericInput("manual_node", "Saisir n° de nœud:", value = 158348, min=1), 
          
          fluidRow(
            column(6, numericInput("api_start_year", "Année début:", value = 1994, min=1994, max=2020)),
            column(6, numericInput("api_end_year", "Année fin:", value = 1998, min=1994, max=2020))
          ),
          selectInput("selected_vars", "Sélectionner les variables:", choices = names(var_groups), multiple = TRUE, selected = c("Hauteur significative des vagues (m)", "Période de pic (s)", "Direction du pic des vagues (°)")),
          actionButton("load_data", "Charger les données"),
          verbatimTextOutput("loading_status")
        )
      )
    ),
    
    tabPanel("Analyse",
      sidebarLayout(
        sidebarPanel(
          sliderInput("year_range", "Période:", min=1995, max=2021, value=c(1995, 1998), step=1, sep=""),
          radioButtons("analysis_type", "Type d'analyse:",
                       choices = c("Descriptive"="descriptive", "GEV"="gev", 
                                   "GPD"="gpd", "Régression quantile"="quantreg", "Régression linéaire"="linreg",
                                  "Contours environnementaux" = "contours")),
          
          conditionalPanel("input.analysis_type == 'descriptive'",
            radioButtons("plot_type", "Type de graphique:", choices=c("Univarié"="Univariate", "Bivarié"="Bivariate", "Analyse fréquentielle"="Frequentielle")),
            
            hr(), strong("Filtre (optionnel)"),
            selectInput("filter_var", "Filtrer par variable :", choices=c("Aucun"="none", names(vars_cat), names(vars_all_num))),
            uiOutput("filter_ui"),
            hr(),
            
            conditionalPanel("input.plot_type == 'Univariate'",
              selectInput("var1", "Variable:", choices=names(vars_all_num)), 
              radioButtons("viz_type", "Visualisation:", choices=c("Histogramme"="hist", "Série temporelle"="ts", "Rose"="rose")),
              
            conditionalPanel("input.viz_type=='hist'", 
              checkboxGroupInput("stats_to_show", "Stats:", 
                                  choices=c("Moyenne"="basic", "Minimum"="min", "Maximum"="max"), 
                                  selected="basic"),
              
              checkboxInput("show_quantiles", "Afficher quantiles personnalisés", FALSE),
              conditionalPanel("input.show_quantiles", textInput("desc_quantiles", "Quantiles (ex: 0.05, 0.95)", "0.05, 0.95"))
            ),  
                        conditionalPanel("input.viz_type=='ts'", 
                               selectInput("ts_agg", "Agrégation:", choices=c("Aucune"="none", "Moyenne mois"="mean_month", "Max an"="max_year")),
                               conditionalPanel("input.ts_agg == 'none'", 
                                                 checkboxInput("show_quantiles_ts", "Afficher courbes quantiles", FALSE),
                                                 conditionalPanel("input.show_quantiles_ts", textInput("desc_quantiles_ts", "Quantiles:", "0.1, 0.5, 0.9"))
                               )
              ),
              conditionalPanel("input.viz_type=='rose'", 
                               selectInput("rose_var", "Variable vitesse:", choices=names(vars_num_pure)),
                               sliderInput("n_speed_cuts", "Classes:", 3, 8, 5)
              )
            ),
            conditionalPanel("input.plot_type == 'Bivariate'",
              selectInput("var_x", "X:", choices=names(vars_all_num)), selectInput("var_y", "Y:", choices=names(vars_all_num)),
              checkboxInput("show_correlation", "Afficher corrélation", FALSE),
              checkboxInput("show_quantiles_bivar", "Afficher Quantiles Y|X", FALSE),
              conditionalPanel("input.show_quantiles_bivar", textInput("desc_quantiles_bivar", "Quantiles:", "0.1, 0.5, 0.9"))
            ),
            conditionalPanel("input.plot_type == 'Frequentielle'",
              selectInput("var_x_freq", "X:", choices=names(vars_all_num)), selectInput("var_y_freq", "Y:", choices=names(vars_all_num)),
              selectInput("statvar_freq", "Variable stats:", choices=names(vars_all_num)), numericInput("threshold_freq", "Seuil:", 3)
            ),
            selectInput("facet_var", "Grouper (facets):", choices=c("Aucun"="none", vars_cat))
          ),
          
          conditionalPanel("input.analysis_type == 'gev'",
            selectInput("gev_response", "Variable:", choices=names(vars_all_num)),
            selectInput("gev_covariates", "Covariables:", choices=names(vars_all_num), multiple=TRUE),
            uiOutput("gev_spline_controls"),
            selectInput("gev_separation", "Séparation:", choices=c("Aucune"="none", vars_cat)),
            actionButton("fit_gev", "Ajuster GEV")
          ),
          
          conditionalPanel("input.analysis_type == 'gpd'",
             selectInput("gpd_response", "Variable:", choices=names(vars_all_num)),
             radioButtons("threshold_type", "Seuil:", choices=c("Stationnaire"="stationary", "Non stationnaire"="nonstationary")),
             conditionalPanel("input.threshold_type=='stationary'", numericInput("threshold", "Valeur seuil:", 2)),
             conditionalPanel("input.threshold_type=='nonstationary'", selectInput("var_for_ns", "Variable seuil :", choices=names(vars_all_num), selected="Cycle annuel"), sliderInput("quantile_threshold", "Quantile:", 0.5, 0.99, 0.95)),
             sliderInput("decluster_r", "Dégroupage (deux tempêtes sont séparées d'au moins X heures):", min=1, max=100, value=24),
             
             actionButton("check_threshold", "Visualiser dépassements de seuil", class="btn-info"),
             br(), br(),
             
             conditionalPanel("output.gpd_data_ready",
               h4("Paramètres du modèle"),
               selectInput("gpd_covariates", "Covariables:", choices=names(vars_all_num), multiple=TRUE),
               uiOutput("gpd_spline_controls"),
               selectInput("gpd_separation", "Séparation:", choices=c("Aucune"="none", vars_cat)),
               actionButton("fit_gpd", "Ajuster modèle GPD", class="btn-success")
             )
          ),
          
          conditionalPanel("input.analysis_type == 'quantreg'",
            selectInput("qr_y", "Y:", choices=names(vars_all_num)), selectInput("qr_x", "X:", choices=names(vars_all_num)),
            textInput("quantile_values", "Quantiles:", "0.1, 0.5, 0.9"), 
            selectInput("qr_separation", "Séparation:", choices=c("Aucune"="none", vars_cat)),
            actionButton("fit_quantreg", "Ajuster")
          ),
          conditionalPanel("input.analysis_type == 'linreg'",
            selectInput("linreg_response", "Y:", choices=names(vars_all_num)), selectInput("linreg_predictors", "X:", choices=names(vars_all_num), multiple=T),
            actionButton("fit_linreg", "Ajuster")
          ),

conditionalPanel("input.analysis_type == 'contours'",
    h4("Contours Environnementaux"),
    
    selectInput("contour_method", "Méthode :", 
                choices = c("Direct IFORM" = "direct",
                            "IFORM (Weibull-Lognormal)" = "wl",
                            "IFORM (Heffernan & Tawn)" = "ht04")),
    
    selectInput("contour_var_x", "Variable X (ex: Hs):", choices=names(vars_all_num), selected="Hauteur significative des vagues (m)"),
    selectInput("contour_var_y", "Variable Y (ex: Tp):", choices=names(vars_all_num), selected="Période de pic (s)"),
    
    textInput("contour_periods", "Périodes de retour (années):", "10, 50, 100"),
    
    # --- AJOUT: PARAMÈTRES BOOTSTRAP ---
    hr(),
    checkboxInput("contour_boot_active", "Incertitude d'échantillonnage (bootstrap)", value = FALSE),
    conditionalPanel("input.contour_boot_active",
       numericInput("contour_boot_n", "Nombre de réplications (B):", value = 50, min = 10, max = 1000),
       helpText("Attention : Le calcul peut être long.")
    ),
    hr(),
    # -----------------------------------

    # Paramètres spécifiques DIRECT
    conditionalPanel("input.contour_method == 'direct'",
      sliderInput("contour_spacing", "Espacement (précision):", min=0.05, max=0.5, value=0.1, step=0.05),
      sliderInput("contour_prob", "Seuil Probabilité:", min=0.01, max=0.2, value=0.1, step=0.01)
    ),
    
    # Paramètres spécifiques HT04
    conditionalPanel("input.contour_method == 'ht04'",
       sliderInput("ht04_threshold", "Seuil quantile (u):", min=0.5, max=0.95, value=0.8, step=0.05),
       helpText("Seuil pour le modèle de dépendance conditionnelle.")
    ),
    
    # Paramètres communs
    numericInput("contour_declust", "Dégroupage (D-IFORM) / durée état de mer (IFORM):", value=1, min=1),
    
    actionButton("fit_contour_btn", "Calculer les contours", class="btn-primary")
),
                   hr(),
          sliderInput("plotHeight", "Hauteur plot:", 300, 1500, 600)
        ),
        
        mainPanel(
          conditionalPanel("input.analysis_type == 'descriptive'", uiOutput("dynamic_descriptive_plot"), verbatimTextOutput("stats_text"), verbatimTextOutput("corr_text")),
          
          conditionalPanel("input.analysis_type == 'gev'", 
                           uiOutput("gev_results_ui")
          ),
          
          conditionalPanel("input.analysis_type == 'gpd'", 
             h4("Étape 1: Dépassements de seuil et clusters"), 
             uiOutput("gpd_exceedances_2"),
             uiOutput("gpd_results_ui")
          ),
          
          conditionalPanel("input.analysis_type == 'quantreg'", uiOutput("quantreg_plot_X")),
          conditionalPanel("input.analysis_type == 'linreg'", plotOutput("linreg_plot"), verbatimTextOutput("linreg_summary")),
          conditionalPanel("input.analysis_type == 'contours'", 
              plotlyOutput("contour_plot", height="600px"),
              verbatimTextOutput("contour_info")
          ),
        )
      )
    ),

    tabPanel("Analyse spectrale",
       sidebarLayout(
         sidebarPanel(
           width = 3,
           h3("Spectres de vagues"),
           helpText("1. Sélectionnez un point bleu sur la carte ci-contre."),
           helpText("2. Choisissez la date et le type."),
           
           hr(),
           htmlOutput("spec_selected_point_ui"),
           
           hr(),
           radioButtons("spec_type", "Type de spectre :", 
                        choices = c("1D (Fréquence)" = "1d", "2D (Fréq-Direction)" = "2d"), 
                        inline = TRUE),
           
           textInput("spec_date", "Date (YYYY-MM-DD HH:00:00)", value = "1994-02-10 06:00:00"),
           helpText("Format: 1994-02-10 06:00:00. Pour minuit : 1994-02-10"),
           
           actionButton("plot_spectrum_btn", "Visualiser le spectre", class = "btn-primary", width = "100%")
         ),
         mainPanel(
           width = 9,
           leafletOutput("map_spectral", height = "350px"),
           br(),
           plotOutput("plot_spectrum_viz", height = "500px")
         )
       )
    )
  )
)

server <- function(input, output, session) {
  
  output$map <- renderLeaflet({ 
    leaflet() %>% 
      addTiles() %>% 
      setView(-4, 48.5, 7)
  })
  
  observe({
    proxy <- leafletProxy("map")
    
    if (input$show_grid_points) {
      proxy %>% 
        addCircleMarkers(data=rscd_field, 
                         lng=~longitude, lat=~latitude, 
                         layerId=~node, 
                         radius=3, color="blue", opacity = 0.8, fillOpacity = 0.5,
                         popup=~paste("Node:", node),
                         group = "grid_points") 
    } else {
      proxy %>% clearGroup("grid_points")
    }
  })
  
  selected_node <- reactiveVal(158348)
  
  observeEvent(input$map_marker_click, { 
    selected_node(input$map_marker_click$id) 
    updateNumericInput(session, "manual_node", value = input$map_marker_click$id)
  })
  
  observeEvent(input$manual_node, {
    req(input$manual_node)
    if (input$manual_node %in% rscd_field$node) {
      selected_node(input$manual_node)
    } else {
      showNotification("Numéro de nœud invalide.", type = "warning", duration = 3)
    }
  })
  
  observe({ 
    req(selected_node()); 
    if (selected_node() %in% rscd_field$node) {
      leafletProxy("map") %>% 
        clearGroup("hl") %>% 
        addCircleMarkers(data=rscd_field[rscd_field$node==selected_node(),], 
                         lng=~longitude, lat=~latitude, group="hl", 
                         radius=8, color="yellow", fillColor="yellow", 
                         fillOpacity=1, opacity=1, weight=1.5)
    }
  })
  
  loaded_data <- reactiveVal(NULL)
  
  load_parameters_logic <- function() {
    output$loading_status <- renderText("Chargement en cours...")
    
      sel <- input$selected_vars; if(is.null(sel)) sel <- c("Hauteur significative des vagues (m)", "Période de pic (s)", "Direction du pic des vagues (°)")
      req_vars <- get_selected_variables(sel)
      final_params <- intersect(req_vars, c("tp", resourcecodedata::rscd_variables$name))
      
      y_start <- if(is.null(input$api_start_year)) 1994 else input$api_start_year
      y_end <- if(is.null(input$api_end_year)) 1998 else input$api_end_year
      start_dt <- as.POSIXct(paste0(y_start,"-01-01 00:00:00"), tz="UTC")
      end_dt <- as.POSIXct(paste0(y_end,"-12-31 23:00:00"), tz="UTC")

      withProgress(message = 'Chargement', value = 0, {
        incProgress(0.3, detail = "API")
        raw <- get_parameters(node = selected_node(), start=start_dt, end=end_dt, parameters=final_params)
        incProgress(0.3, detail = "Transform")
        dat <- transform_dataset(raw, final_params)
      })
      
      loaded_data(dat)
      output$loading_status <- renderText(paste("Données chargées : Node", selected_node()))
  }
  
  observeEvent(input$load_data, { load_parameters_logic() })
  observe({ isolate({ if(is.null(loaded_data())) load_parameters_logic() }) }, priority=10)
  
  output$filter_ui <- renderUI({
    req(input$filter_var != "none", loaded_data())
    v <- input$filter_var 
    data <- loaded_data()
    
    if (is.numeric(data[[v]])) {
      rng <- range(data[[v]], na.rm=TRUE)
      sliderInput("filter_range", paste("Plage", names(var_labels)[var_labels==v]), min=floor(rng[1]), max=ceiling(rng[2]), value=c(floor(rng[1]), ceiling(rng[2])))
    } else {
      levs <- sort(unique(as.character(data[[v]])))
      checkboxGroupInput("filter_cats", paste("Catégories", names(var_labels)[var_labels==v]), choices=levs, selected=levs)
    }
  })
  
  current_data <- reactive({
    req(loaded_data())
    df <- loaded_data()[loaded_data()$annee >= input$year_range[1] & loaded_data()$annee <= input$year_range[2], ]
    
    if(input$filter_var != "none") {
      v <- input$filter_var
      if(is.numeric(df[[v]])) {
        req(input$filter_range)
        df <- df[df[[v]] >= input$filter_range[1] & df[[v]] <= input$filter_range[2], ]
      } else {
        req(input$filter_cats)
        df <- df[df[[v]] %in% input$filter_cats, ]
      }
    }
    df
  })
  
  observeEvent(input$viz_type, {
    req(input$plot_type == "Univariate")
    if(input$viz_type == "rose") {
      updateSelectInput(session, "var1", label="Direction:", choices=names(vars_dir))
    } else {
      updateSelectInput(session, "var1", label="Variable:", choices=names(vars_all_num))
    }
  })
  
  output$facet_selector <- renderUI({ 
    ch <- if(input$viz_type=="ts") get_facet_labels(input$ts_agg) else get_facet_labels() 
    selectInput("facet_var", "Grouper:", choices=c("Aucun"="none", vars_cat)) 
  })
  
  output$main_plot <- renderPlotly({
    req(current_data())
    df <- current_data()
    p <- NULL
    
    if(input$plot_type == "Univariate") {
      v <- var_labels[input$var1]
      if(input$viz_type == "hist") {
        p <- create_histogram_plot(df, v, input$stats_to_show, input$facet_var, show_q=input$show_quantiles, q_vals=input$desc_quantiles)
      } else if(input$viz_type == "ts") {
        p <- create_timeseries_plot(df, v, input$ts_agg, NULL, input$facet_var, show_q=input$show_quantiles_ts, q_vals=input$desc_quantiles_ts)
      } else if(input$viz_type == "rose") {
        p <- create_wind_rose(df, v, var_labels[input$rose_var], input$n_speed_cuts, facet_var=input$facet_var)
      }
    } else if(input$plot_type == "Bivariate") {
      p <- create_bivariate_plot(df, var_labels[input$var_x], var_labels[input$var_y], input$facet_var, show_q=input$show_quantiles_bivar, q_vals=input$desc_quantiles_bivar)
    } else {
      p <- heatmap_with_progress(df, var_labels[input$var_x_freq], var_labels[input$var_y_freq], var_labels[input$statvar_freq], input$threshold_freq)
    }
    if(!is.null(p) && !inherits(p, "plotly")) ggplotly(p + theme_bw()) else p
  })
  
  output$dynamic_descriptive_plot <- renderUI({ plotlyOutput("main_plot", height=paste0(input$plotHeight, "px")) })
  output$corr_text <- renderText({ req(input$show_correlation, input$var_x, input$var_y); paste("Corrélation:", round(cor(current_data()[[var_labels[input$var_x]]], current_data()[[var_labels[input$var_y]]], use="complete.obs"), 3)) })

# --- DANS LE SERVER ---

output$stats_text <- renderPrint({
    req(current_data())
    
    if (input$plot_type == "Univariate" && input$viz_type == "hist") {
      df <- current_data()
      var_code <- var_labels[input$var1]
      var_name <- names(var_labels)[var_labels == input$var1]
      
      # Récupération des choix
      show_mean <- "basic" %in% input$stats_to_show
      show_min <- "min" %in% input$stats_to_show   # AJOUT
      show_max <- "max" %in% input$stats_to_show   # AJOUT
      
      show_quant <- input$show_quantiles && input$desc_quantiles != ""
      qs <- if(show_quant) as.numeric(trimws(strsplit(input$desc_quantiles, ",")[[1]])) else NULL
      
      if(!show_mean && !show_quant && !show_min && !show_max) {
        cat("Aucune statistique sélectionnée.")
        return()
      }
      
      cat("STATISTIQUES DÉTAILLÉES POUR :", var_name, "\n")
      cat("=================================================\n")
      
      print_stats <- function(sub_d, prefix="") {
        if(show_min) cat(prefix, "Minimum :", round(min(sub_d[[var_code]], na.rm=TRUE), 3), "\n") # AJOUT
        if(show_mean) cat(prefix, "Moyenne :", round(mean(sub_d[[var_code]], na.rm=TRUE), 3), "\n")
        if(show_max) cat(prefix, "Maximum :", round(max(sub_d[[var_code]], na.rm=TRUE), 3), "\n") # AJOUT
        
        if(show_quant && !any(is.na(qs))) {
          q_vals <- quantile(sub_d[[var_code]], probs=qs, na.rm=TRUE)
          q_str <- paste(paste0(names(q_vals), ": ", round(q_vals, 3)), collapse=" | ")
          cat(prefix, "Quantiles :", q_str, "\n")
        }
      }
      
      if (input$facet_var == "none") {
        print_stats(df)
      } else {
        facet_col <- input$facet_var
        groups <- sort(unique(as.character(df[[facet_col]])))
        for(g in groups) {
          cat("\n[ GROUPE :", g, "]\n")
          cat("--------------------\n")
          sub_df <- df[df[[facet_col]] == g, ]
          print_stats(sub_df, prefix="  ")
        }
      }
    }
  })

  monthly_maxima <- reactive({
    req(current_data(), input$gev_response)
    v <- var_labels[input$gev_response]
    current_data() %>% group_by(annee, mois) %>% slice_max(!!sym(v), with_ties=F) %>% ungroup() %>% mutate(time=as.Date(sprintf("%d-%02d-01", annee, mois)))
  })
  
  output$gev_spline_controls <- renderUI({
    req(input$gev_covariates)
    lapply(input$gev_covariates, function(v) {
      varn <- var_labels[v]
      fluidRow(
        column(4, strong(v)),
        column(4, selectInput(paste0("gev_bs_", varn), "Base", choices=c('cc','cp','cs','ps','cr','tp'), selected='cs')),
        column(4, numericInput(paste0("gev_k_", varn), "k", value=5, min=3))
      )
    })
  })
  
  gev_fitted_model <- reactiveVal(NULL)
  
  observeEvent(input$fit_gev, {
    req(monthly_maxima())
    
    spline_settings <- list()
    if(!is.null(input$gev_covariates)) {
      for(v in input$gev_covariates) {
        varn <- var_labels[v]
        spline_settings[[varn]] <- list(
          bs = input[[paste0("gev_bs_", varn)]],
          k = input[[paste0("gev_k_", varn)]]
        )
      }
    }
    
    m <- fit_gev_model(monthly_maxima(), var_labels[input$gev_response], var_labels[input$gev_covariates], input$gev_separation, spline_params = spline_settings)
    gev_fitted_model(m)
  })
  
  output$gev_results_ui <- renderUI({
    if(is.null(gev_fitted_model())) return(NULL) 
    
    m <- gev_fitted_model()
    is_stationary <- (length(input$gev_covariates) == 0 && input$gev_separation == "none")
    
    list(
      h4("Maxima mensuels"), plotOutput("gev_monthly_maxima"),
      h4("Diagnostic du modèle"), plotOutput("gev_diagnostics"),
      verbatimTextOutput("gev_summary"),
      
      if(!is_stationary) {
        tagList(h4("Splines (Covariables)"), plotOutput("gev_splines_plot"))
      } else {
        tagList(
          h4("Niveaux de retour (modèle stationnaire)"),
          plotOutput("gev_rl_plot"),
          tableOutput("gev_rl_table")
        )
      }
    )
  })

  output$gev_summary <- renderPrint({ req(gev_fitted_model()); summary(gev_fitted_model()) })

  output$gev_diagnostics <- renderPlot({
    req(gev_fitted_model(), monthly_maxima())
    m <- gev_fitted_model()
    
    data_check <- monthly_maxima()
    obs <- data_check[[var_labels[input$gev_response]]]
    
    preds <- predict(m, newdata=data_check, type="response")
    mu <- preds$location
    sig <- preds$scale
    xi <- preds$shape
    
    term <- 1 + xi * (obs - mu) / sig
    
    res <- ifelse(abs(xi) < 1e-6, 
                  (obs - mu) / sig, 
                  (1 / xi) * log(term))
    
    res <- res[!is.na(res) & is.finite(res)]
    
    n <- length(res)
    p_pos <- (1:n) / (n + 1)
    theo_q <- -log(-log(p_pos)) 
    emp_q <- sort(res)          
    
    df_qq <- data.frame(Theorique = theo_q, Empirique = emp_q)
    
    ggplot(df_qq, aes(x=Theorique, y=Empirique)) +
      geom_point(color="#0073C2", alpha=0.6, size=2) +
      geom_abline(intercept=0, slope=1, color="red", size=1, linetype="dashed") +
      theme_bw() +
      labs(title = "Diagnostic : QQ-Plot des résidus",
           subtitle = "Les points doivent s'aligner sur la ligne rouge",
           x = "Quantiles théoriques",
           y = "Résidus du modèle")
  })
    output$gev_monthly_maxima <- renderPlot({ req(monthly_maxima()); ggplot(monthly_maxima(), aes(x=time, y=!!sym(var_labels[input$gev_response]))) + geom_line() + theme_bw() })
  
  output$gev_splines_plot <- renderPlot({
    req(gev_fitted_model())
    if(length(input$gev_covariates)>0) plot(gev_fitted_model())
  })
  
output$gev_rl_table <- renderTable({
    req(gev_fitted_model(), input$gev_separation)
    m <- gev_fitted_model()
    
    if(length(input$gev_covariates)==0 && input$gev_separation=="none") {
       
       dummy_data <- monthly_maxima()[1, ] 
       preds <- predict(m, newdata=dummy_data, type="response")
       
       mu <- preds$location[1]
       sig <- preds$scale[1]
       xi <- preds$shape[1]  
       
       T_periods <- c(2, 5, 10, 20, 50, 100)
       rls <- sapply(T_periods, function(T_val) { 
         p <- 1 - 1/(T_val * 12) 
         
         if(!is.na(xi) && abs(xi) < 1e-6) {
           mu - sig * log(-log(p))
         } else {
           mu - (sig/xi) * (1 - (-log(p))^(-xi)) 
         }
       })
       data.frame(Retour_annees=as.integer(T_periods), Niveau=rls)
    }
  })

output$gev_rl_plot <- renderPlot({
    req(gev_fitted_model(), input$gev_separation)
    m <- gev_fitted_model()
    
    if(length(input$gev_covariates)==0 && input$gev_separation=="none") {
      
      dummy_data <- monthly_maxima()[1, ]
      preds <- predict(m, newdata=dummy_data, type="response")
      
      mu <- preds$location[1] 
      sig <- preds$scale[1] 
      xi <- preds$shape[1] 
      
      return_periods <- c(2, 5, 10, 20, 50, 100)
      
      rls <- sapply(return_periods, function(T_val) { 
         p <- 1 - 1/(T_val * 12)
         if(!is.na(xi) && abs(xi)<1e-6) mu - sig*log(-log(p)) else mu - (sig/xi)*(1-(-log(p))^(-xi)) 
      })
      df_rl <- data.frame(Period = return_periods, Level = rls)
      
      seq_T <- seq(1.1, 200, length.out=200)
      seq_rl <- sapply(seq_T, function(T_val) {
        p <- 1 - 1/(T_val * 12)
        if(!is.na(xi) && abs(xi)<1e-6) mu - sig*log(-log(p)) else mu - (sig/xi)*(1-(-log(p))^(-xi))
      })
      df_curve <- data.frame(Period = seq_T, Level = seq_rl)
      
      ggplot() +
        geom_line(data=df_curve, aes(x=Period, y=Level), color="#0073C2", size=1) +
        geom_point(data=df_rl, aes(x=Period, y=Level), color="red", size=3) +
        scale_x_log10(breaks=return_periods, labels=return_periods) +
        theme_bw() +
        labs(x="Période de retour (années) [Log]", y="Niveau (m)", title="GEV - Niveaux de retour (Stationnaire / evgam)")
    }
  })

  gpd_prep_data <- reactiveVal(NULL)
  gpd_fitted_model <- reactiveVal(NULL)
  output$gpd_data_ready <- reactive({ !is.null(gpd_prep_data()) })
  outputOptions(output, "gpd_data_ready", suspendWhenHidden=FALSE)
  
  observeEvent(input$check_threshold, {
    req(current_data(), input$gpd_response)
    v <- var_labels[input$gpd_response]
    dat <- current_data()
    
    thr <- input$threshold
    if(input$threshold_type == "nonstationary") {
      dat <- compute_nonstationary_threshold(dat, v, input$quantile_threshold, var_labels[input$var_for_ns])
      thr <- dat$threshold
    } else { dat$threshold <- thr }
    
    exceed_indices <- which(dat[[v]] > dat$threshold)
    n_exc <- length(exceed_indices)
    
    if (n_exc < 10) {
      showNotification(paste0("Trop peu de dépassements (", n_exc, "). Abaissez le seuil."), type = "error", duration = 5)
      gpd_prep_data(NULL); gpd_fitted_model(NULL)
      return()
    }
    
    exc_df <- dat[exceed_indices, ]; exc_df <- exc_df[order(exc_df$time), ]
    time_diffs <- c(Inf, as.numeric(difftime(exc_df$time[-1], exc_df$time[-nrow(exc_df)], units="hours")))
    exc_df$cluster_id <- cumsum(time_diffs > input$decluster_r)
    peaks_df <- exc_df %>% group_by(cluster_id) %>% slice_max(!!sym(v), with_ties=FALSE) %>% ungroup()
    peaks_df$excess <- peaks_df[[v]] - peaks_df$threshold
    
    gpd_prep_data(list(full_data = dat, peaks = peaks_df))
    gpd_fitted_model(NULL)
    showNotification(paste("Analyse préliminaire terminée :", nrow(peaks_df), "clusters identifiés."), type="message")
  })
  
  output$gpd_exc_plot <- renderPlot({
    req(gpd_prep_data(), input$gpd_response)
    v <- var_labels[input$gpd_response]
    L <- gpd_prep_data()
    p <- ggplot(L$full_data, aes(x=time, y=!!sym(v))) + geom_line(color="grey") + theme_bw()
    if(input$threshold_type=="stationary") p <- p + geom_hline(yintercept=input$threshold, color="blue", linetype="dashed", size=1)
    else p <- p + geom_line(aes(y=threshold), color="blue", linetype="dashed", size=1)
    p <- p + geom_point(data=L$peaks, color="red", size=2) + ggtitle(paste("Série temporelle avec seuil et", nrow(L$peaks), "clusters"))
    p
  })
  output$gpd_exceedances_2 <- renderUI({ plotOutput("gpd_exc_plot", height=paste0(input$plotHeight, "px")) })
  
  output$gpd_spline_controls <- renderUI({
      covs <- input$gpd_covariates
      if(input$threshold_type == "nonstationary") covs <- c(covs, input$var_for_ns)
      covs <- unique(covs)
      
      if(length(covs)==0) return(NULL)
      
      lapply(covs, function(v) {
          varn <- var_labels[v]
          fluidRow(
              column(4, strong(v)),
              column(4, selectInput(paste0("gpd_bs_", varn), "Base", choices=c('cc','cp','cs','ps','cr','tp'), selected='cs')),
              column(4, numericInput(paste0("gpd_k_", varn), "k", value=5, min=3))
          )
      })
  })

  observeEvent(input$fit_gpd, {
    req(gpd_prep_data())
    peaks <- gpd_prep_data()$peaks
    
    spline_settings <- list()
    
    covs_to_check <- c(input$gpd_covariates)
    if(input$threshold_type == "nonstationary") covs_to_check <- c(covs_to_check, input$var_for_ns)
    covs_to_check <- unique(covs_to_check)
    
    if(length(covs_to_check) > 0) {
        for(v in covs_to_check) {
            varn <- var_labels[v]
            spline_settings[[varn]] <- list(
                bs = input[[paste0("gpd_bs_", varn)]],
                k = input[[paste0("gpd_k_", varn)]]
            )
        }
    }
    
    m <- fit_gpd_model(peaks, var_labels[input$gpd_covariates], input$threshold, input$gpd_separation, stationary=input$threshold_type, threshold_var=var_labels[input$var_for_ns], spline_params = spline_settings)
    gpd_fitted_model(m)
  })
  
output$gpd_results_ui <- renderUI({
    if(is.null(gpd_fitted_model())) return(NULL)
    
    is_stationary <- (length(input$gpd_covariates) == 0 && input$gpd_separation == "none" && input$threshold_type == "stationary")
    
    list(
      hr(),
      h4("Étape 2: Résultats du modèle GPD"),
      
      h4("Diagnostic du modèle"),
      plotOutput("gpd_diagnostics"),
      
      verbatimTextOutput("gpd_summary"),
      
      if(!is_stationary) {
        tagList(h4("Splines du modèle (diagnostic)"), plotOutput("gpd_splines_plot"))
      } else {
        tagList(
          h4("Niveaux de retour (modèle stationnaire)"),
          plotOutput("gpd_rl_plot"),
          tableOutput("gpd_rl_table")
        )
      }
    )
  })  
  output$gpd_summary <- renderPrint({ req(gpd_fitted_model()); summary(gpd_fitted_model()) })

  output$gpd_diagnostics <- renderPlot({
    req(gpd_fitted_model(), gpd_prep_data())
    m <- gpd_fitted_model()
    peaks <- gpd_prep_data()$peaks # Les données utilisées pour le fit
    
    # 1. Prédiction des paramètres (scale sigma et shape xi) pour chaque point
    # Note: Dans un modèle non-stationnaire, sigma varie pour chaque point.
    preds <- predict(m, newdata=peaks, type="response")
    sig <- preds$scale
    xi <- preds$shape
    
    # 2. Récupération des excès (y = Valeur - Seuil)
    excess <- peaks$excess
    
    # 3. Transformation en résidus exponentiels standard
    # Formule : R = (1/xi) * log(1 + xi * (excess / sig))
    term <- 1 + xi * (excess / sig)
    
    # Gestion du cas xi proche de 0 et des valeurs invalides
    res <- ifelse(abs(xi) < 1e-6, 
                  excess / sig, 
                  (1 / xi) * log(term))
    
    # Nettoyage
    res <- res[!is.na(res) & is.finite(res)]
    
    n <- length(res)
    p_pos <- (1:n) / (n + 1)
    
    theo_q <- -log(1 - p_pos) 
    emp_q <- sort(res)
    
    df_qq <- data.frame(Theorique = theo_q, Empirique = emp_q)
    
    ggplot(df_qq, aes(x=Theorique, y=Empirique)) +
      geom_point(color="#009E73", alpha=0.6, size=2) +
      geom_abline(intercept=0, slope=1, color="red", size=1, linetype="dashed") +
      theme_bw() +
      labs(title = "Diagnostic GPD : QQ-Plot des résidus",
           subtitle = "Les points doivent suivre la ligne rouge",
           x = "Quantiles théoriques",
           y = "Quantiles empiriques")
  })

  output$gpd_splines_plot <- renderPlot({ 
    req(gpd_fitted_model())
    if(length(input$gpd_covariates)>0 || input$threshold_type=="nonstationary") plot(gpd_fitted_model()) 
  })
  
  output$gpd_rl_table <- renderTable({
    req(gpd_fitted_model(), gpd_prep_data(), input$gpd_separation)
    m <- gpd_fitted_model()
    peaks <- gpd_prep_data()$peaks
    
    if(length(input$gpd_covariates) == 0 && input$gpd_separation=="none" && input$threshold_type=="stationary") {
      
      total_yrs <- as.numeric(difftime(max(current_data()$time), min(current_data()$time), units="days")) / 365.25
      lambda <- nrow(peaks) / total_yrs
      
      dummy_data <- peaks[1, ]
      preds <- predict(m, newdata=dummy_data, type="response")
      scale_val <- preds$scale[1]
      shape_val <- preds$shape[1]
      u <- input$threshold
      
      T_periods <- c(2, 5, 10, 20, 50, 100)
      rls <- sapply(T_periods, function(T_val) { 
        if(!is.na(shape_val) && abs(shape_val)<1e-6) {
            u + scale_val * log(T_val * lambda) 
        } else {
            u + (scale_val/shape_val) * ( (T_val * lambda)^shape_val - 1 ) 
        }
      })
      data.frame(Retour_annees=as.integer(T_periods), Niveau=rls)
    }
  })
  
  output$gpd_rl_plot <- renderPlot({
    req(gpd_fitted_model(), gpd_prep_data(), input$gpd_separation)
    m <- gpd_fitted_model()
    peaks <- gpd_prep_data()$peaks
    
    if(length(input$gpd_covariates) == 0 && input$gpd_separation=="none" && input$threshold_type=="stationary") {
       
       total_yrs <- as.numeric(difftime(max(current_data()$time), min(current_data()$time), units="days")) / 365.25
       lambda <- nrow(peaks) / total_yrs
       
       dummy_data <- peaks[1, ]
       preds <- predict(m, newdata=dummy_data, type="response")
       scale_val <- preds$scale[1]
       shape_val <- preds$shape[1]
       u <- input$threshold
       
       return_periods <- c(2, 5, 10, 20, 50, 100)
       
       rls <- sapply(return_periods, function(T_val) { 
         if(!is.na(shape_val) && abs(shape_val)<1e-6) u + scale_val * log(T_val * lambda) else u + (scale_val/shape_val) * ( (T_val * lambda)^shape_val - 1 ) 
       })
       df_rl <- data.frame(Period = return_periods, Level = rls)
       
       seq_T <- seq(1.1, 200, length.out=200)
       seq_rl <- sapply(seq_T, function(T_val) {
         if(!is.na(shape_val) && abs(shape_val)<1e-6) u + scale_val * log(T_val * lambda) else u + (scale_val/shape_val) * ( (T_val * lambda)^shape_val - 1 )
       })
       df_curve <- data.frame(Period = seq_T, Level = seq_rl)
       
       ggplot() +
        geom_line(data=df_curve, aes(x=Period, y=Level), color="darkgreen", size=1) +
        geom_point(data=df_rl, aes(x=Period, y=Level), color="orange", size=3) +
        scale_x_log10(breaks=return_periods, labels=return_periods) +
        theme_bw() +
        labs(x="Période de retour (années) [log]", y="Niveau (m)", title="GPD - Niveaux de retour (stationnaire)")
    }
  })
  
  observeEvent(input$fit_quantreg, {
    req(current_data(), input$qr_separation)
    qs <- as.numeric(strsplit(input$quantile_values, ",")[[1]])
    output$quantreg_plot_X <- renderUI({ plotOutput("qr_plot", height=paste0(input$plotHeight,"px")) })
    output$qr_plot <- renderPlot({ create_quantile_plot(current_data(), var_labels[input$qr_x], var_labels[input$qr_y], qs, separation=input$qr_separation, bins=input$bins2) })
  })
  
  observeEvent(input$fit_linreg, {
    req(current_data())
    f <- as.formula(paste(var_labels[input$linreg_response], "~", paste(var_labels[input$linreg_predictors], collapse="+")))
    m <- lm(f, data=current_data())
    output$linreg_summary <- renderPrint({ summary(m) })
    output$linreg_plot <- renderPlot({ par(mfrow=c(2,2)); plot(m) })
  })

  # --- LOGIQUE CONTOURS ENVIRONNEMENTAUX ---

  
# Remplacez la variable reactive contour_results par celle-ci (pour stocker plus d'infos)
contour_results <- reactiveVal(NULL)

observeEvent(input$fit_contour_btn, {
  req(current_data(), input$contour_var_x, input$contour_var_y)
  
  v_x <- var_labels[input$contour_var_x]
  v_y <- var_labels[input$contour_var_y]
  df_curr <- current_data()
  
  periods <- as.numeric(trimws(strsplit(input$contour_periods, ",")[[1]]))
  periods <- periods[!is.na(periods)]
  if(length(periods) == 0) periods <- c(50)
  
  method <- input$contour_method
  sea_state_h <- input$contour_declust
  
  # Paramètres regroupés pour simplifier les appels
  params_list <- list(
    spacing = input$contour_spacing,
    prob = input$contour_prob,
    sea_state = sea_state_h,
    ht_threshold = input$ht04_threshold
  )
  
  # Variables de stockage
  main_contours <- data.frame()
  ci_bands <- data.frame() # Stockera les polygones d'incertitude
  
  do_boot <- input$contour_boot_active
  n_boot <- if(do_boot) input$contour_boot_n else 0
  
  # Centre des données pour l'interpolation polaire (Moyenne ou Médiane)
  center_point <- c(mean(df_curr[[v_x]], na.rm=TRUE), mean(df_curr[[v_y]], na.rm=TRUE))
  
  # Début du calcul avec barre de progression
  withProgress(message = 'Calcul des contours', value = 0, {
    
    # 1. CALCUL DU CONTOUR PRINCIPAL (ESTIMATION PONCTUELLE)
    incProgress(0.1, detail = "Estimation principale...")
    
    for(rp in periods) {
      res <- run_contour_iteration(df_curr, method, v_x, v_y, rp, params_list)
      if(!is.null(res)) {
        res$Period <- factor(rp)
        main_contours <- rbind(main_contours, res)
      }
    }
    
    # 2. BOUCLE BOOTSTRAP (Si activée)
    if(do_boot && nrow(main_contours) > 0) {
      
      # Pour chaque période demandée, on fait le bootstrap
      for(rp_idx in seq_along(periods)) {
        rp <- periods[rp_idx]
        incProgress(0, detail = paste("Bootstrap RP:", rp, "ans"))
        
        # Matrice temporaire pour stocker les rayons interpolés : lignes=angles, cols=bootstraps
        n_angles <- 360
        r_matrix <- matrix(NA, nrow = n_angles, ncol = n_boot)
        
        valid_boots <- 0
        
        for(b in 1:n_boot) {
          # Update progress bar
          prog_val <- 0.2 + (0.8 * ((rp_idx - 1)/length(periods))) + (0.8/length(periods) * (b/n_boot))
          setProgress(prog_val)
          
          # A. Ré-échantillonnage
          idx_resample <- sample(nrow(df_curr), nrow(df_curr), replace = TRUE)
          df_boot <- df_curr[idx_resample, ]
          
          # B. Calcul Contour
          res_boot <- run_contour_iteration(df_boot, method, v_x, v_y, rp, params_list)
          
          # C. Interpolation Polaire (si succès)
          if(!is.null(res_boot) && nrow(res_boot) > 10) {
             interp <- interpolate_contour_polar(res_boot, center_point, n_angles = n_angles)
             r_matrix[, b] <- sqrt((interp$x - center_point[1])^2 + (interp$y - center_point[2])^2)
             valid_boots <- valid_boots + 1
          }
        }
        
        # D. Calcul des Quantiles (2.5% et 97.5%)
        # On calcule rang par rang (angle par angle)
        if(valid_boots > 10) {
          q_lower <- apply(r_matrix, 1, quantile, probs=0.025, na.rm=TRUE)
          q_upper <- apply(r_matrix, 1, quantile, probs=0.975, na.rm=TRUE)
          
          # Reconstruction coordonnées
          target_theta <- seq(-pi, pi, length.out = n_angles)
          
          df_ci <- data.frame(
            x_min = center_point[1] + q_lower * cos(target_theta),
            y_min = center_point[2] + q_lower * sin(target_theta),
            x_max = center_point[1] + q_upper * cos(target_theta),
            y_max = center_point[2] + q_upper * sin(target_theta),
            Period = factor(rp)
          )
          
          ci_bands <- rbind(ci_bands, df_ci)
        }
      } # Fin loop périodes
    } # Fin bloc bootstrap
    
  }) # Fin withProgress
  
  if(nrow(main_contours) > 0) {
    contour_results(list(
      data = df_curr, 
      contours = main_contours, 
      ci = ci_bands, 
      vars = c(v_x, v_y), 
      method = method,
      boot_info = list(active = do_boot, B = n_boot)
    ))
    showNotification("Calcul terminé.", type="message")
  } else {
    showNotification("Échec du calcul (pas assez de données ou erreur fit).", type="error")
  }
})

output$contour_plot <- renderPlotly({
  req(contour_results())
  res <- contour_results()
  
  df_raw <- res$data
  df_cont <- res$contours
  df_ci <- res$ci
  
  vx <- res$vars[1]
  vy <- res$vars[2]
  
  name_x <- names(var_labels)[var_labels == vx]
  name_y <- names(var_labels)[var_labels == vy]
  
  title_str <- switch(res$method,
                      "direct" = "Direct IFORM",
                      "wl" = "IFORM (Weibull-Lognormal)",
                      "ht04" = "IFORM (Heffernan & Tawn)")
  
  # Construction du graphique
  p <- ggplot()
  
  # 1. Données brutes
  p <- p + geom_point(data = df_raw, aes_string(x = vx, y = vy), 
                      color = "black", alpha = 0.2, size = 0.8)
  
  # 2. Bandes d'incertitude (Bootstrap)
  if(!is.null(df_ci) && nrow(df_ci) > 0) {
    # Astuce pour geom_polygon avec trou ou ribbon fermé :
    # Ici on trace un ruban entre min et max.
    # Pour plotly propre, on peut utiliser geom_ribbon si on trie par angle, 
    # mais geom_polygon est plus générique pour les contours fermés.
    # On va faire simple : un polygone ombré pour l'intervalle.
    
    # On structure les données pour faire un polygone fermé : (x_max -> x_min inverse)
    # Mais geom_ribbon demande x commun.
    # On utilise l'astuce geom_polygon avec un groupe par Période.
    
    # Transformation des données CI pour format polygone
    # On concatène le contour extérieur et le contour intérieur inversé
    
    poly_list <- list()
    for(lev in levels(df_ci$Period)) {
      sub_ci <- df_ci[df_ci$Period == lev, ]
      if(nrow(sub_ci)==0) next
      
      # Contour extérieur
      outer <- data.frame(x = sub_ci$x_max, y = sub_ci$y_max)
      # Contour intérieur (ordre inversé pour fermer la boucle proprement)
      inner <- data.frame(x = rev(sub_ci$x_min), y = rev(sub_ci$y_min))
      
      full_poly <- rbind(outer, inner)
      full_poly$Period <- lev
      poly_list[[lev]] <- full_poly
    }
    df_poly_ci <- do.call(rbind, poly_list)
    
    if(!is.null(df_poly_ci)) {
      p <- p + geom_polygon(data = df_poly_ci, aes(x = x, y = y, fill = Period), 
                            alpha = 0.2, color = NA)
    }
  }
  
  # 3. Contour principal (Ligne)
  p <- p + geom_path(data = df_cont, aes_string(x = vx, y = vy, color = "Period"), size = 1.2)
  
  # Styles
  p <- p + 
    scale_color_viridis_d(option = "plasma", name = "Période (ans)", end = 0.9) +
    scale_fill_viridis_d(option = "plasma", name = "IC 95%", end = 0.9) +
    theme_bw() +
    labs(
      title = paste("Contours :", title_str),
      subtitle = if(res$boot_info$active) paste("Avec intervalle de confiance 95% (B =", res$boot_info$B, ")") else "Estimation simple",
      x = name_x, 
      y = name_y
    )
  
  ggplotly(p)
})

output$contour_info <- renderPrint({
  req(contour_results())
  res <- contour_results()
  cat("Méthode :", switch(res$method, "direct"="Direct IFORM", "wl"="Weibull-Lognormal", "ht04"="Heffernan & Tawn"), "\n")
  cat("Variables :", res$vars[1], "&", res$vars[2], "\n")
  cat("Périodes :", paste(unique(res$contours$Period), collapse=", "), "ans.\n")
  if(res$boot_info$active) {
    cat("Bootstrap : Activé (", res$boot_info$B, " réplications)\n")
    cat("Intervalle : 95% (Quantiles 2.5% - 97.5%)\n")
  } else {
    cat("Bootstrap : Désactivé\n")
  }
})

  spec_node_selected <- reactiveVal(NULL)

  output$map_spectral <- renderLeaflet({
    leaflet(data = rscd_spectral) %>%
      addTiles() %>%
      setView(lng = -4, lat = 48, zoom = 5) %>%
      addCircleMarkers(
        lng = ~longitude, lat = ~latitude, layerId = ~name,
        radius = 6, color = "blue", fillOpacity = 0.6, stroke = TRUE, weight = 1,
        label = ~name, group = "base_points"
      )
  })
  
  observeEvent(input$map_spectral_marker_click, {
    click <- input$map_spectral_marker_click
    spec_node_selected(click$id)
  })
  
  observeEvent(spec_node_selected(), {
    req(spec_node_selected())
    selected_point_data <- rscd_spectral %>% filter(name == spec_node_selected())
    
    leafletProxy("map_spectral") %>%
      clearGroup("highlight_aura") %>% 
      addCircleMarkers(
        data = selected_point_data,
        lng = ~longitude, lat = ~latitude,
        radius = 12, 
        color = "#FFD700",
        fillColor = "#FFD700",
        fillOpacity = 0.5,
        stroke = FALSE,
        group = "highlight_aura"
      )
  })
  
  output$spec_selected_point_ui <- renderUI({
    if(is.null(spec_node_selected())) {
      HTML("<b style='color:red;'>Aucun point sélectionné.</b>")
    } else {
      HTML(paste0("<b>Point :</b> ", spec_node_selected()))
    }
  })
  
  output$plot_spectrum_viz <- renderPlot({
    input$plot_spectrum_btn
    
    isolate({
      req(spec_node_selected(), input$spec_date)
      
      target_date_str <- input$spec_date
      if(nchar(target_date_str) <= 10) target_date_str <- paste0(target_date_str, " 00:00:00")
      
      target_time <- as.POSIXct(target_date_str, tz = "UTC")
      if(is.na(target_time)) return(ggplot() + annotate("text", x=0, y=0, label="Format de date invalide"))
      
      day_str <- format(target_time, "%Y-%m-%d")
      
      if(input$spec_type == "1d") {
        spec_data <- get_1d_spectrum(point = spec_node_selected(), start = day_str, end = day_str)
        if(is.null(spec_data)) return(ggplot() + annotate("text", x=0, y=0, label="Pas de données"))
        idx <- which(spec_data$forcings$time == target_time)
        if(length(idx) == 0) return(ggplot() + annotate("text", x=0, y=0, label="Date non trouvée"))
        
        df_viz <- data.frame(Frequency = spec_data$freq, Energy = spec_data$ef[, idx])
        
        ggplot(df_viz, aes(x=Frequency, y=Energy)) +
          geom_line(color = "dodgerblue", size=1) + geom_area(fill = "dodgerblue", alpha=0.3) +
          theme_minimal() + labs(title = paste("Spectre 1D -", spec_node_selected()), subtitle = target_time, x = "Fréquence (Hz)", y = "Densité spectrale (m²/Hz)")
        
      } else {
        spec_data <- get_2d_spectrum(point = spec_node_selected(), start = day_str, end = day_str)
        if(is.null(spec_data)) return(ggplot() + annotate("text", x=0, y=0, label="Pas de données"))
        idx <- which(spec_data$forcings$time == target_time)
        if(length(idx) == 0) return(ggplot() + annotate("text", x=0, y=0, label="Date non trouvée"))
        
        matrix_efth <- spec_data$efth[, , idx]
        dirs <- spec_data$dir
        freqs <- spec_data$freq
        df_viz <- expand.grid(Direction = dirs, Frequency = freqs)
        df_viz$Energy <- as.vector(matrix_efth)
        df_viz$Direction <- (360 - df_viz$Direction) %% 360
        
        u_freq <- sort(unique(freqs))
        diffs <- diff(u_freq)
        y_mins <- c(u_freq[1] - diffs[1]/2, u_freq[-length(u_freq)] + diffs/2)
        y_maxs <- c(u_freq[-length(u_freq)] + diffs/2, u_freq[length(u_freq)] + diffs[length(diffs)]/2)
        
        bounds_df <- data.frame(Frequency = u_freq, ymin = y_mins, ymax = y_maxs)
        df_viz <- df_viz %>% left_join(bounds_df, by="Frequency")
        
        dir_step <- 360 / length(unique(dirs))
        df_viz$xmin <- df_viz$Direction - dir_step/2
        df_viz$xmax <- df_viz$Direction + dir_step/2
        
        ggplot(df_viz) +
          geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=Energy)) +
          scale_fill_viridis_c(option="H", name = "E (m²s/°)") + 
          theme_minimal() +
          scale_x_continuous(breaks = seq(0, 360, 45), limits = c(0, 360), expand = c(0,0)) +
          scale_y_continuous(expand = c(0,0)) +
          labs(
            title = paste("Spectre 2D -", spec_node_selected()), 
            subtitle = target_time, 
            x = "Direction de provenance (°)", 
            y = "Fréquence (Hz)"
          )
      }
    })
  })
}

shinyApp(ui, server)