library(data.table)
library(stats)

# ==============================================================================
# 1. HELPER: Generate Direction Vectors
# ==============================================================================
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

# ==============================================================================
# 2. HELPER: Constrained GPD Fit
# ==============================================================================
fit_constrained_gpd <- function(x, threshold) {
  nll <- function(theta) {
    scale <- exp(theta[1])
    xi <- theta[2]
    
    if (xi >= -1e-4 || xi < -1) return(1e9)
    limit <- threshold - scale/xi
    if (max(x + threshold) > limit) return(1e9)
    
    term <- 1 + xi * (x / scale)
    if (any(term <= 0)) return(1e9)
    
    val <- length(x) * log(scale) + (1/xi + 1) * sum(log(term))
    return(val)
  }
  
  x_bar <- mean(x)
  s2 <- var(x)
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

# ==============================================================================
# 3. HELPER: Reconstruct 2D Contour
# ==============================================================================
reconstruct_contour_2d <- function(U, R, center_mean, center_sd) {
  n <- nrow(U)
  vertices <- matrix(NA, nrow = n, ncol = 2)
  
  for (i in 1:n) {
    idx1 <- i
    idx2 <- if (i == n) 1 else i + 1 
    
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
  vertices[,1] <- vertices[,1] * center_sd[1] + center_mean[1]
  vertices[,2] <- vertices[,2] * center_sd[2] + center_mean[2]
  
  vertices <- rbind(vertices, vertices[1,])
  return(vertices)
}

# ==============================================================================
# 4. MAIN FUNCTION: Modified to Return Peak Indices
# ==============================================================================
direct_iform <- function(data, 
                         vars = c("hs", "tp"), 
                         return_period = 50, 
                         peak_separation_hours = 72, 
                         threshold_prob = 0.1,
                         spacing = 0.1) {
  
  if (!all(vars %in% names(data))) stop("Variables missing in data")
  if (!"time" %in% names(data)) stop("Time column missing")
  
  dt <- data.table(data)[, c("time", vars), with = FALSE]
  
  # 1. Normalization
  stats_list <- lapply(vars, function(v) {
    list(med = median(dt[[v]], na.rm=TRUE), sd = sd(dt[[v]], na.rm=TRUE))
  })
  names(stats_list) <- vars
  
  Y <- matrix(NA, nrow = nrow(dt), ncol = length(vars))
  for (i in seq_along(vars)) {
    Y[, i] <- (dt[[vars[i]]] - stats_list[[vars[i]]]$med) / stats_list[[vars[i]]]$sd
  }
  
  # 2. Directions
  U <- generate_directions(length(vars), spacing)
  
  # 3. Projection & Extremes
  return_values <- numeric(nrow(U))
  total_years <- as.numeric(difftime(max(dt$time), min(dt$time), units="days"))/365.25
  
  # Calculate window size
  time_diffs <- as.numeric(diff(dt$time), units="hours")
  dt_hours <- median(time_diffs, na.rm=TRUE)
  if(is.na(dt_hours) || dt_hours == 0) dt_hours <- 1
  w_size_steps <- round(peak_separation_hours / dt_hours)
  if (w_size_steps %% 2 == 0) w_size_steps <- w_size_steps + 1
  
  # LIST TO STORE PEAK INDICES
  all_peak_indices <- list()
  
  cat(sprintf("Calculating extremes for %d directions...\n", nrow(U)))
  
  for (k in 1:nrow(U)) {
    # Projection
    r_proj <- Y %*% U[k,]
    
    # Fast Declustering
    roll_max <- frollapply(as.numeric(r_proj), n=w_size_steps, FUN=max, align="center", fill=NA)
    
    # Identify peaks
    is_peak <- (r_proj == roll_max) & !is.na(r_proj) & !is.na(roll_max)
    
    # STORE INDICES OF DECLUSTERED PEAKS
    # These are the candidate storms for this direction
    all_peak_indices[[k]] <- which(is_peak)
    
    peak_values <- r_proj[is_peak]
    
    if (length(peak_values) < 20) {
      return_values[k] <- NA
      next
    }
    
    # Threshold
    u_thresh <- quantile(peak_values, probs = 1 - threshold_prob)
    excesses <- peak_values[peak_values > u_thresh]
    
    # Fit GPD
    fit <- fit_constrained_gpd(excesses - u_thresh, u_thresh)
    
    # Return Level
    lambda <- length(excesses) / total_years
    xi <- fit$shape
    sigma <- fit$scale
    
    if (abs(xi) < 1e-6) {
      rv <- u_thresh + sigma * log(lambda * return_period)
    } else {
      rv <- u_thresh + (sigma/xi) * ((lambda * return_period)^xi - 1)
    }
    return_values[k] <- rv
  }
  
  # 4. Construct Contour
  contour_coords <- NULL
  if (length(vars) == 2) {
    means <- c(stats_list[[1]]$med, stats_list[[2]]$med)
    sds <- c(stats_list[[1]]$sd, stats_list[[2]]$sd)
    contour_coords <- reconstruct_contour_2d(U, return_values, means, sds)
    
    hull_idx <- chull(contour_coords)
    hull_idx <- c(hull_idx, hull_idx[1]) 
    contour_coords <- contour_coords[hull_idx, ]
    colnames(contour_coords) <- vars
  } else {
    contour_coords <- matrix(NA, nrow=nrow(U), ncol=length(vars))
    for(k in 1:nrow(U)) {
      if(!is.na(return_values[k])) {
        p_norm <- return_values[k] * U[k,]
        for(i in 1:length(vars)) {
          contour_coords[k,i] <- p_norm[i] * stats_list[[vars[i]]]$sd + stats_list[[vars[i]]]$med
        }
      }
    }
  }
  
  # 5. Aggregate Unique Peak Indices from all directions
  unique_peaks <- unique(unlist(all_peak_indices))
  
  return(list(contour = contour_coords, 
              directions = U, 
              return_values = return_values,
              peak_indices = unique_peaks)) # <--- New Return Item
}

# ==============================================================================
# EXAMPLE USAGE
# ==============================================================================

library(ggplot2)
library(dplyr)
library(viridis)

#' Fonction de tracé des contours environnementaux
#'
#' @param contours_list Liste nommée contenant les matrices/dataframes des contours (ex: list("1 an" = ..., "100 ans" = ...))
#' @param data_obs Dataframe des observations (points gris en fond)
#' @param x_var Nom de la colonne pour l'axe X (ex: "hs")
#' @param y_var Nom de la colonne pour l'axe Y (ex: "tp")
#' @param var_for_max Nom de la variable à maximiser pour trouver le point de conception (ex: "hs")
#' @param design_period Nom exact de la période à étiqueter (ex: "100 ans")
#' @param x_lab Label de l'axe X
#' @param y_lab Label de l'axe Y
#' @param title Titre du graphique
#' @param nudge_x Décalage horizontal du label (pour éviter de masquer le point)
#' @param nudge_y Décalage vertical du label
#'
#' @return Un objet ggplot
plot_contours_env <- function(contours_list, 
                              data_obs, 
                              x_var = "hs",
                              y_var = "tp", 
                              var_for_max = "hs", 
                              design_period = "100", # Adjusted default to match numeric/string usually found in data
                              x_lab = "Hs [m]", 
                              y_lab = "Tp [s]",
                              title = "Environmental Contours",
                              nudge_x = 0,
                              nudge_y = 2) { # Adjusted default nudge for visibility
  
  # 1. Prepare Contour Data
  # Transform list to single dataframe
  contours_df <- do.call(rbind, lapply(names(contours_list), function(nm) {
    df <- as.data.frame(contours_list[[nm]])
    # Rename columns to match selected variables if needed
    if(ncol(df) >= 2) {
      colnames(df)[1:2] <- c(x_var, y_var)
    }
    df$periode <- nm
    df
  }))
  
  # Ensure factor order matches input list order
  contours_df$periode <- factor(contours_df$periode, levels = names(contours_list))
  
  # 2. Calculate Red Points (Max per contour)
  # Logic: Sort by var_for_max (desc) THEN by other variable (desc)
  other_var <- if(var_for_max == x_var) y_var else x_var
  
  points_max <- contours_df %>%
    group_by(periode) %>%
    arrange(desc(.data[[var_for_max]]), desc(.data[[other_var]])) %>%
    slice(1) %>%
    ungroup()
  
  # 3. Prepare Design Point Label
  design_pt <- points_max %>% filter(periode == design_period)
  
  if(nrow(design_pt) == 0) {
    # Fallback if specific string match fails, try generic search or warn
    warning(paste("Period '", design_period, "' not found in contours. Label skipped."))
    design_label <- ""
  } else {
    val_x <- round(design_pt[[x_var]], 2)
    val_y <- round(design_pt[[y_var]], 2)
    
    design_label <- paste0("Design point\n(", design_period, " years):\n",
                           x_lab, " = ", val_x, "\n",
                           y_lab, " = ", val_y)
  }
  
  # 4. Create Plot
  p <- ggplot() + 
    # Raw Observations (Matched style: black, alpha 0.5, size 0.8)
    geom_point(data = data_obs, aes(x = .data[[x_var]], y = .data[[y_var]]),
               color = "black", alpha = 0.5, size = 0.8) +
    
    # Contour Lines (Using geom_path for smooth lines, but matching color palette)
    geom_path(data = contours_df,
              aes(x = .data[[x_var]], y = .data[[y_var]], color = periode),
              linewidth = 1) + 
    
    # Max Points (Matched style: red, size 3)
    geom_point(data = points_max,
               aes(x = .data[[x_var]], y = .data[[y_var]]),
               color = "red", size = 3, shape = 19) +
    
    # Design Point Label (Matched style: bold, white fill, alpha 0.9)
    geom_label(data = design_pt,
               aes(x = .data[[x_var]], y = .data[[y_var]], label = design_label),
               fill = "white", color = "black", fontface = "bold", size = 3.5,
               hjust = 0.5, vjust = 0, alpha = 0.9,
               nudge_x = nudge_x,
               nudge_y = nudge_y,
               show.legend = FALSE) +
    
    # Scales and Theme (Matched style: viridis plasma end 0.9, theme_bw)
    scale_color_viridis_d(name = "Return period\n(years)", option = "plasma", end = 0.9) +
    labs(title = title, x = x_lab, y = y_lab) +
    theme_bw() +
    theme(legend.position = "right")
  
  return(p)
}




plot_contours_env_noetiquette <- function(contours_list, 
                                          data_obs, 
                                          x_var = "hs",
                                          y_var = "tp", 
                                          var_for_max = "hs", 
                                          design_period = "100 ans",
                                          x_lab = "Hs [m]", 
                                          y_lab = "Tp [s]",
                                          title = "Contours D-IFORM",
                                          subtitle="") {
  
  # 1. Préparation des données contours
  contours_df <- do.call(rbind, lapply(names(contours_list), function(nm) {
    df <- as.data.frame(contours_list[[nm]])
    if(ncol(df) >= 2) {
      colnames(df)[1:2] <- c(x_var, y_var)
    }
    df$periode <- nm
    df
  }))
  
  # Forcer l'ordre des facteurs
  contours_df$periode <- factor(contours_df$periode, levels = names(contours_list))
  
  # 2. Calcul des points rouges (Max par contour)
  other_var <- if(var_for_max == x_var) y_var else x_var
  
  points_max <- contours_df %>%
    group_by(periode) %>%
    arrange(desc(.data[[var_for_max]]), desc(.data[[other_var]])) %>%
    slice(1) %>%
    ungroup()
  
  # 3. Préparation du label (optionnel, utilisé ici juste pour filtrer si besoin)
  design_pt <- points_max %>% filter(periode == design_period)
  
  # 4. Création du graphique
  p <- ggplot() + 
    # Points d'observations (Style matched: black, alpha 0.5, size 0.8)
    geom_point(data = data_obs, aes(x = .data[[x_var]], y = .data[[y_var]]),
               color = "black", alpha = 0.5, size = 0.8) +
    
    # Lignes des contours (Style matched: linewidth 1)
    geom_path(data = contours_df,
              aes(x = .data[[x_var]], y = .data[[y_var]], color = periode),
              linewidth = 1) +
    
    # Points rouges (Style matched: red, size 3)
    geom_point(data = points_max,
               aes(x = .data[[x_var]], y = .data[[y_var]]),
               color = "red", size = 3, shape = 19) + 
    
    # Esthétique (Style matched: plasma end 0.9, theme_bw)
    scale_color_viridis_d(name = "Return period\n(years)", option = "plasma", end = 0.9) +
    labs(x = x_lab, y = y_lab, title = title, subtitle=subtitle) +
    theme_bw() +
    theme(legend.position = "right")
  
  return(p)
}


data <- fread("C:/these_docs/data2024.csv")
#data <- fread("C:/these_docs/presentation-0711/data2020fromveur.csv")
#data <- fread("C:/these_docs/presentation-0711/data2020brehat.csv")
#data <- fread("C:/these_docs/presentation-0711/data2020barfleur.csv")

data[, time := as.POSIXct(time, format="%Y-%m-%d %H:%M:%S")]


res1 <- direct_iform(data, vars = c("hs", "tp"), return_period = 1, spacing = 0.2)
res2 <- direct_iform(data, vars = c("hs", "tp"), return_period = 2, spacing = 0.2)
res5 <- direct_iform(data, vars = c("hs", "tp"), return_period = 5, spacing = 0.2)
res10 <- direct_iform(data, vars = c("hs", "tp"), return_period = 10, spacing = 0.2)
res20 <- direct_iform(data, vars = c("hs", "tp"), return_period = 20, spacing = 0.2)
res50 <- direct_iform(data, vars = c("hs", "tp"), return_period = 50, spacing = 0.2)
res100 <- direct_iform(data, vars = c("hs", "tp"), return_period = 100, spacing = 0.2)


# 1. Créer la liste propre des contours
my_contours <- list(
  "1 an"    = res1$contour,
  "2 ans"   = res2$contour,
  "5 ans"   = res5$contour,
  "10 ans"  = res10$contour,
  "20 ans"  = res20$contour,
  "50 ans"  = res50$contour,
  "100 ans" = res100$contour
)

my_contours <- lapply(my_contours, function(x) {
  x[x < 0] <- 0
  return(x)
})

# 2. Appel de la fonction
p1 <- plot_contours_env_noetiquette(
  contours_list = my_contours,
  data_obs      = data,       # Votre dataframe 'data' original
  x_var         = "hs",       # Nom de la colonne Hs
  y_var         = "tp",       # Nom de la colonne Tp
  var_for_max   = "hs",       # On cherche le Hs max
  design_period = "100 ans",  # Le contour à étiqueter
  title         = "D-IFORM contours (Hs, Tp)",
  subtitle = "Projection on 20 directions"
)
p1

########################################################################################
########################################################################################

res1 <- direct_iform(data, vars = c("Cspd", "hs"), return_period = 1, spacing = 0.2)
res2 <- direct_iform(data, vars = c("Cspd", "hs"), return_period = 2, spacing = 0.2)
res5 <- direct_iform(data, vars = c("Cspd", "hs"), return_period = 5, spacing = 0.2)
res10 <- direct_iform(data, vars = c("Cspd", "hs"), return_period = 10, spacing = 0.2)
res20 <- direct_iform(data, vars = c("Cspd", "hs"), return_period = 20, spacing = 0.2)
res50 <- direct_iform(data, vars = c("Cspd", "hs"), return_period = 50, spacing = 0.2)
res100 <- direct_iform(data, vars = c("Cspd", "hs"), return_period = 100, spacing = 0.2)

my_contours <- list(
  "1 an"    = res1$contour,
  "2 ans"   = res2$contour,
  "5 ans"   = res5$contour,
  "10 ans"  = res10$contour,
  "20 ans"  = res20$contour,
  "50 ans"  = res50$contour,
  "100 ans" = res100$contour
)

# 2. Appel de la fonction
p2 <- plot_contours_env_noetiquette(
  contours_list = my_contours,
  data_obs      = data,       # Votre dataframe 'data' original
  x_var         = "Cspd",       # Nom de la colonne Hs
  y_var         = "hs",       # Nom de la colonne Tp
  var_for_max   = "hs",       # On cherche le Hs max
  design_period = "100 ans",  # Le contour à étiqueter
  title         = "D-IFORM contours (Cs, Hs)",
  subtitle = "Projection on 20 directions",
  x_lab = "Cs [m.s-1]",
  y_lab= "Hs [m]"
)
p2

########################################################################################
########################################################################################

res1 <- direct_iform(data, vars = c("angleDiff", "hs"), return_period = 1, spacing = 0.2)
res2 <- direct_iform(data, vars = c("angleDiff", "hs"), return_period = 2, spacing = 0.2)
res5 <- direct_iform(data, vars = c("angleDiff", "hs"), return_period = 5, spacing = 0.2)
res10 <- direct_iform(data, vars = c("angleDiff", "hs"), return_period = 10, spacing = 0.2)
res20 <- direct_iform(data, vars = c("angleDiff", "hs"), return_period = 20, spacing = 0.2)
res50 <- direct_iform(data, vars = c("angleDiff", "hs"), return_period = 50, spacing = 0.2)
res100 <- direct_iform(data, vars = c("angleDiff", "hs"), return_period = 100, spacing = 0.2)

my_contours <- list(
  "1 an"    = res1$contour,
  "2 ans"   = res2$contour,
  "5 ans"   = res5$contour,
  "10 ans"  = res10$contour,
  "20 ans"  = res20$contour,
  "50 ans"  = res50$contour,
  "100 ans" = res100$contour
)

my_contours <- lapply(my_contours, function(x) {
  x[x < 0] <- 0 
  x[x > 360] <- 360
  return(x)
})

# 2. Appel de la fonction
p3 <- plot_contours_env_noetiquette(
  contours_list = my_contours,
  data_obs      = data,       # Votre dataframe 'data' original
  x_var         = "angleDiff",       # Nom de la colonne Hs
  y_var         = "hs",       # Nom de la colonne Tp
  var_for_max   = "hs",       # On cherche le Hs max
  design_period = "100 ans",  # Le contour à étiqueter
  title         = "D-IFORM contours (Hs, angular difference)",
  subtitle = "Projection on 20 directions",
  x_lab = "Angular difference (°)",
  y_lab = "Hs (m)"
)
p3

########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################

data <- fread("C:/these_docs/data2024.csv")
data <- fread("C:/these_docs/presentation-0711/data2020fromveur.csv")
data <- fread("C:/these_docs/presentation-0711/data2020brehat.csv")
data <- fread("C:/these_docs/presentation-0711/data2020barfleur.csv")

res1 <- direct_iform(data, vars = c("OrbSpeed_current", "tp"), return_period = 1, spacing = 0.2)
res2 <- direct_iform(data, vars = c("OrbSpeed_current", "tp"), return_period = 2, spacing = 0.2)
res5 <- direct_iform(data, vars = c("OrbSpeed_current", "tp"), return_period = 5, spacing = 0.2)
res10 <- direct_iform(data, vars = c("OrbSpeed_current", "tp"), return_period = 10, spacing = 0.2)
res20 <- direct_iform(data, vars = c("OrbSpeed_current", "tp"), return_period = 20, spacing = 0.2)
res50 <- direct_iform(data, vars = c("OrbSpeed_current", "tp"), return_period = 50, spacing = 0.2)
res100 <- direct_iform(data, vars = c("OrbSpeed_current", "tp"), return_period = 100, spacing = 0.2)


# 1. Créer la liste propre des contours
my_contours <- list(
  "1 an"    = res1$contour,
  "2 ans"   = res2$contour,
  "5 ans"   = res5$contour,
  "10 ans"  = res10$contour,
  "20 ans"  = res20$contour,
  "50 ans"  = res50$contour,
  "100 ans" = res100$contour
)

# 2. Appel de la fonction
plot_contours_env(
  contours_list = my_contours,
  data_obs      = data,       # Votre dataframe 'data' original
  x_var         = "OrbSpeed_current",       # Nom de la colonne Hs
  y_var         = "tp",       # Nom de la colonne Tp
  var_for_max   = "OrbSpeed_current",       # On cherche le Hs max
  design_period = "100 ans",  # Le contour à étiqueter
  title         = "Contours D-IFORM (vitesse orbitale, Tp) au Raz de Barfleur",
  nudge_x       = -0.1,       # Ajustement horizontal du label (vers la gauche)
  nudge_y       = -9.2,       # Ajustement vertical
  x_lab = "Vit. orbitale à 10m du fond (m/s)",
  y_lab = "Tp (s)"
)

########################################################################################
########################################################################################

res1 <- direct_iform(data, vars = c("Cspd", "OrbSpeed_current"), return_period = 1, spacing = 0.2)
res2 <- direct_iform(data, vars = c("Cspd", "OrbSpeed_current"), return_period = 2, spacing = 0.2)
res5 <- direct_iform(data, vars = c("Cspd", "OrbSpeed_current"), return_period = 5, spacing = 0.2)
res10 <- direct_iform(data, vars = c("Cspd", "OrbSpeed_current"), return_period = 10, spacing = 0.2)
res20 <- direct_iform(data, vars = c("Cspd", "OrbSpeed_current"), return_period = 20, spacing = 0.2)
res50 <- direct_iform(data, vars = c("Cspd", "OrbSpeed_current"), return_period = 50, spacing = 0.2)
res100 <- direct_iform(data, vars = c("Cspd", "OrbSpeed_current"), return_period = 100, spacing = 0.2)

my_contours <- list(
  "1 an"    = res1$contour,
  "2 ans"   = res2$contour,
  "5 ans"   = res5$contour,
  "10 ans"  = res10$contour,
  "20 ans"  = res20$contour,
  "50 ans"  = res50$contour,
  "100 ans" = res100$contour
)

# 2. Appel de la fonction
plot_contours_env(
  contours_list = my_contours,
  data_obs      = data,       # Votre dataframe 'data' original
  x_var         = "Cspd",       # Nom de la colonne Hs
  y_var         = "OrbSpeed_current",       # Nom de la colonne Tp
  var_for_max   = "OrbSpeed_current",       # On cherche le Hs max
  design_period = "100 ans",  # Le contour à étiqueter
  title         = "Contours D-IFORM (vitesse orbitale, vitesse du courant)\nau Raz de Barfleur",
  nudge_x       = -0.2,       # Ajustement horizontal du label (vers la gauche)
  nudge_y       = -0.85,       # Ajustement vertical
  x_lab = "Vit. du courant (m/s)",
  y_lab = "Vit. orbitale à 10m du fond (m/s)"
)

########################################################################################
########################################################################################

res1 <- direct_iform(data, vars = c("angleDiff", "OrbSpeed_current"), return_period = 1, spacing = 0.3)
res2 <- direct_iform(data, vars = c("angleDiff", "OrbSpeed_current"), return_period = 2, spacing = 0.3)
res5 <- direct_iform(data, vars = c("angleDiff", "OrbSpeed_current"), return_period = 5, spacing = 0.3)
res10 <- direct_iform(data, vars = c("angleDiff", "OrbSpeed_current"), return_period = 10, spacing = 0.3)
res20 <- direct_iform(data, vars = c("angleDiff", "OrbSpeed_current"), return_period = 20, spacing = 0.3)
res50 <- direct_iform(data, vars = c("angleDiff", "OrbSpeed_current"), return_period = 50, spacing = 0.3)
res100 <- direct_iform(data, vars = c("angleDiff", "OrbSpeed_current"), return_period = 100, spacing = 0.3)

my_contours <- list(
  "1 an"    = res1$contour,
  "2 ans"   = res2$contour,
  "5 ans"   = res5$contour,
  "10 ans"  = res10$contour,
  "20 ans"  = res20$contour,
  "50 ans"  = res50$contour,
  "100 ans" = res100$contour
)

# 2. Appel de la fonction
plot_contours_env(
  contours_list = my_contours,
  data_obs      = data,       # Votre dataframe 'data' original
  x_var         = "angleDiff",       # Nom de la colonne Hs
  y_var         = "OrbSpeed_current",       # Nom de la colonne Tp
  var_for_max   = "OrbSpeed_current",       # On cherche le Hs max
  design_period = "100 ans",  # Le contour à étiqueter
  title         = "Contours D-IFORM (vitesse orbitale, écart angulaire)\nau Raz de Barfleur",
  nudge_x       = 60,       # Ajustement horizontal du label (vers la gauche)
  nudge_y       = -0.5,       # Ajustement vertical
  x_lab = "Ecart angulaire (°)",
  y_lab = "Vit. orbitale à 10m du fond (m/s)"
)




########################################################################################
##################################BOOTSTRAP##################################################
########################################################################################

library(ggplot2)

# ==============================================================================
# 5. BOOTSTRAP FUNCTION (Unchanged)
# ==============================================================================
bootstrap_iform <- function(original_data, n_boot = 20, vars = c("hs", "tp"), ...) {
  original_data[, year := year(time)]
  unique_years <- unique(original_data$year)
  n_years <- length(unique_years)
  
  boot_contours <- list()
  
  cat("Bootstrapping")
  for (b in 1:n_boot) {
    if(b %% 10 == 0) cat(".")
    sampled_years <- sample(unique_years, n_years, replace = TRUE)
    
    resampled_list <- list()
    for (y in sampled_years) {
      resampled_list[[length(resampled_list) + 1]] <- original_data[year == y]
    }
    boot_data <- rbindlist(resampled_list)
    
    tryCatch({
      capture.output({
        res <- direct_iform(boot_data, vars = vars, ...)
      })
      boot_contours[[b]] <- res$contour
    }, error = function(e) {})
  }
  cat(" Done.\n")
  boot_contours <- boot_contours[!sapply(boot_contours, is.null)]
  return(boot_contours)
}

# ==============================================================================
# 6. CI CALCULATOR (With Negative Check + Polygon Formatting)
# ==============================================================================
# ==============================================================================
# 6. UPDATED CI CALCULATOR (Forces Convex Hull for shape consistency)
# ==============================================================================
calc_contour_ci_poly <- function(boot_contours, center_point, probs = c(0.05, 0.95)) {
  # 1. Interpolate Radii
  angles_rad <- seq(0, 2*pi, length.out = 360)
  R_matrix <- matrix(NA, nrow = length(angles_rad), ncol = length(boot_contours))
  
  for (i in seq_along(boot_contours)) {
    ctr <- boot_contours[[i]]
    x <- ctr[,1] - center_point[1]
    y <- ctr[,2] - center_point[2]
    r_pts <- sqrt(x^2 + y^2)
    theta_pts <- atan2(y, x) 
    
    ord <- order(theta_pts)
    theta_pts <- theta_pts[ord]; r_pts <- r_pts[ord]
    
    theta_full <- c(theta_pts - 2*pi, theta_pts, theta_pts + 2*pi)
    r_full <- c(r_pts, r_pts, r_pts)
    
    approx_fun <- approxfun(theta_full, r_full)
    R_matrix[, i] <- approx_fun(angles_rad)
  }
  
  # 2. Calculate Quantiles
  ci_lower <- apply(R_matrix, 1, quantile, probs = probs[1], na.rm = TRUE)
  ci_upper <- apply(R_matrix, 1, quantile, probs = probs[2], na.rm = TRUE)
  
  # 3. Convert to Cartesian & Prevent Negatives
  get_coords <- function(r) {
    xx <- center_point[1] + r * cos(angles_rad)
    yy <- center_point[2] + r * sin(angles_rad)
    xx <- pmax(xx, 0)
    yy <- pmax(yy, 0)
    return(data.frame(x=xx, y=yy))
  }
  
  up <- get_coords(ci_upper)
  lo <- get_coords(ci_lower)

  # --- FIX: APPLY CONVEX HULL TO THE CI ---
  # This forces the CI to have the same "straight edge" shape as the original
  
  # Apply Hull to Upper Bound
  h_up <- chull(up)
  h_up <- c(h_up, h_up[1]) # Close the loop
  up <- up[h_up, ]
  
  # Apply Hull to Lower Bound
  h_lo <- chull(lo)
  h_lo <- c(h_lo, h_lo[1]) # Close the loop
  lo <- lo[h_lo, ]
  
  # ----------------------------------------
  
  # 4. Create Closed Polygon for ggplot (Upper Forward -> Lower Backward)
  poly_df <- data.frame(
    hs = c(up$x, rev(lo$x)),
    tp = c(up$y, rev(lo$y))
  )
  
  return(poly_df)
}
# ==============================================================================
# 7. WORKFLOW & GGPLOT
# ==============================================================================
# --- Settings ---
rps <- c(1, 2, 5, 10, 20, 50, 100)
# Define colors mapped to the RP labels
custom_colors <- c(
  "1 Year"   = "forestgreen", 
  "2 Year"   = "green", 
  "5 Year"   = "yellow",
  "10 Year"  = "orange",
  "20 Year"  = "red",
  "50 Year"  = "purple",
  "100 Year" = "black"
)

center_pt <- c(mean(data$hs, na.rm=TRUE), mean(data$tp, na.rm=TRUE))

# Data containers for plotting
df_contours_all <- data.frame()
df_ci_all <- data.frame()

cat("Starting Calculations...\n")

for(rp in rps) {
  rp_label <- paste(rp, "Year")
  cat(sprintf("Processing %s...\n", rp_label))
  
  # A. Main Estimate
  est <- direct_iform(data, vars = c("hs", "tp"), return_period = rp, spacing = 0.1)
  est_df <- as.data.frame(est$contour)
  est_df$RP <- rp_label # Add label
  df_contours_all <- rbind(df_contours_all, est_df)
  
  # B. Bootstrap Uncertainty
  # Note: Increase n_boot (e.g., 100) for final publication quality
  boots <- bootstrap_iform(data, n_boot = 5, vars = c("hs", "tp"), return_period = rp, spacing = 0.1)
  
  # C. CI Polygon
  ci_poly <- calc_contour_ci_poly(boots, center_point = center_pt, probs = c(0.05, 0.95))
  ci_poly$RP <- rp_label
  df_ci_all <- rbind(df_ci_all, ci_poly)
}

# Ensure Factor Order for Legend (so 1 Year is first, 100 Year is last)
rp_levels <- paste(rps, "Year")
df_contours_all$RP <- factor(df_contours_all$RP, levels = rp_levels)
df_ci_all$RP <- factor(df_ci_all$RP, levels = rp_levels)

# --- PLOT ---
cat("Generating Plot...\n")

p <- ggplot() +
  # 1. Scatter Plot of Raw Data
  geom_point(data = data, aes(x = hs, y = tp), color = "gray90", alpha = 0.5) +
  
  # 2. Confidence Intervals (Shaded Bands)
  # We use fill=RP to get the legend, alpha=0.2 for transparency
  geom_polygon(data = df_ci_all, aes(x = hs, y = tp, fill = RP, group = RP), alpha = 0.2) +
  
  # 3. Main Contour Lines
  # We use color=RP to match the fill legend
  geom_path(data = df_contours_all, aes(x = hs, y = tp, color = RP, group = RP), size = 1) +
  
  # 4. Custom Colors
  scale_color_manual(values = custom_colors, name = "Return Period\n(with 90% CI)") +
  scale_fill_manual(values = custom_colors, name = "Return Period\n(with 90% CI)") +
  
  # 5. Formatting
  coord_cartesian(xlim = c(0, 9), ylim = c(0, 30), expand = FALSE) +
  labs(title = "Direct IFORM Contours",
       subtitle = "Shaded areas represent 90% Confidence Intervals via Bootstrap",
       x = "Hs (m)", 
       y = "Tp (s)") +
  theme_minimal() +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "gray"),
    plot.title = element_text(face = "bold")
  )

print(p)