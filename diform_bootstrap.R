library(data.table)
library(stats)
library(ggplot2)
library(viridis)

# ==============================================================================
# 1. FONCTIONS AUXILIAIRES (Génération de directions & Ajustement GPD)
# ==============================================================================
get_convex_hull_coords <- function(coords) {

  # chull retourne les INDICES des points constituant l'enveloppe convexe,
  # triés dans le sens des aiguilles d'une montre.
  hull_idx <- chull(coords)
  
  # On extrait les points et on ferme la boucle (le dernier point = le premier)
  hull_coords <- coords[c(hull_idx, hull_idx[1]), ]
  
  return(hull_coords)
}

generate_directions <- function(ndim, spacing = 0.1) {
  # Génère des vecteurs unitaires uniformément répartis sur l'hypersphère
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

fit_constrained_gpd <- function(x, threshold) {
  # Ajuste une loi GPD avec contraintes physiques sur le paramètre de forme (xi)
  nll <- function(theta) {
    scale <- exp(theta[1])
    xi <- theta[2]
    
    # Contraintes: xi doit être réaliste pour l'environnement (-0.5 à 0 typiquement)
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
  # Estimateurs par la méthode des moments pour l'initialisation
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

reconstruct_contour_2d <- function(U, R, center_mean, center_sd) {
  # Reconstruit les coordonnées cartésiennes (x,y) à partir des rayons R et directions U
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
  
  # Dénormalisation
  vertices <- vertices[complete.cases(vertices), , drop=FALSE]
  if(nrow(vertices) > 0){
    vertices[,1] <- vertices[,1] * center_sd[1] + center_mean[1]
    vertices[,2] <- vertices[,2] * center_sd[2] + center_mean[2]
    # Fermeture du contour pour le tracé
    vertices <- rbind(vertices, vertices[1,])
  }
  return(vertices)
}

# ==============================================================================
# 2. MOTEUR DE CALCUL IFORM (Vectorisé sur les périodes de retour)
# ==============================================================================

compute_iform_radii <- function(dt, vars, U, return_periods, 
                                stats_list, peak_separation_hours, threshold_prob) {
  
  # Normalisation
  Y <- matrix(NA, nrow = nrow(dt), ncol = length(vars))
  for (i in seq_along(vars)) {
    Y[, i] <- (dt[[vars[i]]] - stats_list[[vars[i]]]$med) / stats_list[[vars[i]]]$sd
  }
  
  total_years <- as.numeric(difftime(max(dt$time), min(dt$time), units="days"))/365.25
  if(total_years < 1) total_years <- 1 # Sécurité pour petits jeux de données
  
  # Paramètres de fenêtre glissante
  time_diffs <- as.numeric(diff(dt$time), units="hours")
  dt_hours <- median(time_diffs, na.rm=TRUE)
  if(is.na(dt_hours) || dt_hours == 0) dt_hours <- 1
  w_size_steps <- round(peak_separation_hours / dt_hours)
  if (w_size_steps %% 2 == 0) w_size_steps <- w_size_steps + 1
  
  # Matrice de stockage : Lignes = Directions, Colonnes = Périodes de retour
  radii_matrix <- matrix(NA, nrow = nrow(U), ncol = length(return_periods))
  
  for (k in 1:nrow(U)) {
    # 1. Projection sur la direction k
    r_proj <- Y %*% U[k,]
    
    # 2. Dégroupage (Declustering) pour trouver les pics indépendants
    roll_max <- frollapply(as.numeric(r_proj), n=w_size_steps, FUN=max, align="center", fill=NA)
    is_peak <- (r_proj == roll_max) & !is.na(r_proj) & !is.na(roll_max)
    peak_values <- r_proj[is_peak]
    
    if (length(peak_values) < 20) next 
    
    # 3. Sélection des excès au-dessus du seuil
    u_thresh <- quantile(peak_values, probs = 1 - threshold_prob)
    excesses <- peak_values[peak_values > u_thresh]
    
    if (length(excesses) < 5) next
    
    # 4. Ajustement GPD (Une seule fois par direction)
    fit <- fit_constrained_gpd(excesses - u_thresh, u_thresh)
    
    # 5. Calcul des niveaux de retour pour TOUTES les périodes demandées
    lambda <- length(excesses) / total_years
    xi <- fit$shape
    sigma <- fit$scale
    
    for(j in seq_along(return_periods)) {
      rp <- return_periods[j]
      # Formule inverse de la GPD
      if (abs(xi) < 1e-6) {
        rv <- u_thresh + sigma * log(lambda * rp)
      } else {
        rv <- u_thresh + (sigma/xi) * ((lambda * rp)^xi - 1)
      }
      radii_matrix[k, j] <- rv
    }
  }
  
  return(radii_matrix)
}

# ==============================================================================
# 3. PROCÉDURE BOOTSTRAP PRINCIPALE
# ==============================================================================

bootstrap_iform_contours <- function(data, 
                                     vars = c("hs", "tp"), 
                                     return_periods = c(1, 10, 50, 100),
                                     n_bootstrap = 50, # Nombre de rééchantillonnages
                                     peak_separation_hours = 72, 
                                     threshold_prob = 0.1,
                                     spacing = 0.1) {
  
  if (!all(vars %in% names(data))) stop("Variables manquantes dans les données")
  
  # Préparation des données
  dt_orig <- data.table(data)[, c("time", vars), with = FALSE]
  dt_orig[, year := year(time)]
  unique_years <- unique(dt_orig$year)
  
  # Statistiques globales (fixées pour tout le processus pour cohérence de l'espace normalisé)
  stats_list <- lapply(vars, function(v) {
    list(med = median(dt_orig[[v]], na.rm=TRUE), sd = sd(dt_orig[[v]], na.rm=TRUE))
  })
  names(stats_list) <- vars
  
  # Directions (fixées)
  U <- generate_directions(length(vars), spacing)
  
  # ----------------------------------------------------------------------------
  # A. Estimation sur les données originales (Best Estimate)
  # ----------------------------------------------------------------------------
  cat("Calcul des contours originaux...\n")
  radii_orig <- compute_iform_radii(dt_orig, vars, U, return_periods, 
                                    stats_list, peak_separation_hours, threshold_prob)
  
  # ----------------------------------------------------------------------------
  # B. Boucle Bootstrap (Block Bootstrap par année)
  # ----------------------------------------------------------------------------
  cat(sprintf("Lancement du Bootstrap (%d itérations)...\n", n_bootstrap))
  
  # Stockage : Liste de matrices [N_dir x N_return]
  boot_results <- list()
  
  for (b in 1:n_bootstrap) {
    if(b %% 10 == 0) cat(sprintf("  Iteration %d/%d\n", b, n_bootstrap))
    
    # 1. Rééchantillonnage des années avec remise
    sampled_years <- sample(unique_years, length(unique_years), replace = TRUE)
    
    # 2. Reconstruction du dataset temporel
    dt_boot_list <- list()
    fake_year_start <- 2000 # Année fictive pour maintenir la continuité chronologique relative
    
    for (y in sampled_years) {
      subset_data <- dt_orig[year == y]
      # On décale le temps pour éviter les trous ou chevauchements bizarres
      # On crée une séquence continue artificielle basée sur l'ordre de tirage
      subset_data[, temp_id := 1:.N]
      dt_boot_list[[length(dt_boot_list) + 1]] <- subset_data
    }
    
    # Fusionner et recréer un axe temporel continu approximatif (suffisant pour frollapply)
    dt_boot <- rbindlist(dt_boot_list)
    # Recréer un vecteur temps artificiel continu pour le calcul des écarts
    dt_boot[, time := seq(from=as.POSIXct("2000-01-01 00:00:00"), 
                          by="1 hour", 
                          length.out=.N)]
    
    # 3. Calcul IFORM sur le jeu rééchantillonné
    # Note : On utilise les stats de normalisation d'origine pour projeter dans le même espace
    radii_b <- compute_iform_radii(dt_boot, vars, U, return_periods, 
                                   stats_list, peak_separation_hours, threshold_prob)
    
    boot_results[[b]] <- radii_b
  }
  
# ----------------------------------------------------------------------------
  # C. Agrégation et Calcul des Intervalles de Confiance (CORRIGÉ AVEC CHULL)
  # ----------------------------------------------------------------------------
  arr_boot <- array(unlist(boot_results), dim = c(nrow(U), length(return_periods), n_bootstrap))
  
  final_contours <- list()
  
  means <- c(stats_list[[1]]$med, stats_list[[2]]$med)
  sds <- c(stats_list[[1]]$sd, stats_list[[2]]$sd)
  
  for (j in seq_along(return_periods)) {
    rp <- return_periods[j]
    
    # Extraction des rayons
    R_est <- radii_orig[, j]
    R_boot_mat <- arr_boot[, j, ]
    R_lower <- apply(R_boot_mat, 1, quantile, probs = 0.025, na.rm=TRUE)
    R_upper <- apply(R_boot_mat, 1, quantile, probs = 0.975, na.rm=TRUE)
    
    # 1. Reconstruction brute des coordonnées 2D
    coords_est <- reconstruct_contour_2d(U, R_est, means, sds)
    coords_low <- reconstruct_contour_2d(U, R_lower, means, sds)
    coords_upp <- reconstruct_contour_2d(U, R_upper, means, sds)
    
    # 2. APPLICATION DE LA CORRECTION CONVEX HULL
    # On applique chull() pour éliminer les concavités et boucles
    coords_est <- get_convex_hull_coords(coords_est)
    coords_low <- get_convex_hull_coords(coords_low)
    coords_upp <- get_convex_hull_coords(coords_upp)
    
    # 3. Formatage pour ggplot
    df_est <- as.data.frame(coords_est); df_est$type <- "Estimate"
    df_low <- as.data.frame(coords_low); df_low$type <- "Lower"
    df_upp <- as.data.frame(coords_upp); df_upp$type <- "Upper"
    
    df_combined <- rbind(df_est, df_low, df_upp)
    df_combined$ReturnPeriod <- as.factor(rp)
    
    final_contours[[j]] <- df_combined
  }
  
  all_data <- do.call(rbind, final_contours)
  colnames(all_data)[1:2] <- vars 
  
  return(all_data)
}

# ==============================================================================
# 4. EXEMPLE D'UTILISATION ET VISUALISATION GGPLOT
# ==============================================================================

set.seed(123)
data <- fread("C:/these_docs/data2024.csv")

results_df <- bootstrap_iform_contours(
  data, 
  vars = c("hs", "tp"), 
  return_periods = c(1, 10, 50, 100), 
  n_bootstrap = 50, 
  spacing = 0.2
)

df_lines <- results_df[results_df$type == "Estimate", ]
df_lower <- results_df[results_df$type == "Lower", ]
df_upper <- results_df[results_df$type == "Upper", ]

ci_polygons <- list()
for(rp in unique(results_df$ReturnPeriod)) {
  low <- df_lower[df_lower$ReturnPeriod == rp, ]
  upp <- df_upper[df_upper$ReturnPeriod == rp, ]
  
  poly_df <- rbind(upp, low[nrow(low):1, ])
  poly_df$ReturnPeriod <- rp
  ci_polygons[[as.character(rp)]] <- poly_df
}
df_ribbons <- do.call(rbind, ci_polygons)

sorted_levels <- sort(unique(results_df$ReturnPeriod))
df_ribbons$ReturnPeriod <- factor(df_ribbons$ReturnPeriod, levels = sorted_levels)
df_lines$ReturnPeriod <- factor(df_lines$ReturnPeriod, levels = sorted_levels)

df_lines$hs <- pmax(df_lines$hs, 0)
df_ribbons$hs <- pmax(df_ribbons$hs, 0)

p <- ggplot() +
  geom_point(data = data, aes(x = hs, y = tp),alpha=0.2) +
  
  geom_polygon(data = df_ribbons, 
               aes(x = hs, y = tp, fill = ReturnPeriod, group = ReturnPeriod), 
               alpha = 0.2) +
  
  geom_path(data = df_lines, 
            aes(x = hs, y = tp, color = ReturnPeriod, group = ReturnPeriod), 
            linewidth = 1) +
  
  scale_color_viridis_d(option = "plasma", name = "Période de\nretour (ans)") +
  scale_fill_viridis_d(option = "plasma", name = "Période de\nretour (ans)") +
  labs(title = "Contours D-IFORM",
       subtitle = "IC 95% par échantillonnage par blocs avec remise",
       x = "Hauteur significative (m)",
       y = "Période pic (s)") +
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())

print(p)
