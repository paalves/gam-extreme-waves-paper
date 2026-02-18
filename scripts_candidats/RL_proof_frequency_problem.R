# Chargement des bibliothèques nécessaires
library(evgam)
library(dplyr)
library(ggplot2)
library(viridis) # Pour de belles palettes de couleurs

# ==============================================================================
# 1. Préparation et Simulation (À remplacer par vos données réelles)
# ==============================================================================
set.seed(42)

# Simulation de 31 ans de données horaires (approximatif)
n_hours <- 31 * 365 * 24
data <- fread("data2024.csv")

library(dplyr)
library(extRemes) 
u <- quantile(data$hs, 0.96)

declust_obj <- decluster(data$hs, threshold = u, r = 72, verbose = FALSE)
peak_indices <- which(declust_obj > u)

peaks <- data[peak_indices, ] %>%
  mutate(excess = hs - u)

# --- Vérification Rapide ---
cat("Nombre d'excès bruts : ", sum(data$hs > u), "\n")
cat("Nombre de peaks indépendants (après declustering) : ", nrow(peaks), "\n")
cat("Ratio de réduction : ", round(nrow(peaks) / sum(data$hs > u) * 100, 1), "%\n")

# ==============================================================================
# 2. Ajustement du Modèle evgam
# ==============================================================================

formula_evgam_aggressive <- list(
  excess ~ te(Cspd, angleDiff, bs = c("cr", "cc"), k = c(12, 20)),
  ~ 1
)
fit <- evgam(formula_evgam_aggressive, data = peaks, family = "gpd")
plot(fit)
cat("Résumé du modèle :\n")
print(summary(fit))


# ==============================================================================
# 3. Discrétisation de l'Espace (CORRIGÉ)
# ==============================================================================

n_bins_cspd <- 30
n_bins_angle <- 50

# Définition des bornes
breaks_cspd <- seq(min(data$Cspd, na.rm=TRUE), max(data$Cspd, na.rm=TRUE), length.out = n_bins_cspd + 1)
breaks_angle <- seq(0, 360, length.out = n_bins_angle + 1)

# Fonction d'assignation
assign_bins <- function(df) {
  df %>%
    mutate(
      bin_cspd = cut(Cspd, breaks = breaks_cspd, include.lowest = TRUE, labels = FALSE),
      bin_angle = cut(angleDiff, breaks = breaks_angle, include.lowest = TRUE, labels = FALSE)
    )
}

data_binned <- assign_bins(data)
peaks_binned <- assign_bins(peaks)

# Calcul statistiques (Identique à votre logique)
nyears <- 31
grid_stats <- data_binned %>%
  group_by(bin_cspd, bin_angle) %>%
  summarise(
    nb_obs_total = n(),
    Cspd_center = mean(Cspd, na.rm=TRUE),
    angleDiff_center = mean(angleDiff, na.rm=TRUE),
    .groups = 'drop'
  ) %>%
  left_join(
    peaks_binned %>% 
      group_by(bin_cspd, bin_angle) %>% 
      summarise(nb_peaks = n(), .groups = 'drop'),
    by = c("bin_cspd", "bin_angle")
  ) %>%
  mutate(
    nb_peaks = coalesce(nb_peaks, 0),
    nb_obs_moyen_an_bin = nb_obs_total / nyears,
    lambda_bin = (nb_peaks / nb_obs_total) * nb_obs_moyen_an_bin
  ) %>%
  filter(nb_obs_total > 0)
  T_return <- 100


# ==============================================================================
# CORRECTION ET VISUALISATION OPTIMISÉE
# ==============================================================================

# 1. Définition stricte des centres de bins pour le graphique
# On calcule le milieu géométrique de chaque intervalle défini par 'breaks'
mid_cspd <- (breaks_cspd[-1] + breaks_cspd[-length(breaks_cspd)]) / 2
mid_angle <- (breaks_angle[-1] + breaks_angle[-length(breaks_angle)]) / 2

# Création d'une grille de référence parfaite pour ggplot
grid_ref <- expand.grid(
  bin_cspd = 1:n_bins_cspd,
  bin_angle = 1:n_bins_angle
)
grid_ref$plot_x_angle <- mid_angle[grid_ref$bin_angle]
grid_ref$plot_y_cspd <- mid_cspd[grid_ref$bin_cspd]

# 2. Jointure avec vos statistiques calculées
# On conserve vos calculs (Cspd_center) pour la prédiction, mais on ajoute les coords de plot
grid_visualization <- grid_stats %>%
  right_join(grid_ref, by = c("bin_cspd", "bin_angle"))

# 3. Prédiction (Sur les données observées, c'est plus précis)
# Attention : Si une case de la grille est vide de données, Cspd_center sera NA.
# On doit gérer ces cas.
newdata <- grid_visualization %>% 
  dplyr::select(Cspd = Cspd_center, angleDiff = angleDiff_center)

# On ne prédit que là où on a des données pour éviter les erreurs
valid_rows <- !is.na(newdata$Cspd) & !is.na(newdata$angleDiff)
pred_params <- matrix(NA, nrow = nrow(newdata), ncol = 2) # Matrice vide
colnames(pred_params) <- c("scale", "shape")

# Prédiction uniquement sur les lignes valides
if(sum(valid_rows) > 0){
  preds <- predict(fit, newdata[valid_rows, ], type = "response")
  pred_params[valid_rows, ] <- as.matrix(preds)
}

# 4. Calcul du Niveau de Retour
grid_final <- grid_visualization %>%
  bind_cols(as.data.frame(pred_params)) %>%
  rename(scale_pred = scale, shape_pred = shape) %>%
  mutate(
    # Gestion du taux d'occurrence nul
    lambda_bin = coalesce(lambda_bin, 0),
    
    # Calcul du RL (Return Level)
    # Si lambda est 0, le terme (T*lambda)^shape devient Infini (si shape < 0).
    # On force à NA proprement.
    rl_term = ifelse(lambda_bin > 0, (T_return * lambda_bin)^shape_pred - 1, NA),
    rl_50 = u + (scale_pred / shape_pred) * rl_term
  )

# Nettoyage final pour le plot : Si rl_50 est NA, c'est qu'il n'y a pas de risque
# calculable (pas d'excès).

# ==============================================================================
# 5. Le Plot (Fonctionnel)
# ==============================================================================

ggplot(grid_final, aes(x = plot_y_cspd, y = plot_x_angle, fill = rl_50)) +
  geom_tile() + 
  scale_fill_viridis(
    name = paste("Hs", T_return, "ans (m)"), 
    option = "plasma", 
    na.value = "gray10" # Les zones NA (lambda=0 ou pas de données) en gris foncé
  ) +
  labs(
    title = paste("Niveaux de retour Hs", T_return, "ans (Modèle evgam)"),
    subtitle = "Zones grises : Pas d'excès observés (Lambda = 0) ou données insuffisantes",
    x = "Direction (Angle °)",
    y = "Vitesse du Courant (m/s)"
  ) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0, 360, 45), limits = c(0, 360), expand = c(0,0)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), # Pas de grille, geom_tile remplit tout
    panel.background = element_rect(fill = "gray10") # Fond noir pour les trous
  )

mean(grid_final$lambda_bin)
