# Installation des paquets si nécessaire
# install.packages("ggplot2")
# install.packages("patchwork")

library(ggplot2)
library(patchwork)

# --- 1. Génération de données fictives pour la forme de la courbe ---
# On crée une fonction qui ressemble à votre dessin (une bosse asymétrique)
set.seed(123)
x <- seq(0, 10, length.out = 200)

# Fonction de base pour la forme de la courbe (un mélange de gaussiennes)
base_shape <- function(x) {
  # Une bosse principale vers x=6 et une petite ondulation avant
  y <- 2 + 0.5 * sin(x * 1.5) + 4 * exp(-0.2 * (x - 6.5)^2)
  return(y)
}

y_vals <- base_shape(x)

# --- 2. Création du Graphique de GAUCHE (k=3, Niveau ~ 8) ---
# Ici, le niveau est haut car l'espace est peu discrétisé (grand lambda)
df1 <- data.frame(x = x, y = y_vals * 1.2) # On scale pour atteindre ~8
max_y1 <- max(df1$y) # Le max est environ 8

# Positions des barres bleues (discrétisation grossière)
knots1 <- c(3.5, 7.5)

p1 <- ggplot(df1, aes(x, y)) +
  # La ligne pointillée du maximum
  geom_segment(aes(x = 0, xend = x[which.max(y)], y = max_y1, yend = max_y1), 
               color = "red", linetype = "dashed", alpha = 0.6) +
  # La courbe rouge
  geom_line(color = "red", linewidth = 1.2) +
  # Les barres bleues (discrétisation)
  geom_segment(data = data.frame(k = knots1), 
               aes(x = k, xend = k, y = 0, yend = 8.5), 
               color = "royalblue", linewidth = 1.5) +
  # Annotation du niveau 8
  annotate("text", x = 0.5, y = max_y1, label = round(max_y1, 0), 
           color = "red", size = 8, fontface = "bold") +
  # Titres et axes
  labs(title = "Retour calculé GAM-GPD (100 ans)\nDiscrétisation k = 3",
       y = "Niveau de retour",
       x = "Covariable X") +
  theme_classic() +
  theme(axis.line = element_line(arrow = arrow(length = unit(0.3, "cm")))) +
  scale_y_continuous(limits = c(0, 10), breaks = NULL) +
  scale_x_continuous(limits = c(0, 11), breaks = NULL)

# --- 3. Création du Graphique de DROITE (k=5, Niveau ~ 6) ---
# Ici, le niveau baisse car la discrétisation est fine (petit lambda -> dilution)
df2 <- data.frame(x = x, y = y_vals * 0.9) # On scale pour baisser vers 6
max_y2 <- max(df2$y) 

# Positions des barres bleues (discrétisation fine)
knots2 <- c(2, 4, 6, 8)

p2 <- ggplot(df2, aes(x, y)) +
  # La ligne pointillée du maximum
  geom_segment(aes(x = 0, xend = x[which.max(y)], y = max_y2, yend = max_y2), 
               color = "red", linetype = "dashed", alpha = 0.6) +
  # La courbe rouge (plus basse)
  geom_line(color = "red", linewidth = 1.2) +
  # Les barres bleues (plus nombreuses)
  geom_segment(data = data.frame(k = knots2), 
               aes(x = k, xend = k, y = 0, yend = 8.5), 
               color = "royalblue", linewidth = 1.5) +
  # Annotation du niveau 6
  annotate("text", x = 0.5, y = max_y2, label = round(max_y2, 0), 
           color = "red", size = 8, fontface = "bold") +
  # Titres et axes
  labs(title = "Pareil avec k = 5\n(Effet de dilution du lambda)",
       y = "Niveau de retour",
       x = "Covariable X") +
  theme_classic() +
  theme(axis.line = element_line(arrow = arrow(length = unit(0.3, "cm")))) +
  scale_y_continuous(limits = c(0, 10), breaks = NULL) +
  scale_x_continuous(limits = c(0, 11), breaks = NULL)

# --- 4. Assemblage ---
combined_plot <- p1 + p2
print(combined_plot)