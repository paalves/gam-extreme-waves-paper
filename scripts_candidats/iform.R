# Load necessary library
library(ggplot2)

# --- 1. Simulation Parameters ---
set.seed(42)
n_samples <- 2000
# Hs distribution (Weibull) parameters
shape_hs <- 1.5
scale_hs <- 2.5
# Tp | Hs distribution (Lognormal) parameters
sigma_ln <- 0.2

# Function for conditional mu of Tp
get_mu_ln <- function(hs) {
  return(1.5 + 0.1 * hs)
}

# --- 2. Simulate Data (Physical Space) ---
# Generate uniform random variables
u1_sim <- runif(n_samples)
u2_sim <- runif(n_samples)

# Inverse sampling
hs_sim <- qweibull(u1_sim, shape = shape_hs, scale = scale_hs)
tp_sim <- numeric(n_samples)
for(i in 1:n_samples) {
  mu <- get_mu_ln(hs_sim[i])
  tp_sim[i] <- qlnorm(u2_sim[i], meanlog = mu, sdlog = sigma_ln)
}

# --- 3. Transform to U-space (Rosenblatt) ---
U1_sim <- qnorm(pweibull(hs_sim, shape = shape_hs, scale = scale_hs))
U2_sim <- numeric(n_samples)
for(i in 1:n_samples) {
  mu <- get_mu_ln(hs_sim[i])
  # p_val is the CDF of Tp given Hs
  p_val <- plnorm(tp_sim[i], meanlog = mu, sdlog = sigma_ln)
  U2_sim[i] <- qnorm(p_val)
}

# --- 4. Define IFORM Contour ---
# Return period T = 25 years, 3-hour sea states
n_states <- 25 * 365 * 8
beta <- qnorm(1 - 1/n_states)

# Circle in U-space
theta <- seq(0, 2*pi, length.out = 300)
u1_cont <- beta * cos(theta)
u2_cont <- beta * sin(theta)

# Transform Contour back to Physical Space (Inverse Rosenblatt)
hs_cont <- qweibull(pnorm(u1_cont), shape = shape_hs, scale = scale_hs)
tp_cont <- numeric(length(hs_cont))
for(i in 1:length(hs_cont)) {
  mu <- get_mu_ln(hs_cont[i])
  tp_cont[i] <- qlnorm(pnorm(u2_cont[i]), meanlog = mu, sdlog = sigma_ln)
}

# --- 5. Define Limit State Surface (LSS) ---
# Tangent plane at 45 degrees in U-space
alpha <- pi / 4
# Equation in U-space: u1*cos(alpha) + u2*sin(alpha) = beta
# Let's plot a line segment near the tangency point
u1_lss <- seq(0, beta + 2, length.out = 100)
u2_lss <- (beta - u1_lss * cos(alpha)) / sin(alpha)

# Transform LSS back to Physical Space
hs_lss <- qweibull(pnorm(u1_lss), shape = shape_hs, scale = scale_hs)
tp_lss <- numeric(length(hs_lss))
for(i in 1:length(hs_lss)) {
  mu <- get_mu_ln(hs_lss[i])
  tp_lss[i] <- qlnorm(pnorm(u2_lss[i]), meanlog = mu, sdlog = sigma_ln)
}

# --- 6. Plotting ---

# Prepare Data Frames
df_sim <- data.frame(Hs = hs_sim, Tp = tp_sim, U1 = U1_sim, U2 = U2_sim)
df_cont <- data.frame(Hs = hs_cont, Tp = tp_cont, U1 = u1_cont, U2 = u2_cont)
df_lss <- data.frame(Hs = hs_lss, Tp = tp_lss, U1 = u1_lss, U2 = u2_lss)

# Plot Physical Space
p1 <- ggplot() +
  geom_point(data = df_sim, aes(x = Hs, y = Tp), alpha = 0.3, color = "grey") +
  geom_path(data = df_cont, aes(x = Hs, y = Tp), color = "red", size = 1.2) +
  geom_path(data = df_lss, aes(x = Hs, y = Tp), color = "black", linetype = "dashed", size = 1) +
  labs(title = "Physical Space (X-space)", x = "Hs [m]", y = "Tp [s]") +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 15), ylim = c(0, 25))

# Plot U-Space
p2 <- ggplot() +
  geom_point(data = df_sim, aes(x = U1, y = U2), alpha = 0.3, color = "grey") +
  geom_path(data = df_cont, aes(x = U1, y = U2), color = "red", size = 1.2) +
  geom_path(data = df_lss, aes(x = U1, y = U2), color = "black", linetype = "dashed", size = 1) +
  labs(title = "Standard Normal Space (U-space)", x = "u1", y = "u2") +
  theme_minimal() +
  coord_fixed(ratio = 1, xlim = c(-4, 6), ylim = c(-4, 6))

# Print plots
print(p1)
print(p2)

















# --- MODIFIED 5. Define Limit State Surface for max(u1) ---

# The design point is now at alpha = 0 (Rightmost point)
# The LSS is a vertical line where u1 = beta

# Define the line in U-space
# u1 is constant (equal to beta)
u1_lss <- rep(beta, 100) 
# u2 varies freely (from bottom to top of the graph)
u2_lss <- seq(-4, 4, length.out = 100)

# Transform LSS back to Physical Space (Inverse Rosenblatt)
hs_lss <- qweibull(pnorm(u1_lss), shape = shape_hs, scale = scale_hs)
tp_lss <- numeric(length(hs_lss))

for(i in 1:length(hs_lss)) {
  mu <- get_mu_ln(hs_lss[i])
  tp_lss[i] <- qlnorm(pnorm(u2_lss[i]), meanlog = mu, sdlog = sigma_ln)
}

# --- Plotting (Same as before) ---
df_lss <- data.frame(Hs = hs_lss, Tp = tp_lss, U1 = u1_lss, U2 = u2_lss)

# You can run the plotting section from the previous code using this new df_lss


# Plot Physical Space
p1 <- ggplot() +
  geom_point(data = df_sim, aes(x = Hs, y = Tp), alpha = 0.3, color = "grey") +
  geom_path(data = df_cont, aes(x = Hs, y = Tp), color = "red", size = 1.2) +
  geom_path(data = df_lss, aes(x = Hs, y = Tp), color = "black", linetype = "dashed", size = 1) +
  labs(title = "Physical Space (X-space)", x = "Hs [m]", y = "Tp [s]") +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 15), ylim = c(0, 25))

# Plot U-Space
p2 <- ggplot() +
  geom_point(data = df_sim, aes(x = U1, y = U2), alpha = 0.3, color = "grey") +
  geom_path(data = df_cont, aes(x = U1, y = U2), color = "red", size = 1.2) +
  geom_path(data = df_lss, aes(x = U1, y = U2), color = "black", linetype = "dashed", size = 1) +
  labs(title = "Standard Normal Space (U-space)", x = "u1", y = "u2") +
  theme_minimal() +
  coord_fixed(ratio = 1, xlim = c(-4, 6), ylim = c(-4, 6))

# Print plots
print(p1)
print(p2)
