library(data.table)
library(evgam)
library(extRemes)

# -----------------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------------
#data <- fread("C:/these_docs/data2024.csv")
data <- fread("C:/Users/minim/Documents/teletravail/data2024.csv")

data$angularDifference <- data$angleDiff

# -----------------------------------------------------------------------------
# 2. FIT DYNAMIC THRESHOLD (ALD)
# -----------------------------------------------------------------------------
fit_ald <- evgam(
  formula = list(
    hs ~ s(Cspd, bs = "tp", k = 5) + s(angularDifference, bs = "cc", k = 5, m = 1),
    ~ 1
  ),
  data = data,
  family = "ald",
  ald.args = list(tau = 0.96)
)

plot(fit_ald, scheme = 2)
predict(fit_ald, newdata = data.frame(Cspd = 3.01, angularDifference = 102), type = "response")

# -----------------------------------------------------------------------------
# 3. PREDICT DYNAMIC THRESHOLD
# -----------------------------------------------------------------------------
thresh_pred <- predict(fit_ald, newdata = data, type = "response")
data$threshold <- thresh_pred[, 1]

# -----------------------------------------------------------------------------
# 4. DECLUSTER
# -----------------------------------------------------------------------------
decl <- decluster(data$hs, threshold = data$threshold, r = 72)

data$excess <- pmax(decl - data$threshold, 0)
excess_df <- subset(data, excess > 0)

# -----------------------------------------------------------------------------
# 5. FIT GPD ON EXCESSES
# -----------------------------------------------------------------------------
fit_gpd <- evgam(
  formula = list(
    excess ~ s(Cspd, bs = "tp", k = 5, m = 1) +
      s(angularDifference, bs = "cc", k = 5, m = 2) +
      ti(Cspd, angularDifference, bs = c("tp", "cc"), k = c(5, 5), m = c(1, 2)),
    ~ 1
  ),
  family = "gpd",
  data = excess_df
)

AIC(fit_gpd); BIC(fit_gpd)
summary(fit_gpd)
plot(fit_gpd, scheme = 2, shade = FALSE)

par(mfrow = c(1, 1))
predict(fit_gpd, newdata = excess_df, type = "qqplot")