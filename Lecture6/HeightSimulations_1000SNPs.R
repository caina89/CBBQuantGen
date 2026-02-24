set.seed(42)
n <- 5000 # Increased n slightly to detect the tiny effect
m_total <- 1000
h2_total <- 0.5
maf <- 0.5

# 1. Simulate 1000 SNPs
# Genotype matrix: n x m
G <- matrix(rbinom(n * m_total, 2, maf), nrow = n, ncol = m_total)

# 2. Assign tiny effects to all SNPs to reach total h2 = 0.5
# Genetic Variance Vg = sum(2pq * a^2). If all 'a' are same: a = sqrt(Vg / (m * 2pq))
p <- maf; q <- 1-maf
a_each <- sqrt(h2_total / (m_total * 2 * p * q))
effects <- rep(a_each, m_total)

# 3. Calculate Genetic Component and Phenotype
g_component <- G %*% effects
height <- g_component + rnorm(n, sd = sqrt(1 - h2_total)) + 170

# 4. Focal SNP (SNP #1) and Binarization (75th percentile)
focal_snp <- G[, 1]
z_focal   <- as.numeric(scale(focal_snp))
cutoff    <- quantile(height, 0.75)
case_control <- ifelse(height >= cutoff, 1, 0)

# 5. Run Models
mod_lin  <- lm(scale(height) ~ z_focal)
mod_log  <- glm(case_control ~ z_focal, family = binomial(link = "logit"))
mod_prob <- glm(case_control ~ z_focal, family = binomial(link = "probit"))

# 6. Extract Z-scores
z_scores <- data.frame(
  Model = c("Linear", "Logistic", "Probit"),
  Z_stat = c(summary(mod_lin)$coefficients[2,3], 
             summary(mod_log)$coefficients[2,3], 
             summary(mod_prob)$coefficients[2,3]),
  Beta = c(coef(mod_lin)[2], coef(mod_log)[2], coef(mod_prob)[2])
)

print(z_scores)
