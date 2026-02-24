# 1. Setup Data
set.seed(42)
n <- 1000
maf <- 0.5
h2 <- 0.5
mu <- 170

# Calculate additive effect
p <- maf; q <- 1-maf
a <- sqrt(h2 / (2 * p * q))

# Simulate Genotype and Height
genotype <- rbinom(n, 2, maf)
z_geno   <- as.numeric(scale(genotype)) # Pre-scale to avoid naming errors
liability <- (z_geno * sqrt(h2)) + rnorm(n, sd = sqrt(1 - h2))
height    <- (liability * 1) + mu # Using SD=1 for simplicity

# Threshold at 75th percentile
cutoff <- quantile(height, 0.75)
case_control <- ifelse(height >= cutoff, 1, 0)

# 2. Perform Regressions
mod_lin  <- lm(scale(height) ~ z_geno)
mod_log  <- glm(case_control ~ z_geno, family = binomial(link = "logit"))
mod_prob <- glm(case_control ~ z_geno, family = binomial(link = "probit"))

# 3. Results Table
results <- data.frame(
  Model     = c("Linear", "Logistic", "Probit"),
  Beta      = c(coef(mod_lin)[2], coef(mod_log)[2], coef(mod_prob)[2]),
  Error_Var = c(summary(mod_lin)$sigma^2, pi^2/3, 1.0),
  P_Value   = c(summary(mod_lin)$coefficients[2,4], 
                summary(mod_log)$coefficients[2,4], 
                summary(mod_prob)$coefficients[2,4])
)
print(results)

# 4. Plots
par(mfrow=c(2,2))

# (a) Linear: Height Distribution
hist(height, breaks=30, col="lightblue", border="white", main="Height & 75th Cutoff")
abline(v=cutoff, col="red", lwd=2, lty=2)

# (b) Linear: Genotype vs Phenotype
plot(jitter(genotype), height, pch=16, col=rgb(0,0,0,0.2), main="Linear Fit")
abline(lm(height ~ genotype), col="blue", lwd=2)

# (c) Logistic & Probit Curves
plot(z_geno, case_control, pch=16, col=rgb(0,0,0,0.1), 
     main="S-Curves (Logit vs Probit)", xlab="Std Genotype", ylab="P(Case)")

# Create smooth range for curve plotting
x_seq <- seq(min(z_geno), max(z_geno), length.out=100)
y_log  <- predict(mod_log, newdata=data.frame(z_geno=x_seq), type="response")
y_prob <- predict(mod_prob, newdata=data.frame(z_geno=x_seq), type="response")

lines(x_seq, y_log, col="firebrick", lwd=3)
lines(x_seq, y_prob, col="royalblue", lwd=3, lty=2)
legend("topleft", legend=c("Logit", "Probit"), col=c("firebrick", "royalblue"), lty=1:2, lwd=2)
