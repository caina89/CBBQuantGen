library(glmnet)

# 1. Simulate Genotypes (100 people, 10 SNPs in high LD)
set.seed(42)
base_snp <- rnorm(100)
X <- matrix(NA, nrow=100, ncol=10)
for(i in 1:10) { X[,i] <- base_snp + rnorm(100, sd=0.2) } # High correlation (LD)

# 2. Simulate Phenotype (Only the first SNP is actually causal)
y <- 2 * X[,1] + rnorm(100)

# 3. Fit the three models
# alpha = 0 (Ridge), alpha = 1 (LASSO), alpha = 0.5 (Elastic Net)
fit_ridge <- glmnet(X, y, alpha = 0)
fit_lasso <- glmnet(X, y, alpha = 1)
fit_enet  <- glmnet(X, y, alpha = 0.5)

# 4. Compare the coefficients at a specific lambda
print("Ridge Coefficients (Spreads the effect):")
print(as.vector(coef(fit_ridge, s = 0.1)[-1]))

print("LASSO Coefficients (Picks one, zeros others):")
print(as.vector(coef(fit_lasso, s = 0.1)[-1]))

print("Elastic Net Coefficients (The Balance):")
print(as.vector(coef(fit_enet, s = 0.1)[-1]))
