# Load the necessary library for Genomic REML
if (!require("rrBLUP")) install.packages("rrBLUP")
library(rrBLUP)

# 1. Reuse the Simulated Genotypes from the previous step
set.seed(42)
base_snp <- rnorm(100)
X <- matrix(NA, nrow=100, ncol=10)
for(i in 1:10) { X[,i] <- base_snp + rnorm(100, sd=0.2) }

# 2. Simulate Phenotype (Heritability ~ 0.5)
# Genetic effect + Environmental noise
g_effect <- 2 * X[,1] 
e_noise <- rnorm(100, sd = sd(g_effect)) 
y <- g_effect + e_noise

# 3. Estimate Variance Components using REML
# mixed.solve is the standard rrBLUP function for REML
reml_fit <- mixed.solve(y = y, Z = X)

# 4. Extract Variances
sigma2_g_marker <- reml_fit$Vu  # Variance per marker
sigma2_e <- reml_fit$Ve         # Residual (Environmental) variance

# 5. Calculate the REML-derived Lambda
lambda_reml <- sigma2_e / sigma2_g_marker

# --- Results Output ---
cat("--- REML Estimates ---\n")
cat("Estimated Genetic Variance (per marker):", sigma2_g_marker, "\n")
cat("Estimated Residual Variance:", sigma2_e, "\n")
cat("REML-derived Lambda (Penalty):", lambda_reml, "\n")

# Compare with the Ridge solution
# Note: rrBLUP's 'u' are the actual BLUPs (Ridge coefficients)
cat("\nTop 3 Marker Effects (BLUPs):\n")
print(head(reml_fit$u, 3))
