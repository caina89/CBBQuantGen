library(glmnet)
library(ggplot2)

simulate_and_correlate <- function(N, M = 500, structure_strength = 2.0) {
  # 1. Create a "True" ancestry gradient for the individuals
  ancestry <- seq(-1, 1, length.out = N)
  
  # 2. Assign SNP 'loadings' (how much each SNP reflects ancestry)
  # Most SNPs are noise, but some carry the ancestry signal
  loadings <- rnorm(M, 0, 1)
  
  # 3. Generate Genotype Matrix X: (Ancestry * Loadings) + Random Noise
  # This creates a matrix where PC1 is a clear "Ancestry" signal
  X <- outer(ancestry, loadings) * structure_strength + matrix(rnorm(N * M), N, M)
  X <- scale(X) # Standardize
  
  # 4. Extract the top Principal Component (PC1)
  pc_res <- prcomp(X, rank. = 1)
  pc1 <- pc_res$x[, 1]
  
  # 5. Generate a NULL phenotype (pure random noise, NO genetic effect)
  y <- rnorm(N)
  
  # 6. Fit Ridge Regression (alpha = 0)
  # Use cross-validation to find the "best" lambda
  cv_fit <- cv.glmnet(X, y, alpha = 0)
  
  # 7. Calculate the PGS for these individuals
  pgs <- as.vector(predict(cv_fit, X, s = "lambda.min"))
  
  # 8. Return the absolute correlation between the Null PGS and PC1
  return(abs(cor(pgs, pc1)))
}

# --- Execution ---
n_sizes <- c(100, 250, 500, 1000, 2500)
results <- data.frame()

for (n in n_sizes) {
  message(paste("Running simulation for N =", n))
  # Run 5 iterations per sample size to get an average
  corrs <- replicate(5, simulate_and_correlate(n))
  results <- rbind(results, data.frame(N = n, Mean_Corr = mean(corrs), SD = sd(corrs)))
}

# --- Plotting ---
ggplot(results, aes(x = N, y = Mean_Corr)) +
  geom(line() + 
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean_Corr - SD, ymax = Mean_Corr + SD), width = 50) +
  labs(title = "Correlation of Null PGS with Ancestry (PC1)",
       subtitle = "As Sample Size (N) increases, the 'Spurious' Ancestry signal decreases",
       x = "Sample Size (N)",
       y = "Abs(Correlation: PGS vs PC1)") +
  theme_minimal()
