set.seed(42)
N_pop <- 100000
maf <- 0.5
h2 <- 0.5

# 1. Generate Population
genotype <- rbinom(N_pop, 2, maf)
additive_effect <- sqrt(h2 / (2 * maf * (1 - maf)))
liability <- (genotype * additive_effect) + rnorm(N_pop, sd = sqrt(1 - h2))
height <- liability + 170

# 2. Thresholds for "Case" status (Median, Top 25%, Top 10%)
thresholds <- c(0.5, 0.75, 0.9)
results_comp <- data.frame()

for(k in thresholds) {
  T_val <- quantile(height, k)
  pop_cases <- ifelse(height >= T_val, 1, 0)
  
  # --- A. Population Model (The Truth) ---
  # Prevalence is (1-k)
  mod_pop <- glm(pop_cases ~ genotype, family = binomial)
  beta_pop <- coef(mod_pop)[2]
  int_pop  <- coef(mod_pop)[1]
  
  # --- B. Ascertained Model (The Biased Sample) ---
  # Sample 500 cases and 500 controls
  cases_idx <- sample(which(pop_cases == 1), 500)
  ctrls_idx <- sample(which(pop_cases == 0), 500)
  
  y_study <- c(rep(1, 500), rep(0, 500))
  x_study <- genotype[c(cases_idx, ctrls_idx)]
  
  mod_study <- glm(y_study ~ x_study, family = binomial)
  beta_study <- coef(mod_study)[2]
  int_study  <- coef(mod_study)[1]
  
  # Store for comparison
  results_comp <- rbind(results_comp, data.frame(
    Threshold = k,
    Prevalence = 1 - k,
    Beta_Pop = beta_pop,
    Beta_Study = beta_study,
    Int_Pop = int_pop,
    Int_Study = int_study
  ))
}

print(results_comp)
