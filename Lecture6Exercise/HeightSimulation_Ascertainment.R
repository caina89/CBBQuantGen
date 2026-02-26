set.seed(42)
N_pop <- 100000
maf <- 0.5
h2 <- 0.5

# 1. GENERATE THE POPULATION
# We simulate a "Giant SNP" (h2=0.5) to make the mathematical behavior 
# extremely clear and easy to see without being buried in noise.
genotype <- rbinom(N_pop, 2, maf)
additive_effect <- sqrt(h2 / (2 * maf * (1 - maf)))

# Liability follows the Standard Normal assumption: Genetic + Environmental Error
# Liability ~ N(mean_genotype, variance_error)
liability <- (genotype * additive_effect) + rnorm(N_pop, sd = sqrt(1 - h2))
height <- liability + 170

# 2. DEFINE THRESHOLDS
# 0.5  = Common trait (50% prevalence)
# 0.9  = Rare trait (10% prevalence)
# 0.99 = Very rare trait (1% prevalence)
thresholds <- c(0.5, 0.9, 0.99)
results_inv <- data.frame()

for(k in thresholds) {
  # Determine the height cutoff for this specific prevalence
  T_val <- quantile(height, k)
  pop_cases <- ifelse(height >= T_val, 1, 0)
  
  # --- STEP A: THE POPULATION MODELS ("THE TRUTH") ---
  # These models look at all 100,000 people. This is the "true" genetic 
  # effect on the Logit and Probit scales in the real world.
  mod_pop_log  <- glm(pop_cases ~ genotype, family = binomial(link = "logit"))
  mod_pop_prob <- glm(pop_cases ~ genotype, family = binomial(link = "probit"))
  
  # --- STEP B: THE ASCERTAINED STUDY (THE "BIASED" SAMPLE) ---
  # We artificially create a 50/50 study by picking 1000 Cases and 1000 Controls.
  # This mimics a typical clinical GWAS where you over-sample sick people.
  case_idx <- sample(which(pop_cases == 1), 1000)
  ctrl_idx <- sample(which(pop_cases == 0), 1000)
  
  y_study  <- c(rep(1, 1000), rep(0, 1000))
  x_study  <- genotype[c(case_idx, ctrl_idx)]
  
  mod_study_log  <- glm(y_study ~ x_study, family = binomial(link = "logit"))
  mod_study_prob <- glm(y_study ~ x_study, family = binomial(link = "probit"))
  
  # --- STEP C: EXTRACT AND COMPARE ---
  results_inv <- rbind(results_inv, data.frame(
    Threshold    = k,
    Prevalence   = 1 - k,
    
    # Logistic Results
    Logit_Int_Pop   = coef(mod_pop_log)[1],
    Logit_Int_Study = coef(mod_study_log)[1],
    Logit_Beta_Pop  = coef(mod_pop_log)[2],
    Logit_Beta_Study= coef(mod_study_log)[2],
    
    # Probit Results
    Probit_Int_Pop   = coef(mod_pop_prob)[1],
    Probit_Int_Study = coef(mod_study_prob)[1],
    Probit_Beta_Pop  = coef(mod_pop_prob)[2],
    Probit_Beta_Study= coef(mod_study_prob)[2]
  ))
}

# View the full breakdown
# Look at how Logit_Beta_Study stays close to Logit_Beta_Pop
# Look at how Probit_Beta_Study starts to drift as the Threshold moves
print(results_inv)
