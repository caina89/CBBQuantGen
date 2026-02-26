set.seed(42)
N_pop <- 100000
maf <- 0.5
h2 <- 0.5

# 1. GENERATE THE POPULATION (The Giant SNP)
genotype <- rbinom(N_pop, 2, maf)
additive_effect <- sqrt(h2 / (2 * maf * (1 - maf)))
liability <- (genotype * additive_effect) + rnorm(N_pop, sd = sqrt(1 - h2))
height <- liability + 170

# 2. DEFINE THRESHOLDS (Defining "Cases" as the Top X%)
# Note: We will use the same threshold to define "Controls" as the Bottom X%
thresholds <- c(0.5, 0.75, 0.90)
results_extreme <- data.frame()

for(k in thresholds) {
  # Population Definitions
  T_upper <- quantile(height, k)     # Cutoff for Cases
  T_lower <- quantile(height, 1 - k) # Cutoff for "Extreme" Controls
  
  pop_cases <- ifelse(height >= T_upper, 1, 0)
  
  # --- STEP A: POPULATION MODEL (The Truth) ---
  # Compares Cases (Top X%) to EVERYONE ELSE in the 100k population
  mod_pop <- glm(pop_cases ~ genotype, family = binomial(link = "logit"))
  
  # --- STEP B: EXTREME TAIL STUDY ---
  # We sample 500 from the Top Tail and 500 from the Bottom Tail
  # We DELETE the middle of the distribution entirely.
  case_idx <- sample(which(height >= T_upper), 500)
  ctrl_idx <- sample(which(height <= T_lower), 500)
  
  y_study <- c(rep(1, 500), rep(0, 500))
  x_study <- genotype[c(case_idx, ctrl_idx)]
  
  mod_study <- glm(y_study ~ x_study, family = binomial(link = "logit"))
  
  # --- STEP C: EXTRACT DATA ---
  results_extreme <- rbind(results_extreme, data.frame(
    Threshold_Tail = 1 - k, # E.g., 0.1 means Top 10% vs Bottom 10%
    Logit_Int_Pop  = coef(mod_pop)[1],
    Logit_Int_Study= coef(mod_study)[1],
    Logit_Beta_Pop = coef(mod_pop)[2],
    Logit_Beta_Study= coef(mod_study)[2],
    Beta_Inflation = (coef(mod_study)[2] / coef(mod_pop)[2])
  ))
}

# View the breakdown
# Notice how Beta_Study is much larger than Beta_Pop
print(results_extreme)
