args <- commandArgs(trailingOnly = TRUE)
causal_file <- args[1]
output_name <- args[2]
h2 <- 0.5

# 1. Generate random effect sizes for the causal SNPs
snps <- read.table(causal_file, header = FALSE)
# Assume additive model: SNP_ID, Effect_Allele (we'll just use 'A'), Beta
# Note: In real data, you'd check the BIM file for the actual A1 allele.
effects <- data.frame(SNP = snps$V1, A1 = "A", Beta = rnorm(nrow(snps)))
write.table(effects, paste0(output_name, ".effects"), quote=F, row.names=F, col.names=F)

# 2. (Run PLINK externally to get the GRS - see Bash step below)

# 3. Scale phenotype (to be run after PLINK)
if (file.exists(paste0(output_name, ".profile"))) {
    grs_data <- read.table(paste0(output_name, ".profile"), header = TRUE)
    genetic_var <- var(grs_data$SCORE)
    
    # Calculate required environmental variance: h2 = Vg / (Vg + Ve)
    # Ve = Vg * (1 - h2) / h2
    env_var <- genetic_var * (1 - h2) / h2
    noise <- rnorm(nrow(grs_data), mean = 0, sd = sqrt(env_var))
    
    phenotype <- grs_data$SCORE + noise
    
    final_pheno <- data.frame(FID = grs_data$FID, IID = grs_data$IID, PHENO = phenotype)
    write.table(final_pheno, paste0(output_name, ".pheno"), quote=F, row.names=F, sep="\t")
}
