# Install if necessary
if (!require("qqman")) install.packages("qqman")
library(qqman)

# List of your GWAS results from the previous step
files <- c("gwas_results_1.PHENO1.glm.linear", 
           "gwas_results_5.PHENO1.glm.linear", 
           "gwas_results_1000.PHENO1.glm.linear")

# Set up a 3x2 plotting grid (3 rows for scenarios, 2 columns for Manhattan/QQ)
par(mfrow=c(3, 2))

for (f in files) {
  # Load PLINK 2 output
  gwas <- read.table(f, header=T)
  
  # PLINK 2 column names: #CHROM, POS, ID, REF, ALT, A1, TEST, OBS_CT, BETA, SE, T_STAT, P
  # We need: CHR, BP, P, SNP
  df <- data.frame(SNP=gwas$ID, CHR=gwas$X.CHROM, BP=gwas$POS, P=gwas$P)
  
  # Remove NAs for plotting
  df <- df[!is.na(df$P), ]
  
  # Title based on the file
  trait_name <- gsub("gwas_results_", "", f)
  
  # 1. Manhattan Plot
  manhattan(df, main=paste("Manhattan:", trait_name), 
            suggestiveline = -log10(1e-5), genomewideline = -log10(5e-8),
            col = c("blue4", "orange3"))
  
  # 2. QQ Plot
  qq(df$P, main=paste("QQ Plot:", trait_name))
}
