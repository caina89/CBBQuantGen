args <- commandArgs(trailingOnly = TRUE)
bim_file <- args[1]     # 1kg_eur.bim
n_causal <- as.numeric(args[2])
output_prefix <- args[3]
h2 <- 0.5

# Load the BIM file
bim <- read.table(bim_file, header = FALSE)
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")

# 1. Randomly sample causal SNPs from the BIM
# We sample from the whole BIM to ensure they exist
set.seed(42) # For reproducibility
causal_indices <- sample(1:nrow(bim), n_causal)
causal_snps <- bim[causal_indices, ]

# 2. Assign effect sizes
# Column 1: SNP ID, Column 2: Effect Allele (A1), Column 3: Beta
effects <- data.frame(
    SNP = causal_snps$SNP,
    A1 = causal_snps$A1,
    Beta = rnorm(n_causal, mean = 0, sd = 1)
)

# Write the effects file for PLINK --score
write.table(effects, paste0(output_prefix, ".effects"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Save the list of causal SNPs for record-keeping
write.table(causal_snps$SNP, paste0(output_prefix, ".causal_list"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
