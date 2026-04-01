args <- commandArgs(trailingOnly = TRUE)
bim_pvar_file <- args[1]
n_causal <- as.numeric(args[2])
output_prefix <- args[3]

# Read the variant file (handles .bim or .pvar)
# PLINK2 pvar files might have a header starting with #
is_pvar <- grepl(".pvar$", bim_pvar_file)
if(is_pvar){
  vars <- read.table(bim_pvar_file, comment.char = "#")
  colnames(vars)[1:5] <- c("CHROM", "POS", "ID", "REF", "ALT")
} else {
  vars <- read.table(bim_pvar_file)
  colnames(vars)[1:6] <- c("CHR", "SNP", "CM", "BP", "A1", "A2")
  colnames(vars)[2] <- "ID" # Standardize column name
}

# 1. Randomly sample causal SNPs
set.seed(42)
# Filter out variants with missing IDs if necessary
vars <- vars[vars$ID != ".", ]
causal_indices <- sample(1:nrow(vars), n_causal)
causal_vars <- vars[causal_indices, ]

# 2. Create the effects file (SNP ID, then Beta)
effects <- data.frame(
    ID = causal_vars$ID,
    Beta = rnorm(n_causal, mean = 0, sd = 1)
)

write.table(effects, paste0(output_prefix, ".effects"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
