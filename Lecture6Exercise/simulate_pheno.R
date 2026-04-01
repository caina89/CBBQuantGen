args <- commandArgs(trailingOnly = TRUE)
bim_file <- args[1]
n_causal <- as.numeric(args[2])
output_prefix <- args[3]

# Read the BIM file
# 1:CHR, 2:ID, 3:CM, 4:BP, 5:A1 (Effect), 6:A2 (Other)
vars <- read.table(bim_file, stringsAsFactors=FALSE)
colnames(vars) <- c("CHR", "ID", "CM", "BP", "A1", "A2")

# 1. Filter out variants with no ID or missing alleles
vars <- vars[vars$ID != "." & vars$A1 != "0", ]

# 2. Handle potential duplicate IDs by creating a unique key if necessary
# But for simplicity, we'll just unique them
vars <- vars[!duplicated(vars$ID), ]

# 3. Randomly sample
set.seed(42)
causal_vars <- vars[sample(nrow(vars), n_causal), ]

# 4. Create the effects file
# Column 1: ID
# Column 2: The EXACT allele PLINK sees as A1
# Column 3: The Weight
effects <- data.frame(
    ID = causal_vars$ID,
    A1 = causal_vars$A1,
    Beta = rnorm(n_causal, mean = 0, sd = 1)
)

write.table(effects, paste0(output_prefix, ".effects"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
