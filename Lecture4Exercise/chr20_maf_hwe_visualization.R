# Load libraries
library(ggplot2)
library(dplyr)
library(data.table)

# 1. Load the data
# Use fread for speed since 1000G files are large
freq <- fread("chr20_stats.afreq") # or .frq
hwe <- fread("chr20_stats.hardy")   # or .hwe
panel <- fread("integrated_call_samples_v3.20130502.ALL.panel")

# 2. Merge Population Info
# Note: Since the VCF-wide HWE is calculated across everyone, 
# strictly speaking, "per population" stats usually require 
# running PLINK with the --within flag. 
# For now, we will plot the global distribution.

# 3. Clean and Calculate -log10
hwe_clean <- hwe %>%
  filter(P > 0) %>% # Remove zeros to avoid infinite values
  mutate(negLogP = -log10(P))

# 4. Plot MAF Distribution
maf_plot <- ggplot(freq, aes(x = ALT_FREQS)) +
  geom_histogram(binwidth = 0.01, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(title = "Minor Allele Frequency (MAF) Distribution: Chr 20",
       x = "MAF (Alternate Allele Frequency)",
       y = "Count")

# 5. Plot -log10(HWE)
# We highlight the standard 1e-6 threshold
hwe_plot <- ggplot(hwe_clean, aes(x = negLogP)) +
  geom_histogram(bins = 50, fill = "firebrick", color = "white") +
  geom_vline(xintercept = 6, linetype = "dashed", color = "blue") +
  theme_minimal() +
  labs(title = "HWE p-value Distribution (-log10)",
       subtitle = "Blue dashed line indicates standard 1e-6 threshold",
       x = "-log10(p-value)",
       y = "Count")

# Print plots in pdf 
pdf("chr20_stats.pdf",height=5, width=5)
print(maf_plot)
print(hwe_plot)
dev.off()
