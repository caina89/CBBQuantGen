library(data.table)
library(ggplot2)
library(dplyr)

# Load data
eur_freq <- fread("eur_stats.afreq") %>% mutate(Pop = "EUR")
afr_freq <- fread("afr_stats.afreq") %>% mutate(Pop = "AFR")
eur_hwe <- fread("eur_stats.hardy") %>% mutate(Pop = "EUR")
afr_hwe <- fread("afr_stats.hardy") %>% mutate(Pop = "AFR")

# Combine datasets
freq_all <- rbind(eur_freq, afr_freq)
hwe_all <- rbind(eur_hwe, afr_hwe) %>% 
  filter(P > 0) %>% 
  mutate(negLogP = -log10(P))

# Start PDF output
pdf("chr20_maf_hwe_compare.pdf", width = 10, height = 7)

# Plot 1: MAF Comparison (Density)
p1 <- ggplot(freq_all, aes(x = ALT_FREQS, fill = Pop)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(title = "MAF Density Comparison: EUR vs AFR (Chr 20)",
       x = "Alternate Allele Frequency", y = "Density") +
  scale_fill_manual(values = c("EUR" = "#377eb8", "AFR" = "#e41a1c"))

print(p1)

# Plot 2: HWE Comparison (-log10 P)
p2 <- ggplot(hwe_all, aes(x = negLogP, fill = Pop)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  facet_wrap(~Pop) +
  geom_vline(xintercept = 6, linetype = "dashed") +
  theme_light() +
  labs(title = "Hardy-Weinberg Equilibrium Comparison",
       subtitle = "Dashed line at -log10(P) = 6",
       x = "-log10(P-value)", y = "SNP Count") +
  scale_fill_manual(values = c("EUR" = "#377eb8", "AFR" = "#e41a1c"))

print(p2)

dev.off()
