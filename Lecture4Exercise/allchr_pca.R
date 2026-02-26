library(data.table)
library(ggplot2)

# 1. Load PCA results (PLINK 2.0 format)
pca <- fread("allchr_unrelated_pruned_pca.eigenvec")
# Remove the '#' from the first column name
setnames(pca, "#FID", "FID")

# 2. Load Population Metadata
panel <- fread("integrated_call_samples_v3.20130502.ALL.panel")

# 3. Merge data on Sample ID
# PCA 'IID' matches Panel 'sample'
df <- merge(pca, panel, by.x = "IID", by.y = "sample")

# 4. Generate PDF
pdf("allchr_pca.pdf", width = 10, height = 8)

ggplot(df, aes(x = PC1, y = PC2, color = pop)) +
  geom_point(alpha = 0.7, size = 1.5) +
  theme_minimal() +
  labs(title = "1000 Genomes PCA: PC1 vs PC2",
       subtitle = "Colored by Population (Sub-population)",
       x = "Principal Component 1",
       y = "Principal Component 2",
       color = "Population") +
  # Optional: use super_pop for broader clustering (AFR, EUR, etc.)
  # aes(color = super_pop) 
  guides(color = guide_legend(ncol = 2)) # Make legend readable

dev.off()
