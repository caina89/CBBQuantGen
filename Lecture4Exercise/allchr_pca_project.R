library(ggplot2)

# Load original PCA and Panel
pca_orig <- read.table("allchr.EUR.biallelicsnps_unrelated_pruned_pca.eigenvec")
colnames(pca_orig)=c("FID", "IID",paste0("PC",c(1:10)))
panel <- read.table("integrated_call_samples_v3.20130502.ALL.panel",header=T)
df_orig <- merge(pca_orig, panel, by.x = "IID", by.y = "sample")

# Load Projected Individual (from .sscore file)
# Note: PLINK 2 score output columns are usually SCORE1_AVG, etc.
projected <- read.table("related_projection.sscore")
colnames(projected)=c("FID","IID","SCORE1_AVG", "SCORE2_AVG",paste0("PC",c(1:10)))

pdf("allchr.EUR.biallelicsnps_pca_project.pdf", width = 10, height = 8)

ggplot() +
  # Plot the background (original unrelated individuals)
  geom_point(data = df_orig, aes(x = PC1, y = PC2, color = super_pop), alpha = 0.3) +
  # Plot the projected individual as a large distinct star/point
  geom_point(data = projected, aes(x = PC1, y = PC2), color = "black", shape = 18, size = 5) +
  geom_text(data = projected, aes(x = PC1, y = PC2, label = "Related Individual"), vjust = -1.5) +
  theme_minimal() +
  labs(title = "PCA Projection: Related Individual onto 1000G Space",
       x = "PC1", y = "PC2")

dev.off()
