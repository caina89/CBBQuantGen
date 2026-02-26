import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

# 1. Load PCA results
pca = pd.read_csv("allchr_unrelated_pruned_pca.eigenvec", sep='\t')
pca.columns = [c.replace('#', '') for c in pca.columns]

# 2. Load Population Metadata
panel = pd.read_csv("integrated_call_samples_v3.20130502.ALL.panel", sep='\t')

# 3. Merge data
df = pd.merge(pca, panel, left_on='IID', right_on='sample')

# 4. Create PDF Report
with PdfPages('allchr_pca.pdf') as pdf:
    plt.figure(figsize=(12, 8))
    
    # Plotting by 'pop' (Specific populations like YRI, CEU, etc.)
    # Change 'hue' to 'super_pop' for Continental clusters
    sns.scatterplot(data=df, x='PC1', y='PC2', hue='pop', 
                    palette='turbo', s=40, alpha=0.8, edgecolor='none')
    
    plt.title("1000 Genomes PCA (PC1 vs PC2)")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    
    # Place legend outside the plot for clarity
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', ncol=2, fontsize='small')
    plt.grid(True, linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    pdf.savefig()
    plt.close()

print("PCA plots saved to PDF.")
