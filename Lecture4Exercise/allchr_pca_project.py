import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

# Load data
pca = pd.read_csv("chr20_pca_results.eigenvec", sep='\t').rename(columns={'#FID': 'FID'})
panel = pd.read_csv("integrated_call_samples_v3.20130502.ALL.panel", sep='\t')
df_orig = pd.merge(pca, panel, left_on='IID', right_on='sample')

# Load projected score (PLINK 2 .sscore)
projected = pd.read_csv("related_projection.sscore", sep='\t')
# Map score columns to PC names
projected = projected.rename(columns={'IID': 'IID', 'SCORE1_AVG': 'PC1', 'SCORE2_AVG': 'PC2'})

with PdfPages('allchr_pca_project.pdf') as pdf:
    plt.figure(figsize=(10, 7))
    
    # Background: Unrelated individuals
    sns.scatterplot(data=df_orig, x='PC1', y='PC2', hue='super_pop', alpha=0.3, palette='Set2', edgecolor=None)
    
    # Foreground: The projected individual
    plt.scatter(projected['PC1'], projected['PC2'], color='black', marker='X', s=200, label='Projected Indiv')
    
    plt.title("Projection of Filtered Individual onto 1000G PCA space")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    pdf.savefig()
    plt.close()
