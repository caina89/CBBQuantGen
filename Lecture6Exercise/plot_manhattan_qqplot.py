import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_gwas(file_path, title):
    # Load PLINK 2 data
    df = pd.read_csv(file_path, sep='\t')
    df['-log10P'] = -np.log10(df['P'])
    df = df.dropna(subset=['P'])

    # 1. Manhattan Plot Logic
    df['ind'] = range(len(df))
    df_grouped = df.groupby('#CHROM')
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
    
    colors = ['#272e38', '#7d7d7d']
    x_labels = []
    x_ticks = []
    
    for i, (name, group) in enumerate(df_grouped):
        group.plot(kind='scatter', x='ind', y='-log10P', color=colors[i % 2], ax=ax1, s=10)
        x_labels.append(name)
        x_ticks.append(group['ind'].iloc[len(group)//2])
    
    ax1.set_xticks(x_ticks)
    ax1.set_xticklabels(x_labels)
    ax1.set_title(f"Manhattan: {title}")
    ax1.axhline(y=-np.log10(5e-8), color='r', linestyle='--') # Genome-wide significance

    # 2. QQ Plot Logic
    observed = sorted(df['P'])
    expected = np.linspace(1/len(observed), 1, len(observed))
    ax2.scatter(-np.log10(expected), -np.log10(observed), s=10)
    ax2.plot([0, max(-np.log10(expected))], [0, max(-np.log10(expected))], color='red', linestyle='--')
    ax2.set_xlabel('Expected -log10(P)')
    ax2.set_ylabel('Observed -log10(P)')
    ax2.set_title(f"QQ Plot: {title}")
    
    plt.tight_layout()
    plt.show()

# Run for the 1000-SNP scenario
plot_gwas('gwas_results_1000.PHENO1.glm.linear', '1000 Small Effects')
