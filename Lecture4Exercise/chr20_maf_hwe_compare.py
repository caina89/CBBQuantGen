import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

# Load datasets
eur_f = pd.read_csv("eur_stats.afreq", sep='\t').assign(Pop='EUR')
afr_f = pd.read_csv("afr_stats.afreq", sep='\t').assign(Pop='AFR')
eur_h = pd.read_csv("eur_stats.hardy", sep='\t').assign(Pop='EUR')
afr_h = pd.read_csv("afr_stats.hardy", sep='\t').assign(Pop='AFR')

# Combine and calculate -log10(P)
freq_all = pd.concat([eur_f, afr_f])
hwe_all = pd.concat([eur_h, afr_h])
hwe_all = hwe_all[hwe_all['P'] > 0].copy()
hwe_all['neg_log10_p'] = -np.log10(hwe_all['P'])

with PdfPages('chr20_maf_hwe_compare.pdf') as pdf:
    # Page 1: MAF Comparison
    plt.figure(figsize=(10, 6))
    sns.kdeplot(data=freq_all, x="ALT_FREQS", hue="Pop", fill=True, palette="Set1", common_norm=False)
    plt.title("Allele Frequency Comparison (EUR vs AFR)")
    plt.xlabel("Alternate Allele Frequency")
    plt.ylabel("Density")
    pdf.savefig()
    plt.close()

    # Page 2: HWE Comparison
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for i, pop in enumerate(['EUR', 'AFR']):
        data = hwe_all[hwe_all['Pop'] == pop]
        sns.histplot(data['neg_log10_p'], bins=50, ax=axes[i], color='crimson' if pop=='AFR' else 'steelblue')
        axes[i].set_title(f"HWE Distribution: {pop}")
        axes[i].axvline(6, color='black', linestyle='--')
        axes[i].set_xlabel("$-log_{10}(P)$")
    
    plt.tight_layout()
    pdf.savefig()
    plt.close()

print("Python PDF Report Generated.")
