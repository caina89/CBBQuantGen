import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

# 1. Load the data (Handles PLINK 1.9 and 2.0 formats)
try:
    freq = pd.read_csv("chr20_stats.afreq", sep='\t')
    hwe = pd.read_csv("chr20_stats.hardy", sep='\t')
except FileNotFoundError:
    freq = pd.read_csv("chr20_stats.frq", sep='\s+')
    hwe = pd.read_csv("chr20_stats.hwe", sep='\s+')

# 2. Data Cleaning & Stats Calculation
hwe = hwe[hwe['P'] > 0].copy()
hwe['neg_log10_p'] = -np.log10(hwe['P'])

# 3. Create PDF and Save Plots
with PdfPages('chr20_qc_report.pdf') as pdf:
    
    # --- Page 1: MAF Distribution ---
    plt.figure(figsize=(10, 6))
    sns.histplot(freq['ALT_FREQS'], bins=50, color='steelblue')
    plt.title('Minor Allele Frequency (MAF) - Chromosome 20')
    plt.xlabel('Alternate Allele Frequency')
    plt.ylabel('SNP Count')
    pdf.savefig()  # Save current figure to PDF
    plt.close()    # Close to free memory

    # --- Page 2: HWE Distribution ---
    plt.figure(figsize=(10, 6))
    sns.histplot(hwe['neg_log10_p'], bins=50, color='firebrick')
    plt.axvline(x=6, color='blue', linestyle='--', label='QC Threshold (1e-6)')
    plt.title('Hardy-Weinberg Equilibrium (-log10 P)')
    plt.xlabel('-log10(p-value)')
    plt.ylabel('SNP Count')
    plt.legend()
    pdf.savefig()  # Save current figure to PDF
    plt.close()

print("Report saved as 'chr20_qc_report.pdf'")
