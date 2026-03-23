import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler

def simulate_and_correlate(N, M=500, structure_strength=2.0):
    # 1. Generate population structure
    # Latent ancestry gradient
    ancestry = np.linspace(-1, 1, N)
    # SNPs have different 'loadings' on this ancestry
    loadings = np.random.normal(0, 1, M)
    # Genotype matrix: Structure + Noise
    # We add structure to about 20% of SNPs to make a clear PC1
    X = np.outer(ancestry, loadings) * structure_strength + np.random.normal(0, 1, (N, M))
    
    # 2. Standardize genotypes
    scaler = StandardScaler()
    X_std = scaler.fit_transform(X)
    
    # 3. Get PC1
    pca = PCA(n_components=1)
    pc1 = pca.fit_transform(X_std).flatten()
    
    # 4. Generate NULL phenotype (random noise, no relation to X or ancestry)
    y = np.random.normal(0, 1, N)
    
    # 5. Fit Ridge Regression (PGS)
    # RidgeCV automatically finds the best alpha (lambda)
    model = RidgeCV(alphas=np.logspace(-2, 5, 10))
    model.fit(X_std, y)
    pgs = model.predict(X_std)
    
    # 6. Calculate absolute correlation
    correlation = np.abs(np.corrcoef(pgs, pc1)[0, 1])
    return correlation

# Range of sample sizes
n_sizes = [100, 250, 500, 1000, 2500, 5000]
results = []
errors = []

for n in n_sizes:
    iterations = 10
    corrs = [simulate_and_correlate(n) for _ in range(iterations)]
    results.append(np.mean(corrs))
    errors.append(np.std(corrs))

# Create the plot
plt.errorbar(n_sizes, results, yerr=errors, fmt='-o', capsize=5, color='darkblue', label='Abs(Corr(PGS, PC1))')
plt.axhline(y=0, color='black', linestyle='--', alpha=0.3)
plt.xlabel('Sample Size (N)')
plt.ylabel('Correlation with PC1')
plt.title('Effect of Sample Size on Null PGS Ancestry Correlation')
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.savefig('null_pgs_sample_size.png')

# Output data for confirmation
results_df = pd.DataFrame({'N': n_sizes, 'Mean_Corr': results, 'Std_Dev': errors})
print(results_df)
