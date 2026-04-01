# Tutorial for Lecture 6 - GWAS
In this exercise we want to simulate three phenotypes: one with 1 major causal effect, one with 5 moderate causal effects, one with 1000 small causal effects, each accounting for heritablity = 0.5, and see how they perform in a GWAS.
To do this we will require: 

## Data 
* Genotypes of 1KGP European unrelated individuals we have used in Lecture 4 and Lecture 5 Excercises `allchr.EUR.biallelicsnps_unrelated`
* PCA of 1KGP European unrelated individuals performed using LD-pruned genotypes `allchr.EUR.biallelicsnps_unrelated_pruned`, in the file `allchr.EUR.biallelicsnps_unrelated_pruned_pca.eigenvec`. This will be used as fixed effect covariates. 

## Software 
We will continue to use `plink` as in previous exercises

## Step 1: Define causal SNPs and simulate phenotypes
We will create three SNP effect files. For the phenotype with 1 major SNP effect, we pick a SNP with a decent MAF so the variance is high. As before, we assume all causal SNPs are not in LD with each other, as such we use the `allchr.EUR.biallelicsnps_unrelated_pruned` genotype file for selecting causal SNPs. For each set of causal variants selectged, use the R script `simulate_pheno.R` to first generate random effect sizes for SNP effects at all causal SNPs. Then, obtain the genotypes of all causal SNPs and use the generated random effect sizes per causal SNP to make phenotypes using `plink`. Finally, use the R script `simulate_pheno.R` again to scale the phenotype such that $h^2$ is exactly 0.5. 

```
for N in 1 5 1000; do
    echo "Processing $N SNPs..."
    
    # 1. Generate effects
    Rscript simulate_pheno.R allchr.EUR.biallelicsnps_unrelated_pruned.bim $N sim_$N
    
    # 2. PLINK 2 Score
    # 'header' is omitted because our R script has no header
    # 1=ID, 2=Allele, 3=Weight
    plink2 --bfile allchr.EUR.biallelicsnps_unrelated_pruned \
           --score sim_$N.effects 1 2 3 \
           --out sim_$N \
    
    # 3. Scale Phenotype
    # PLINK2 sscore header: #IID, (sometimes FID), NMISS_ALL_CTS, SCORE1_AVG
    Rscript -e "
        df <- read.table('sim_$N.sscore', header=T, comment.char='');
        score_col <- grep('SCORE1_AVG', colnames(df));
        vg <- var(df[,score_col]);
        h2 <- 0.5;
        ve <- vg * (1 - h2) / h2;
        df\$PHENO <- df[,score_col] + rnorm(nrow(df), 0, sqrt(ve));
        # Standard 3-column output
        out <- data.frame(FID = df[[2]], IID = df[[2]], PHENO = df\$PHENO);
        write.table(out, 'sim_$N.pheno', quote=F, row.names=F, sep='\t')
    "
done
```
## Step 3: GWAS 
Now we run the GWAS for all three phenotypes using `plink` 
```
for N in 1 5 1000; do
    plink2 --bfile allchr.EUR.biallelicsnps_unrelated \
           --pheno sim_$N.txt \
           --covar allchr.EUR.biallelicsnps_unrelated_pruned_pca.eigenvec \
           --glm \
           --out gwas_results_sim$N
done
```
## Step 4: Look at results 
If we were to plot these (Manhattan plots), here is what the distribution of p-values would reveal:

* 1 Major Effect: We will see a single, massive "skyscraper" (extremely low p-value) on one chromosome. The rest of the genome will be a flat line of noise. This is typical of Mendelian-like traits or strong drug-response loci.
* 5 Moderate Effects: We will see five distinct "towers" scattered across the genome. They will be statistically significant, but the peaks won't be as high as the single major effect.
* 1000 Small Effects: This is the Infinitesimal Model. We might see many small "bumps" or maybe nothing even hits the genome-wide significance line ($5 \times 10^{-8}$) since our sample size is small (like the 1000 Genomes $N \approx 500$).

Plot your results using the `plot_manhattan_qqplots.R` in R or `plot_manhattan_qqplots.py` in python.  
## Step 5: Calculate genomic inflation factor $\lambda_{GC}$
To numerically verify if our GWAS is well-calibrated (meaning the signal is truly genetic and not due to population structure or technical artifacts), we calculate the Genomic Inflation Factor ($\lambda_{GC}$).Mathematically, $\lambda_{GC}$ is the ratio of the median observed $\chi^2$ statistic to the expected median $\chi^2$ statistic under the null hypothesis (which is approximately $0.454$). We can do this in R using the following script: 
```
# Function to calculate Lambda GC
calculate_lambda <- function(p_values) {
  # Convert P-values to Chi-squared statistics (1 degree of freedom)
  chisq <- qchisq(1 - p_values, df = 1)
  
  # Calculate Median Observed / Median Expected (0.454)
  lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
  
  return(lambda)
}

# Apply to your three simulated scenarios
files <- c("gwas_results_1.PHENO1.glm.linear", 
           "gwas_results_5.PHENO1.glm.linear", 
           "gwas_results_1000.PHENO1.glm.linear")

for (f in files) {
  data <- read.table(f, header = TRUE)
  l_val <- calculate_lambda(data$P)
  cat("Scenario:", f, "| Lambda_GC:", round(l_val, 4), "\n")
}
```
 
