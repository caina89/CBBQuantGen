# Simulating a phenotype with heritability of 0.5 
We want to simulate a phenotype with SNP heritability = 0.5 with the European unrelated individuals in 1KGP Phase 3 (which we have used in Lecture 4 Exercise). To do so we need the following:

## Data 
1000 Genomes Phase 3 files in PLINK binary format (.bed, .bim, .fam). We will use the `allchr.EUR.biallelicsnps_unrelated.bed` and related `.bim` and `.fam` files that you have created last week. 

## Softwares 
Download [GCTA](https://yanglab.westlake.edu.cn/software/gcta) using the following
```
cd ~/bin
wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip
unzip gcta-1.94.1-linux-kernel-3-x86_64.zip
# Make the binary executable
chmod +x gcta64
# Test the installation
./gcta64 --help | head -n 15
```
## Steps 
### Step 1: Simulate phenotype and causal SNPs 
To simulate phenotype with 1000 causal variants (SNPs) where each causal SNP is not in LD with another one (i.e. only one causal SNP per LD block), we'd simulate the causal effects using the LD-pruned genotype file we obtained in Lecture 4 Exercises `allchr.EUR.biallelicsnps_unrelated_pruned`
```
# 1. Pick 1000 random SNPs from the .bim file
shuf -n 1000 allchr.EUR.biallelicsnps_unrelated_pruned.bim | awk '{print $2}' > causal_snps.txt

# 2. Generate random effect sizes (Beta ~ Normal(0,1))
awk '{print $1, "A1", rand()}' causal_snps.txt > causal_effects.txt

# 3. Calculate the Raw Genetic Score (G) using PLINK 2 in all individuals
# 'cols=scoresums' gives the unscaled sum of (allele_count * beta)
plink2 --bfile allchr.EUR.biallelicsnps_unrelated_pruned \
       --score causal_effects.txt 1 2 3 cols=scoresums \
       --out genetic_values
```
### Step 2: Adjust for heritability (in R)
Now we take that genetic score and add enough noise to make $h^2 = 0.5$. Since $h^2 = \frac{Var(G)}{Var(G) + Var(E)}$, for $h^2=0.5$, we need $Var(E) = Var(G)$.
In R do: 
```
# Load the PLINK 2 score output
gv <- read.table("genetic_values.sscore", header=T)
G <- gv$SCORE1_SUM

# Calculate environmental noise where Var(E) = Var(G)
set.seed(123)
E <- rnorm(length(G), mean=0, sd=sd(G))

# Final Phenotype Y = G + E
Y <- G + E

# Save in GCTA/PLINK format (FID, IID, Pheno)
final_pheno <- data.frame(FID=gv$#FID, IID=gv$IID, Y=Y)
write.table(final_pheno, "simulated_h2_0.5.phen", row.names=F, col.names=F, quote=F)
```
### Step 3: Estimate heritability using GCTA (unpruned SNPs) 
We now want to estimate heritability of this simulated phenotype using GCTA, using all SNPs in `allchr.EUR.biallelicsnps_unrelated`
```
# Build GRM
gcta64 --bfile allchr.EUR.biallelicsnps_unrelated --make-grm --out grm_all
# Estimate h2
gcta64 --grm grm_all --pheno simulated_h2_0.5.phen --reml --out h2_unpruned
```
### Step 4: Estimate heritability using GCTA (LD-pruned SNPs)
We now want to estimate heritability of this simulated phenotype using GCTA, but this time using LD-pruned SNPs we obtained in Lecture 4 Exercises `allchr.EUR.biallelicsnps_unrelated_pruned`
```
# Build GRM with pruned SNPs
gcta64 --bfile allchr.EUR.biallelicsnps_unrelated_pruned --make-grm --out grm_pruned
# Estimate h2
gcta64 --grm grm_pruned --pheno simulated_h2_0.5.phen --reml --out h2_pruned
```

