# Tutorial for Lecture 9 - Prediction 
In this exercise we use the three simulated phenotypes from exercise in Lecture 6: one with 1 major causal effect, one with 5 moderate causal effects, one with 1000 small causal effects, each accounting for heritablity = 0.5. You'd have already performed their GWAS in lecture 6. 

## Data 
* Genotypes of 1KGP European unrelated individuals we have used in Lecture 4 and Lecture 5 Excercises `allchr.EUR.biallelicsnps_unrelated`
* PCA of 1KGP European unrelated individuals performed using LD-pruned genotypes `allchr.EUR.biallelicsnps_unrelated_pruned`, in the file `allchr.EUR.biallelicsnps_unrelated_pruned_pca.eigenvec`. This will be used as fixed effect covariates. 

## Software 
We will continue to use `plink` as in previous exercises
We will also use `prsice2` (pronounced "precise 2") for P+T PRS calculations
```
cd ~/bin 
wget https://github.com/choishingwan/PRSice/releases/download/2.3.5/PRSice_linux.zip
unzip PRSice_linux.zip
``` 
## Step 1: Splitting the unrelated EUR dataset in 1000G Phase 3 
We first want to split the unrelated EUR data in 1000G Phase 3 into a 80% GWAS discovery set, and a 20% PRS testing set.  
```
## Get the 80/20 split of the 1000G Phase 3 EUR dataset 
Rscript -e "
    fam=read.table("allchr.EUR.biallelicsnps_unrelated.fam");
    fam_dis=fam[sample(nrow(fam),0.8*nrow(fam),replace=F),];
    fam_test=fam[-which(fam$V1%in%fam_dis$V1),];
    write.table(fam_dis, 'allchr.EUR.biallelicsnps_unrelated.discovery.fam', quote=F, row.names=F, sep='\t');
    write.table(fam_test, 'allchr.EUR.biallelicsnps_unrelated.test.fam', quote=F, row.names=F, sep='\t')
"
## Output the 80/20 split genotype files
plink2 --bfile allchr.EUR.biallelicsnps_unrelated \
       --keep allchr.EUR.biallelicsnps_unrelated.discovery.fam \
       --make-bed \
       --out allchr.EUR.biallelicsnps_unrelated.discovery
plink2 --bfile allchr.EUR.biallelicsnps_unrelated \
       --keep allchr.EUR.biallelicsnps_unrelated.test.fam \
       --make-bed \
       --out allchr.EUR.biallelicsnps_unrelated.test
## perform GWAS on all three simulated phenotypes in the discovery set only
for N in 1 5 1000; do
    plink2 --bfile allchr.EUR.biallelicsnps_unrelated.discovery \
           --pheno sim_$N.pheno \
           --covar allchr.EUR.biallelicsnps_unrelated_pruned_pca.eigenvec \
           --glm \
           --out discovery_gwas_results_sim$N
done
## 
```
## Step 2: Getting PRS in the test 20% using P+T method  
Here we want to try the C/P+T method where clumping is used to keep the top SNP (lowest P value) within a LD window (--clump-kb 250kb) where all SNPs with LD r2 > 0.1 are removed (--clump-r2 0.1), and all SNPs below P value of 1 (essentially all SNPs) are eligible to be the kept SNP (--clump-p 1.000000).  Note that the default PRS model is $PRS_j = \sum_i{\frac{S_i\times G_{ij}}{M_j}}$ (with --score avg) while in our simulated situations this may or may not work well for all architectures. We will therefore try the other PRS models $PRS_j = \sum_i{S_i\times G_{ij}}$ (with --score sum) and PRS_j = \frac{\sum_i({S_i\times G_{ij}}) - \text{Mean}(PRS)}{\text{SD}(PRS)} (with --score std). Also note we need to put in PCs as covariates as we have learnt that when most SNPs in the genome and included in PRS have no effect on phenotype (which is the case in our model where max 1000 causal SNPs), PRS defaults to capturing largest axes of variation in the genome (PCs). 
```
## try the default model with --score avg 
for N in 1 5 1000; do
    for model in $(echo avg sum std); do 
        Rscript ~/bin/PRSice.R \
            --prsice ~/bin/PRSice_linux \
            --score $model \
            --base discovery_gwas_results_sim$N \
            --target allchr.EUR.biallelicsnps_unrelated.test \
            --binary-target F \
            --clump-kb 250kb --clump-p 1.000000 --clump-r2 0.100000 \
            --interval 5e-05 --upper 0.5 --lower 5e-08 --num-auto 22 \
            --covar allchr.EUR.biallelicsnps_unrelated_pruned_pca.eigenvec \
            --column-adapter ID:ID,BP:POS,CHR:#CHROM,A1:A1,A2:AX,P:P,OR:OR \
            --out test_gwas_results_sim${N}_model${model}\
    done
done 
```

## Lambda estimation using REML 
REML uses likelihood function to estimate the Genetic Variance ($\sigma_g^2$) and the Residual Variance ($\sigma_e^2$). Their ratio then defines the $\lambda$ for BLUP/Ridge Regression. 
To demonstrate this we use the p3d or lme4 R package, adapted for genomic data in the R package rrBLUP, as shown in `reml_lambda.R`. 

In this script, the mixed.solve function performs Restricted Maximum Likelihood (REML). It looks at the distribution of the phenotypes $y$ and asks: "Given the relationship structure of these 10 SNPs, what values of $\sigma_g^2$ and $\sigma_e^2$ make this specific pattern of $y$ most likely?" It ignores the fixed effects (the mean) to ensure the variance estimates aren't biased.

### The Relationship to $\lambda$
* If the data is very noisy, REML will estimate a large $\sigma_e^2$.
* A large $\lambda$ automatically "shrinks" your marker effects in the BLUP solution.

### Merits of using REML 
REML is more stable because it uses the mathematical properties of the entire dataset at once to find the "biological" $\lambda$. As long as your assumption of a Normal distribution for SNP effects is mostly true, the REML $\lambda$ will be very close to the "optimal" prediction $\lambda$.

## Elastic net 
In Elastic Net the penalty term is: 

$\text{Penalty} = \lambda\left[\alpha\underbrace{\sum|\beta|}_{L_1\text{(LASSO)}} + (1 - \alpha)\underbrace{\frac{1}{2}\sum\beta^2}{L_2\text{ (Ridge)}}\right]$
* When $\alpha = 1$: The $L_2$ part becomes $0$. You have a Pure LASSO model.
* When $\alpha = 0$: The $L_1$ part becomes $0$. You have a Pure Ridge model.
* When $\alpha = 0.5$: You are giving equal "weight" to both penalties.
  
To compare coefficients ($\beta$) given to SNPs in Ridge Regression, LASSO and Elastic Net, we can use the script shown in `compare_coefficients.R` which uses the glmnet, a package for penalized regression in R.

When you run the script `compare_coefficients.R` you will see that: 

Ridge ($\alpha = 0$)
* Result: $[0.21, 0.19, 0.20, 0.18, \dots]$
* Behavior: Because the SNPs are in LD, Ridge refuses to pick a favorite. It sees 10 SNPs that all look like the causal variant and spreads the credit almost equally among them. This captures the heritability of the region but fails to identify the specific causal SNP.

LASSO ($\alpha = 1$)
* Result: $[1.85, 0, 0, 0, \dots]$
* Behavior: LASSO sees the redundancy and "fires" 9 of the SNPs. It picks the one that fits the data slightly better and gives it a large weight, setting the rest to exactly zero. This handles the kurtosis (sparsity) perfectly but can be unstable if that one SNP has a measurement error.

Elastic Net ($\alpha = 0.5$)
* Result: $[0.85, 0.42, 0.38, 0, 0, \dots]$
* Behavior: This is the compromise. It zeroes out the SNPs that are clearly noise (Selection), but for the SNPs that are highly correlated and close to the signal, it groups them together (Grouping). It keeps a few strong candidates in the model, providing a more stable prediction than LASSO.
