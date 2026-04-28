# Tutorial for Lecture 9 - Prediction 
In this exercise we use the three simulated phenotypes from exercise in Lecture 6: one with 1 major causal effect, one with 5 moderate causal effects, one with 1000 small causal effects, each accounting for heritablity = 0.5. You'd have already performed their GWAS in lecture 6. 

## Data 
* Genotypes of 1KGP European unrelated individuals we have used in Lecture 4 and Lecture 5 Excercises `allchr.EUR.biallelicsnps_unrelated`
* PCA of 1KGP European unrelated individuals performed using LD-pruned genotypes `allchr.EUR.biallelicsnps_unrelated_pruned`, in the file `allchr.EUR.biallelicsnps_unrelated_pruned_pca.eigenvec`. This will be used as fixed effect covariates. 

## Software 
We will continue to use `plink` as in previous exercises
We will also use `prsice2` (pronounced "precise 2") for P+T PRS, and `gcta64` for BLUP estimates 
```
cd ~/bin 
wget https://github.com/choishingwan/PRSice/releases/download/2.3.5/PRSice_linux.zip
unzip PRSice_linux.zip
wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.95.1-linux-x86_64.zip
unzip gcta-1.95.1-linux-x86_64.zip
```
Note that gcta64 may throw an error in Euler as shown below:
```
> ~/bin/gcta-1.95.1-linux-x86_64/gcta64
Error: No suitable fusermount binary found on the $PATH
Error: $FUSERMOUNT_PROG not set

Cannot mount AppImage, please check your FUSE setup.
You might still be able to extract the contents of this AppImage 
if you run it with the --appimage-extract option. 
See https://github.com/AppImage/AppImageKit/wiki/FUSE 
for more information
open dir error: No such file or directory
```
If so, run the following 
```
~/bin/gcta-1.95.1-linux-x86_64/gcta64 --appimage-extract
```
This will create a new directory named `squashfs-root`. Inside that new directory, you will find the actual executable. It is usually located at `squashfs-root/AppRun` or `squashfs-root/usr/bin/gcta64`. 

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
           --glm hide-covar \
           --out discovery_gwas_results_sim$N
done
## 
```
## Step 2: Getting PRS in the test 20% using P+T method  
Here we want to try the C/P+T method where 
* clumping is used to keep the top SNP (lowest P value) within a LD window (--clump-kb 250kb) where all SNPs with LD r2 > 0.1 are removed (--clump-r2 0.1)
* all SNPs below P value of 1 (essentially all SNPs) are eligible to be the kept SNP (--clump-p 1.000000)
* Note that the default PRS model is $PRS_j = \sum_i{\frac{S_i\times G_{ij}}{M_j}}$ (with --score avg) while in our simulated situations this may or may not work well for all architectures. We will therefore try the other PRS models $PRS_j = \sum_i{S_i\times G_{ij}}$ (with --score sum) and $PRS_j = \frac{\sum_i({S_i\times G_{ij}}) - \text{Mean}(PRS)}{\text{SD}(PRS)}$ (with --score std).
* Also note we need to put in PCs as covariates as we have learnt that when most SNPs in the genome and included in PRS have no effect on phenotype (which is the case in our model where max 1000 causal SNPs), PRS defaults to capturing largest axes of variation in the genome (PCs). 
```
## try the default model with --score avg 
for N in 1 5 1000; do
    for model in $(echo avg sum std); do 
        ~/bin/PRSice_linux --score $model \
            --base discovery_gwas_results_sim$N.PHENO.glm.linear \
            --a1 A1 --a2 A0 --snp ID --bp POS  --stat BETA --pvalue P --base-maf A1_FREQ:0.05 \
            --target allchr.EUR.biallelicsnps_unrelated.test \
            --binary-target F \
            --pheno sim_$N.pheno \
            --clump-kb 250kb --clump-p 1.000000 --clump-r2 0.100000 \
            --interval 5e-05 --upper 0.5 --lower 5e-08 --num-auto 22 \
            --cov allchr.EUR.biallelicsnps_unrelated_pruned_pca.eigenvec \
            --out test_gwas_results_sim${N}_model${model}
    done
done 
```
We will get the following files: 
* `.summary` files will show us the PRS incremental $R^2$ (full model $R^2$ - null model $R^2$) at the P-value threshold that gives the highest $R^2$.
* `.best` files will give us the PRS values for all individuals in the 20% test dataset at the best P-value threshold that gave the highest $R^2$.
* `.prsice` files will give us the PRS incremental $R^2$ (full model $R^2$ - null model $R^2$) all P-value thresholds for comparison.
For example, the `.summary` files for the avg and sum models for the 5 SNP model look like this:
```
> head test_gwas_results_sim5_modelavg.summary
Phenotype	Set	Threshold	PRS.R2	Full.R2	Null.R2	Prevalence	Coefficient	Standard.Error	P	Num_SNP
-	Base	0.00155005	0.0616424	0.211573	0.15978	-	-90.872	51.7161	0.0840831	1442
> head test_gwas_results_sim5_modelsum.summary
Phenotype	Set	Threshold	PRS.R2	Full.R2	Null.R2	Prevalence	Coefficient	Standard.Error	P	Num_SNP
-	Base	0.00155005	0.0616424	0.211573	0.15978	-	-0.031509	0.0179321	0.0840831	1442
> head test_gwas_results_sim5_modelstd.summary
Phenotype	Set	Threshold	PRS.R2	Full.R2	Null.R2	Prevalence	Coefficient	Standard.Error	P	Num_SNP
-	Base	0.00155005	0.0616424	0.211573	0.15978	-	-0.0172822	0.00983546	0.0840831	1442
```
Ok seems no difference in this instance... maybe that's why defaults are good. :) We haven't touched a fringe situation where this matters it seems. Let's now look at the three different phenotypes that are all supposed to have h2 = 0.5 
```
> head test_gwas_results_sim1_modelavg.summary
Phenotype	Set	Threshold	PRS.R2	Full.R2	Null.R2	Prevalence	Coefficient	Standard.Error	P	Num_SNP
-	Base	0.00130005	0.108167	0.186389	0.0877093	-	167.206	70.0321	0.0201826	1471
> head test_gwas_results_sim5_modelsum.summary
Phenotype	Set	Threshold	PRS.R2	Full.R2	Null.R2	Prevalence	Coefficient	Standard.Error	P	Num_SNP
-	Base	0.00155005	0.0616424	0.211573	0.15978	-	-90.872	51.7161	0.0840831	1442
> head test_gwas_results_sim1000_modelstd.summary
Phenotype	Set	Threshold	PRS.R2	Full.R2	Null.R2	Prevalence	Coefficient	Standard.Error	P	Num_SNP
-	Base	0.00105005	0.191741	0.270406	0.0973269	-	104.022	31.1525	0.00146065	903
```
Ok now there are big differences - in terms of the P-value threshold at which C/P+T model decides to draw the line for best $R^2$, as well as the best $R^2$ itself. Notice none of them are close to h2=0.5, though it's pretty great for a small discovery cohort of N=231. Usually we don't get such high $R^2$ with much larger sample sizes - why? Our simulated phenotypes have much larger simulated SNP effects than SNP effects on most complex phenotypes.  
## Step 3: Try getting BLUPs rather than C/P+T PRS predictions 
Here since we have a small discovery and test cohort we can also try to get the BLUPs of this simulated phenotype and ask how well it predicts relatives to C/P+T. To do that we first need to format the summary statistics into what GCTA wants 
```
for N in 1 5 1000; do
    awk 'BEGIN{OFS="\t"; print "SNP","A1","A2","freq","b","se","p","N"} 
        NR>1 {print $3, $7, $4, $9, $12, $13, $15, $11}' \
       discovery_gwas_results_sim$N.PHENO.glm.linear > discovery_gwas_results_sim${N}_for_gcta.ma
done 
```
We then calculate the BLUP weights while adjusting for LD. We need a reference panel (like the genotypes of your discovery set or all of 1000 Genomes EUR as shown) to provide the LD structure. Here gcta64 learns what the shrinkage parameter lambda should be, and then use that to get the BLUP weights at every SNP. The shrinkage value is calculated as $\lambda = m(1/h^2 - 1)$, where $m$ is the number of SNPs and $h^2$ is the heritability. If unknown, --cojo-sblup 0.01 is often used as a starting heuristic.
```
for N in 1 5 1000; do
    ~/bin/gcta-1.95.1-linux-x86_64/squashfs-root/AppRun --bfile allchr.EUR.biallelicsnps_unrelated \
        --cojo-file discovery_gwas_results_sim${N}_for_gcta.ma \
        --cojo-sblup 0.01 \
        --thread-num 1 \
        --out discovery_gwas_results_sim${N}
done 
```
Finally, we have to apply the weights to our test dataset. The output from the previous step will be a file named blup_weights.sblup.cojo.
```
for N in 1 5 1000; do
    ~/bin/gcta-1.95.1-linux-x86_64/squashfs-root/AppRun --bfile allchr.EUR.biallelicsnps_unrelated.test \
        --score discovery_gwas_results_sim${N}.sblup.cojo 1 2 4 \
        --out discovery_gwas_results_sim${N}_blup_prs
```
## Step 4: Compare BLUP vs C/P+T PRS 
From PRSice2 we get incremental R2 for all P value thresholds between 0.5 and 5e-08 at intervals of 5e-05 (--interval 5e-05 --upper 0.5 --lower 5e-08). PRSice2 even tells you which is best in the .best file output. But we need to calculate this for BLUPs to compare them with the C/P+T, using the following R script: 
```
for (N in c(1,5,1000)){
    # Load data
    data <- read.table(paste0("discovery_gwas_results_sim",N,"_blup_prs.profile"), header=T)
    pheno <- read.table(paste0("sim_",N,".pheno"), header=F)  
    colnames(pheno)=c("FID","IID","pheno")
    covar <- read.table("allchr.EUR.biallelicsnps_unrelated_pruned_pca.eigenvec", header=F)
    colnames(covar)=c("FID","IID",paste0("PC",1:10))
    merged <- merge(data, pheno, by="IID")
    
    # 1. Null Model (Covariates only)
    null_model <- lm(pheno ~ PC1 + PC2 + PC3, data=merged)
    r2_null <- summary(null_model)$r.squared
    
    # 2. Full Model (Covariates + PRS)
    full_model <- lm(pheno ~ Age + Sex + PC1 + PC2 + PC3 + SCORE, data=merged)
    r2_full <- summary(full_model)$r.squared
    
    # 3. Calculate Incremental R2
    incremental_r2 <- r2_full - r2_null
    print(paste("The PRS explains", round(incremental_r2 * 100, 2), "% of the phenotypic variance."))
}
``` 

# Extra stuff 
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
