# Tutorial for Lecture 4 - Relatedness and population structure 
# Tools
We will need [bcftools](https://samtools.github.io/bcftools/bcftools.html) a set of utilities that manipulate variant calls in the [Variant Call Format (VCF)](https://samtools.github.io/hts-specs/VCFv4.1.pdf) and its binary counterpart BCF, and [plink2](https://www.cog-genomics.org/plink/2.0/) a free, open-source whole genome association analysis toolset, designed to perform a range of basic, large-scale analyses in a computationally efficient manner. for processing If you haven't already install them in the past weeks, we'd now install them using conda 
```
conda install -c bioconda bcftools plink2 
``` 
We will also be needing KING which is the software we will use to determine relatedness between individuals in the 1000 Genomes Project. 
``` 
# Download the Linux version
wget https://www.kingrelatedness.com/Linux-king.tar.gz
tar -xzvf Linux-king.tar.gz
# Move to your local bin
chmod +x king
sudo mv king /usr/local/bin/
```
# Data 
In terms of data, we will be using all variants on chr20 of 2504 individuals in the [1000 Genomes Project Phase 3](https://www.internationalgenome.org/) release to demonstrate how phased variant calls are filtered, and a filtered set of xxx variants from across all autosomes in all 2504 individuals in the 1000 Genomes Project Phase 3 (that I've performed beforehand) to demonstrate how relatedness between them may be calculated, and how population structure between them can be obtained and visualized. 
We will first download the chr20 phased variants of all individuals in the 1000 Genomes Project. Note that while Chr20 is one of the smaller chromosomes, but this file is still large. Ensure you have 10-15 GB of free space.
```
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/phase3_liftover_nygc_dir/phase3.chr20.GRCh38.GT.crossmap.vcf.gz; 
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/phase3_liftover_nygc_dir/phase3.chr20.GRCh38.GT.crossmap.vcf.gz.tbi; 
``` 
We will then download the [PLINK files in .bed, .bim and .fam format](https://www.cog-genomics.org/plink/1.9/formats#bed) for xxx variants (filtered) across all individuals in the 1000 Genomes Project [here](https://zenodo.org/records/19068199). 
```
wget https://zenodo.org/records/19068199/files/allchr.EUR.biallelicsnps.bed
wget https://zenodo.org/records/19068199/files/allchr.EUR.biallelicsnps.bim
wget https://zenodo.org/records/19068199/files/allchr.EUR.biallelicsnps.fam 
```
Finally we need the metadata of the individuals in the 1000 Genomes Project, downloadable directly from [here](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel) or using
```
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
```
# Filtering the chr20 data 
We use plink2 to identify all chr20 variants that are SNPs, biallelic, with MAF ($> 0.01$) and P-value for violation of Hardy-Weinberg Equilibrium (HWE) ($> 10^{-6}$):
* Note: we are performing the MAF and HWE P-value calculations and filtering on all individuals across many populations in 1000 Genomes Project here. Is this the right thing to do? How differently would you do it this if you were interested in a particular population?
```
vcf="phase3.chr20.GRCh38.GT.crossmap.vcf.gz"
plink2 --vcf $vcf --snps-only just-acgt \
--max-alleles 2 --min-alleles 2 --maf 0.01 --geno 0.01 --hwe 1e-6 \
--rm-dup exclude-all --set-missing-var-ids '@:#:$r:$a' \
--double-id --new-id-max-allele-len 10 missing --make-bed --threads 1 \
--out chr20_snps
``` 
We can also ask plink2 to output the calculated MAF and 
```
# Generate Allele Frequency report (.afreq) and Hardy-Weinberg Equilibrium report (.hardy)
plink2 --bfile chr20_snps --freq --hardy --out chr20_snps_stats
```
There are several things to note about the output files:
* In earlier versions of plink (e.g. plink1.9), the output of `--maf` is a `.frq` file focusing on the Minor Allele (MAF). In plink2, the output is an `.afreq` file focusing on the Alternate Allele (ALT).
* plink2 can calculate much more precise p-values ($< 10^{-300}$) for HWE violation. To avoid scientific notation, we can use the log10 modifier: `--hardy log10`
###  Inspecting the MAFs and HWE P values 
We then want to inspect the MAFs and HWE violation P values across all SNPs and individuals in 1000 Genomes Project, grouping individuals by populations. To visualize this in R, use script `chr20_maf_hwe.R` and to do this in python, use script `chr20_maf_hwe.py`. 
### Getting MAFs and HWE P values per population 
Remember, previously we are performed the MAF and HWE P-value calculations and filtering on all individuals across many populations in 1000 Genomes Project here. To obtain the MAFs and HWE violation P values of all chr20 SNPs per population (E.g. EUR and AFR), do the following: 
```
# Get EUR IDs
awk '$3=="EUR" {print $1, $1}' integrated_call_samples_v3.20130502.ALL.panel > eur_ids.txt
# Get AFR IDs
awk '$3=="AFR" {print $1, $1}' integrated_call_samples_v3.20130502.ALL.panel > afr_ids.txt
# Get stats for each
plink2 --bfile chr20_snps --keep eur_ids.txt --hardy --freq --out eur_stats
plink2 --bfile chr20_snps --keep afr_ids.txt --hardy --freq --out afr_stats
```
### Inspecting the difference between getting MAF and HWE P values between populations
We are now better able to compare the MAF and HWE P values at SNPs between populations. To visualize this in R use script `chr20_maf_hwe_compare.R`, and to do this in python, use script `chr20_maf_hwe_compare.py`. Note that we have previously already filtered the SNPs such that they all have MAF ($> 0.01$) and P-value for violation of Hardy-Weinberg Equilibrium (HWE) ($> 10^{-6}$). 
Key differences you might still observe: 
* MAF Distribution: You will likely notice that the AFR population has different distribution of MAF compared to the EUR population, reflecting the greater genetic diversity found in African populations.
* HWE Outliers: If one population has more high $-log_{10}(P)$ values than the other, it may indicate a population-specific technical artifact or a region under intense natural selection in that specific ancestry group.
# Relatedness 
Now let's use the filtered bi-allelic SNP data (MAF ($< 0.05$) and HWE P ($< 10^{-6}$)) from all 1000 Genomes individuals that I've prepared. Let's first calculate the relatedness between all pairs of individuals using KING. Because the 1000 Genomes dataset contains individuals from different populations, the KING-robust algorithm is the best choice because it is specifically designed to handle population structure without needing pre-defined allele frequencies. 
First, to get the kinship coefficient for all pairs of individuals we can use: 
```
# -b: input the PLINK .bed file (use the prefix of your files)
# --kinship: calculate the kinship coefficient (KING-robust)
# --prefix: name your output files
king -b allchr.EUR.biallelicsnps.bed --kinship --prefix allchr_kinship --cpus 1 
``` 
Second, if we only want to see pairs that are actually related (e.g., 1st, 2nd, or 3rd degree) and ignore the unrelated pairs, use the --related flag. This is much more efficient for large datasets like the 1000 Genomes.
```
# --related: specifically identifies and classifies relative pairs
king -b allchr.EUR.biallelicsnps.bed --related --prefix allchr_relatives --cpus 1 
```
### Filtering out related individuals 
To filter out relatives (those closer than 3rd degree relatives, who have relatedness = 0.044) and keep a set of unrelated individuals, we will use a two-step process: first, let KING identify which individuals to remove based on your specific threshold, and then use plink2 to create the new, "clean" dataset.
To do the first step, KING has a built-in "unrelated" command that uses a greedy algorithm to find the largest possible subset of unrelated individuals. This ensures we don't filter out people unnecessarily: If Person A is related to Person B and Person C, it will remove only Person A to save the other two, rather than removing all three.
```
# --unrelated: identifies a subset of unrelated individuals
# --degree 3: removes up to 3rd-degree relatives (kinship > 0.044)
king -b allchr.EUR.biallelicsnps.bed --unrelated --degree 3 --prefix allchr.EUR.biallelicsnps_unrelated
```
Once KING finishes, it will produce a file ending in .unrelated.id. This file contains the list of people we should keep.
```
# --vcf or --bfile: your input
# --keep: tells PLINK to only retain the IDs in the KING output file
# --make-bed: saves the new unrelated dataset
plink --bfile allchr \
      --keep allchr.EUR.biallelicsnps_unrelated.unrelated.id \
      --make-bed \
      --out allchr.EUR.biallelicsnps_unrelated 
```
# Population structure 
In population studies (like PCA) or GWAS, having relatives in the data can cause "inflation" ($\lambda > 1$), making our results look more significant than they actually are because the allele frequencies are skewed by family clusters. Now that we have filtered out related individuals, we can prepare for performing PCA, generally following a three-step process: Linkage Disequilibrium (LD) Pruning, the PCA calculation, and visualization.
### LD Pruning 
PCA is highly sensitive to "clumps" of variants that are inherited together (LD). If we don't prune your SNPs first, the PCA will reflect local LD patterns (like the MHC region) rather than global ancestry. We use PLINK to identify a subset of SNPs that are independent of one another.
```
# --indep-pairwise <window size> <step size> <r^2 threshold>
# This looks at 50kb windows, slides by 5 SNPs, and flags SNPs with r^2 > 0.2
plink2 --bfile allchr.EUR.biallelicsnps_unrelated \
      --indep-pairwise 50 5 0.2 \
      --out allchr.EUR.biallelicsnps_unrelated_prune
# Now create a new pruned dataset using the 'prune.in' list
plink2 --bfile allchr.EUR.biallelicsnps_unrelated \
      --extract allchr.EUR.biallelicsnps_unrelated_prune.prune.in \
      --make-bed \
      --out allchr.EUR.biallelicsnps_unrelated_pruned
``` 
### PCA 
Now that we have a set of independent SNPs, we run the actual Principal Component Analysis. By default, most researchers extract the first 10 or 20 PCs.
```
# --pca: calculate principal components
# Default is 20 PCs, but you can specify a number (e.g., --pca 10)
plink2 --bfile allchr.EUR.biallelicsnps_unrelated_pruned \
      --pca 10 \
      --out allchr.EUR.biallelicsnps_unrelated_pruned_pca
```
`--pca` in plink2 will produce two main files:
* `allchr.EUR.biallelicsnps_unrelated_pruned_pca.eigenvec`: This contains the actual PC coordinates for every person. Column 1 & 2 are IDs, Column 3 is PC1, Column 4 is PC2, etc.
* `allchr.EUR.biallelicsnps_unrelated_pruned_pca.eigenval`: This contains the eigenvalues. These represent the amount of variance explained by each PC.
### Visualization 
We can now visualize the PCA results, plotting PC1 vs PC2, and colouring each sample by their reported populations. To visualize this in R use script `allchr_pca.R`, and to do this in python, use script `allchr_pca.py`. 
### Projecting new individuals onto PCA 
Remember those individuals we threw out because they were related to individuals we used in the PCA? We can project them onto the PCA already performed, using SNP Loadings (the weight each SNP contributes to each PC). To do this we first have to obtain the SNP loadings for each PC using plink2
```
#Save SNP scores (loadings) of PCA run on unrelated people 
plink2 --bfile allchr.EUR.biallelicsnps_unrelated_pruned \
       --pca vscore \
       --out allchr.EUR.biallelicsnps_unrelated_pruned_pca
# This creates 'allchr.EUR.biallelicsnps_unrelated_pruned_pca.eigenvec.var', which contains the weights for each SNP.
```
Now, we need to generate the list of individuals we filtered out with the greedy algorithm - taking the difference between allchr.fam and allchr_unrelated.unrelated.id. We can do this with AWK: 
```
awk 'NR==FNR {keep[$1,$2]; next} !(($1,$2) in keep) {print $1, $2}' \
    allchr.EUR.biallelicsnps_unrelated.unrelated.id \
    allchr.EUR.biallelicsnps.fam \
    > related_indiv.txt
```
We can now apply the SNP loadings to their genotypes: 
```
# Project the related individuals (not used in PCA) onto the PCs
# --score: Applies the weights from the .eigenvec.var file to the target individual
plink2 --bfile allchr.EUR.biallelicsnps \
       --keep related_indiv.txt \
       --score allchr.EUR.biallelicsnps_unrelated_pruned_pca.eigenvec.var 2 3 header-read no-mean-imputation \
       --out related_projection
``` 
And we can now visualize the where these related individuals lie on the PCA results from unrelated individuals: we'd still be plotting all unrelatedness individuals PC1 vs PC2, and colouring them by their reported populations, but now we have their related individuals projected onto the PC plot as black points. To visualize this in R use script `allchr_pca_project.R`, and to do this in python, use script `allchr_pca_project.py`. 
