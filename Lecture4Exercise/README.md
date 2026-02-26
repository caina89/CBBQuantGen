# Tutorial for Lecture 4 - Relatedness and population structure 
## Tools
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
## Data 
In terms of data, we will be using all variants on chr20 of 2504 individuals in the [1000 Genomes Project Phase 3](https://www.internationalgenome.org/) release to demonstrate how phased variant calls are filtered, and a filtered set of xxx variants from across all autosomes in all 2504 individuals in the 1000 Genomes Project Phase 3 (that I've performed beforehand) to demonstrate how relatedness between them may be calculated, and how population structure between them can be obtained and visualized. 
We will first download the chr20 phased variants of all individuals in the 1000 Genomes Project. Note that while Chr20 is one of the smaller chromosomes, but this file is still large. Ensure you have 10-15 GB of free space.
```
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
``` 
We will then download the [PLINK files in .bed, .bim and .fam format](https://www.cog-genomics.org/plink/1.9/formats#bed) for xxx variants (filtered) across all individuals in the 1000 Genomes Project [here](link). 
```
```
Finally we need the metadata of the individuals in the 1000 Genomes Project, downloadable directly from [here](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel) or using
```
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
```
## Filtering the chr20 data 
We first remove from all chr20 variants those that are non-SNPs, multiallelic, those that didn't pass the initial 1000 Genomes quality checks (indicated in the "FILTER" column in the VCF file). We can do this using bcftools: 
```
# -f PASS: Only keep variants that passed initial QC
# -v snps: Keep only SNPs (removes Indels)
# -m2 -M2: Biallelic sites only
bcftools view -f PASS -v snps -m2 -M2 \
  ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
  -Oz -o chr20_snps_only.vcf.gz
```
## Further filtering the chr20 data on variant statistics  
We then use plink2 to calculate the MAF and P-value for violation of Hardy-Weinberg Equilibrium (HWE) at each SNP, and filter out those with MAF ($< 0.05$) and HWE P ($< 10^{-6}$). 
* Note: we are performing the MAF and HWE P-value calculations and filtering on all individuals across many populations in 1000 Genomes Project here. Is this the right thing to do? How differently would you do it this if you were interested in a particular population? 
```
# --vcf: Input your filtered VCF
# --maf 0.05: Filter for Minor Allele Frequency > 5%
# --hwe 1e-6: Filter for Hardy-Weinberg p-value > 10^-6
# --make-bed: Create the standard .bed/.bim/.fam binary fileset
# --out: The prefix for your final files
plink2 --vcf chr20_snps_only.vcf.gz \
      --maf 0.05 \
      --hwe 1e-6 \
      --make-bed \
      --out chr20_final_cleaned
```
We can also ask plink2 to output the calculated MAF and 
```
# Generate Allele Frequency report (.afreq) and Hardy-Weinberg Equilibrium report (.hardy)
plink2 --bfile chr20_final_cleaned --freq --hardy --out chr20_stats
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
plink2 --bfile chr20_final_cleaned --keep eur_ids.txt ---hardy --freq --out eur_stats
plink2 --bfile chr20_final_cleaned --keep afr_ids.txt ---hardy --freq --out afr_stats
```
### Inspecting the difference between getting MAF and HWE P values between populations
We are now better able to compare the MAF and HWE P values at SNPs between populations. To visualize this in R use script `chr20_maf_hwe_compare.R`, and to do this in python, use script `chr20_maf_hwe_compare.py`. 
Key differences to observe: 
* MAF Distribution: You will likely notice that the AFR population has a higher density of rare variants compared to the EUR population, reflecting the greater genetic diversity found in African populations.
* HWE Outliers: If one population has a massive spike at high $-log_{10}(P)$ values that the other doesn't, it may indicate a population-specific technical artifact or a region under intense natural selection in that specific ancestry group.
## Relatedness 
Now let's use the filtered bi-allelic SNP data (MAF ($< 0.05$) and HWE P ($< 10^{-6}$)) from all 1000 Genomes individuals that I've prepared. Let's first calculate the relatedness between all pairs of individuals using KING. Because the 1000 Genomes dataset contains individuals from different populations, the KING-robust algorithm is the best choice because it is specifically designed to handle population structure without needing pre-defined allele frequencies. 
First, to get the kinship coefficient for all pairs of individuals we can use: 
```
# -b: input the PLINK .bed file (use the prefix of your files)
# --kinship: calculate the kinship coefficient (KING-robust)
# --prefix: name your output files
king -b allchr.bed --kinship --prefix allchr_kinship
``` 
Second, if we only want to see pairs that are actually related (e.g., 1st, 2nd, or 3rd degree) and ignore the unrelated pairs, use the --related flag. This is much more efficient for large datasets like the 1000 Genomes.
```
# --related: specifically identifies and classifies relative pairs
king -b allchr.bed --related --prefix allchr_relatives 
```
### Filtering out related individuals 
To filter out relatives (those closer than 3rd degree relatives, who have relatedness = 0.044) and keep a set of unrelated individuals, we will use a two-step process: first, let KING identify which individuals to remove based on your specific threshold, and then use plink2 to create the new, "clean" dataset.
To do the first step, KING has a built-in "unrelated" command that uses a greedy algorithm to find the largest possible subset of unrelated individuals.
```
# --unrelated: identifies a subset of unrelated individuals
# --degree 3: removes up to 3rd-degree relatives (kinship > 0.044)
king -b allchr.bed --unrelated --degree 3 --prefix allchr_unrelated
```
Once KING finishes, it will produce a file ending in .unrelated.id. This file contains the list of people we should keep.
```
# --vcf or --bfile: your input
# --keep: tells PLINK to only retain the IDs in the KING output file
# --make-bed: saves the new unrelated dataset
plink --bfile allchr \
      --keep allchr_unrelated.unrelated.id \
      --make-bed \
      --out allchr_unrelated 
``` 
