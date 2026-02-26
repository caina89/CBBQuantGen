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
We then use plink2 to calculate the MAF and P-value for violation of Hardy-Weinberg Equilibrium (HWE) at each SNP, and filter out those with MAF < 0.05 and HWE P < 10-6. 
* Note: we are performing the MAF and HWE P-value calculations and filtering on all individuals across many populations in 1000 Genomes Project here. Is this the right thing to do? How differently would you do it this if you were interested in a particular population? 
```
# --vcf: Input your filtered VCF
# --maf 0.05: Filter for Minor Allele Frequency > 5%
# --hwe 1e-6: Filter for Hardy-Weinberg p-value > 10^-6
# --make-bed: Create the standard .bed/.bim/.fam binary fileset
# --out: The prefix for your final files
plink --vcf chr20_snps_only.vcf.gz \
      --maf 0.05 \
      --hwe 1e-6 \
      --make-bed \
      --out chr20_final_cleaned
```
