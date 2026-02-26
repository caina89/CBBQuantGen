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
## File formats 
