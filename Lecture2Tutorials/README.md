# Tutorial for Lecture 2 - Human genetic variation

## Get data 
We will be using  
```bash
bfile= #genotype files, plink format, excluding extension
phenofile= #binary phenotype file, header: FID IID PHENO
covar= #covariate file, header FID IID cov1 cov2 ...
outfile= #full basename of output file to be created
memory_use= #on cluster submission to match with specifications from job submission to avoid out-of-memory job failure

##MarkDuplicates

PCR duplicate reads are identified and marked before variant calling to reduce bias. Mark Duplicates is performed using [Picard tools](https://broadinstitute.github.io/picard/) and detailed description of MarkDuplicates is given [here](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)

```
java -jar picard.jar MarkDuplicates \
      I=input.bam \
      O=marked_duplicates.bam \
      M=marked_dup_metrics.txt
```
