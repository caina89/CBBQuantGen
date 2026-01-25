# Tutorial for Lecture 2 - Human genetic variation

## Get data 
We will be using the hg38 human genome reference for all our exercises throughout the course. See a description of the human genome reference [here]. The [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file for hg38 is already downloaded for you in xxx, but to be able to download it for yourselves see "Genome sequence, primary assembly (GRCh38)" in FASTA section in [Gencode](https://www.gencodegenes.org/human/). 
We will be using reads aligning to the human genome reference hg38 on chr20 from three individuals xxx, xxx and xxx from the [1000 Genomes Project Phase 3](https://www.internationalgenome.org/data-portal/data-collection/phase-3) dataset for this exercise. This is already downloaded for your in xxx, but to be able to download them yourselves see xxx. 
## Index and creasing dictionary for reference genome 
Indexing the reference FASTA file enables enables efficient access to arbitrary regions within those reference sequences.
```
samtools faidx hg38.fa
``` 
Creating a dictionary .dict file for a reference FASTA enables bioinformatics tools, such as [picardtools](https://broadinstitute.github.io/picard/) and [genome analysis toolkit (GATK)](https://gatk.broadinstitute.org/hc/en-us), to check that essential information about the contigs (chromosomes) present in the reference FASTA file, such as their names, lengths, and MD5 checksums. 
```
java -jar picard.jar CreateSequenceDictionary R=hg38.fasta O=hg38.dict
```
## Sorting, indexing and viewing your bams 
We use [samtools](https://www.htslib.org/) and [picardtools](https://broadinstitute.github.io/picard/) for inspecting and processing sequencing files in the [sequence alignment formats SAM/BAM](https://samtools.github.io/hts-specs/SAMv1.pdf). As these files contain a large number of aligned sequences, identifying and accessing sequencing reads aligned to specific regions in the genome (e.g. chr20:500000-600000) is made faster through sorting them by genome coordinates (in the reference genome to which the sequences are aligned) and indexing them for fast random access.  
### Sort coordinates 
```
java -jar picard.jar SortSam I=input.bam O=sorted.bam SORT_ORDER=coordinate
```
### Indexing
```
samtools index input.bam
```
### Looking at a BAM files 
[samtools view](https://www.htslib.org/doc/samtools-view.html) allows you to access aligned reads in stdout or piped to an output file (in SAM or BAM format) for later use. For options see descriptions in the samtools view page.  
This example command shows you the first 10 lines of a sequence alignment BAM file, which contains the first 10 reads aligned to the chr20:500000-600000 region in the human genome reference file xxx.  
```
samtools view input.bam chr20:500000-600000 | head -n 10
```
## MarkDuplicates
PCR duplicate reads are identified and marked before variant calling to reduce bias. Mark Duplicates is performed using [picardtools](https://broadinstitute.github.io/picard/) and detailed description of MarkDuplicates is given [here](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)
```
java -jar picard.jar MarkDuplicates I=input.bam O=marked_duplicates.bam M=marked_dup_metrics.txt
```
