# Tutorial for Lecture 2 - Human genetic variation

## Get data 

We will be using chr20 data from three individuals xxx, xxx and xxx from the [1000 Genomes Project Phase 3](https://www.internationalgenome.org/data-portal/data-collection/phase-3) dataset for this exercise. 

To download the data, 

## Sorting, indexing and viewing your bams 

We use samtools and picardtools for inspecting and processing sequencing files in the [sequence alignment formats SAM/BAM](https://samtools.github.io/hts-specs/SAMv1.pdf). As these files contain a large number of aligned sequences, identifying and accessing sequencing reads aligned to specific regions in the genome (e.g. chr20:500000-600000) is made faster through sorting them by genome coordinates (in the reference genome to which the sequences are aligned) and indexing them for fast random access.  

### Sort coordinates 
```
java -jar picard.jar SortSam \
      I=input.bam \
      O=sorted.bam \
      SORT_ORDER=coordinate
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
java -jar picard.jar MarkDuplicates \
      I=input.bam \
      O=marked_duplicates.bam \
      M=marked_dup_metrics.txt
```
