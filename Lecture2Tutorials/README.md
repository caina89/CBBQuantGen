# Tutorial for Lecture 2 - Human genetic variation

## Get data 

We will be using chr20 data from three individuals xxx, xxx and xxx from the [1000 Genomes Project Phase 3](https://www.internationalgenome.org/data-portal/data-collection/phase-3) dataset for this exercise. 

To download the data, 


## MarkDuplicates

PCR duplicate reads are identified and marked before variant calling to reduce bias. Mark Duplicates is performed using [Picard tools](https://broadinstitute.github.io/picard/) and detailed description of MarkDuplicates is given [here](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)

```
java -jar picard.jar MarkDuplicates \
      I=input.bam \
      O=marked_duplicates.bam \
      M=marked_dup_metrics.txt
```
