# Tutorial for Lecture 2 - Human genetic variation

## Data 
Human genome reference: We will be using the GRCh38 human genome reference for all our exercises throughout the course. This is already downloaded for you in xxx, but to be able to download it for yourselves see "Genome sequence, primary assembly (GRCh38)" in FASTA section in [Gencode](https://www.gencodegenes.org/human/). 
Human sequence data: We will be using reads aligning to the human genome reference GRCh38 on chr20 from three individuals xxx, xxx and xxx from the [1000 Genomes Project Phase 3](https://www.internationalgenome.org/data-portal/data-collection/phase-3) dataset for this exercise. This is already downloaded for your in xxx, but to be able to download them yourselves see xxx. 
Human reference variations: 

## Tools 
[samtools](https://www.htslib.org/)
[picardtools](https://broadinstitute.github.io/picard/)
[Genome analysis toolkit (GATK)](https://gatk.broadinstitute.org/hc/en-us)
## File formats 
[FASTA](https://en.wikipedia.org/wiki/FASTA_format)
[Sequence alignment formats SAM/BAM](https://samtools.github.io/hts-specs/SAMv1.pdf)
[Variant call format VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf)

## Index and creating dictionary for reference genome 
Indexing the reference FASTA file enables enables efficient access to arbitrary regions within those reference sequences.
```
samtools faidx hg38.fa
``` 
Creating a dictionary .dict file for a reference FASTA enables bioinformatics tools, such picardtools and GATK, to check that essential information about the contigs (chromosomes) present in the reference FASTA file, such as their names, lengths, and MD5 checksums. 
```
java -jar picard.jar CreateSequenceDictionary R=hg38.fasta O=hg38.dict
```
## Sorting, indexing and viewing your bams 
We use samtools and picardtools for inspecting and processing sequencing files in the SAM/BAM formats. As these files contain a large number of aligned sequences, identifying and accessing sequencing reads aligned to specific regions in the genome (e.g. chr20:500000-600000) is made faster through sorting them by genome coordinates (in the reference genome to which the sequences are aligned) and indexing them for fast random access.  
### Sort coordinates 
```
java -jar picard.jar SortSam I=NA12878.bam O=NA12878.sorted.bam SORT_ORDER=coordinate
```
### Indexing
```
samtools index $wdir/data/input.bam
```
### Looking at a BAM files 
[samtools view](https://www.htslib.org/doc/samtools-view.html) allows you to access aligned reads in stdout or piped to an output file (in SAM or BAM format) for later use. For options see descriptions in the samtools view page.  
This example command shows you the first 10 lines of a sequence alignment BAM file, which contains the first 10 reads aligned to the chr20:500000-600000 region in the human genome reference file xxx.  
```
samtools view $wdir/data/NA12878.sorted.bam chr20:500000-600000 | head -n 10
```
## Preprocessing before variant calling 
### MarkDuplicates
PCR duplicate reads are identified and marked before variant calling to reduce bias. Mark Duplicates is performed using [MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) in picardtools. 
```
java -jar picard.jar MarkDuplicates I=$wdir/data/NA12878.sorted.bam O=$wdir/data/NA12878.markdup.bam M=$wdir/data/NA12878.markdup.metrics.txt
```
### Base quality score recalibration
Systematic bias can originate from library preparation, sequencing, manufacturing defects in the flowcell chips, sequencer variation, and sequencing chemistry, resulting in over- or underestimation of quality scores. Base quality score recalibration in GATK involves two steps. In step 1, in BaseRecalibrator, an error model is built through comparing the base quality scores at all bases in input file (raw, from sequencers) to those at known variants (previously identified to be true human genetic variations). The error model calibrates the base quality scores such that those at known human variations are more likely to be adjusted higher, and those at novel variations identified in the input sequencing file are likely to be adjusted lower (since they are more likely to be sequencing errors). 
```
gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=2" BaseRecalibrator \
  -I $wdir/data/NA12878.markdup.bam \
  -R $wdir/hg38/hg38.fa \
  -O $wdir/data/NA12878.markdup.bqsr.report \
  --known-sites $wdir/hg38/dbsnp_146.hg38.vcf.gz \
  --known-sites $wdir/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz \
  --known-sites $wdir/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
```
This error model is applied in step 2, using ApplyBQSR. 
```
gatk --java-options "-Xms2G -Xmx2G -XX:ParallelGCThreads=2" ApplyBQSR \
  -I $wdir/data/NA12878.markdup.bam \
  -R $wdir/hg38/hg38.fa \
  --bqsr-recal-file $wdir/data/NA12878.markdup.bqsr.report \
  -O $wdir/data/NA12878.markdup.bqsr.bam
```
## Variant calling 
The [GATK HaplotypeCaller] calls SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. In other words, whenever the program encounters a region showing signs of variation, it discards the existing mapping information and completely reassembles the reads in that region. This allows the HaplotypeCaller to be more accurate when calling regions that are traditionally difficult to call, for example when they contain different types of variants close to each other. For each potentially variant site, the program applies Bayes' rule, using the likelihoods of alleles given the read data to calculate the likelihoods of each genotype per sample given the read data observed for that sample. The most likely genotype is then assigned to the sample.
### Single-sample variant calling 
Single-sample variant calling generates GVCF files, which records variant sites and groups non-variant sites into blocks during the calling process based on genotype quality. This is a way of compressing the VCF file without losing any sites in order to do joint analysis in subsequent steps.
```
gatk --java-options "-Xms20G -Xmx20G -XX:ParallelGCThreads=2" HaplotypeCaller \
  -R $wdir/hg38/hg38.fa \
  -I $wdir/data/NA12878.markdup.bqsr.bam \
  -O $wdir/data/NA12878.markdup.bqsr.g.vcf.gz \
  -ERC GVCF
```
### Generating a variant database 
GVCFs are consolidated into a GenomicsDB datastore in order to improve scalability and speedup the next step: joint genotyping. This step can be done per chromosome, so that the next step can be run in parallel across all chromosomes and therefore speedup the process. 
```
j=20
gatk --java-options "-Xms2G -Xmx2G -XX:ParallelGCThreads=2" GenomicsDBImport \
  --genomicsdb-workspace-path $wdir/data/chr${j}db \
  -R $wdir/hg38/hg38.fa \
  --sample-name-map $wdir/data/sample_map.txt \
```
For building this database on a number of samples, it is recommended to use a sample name map file rather than keying in each sample using input option -V individually. The sample map file looks like this: 
```
  NA12878      NA12878.markdup.bqsr.g.vcf.gz
  NA12892      NA12892.markdup.bqsr.g.vcf.gz
  NA12891      NA12891.markdup.bqsr.g.vcf.gz
``` 
### Joint genotype calling 
GenotypeGVCFs uses the potential variants from the HaplotypeCaller recorded in the GCVFs of all samples in the cohort and does the joint genotyping. It will look at the available information for each site from both variant and non-variant alleles across all samples, and will produce a VCF file containing only the sites that it found to be variant in at least one sample.
```
j=20
gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms2G -Xmx2G -XX:ParallelGCThreads=2" GenotypeGVCFs \
  -R $wdir/hg38/hg38.fa \
  -V gendb://$wdir/data/chr${j}db -O $wdir/data/chr${j}.vcf.gz
``` 
