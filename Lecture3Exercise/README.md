# Tutorial for Lecture 3 - Human genetic variation
## Tools 
[samtools](https://www.htslib.org/)
[picardtools](https://broadinstitute.github.io/picard/)
[Genome analysis toolkit (GATK)](https://gatk.broadinstitute.org/hc/en-us)
You can install all three using conda 
```
conda install -c bioconda samtools picard gatk4
conda install -c conda-forge ncurses
## because of some dependencies samtools need 
```
## File formats 
[FASTA](https://en.wikipedia.org/wiki/FASTA_format)
[Sequence alignment formats SAM/BAM](https://samtools.github.io/hts-specs/SAMv1.pdf)
[Variant call format VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf)
## Data 
Human genome reference: We will be using the GRCh38 human genome reference for all our exercises throughout the course. 
Human sequence data: We will be using reads aligning to the human genome reference GRCh38 on chr20 from three individuals HG00096, HG00100 and NA12878 from the [1000 Genomes Project Phase 3](https://www.internationalgenome.org/data-portal/data-collection/phase-3) dataset for this exercise. 
1. Download the Human Genome Reference hg38 from the UCSC Genome Browser - note this version of the hg38 genome reference uses the chr1, chr2 naming convention:
```
cd
mkdir data
cd data 
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gzip -d hg38.fa.gz
```
2. Download the chr20 bams (actually, crams, even more compact version) from the 1000 Genomes Project FTP
```
# 1. Identify a random samples
SAMPLES=$(echo "HG00100" "HG00109" "HG00132")
# 2. Identify the Base URL for the High Coverage (30x) data
BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR"
# 3. Download
for ID in $SAMPLES; do
    echo "------------------------------------------"
    echo "Downloading Chromosome 20 for $ID..."
    # Construct the URL for the CRAM file
    REMOTE_URL="${BASE_URL}/${ID}/alignment/${ID}.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram"
    # Use samtools to download ONLY chr20 and save it as a BAM locally
    # -b: output in BAM format
    # -h: include the header (essential for the file to be valid)
    samtools view -T hg38.fa -b -h "$REMOTE_URL" chr20 > "${ID}_chr20.bam"
    # Create an index for the new BAM file so you can view it in IGV
    samtools index "${ID}_chr20.bam"
    echo "Done! Saved as ${ID}_chr20.bam and indexed."
done
```
## Index and creating dictionary for reference genome 
Indexing the reference FASTA file enables enables efficient access to arbitrary regions within those reference sequences.
```
samtools faidx hg38.fa
``` 
Creating a dictionary .dict file for a reference FASTA enables bioinformatics tools, such picardtools and GATK, to check that essential information about the contigs (chromosomes) present in the reference FASTA file, such as their names, lengths, and MD5 checksums. 
```
picard CreateSequenceDictionary R=hg38.fa O=hg38.dict
```
## Sorting, indexing and viewing your bams 
We use samtools and picardtools for inspecting and processing sequencing files in the SAM/BAM formats. As these files contain a large number of aligned sequences, identifying and accessing sequencing reads aligned to specific regions in the genome (e.g. chr20:500000-600000) is made faster through sorting them by genome coordinates (in the reference genome to which the sequences are aligned) and indexing them for fast random access.  
### Sort coordinates 
```
SAMPLES=$(echo "HG00100" "HG00109" "HG00132")
for ID in $SAMPLES; do
picard SortSam I=${ID}_chr20.bam O=${ID}_chr20.sorted.bam SORT_ORDER=coordinate
samtools index ${ID}_chr20.sorted.bam
done 
```
### Looking at a BAM files 
[samtools view](https://www.htslib.org/doc/samtools-view.html) allows you to access aligned reads in stdout or piped to an output file (in SAM or BAM format) for later use. For options see descriptions in the samtools view page.  
This example command shows you the first 10 lines of a sequence alignment BAM file, which contains the first 10 reads aligned to the chr20:500000-600000 region in the human genome reference file `hg38.fa`.  
```
samtools view -T hg38.fa HG00100_chr20.sorted.bam chr20:500000-600000 | head -n 10
```
## Preprocessing before variant calling 
### MarkDuplicates
PCR duplicate reads are identified and marked before variant calling to reduce bias. Mark Duplicates is performed using [MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) in picardtools. 
```
picard MarkDuplicates I=$wdir/data/HG00100_chr20.sorted.bam O=$wdir/data/HG00100_chr20.markdup.bam M=$wdir/data/HG00100_chr20.markdup.metrics.txt
```
### Base quality score recalibration
Systematic bias can originate from library preparation, sequencing, manufacturing defects in the flowcell chips, sequencer variation, and sequencing chemistry, resulting in over- or underestimation of quality scores. Base quality score recalibration in GATK involves two steps. In step 1, in BaseRecalibrator, an error model is built through comparing the base quality scores at all bases in input file (raw, from sequencers) to those at known variants (previously identified to be true human genetic variations). The error model calibrates the base quality scores such that those at known human variations are more likely to be adjusted higher, and those at novel variations identified in the input sequencing file are likely to be adjusted lower (since they are more likely to be sequencing errors). 

First, to download the known variants and their index files. In this instance, we'd use the hg38 VCF file for all common variants documented in dbSNP
```
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/GATK/common_all_20180418.vcf.gz
```
Then perform step1, the BaseRecalibrator: 
```
SAMPLES=$(echo "HG00100" "HG00109" "HG00132")
for ID in $SAMPLES; do
gatk --java-options "-Xms4G -Xmx4G" BaseRecalibrator \
  -I ${ID}_chr20.markdup.bam \
  -R hg38.fa \
  -O ${ID}_chr20.markdup.bqsr.report \
  --known-sites common_all_20180418.vcf.gz
done
```
This error model is applied in step 2, using ApplyBQSR: 
```
SAMPLES=$(echo "HG00100" "HG00109" "HG00132")
for ID in $SAMPLES; do
gatk --java-options "-Xms4G -Xms4G" ApplyBQSR \
  -I $wdir/data/${ID}_chr20.markdup.bam \
  -R $wdir/hg38/hg38.fa \
  --bqsr-recal-file $wdir/data/${ID}_chr20.markdup.bqsr.report \
  -O $wdir/data/${ID}_chr20.markdup.bqsr.bam
done 
```
## Variant calling 
The [GATK HaplotypeCaller] calls SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. In other words, whenever the program encounters a region showing signs of variation, it discards the existing mapping information and completely reassembles the reads in that region. This allows the HaplotypeCaller to be more accurate when calling regions that are traditionally difficult to call, for example when they contain different types of variants close to each other. For each potentially variant site, the program applies Bayes' rule, using the likelihoods of alleles given the read data to calculate the likelihoods of each genotype per sample given the read data observed for that sample. The most likely genotype is then assigned to the sample.
### Single-sample variant calling 
Single-sample variant calling generates GVCF files, which records variant sites and groups non-variant sites into blocks during the calling process based on genotype quality. This is a way of compressing the VCF file without losing any sites in order to do joint analysis in subsequent steps.
```
SAMPLES=$(echo "HG00100" "HG00109" "HG00132")
for ID in $SAMPLES; do
gatk --java-options "-Xms4G -Xmx4G" HaplotypeCaller \
  -R $wdir/hg38/hg38.fa \
  -I $wdir/data/${ID}_chr20.markdup.bqsr.bam \
  -O $wdir/data/${ID}_chr20.markdup.bqsr.g.vcf.gz \
  -ERC GVCF
done 
```
### Generating a variant database 
GVCFs are consolidated into a GenomicsDB datastore in order to improve scalability and speedup the next step: joint genotyping. This step can be done per chromosome, so that the next step can be run in parallel across all chromosomes and therefore speedup the process. 
```
j=20
gatk --java-options "-Xms4G -Xmx4G" GenomicsDBImport \
  --genomicsdb-workspace-path $wdir/data/chr${j}db \
  -R $wdir/hg38/hg38.fa \
  --sample-name-map $wdir/data/sample_map.txt \
```
For building this database on a number of samples, it is recommended to use a sample name map file rather than keying in each sample using input option -V individually. The tab-delimited sample map file looks like this: 
```
  HG00100      HG00100_chr20.markdup.bqsr.g.vcf.gz
  HG00109      HG00109_chr20.markdup.bqsr.g.vcf.gz
  HG00132      HG00132_chr20.markdup.bqsr.g.vcf.gz
``` 
### Joint genotype calling 
GenotypeGVCFs uses the potential variants from the HaplotypeCaller recorded in the GCVFs of all samples in the cohort and does the joint genotyping. It will look at the available information for each site from both variant and non-variant alleles across all samples, and will produce a VCF file containing only the sites that it found to be variant in at least one sample.
```
j=20
gatk --java-options "-Xms4G -Xmx4G" GenotypeGVCFs \
  -R hg38.fa \
  -V gendb://chr${j}db -O $wdir/data/chr${j}.vcf.gz
``` 
## Variant call quality control 
### Variant annotations 
Raw variant calls can include many artifacts, and this is captured in many metrics (called "annotations") that are output in the joint calling step. These include: QD (Variant confidence normalized by unfiltered depth of variant samples), MQ (Root mean square of the mapping quality of reads across all samples), MQRankSum (Rank sum test for mapping qualities of REF versus ALT reads), ReadPosRankSum (Rank sum test for relative positioning of REF versus ALT alleles within reads), FS (Strand bias estimated using Fisher's exact test) and SOR (Strand bias estimated by the symmetric odds ratio test). A comprehensive list of all available annotations can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/30331989211419--Tool-Documentation-Index#VariantAnnotations) though what to use is dependent on the project and the sequencing data input. 
### Variant quality score recalibration
[Variant quality score recalibration (VQSR) in GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR) works by builing a Gaussian Mixture Model of these annotations for all variants identified at the joint-calling step, and asking if they cluster with those of known human variations (training sets - can supply multiple) previously identified in other projects (note, this can be the same ones as used in BQSR step earlier). 
We would also supply a range of "tranches", which are percentage specificity to known variants in the "positive" set we'd identify as true variants in the analysed dataset i.e. 100.0 would contain only known variants supplied in the training sets, 90.0 would mean the positive set contains 90% known variants. Usually, a lot of tranches in the 90-100 range are given, since we expect a new variant call set to be majority known (and contain very few novel, rare variants). Through doing this, VQSR generates a new quality score called the VQSLOD, a log-ratio of the variant’s probabilities belonging to the positive and negative sets given each of the tranches (selectivities). The purpose of this new score is to enable variant filtering in a way that allows analysts to balance sensitivity (trying to discover all the real variants) and specificity (trying to limit the false positives that creep in when filters get too lenient) as finely as possible. 
Also, notably, this step is usually better when performed on as many variants as possible. As the earlier step is done per chromosome to parallelize the process, usually VCFs from all chromosomes are merged before the VQSR step. 
```
gatk --java-options "-Xms4G -Xmx4G" VariantRecalibrator \
  -tranche 100.0 -tranche 99.95 -tranche 99.9 \
  -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
  -tranche 95.0 -tranche 94.0 \
  -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
  -R hg38.fa \
  -V allchr.vcf.gz \
  --resource:dbsnp,known=true,training=true,truth=true,prior=15.0 common_all_20180418.vcf.gz  \
  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR  \
  -mode SNP -O allchr.recal --tranches-file allchr.tranches \
  --rscript-file allchr.plots.R
```
The output file ```$wdir/data/allchr.plots.R``` is inspected to identify the right "tranche" to select for SNPs, usually through assessing if the SNP Ti/Tv ratios of all variants in the tranche, including known and novel, are close to 2. The VQSLOD scores corresponding to the tranche the user finds ideal, in the ```$wdir/data/allchr.tranches``` file, are then applied in the ApplyVQSR step to determine if a variant has a PASS or FAIL in the QUAL column of the final, recalibrated, VCF file. As the Ti/Tv ratio metric can only be used on SNPs, usually the recalibration tranche is determined using SNPs, and the same tranche (and therefore VQSLOD cutoff) is applied to call INDELs.  
```
gatk --java-options "-Xms4G -Xmx4G" ApplyVQSR \
  -V allchr.vcf.gz \
  --recal-file allchr.recal \
  -mode SNP \
  --tranches-file allchr.tranches \
  --truth-sensitivity-filter-level 99.9 \
  --create-output-variant-index true \
  -O allchr.recalibrated_99.9.vcf.gz
```

