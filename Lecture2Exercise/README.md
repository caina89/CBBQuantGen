# Tutorial for Lecture 2 - Compute set up 
## Conda environment 
Usually for computational biology and bioinformatics projects Miniconda (lightweight) is sufficient. This installer is small and only includes conda and Python. You then install only the packages you need.
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
## Installation Steps (Quick Guide)
Use the bash commands above to download Miniconda and install it.
Review License: Press Enter to view the terms, and then type `yes` to accept.
Confirm Path: It will ask to install in `~/miniconda3` or `~/anaconda3`. Press `Enter` to accept the default.
Initialize: At the end, it will ask if you want to `run conda init`. Type `yes`. This adds conda to your shell so it works every time you open the terminal.
Restart Terminal: Close your terminal window and reopen it (or run source `~/.bashrc`). You should see (`base`) next to your cursor.

## Tools 
[samtools](https://www.htslib.org/)
[picardtools](https://broadinstitute.github.io/picard/)
[Genome analysis toolkit (GATK)](https://gatk.broadinstitute.org/hc/en-us)
You can install all three using conda 
```
conda install -c bioconda samtools picard gatk4
```
## File formats 
[FASTA](https://en.wikipedia.org/wiki/FASTA_format)
[Sequence alignment formats SAM/BAM](https://samtools.github.io/hts-specs/SAMv1.pdf)
[Variant call format VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf)
