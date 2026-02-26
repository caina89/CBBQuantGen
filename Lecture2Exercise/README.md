# Tutorial for Lecture 2 - Compute set up 
## Euler 
You should all have access to [Euler](https://sis.id.ethz.ch/services/hpc/#the-euler-cluster), the ETH High-Performance-Cluster (HPC) available to all ETH researchers and students. By default everyone gets 45G of space on Euler for study uses. The current system, called “Euler”, was introduced in early 2014 and has been expanded regularly since then. The latest “Euler IX” expansion in 2023 brought the cluster’s overall computing capacity to over 200,000 CPU cores and 1,300 GPUs. To access Euler you need to first [create an account as shown here](https://docs.hpc.ethz.ch/tutorials/getting-started/). You can then log into Euler using SSH. 
```
ssh <USER>@euler.ethz.ch #(e.g. ssh f.muster@euler.ethz.ch)
```
Euler is a linux machine that enables running of softwares on the linux command line. Most bioinformatics and quantitative genetics tools have a command line version, or only exist as a command line executable. For basic linux commands [see](https://docs.hpc.ethz.ch/tutorials/linux-command-line/).  
## Working on Euler
* Login to a login-node via your computer.
* Submit your jobs on the Login nodes. Jobs get sent to the Compute nodes by the scheduler.
* Results will becomes available as defined by job script, and results will be saved in your your working directory (Scratch) as specified by job script.
* Collect the results (select, filter, compress etc.).

## Conda environment 
Usually for computational biology and bioinformatics projects Miniconda (lightweight) is sufficient. This installer is small and only includes conda and Python. You then install only the packages you need.
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
## Installation Steps (Quick Guide)
Use the bash commands above to download Miniconda and install it.
* Review License: Press Enter to view the terms, and then type `yes` to accept.
* Confirm Path: It will ask to install in `~/miniconda3` or `~/anaconda3`. Press `Enter` to accept the default.
* Initialize: At the end, it will ask if you want to `run conda init`. Type `yes`. This adds conda to your shell so it works every time you open the terminal.
* Restart Terminal: Close your terminal window and reopen it (or run source `~/.bashrc`). You should see (`base`) next to your cursor.
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
