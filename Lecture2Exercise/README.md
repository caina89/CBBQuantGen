# Exercise for Lecture 2 - Compute set up 
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
## Organization of directories (usually)
* After logging into Euler you'd land in the `/cluster/home/$username/` directory. We will be using this directory for the whole course. 
* Usually, if you're starting work on a new HPC, you will be given a place in `/home` that usually has very little storage space but is well backed-up. `/home` is usually abbreviated `~`. Usually it is not a great idea to put a lot of data in this space. Rather, this is where you set up your conda, and directories for your code `~/src` and the softwares you will use `~/bin`.
* Usually we'd store the really important and "hot" (meaning it's accessed a lot) data in a `/data` partition that has limited storage but is well backed-up. This is where all the raw data of important, on-going projects is kept. We usually do not perform analysis directly on data in this directory as REALLY BAD THINGS HAPPEN if you overwrite it by accident. We therefore usually create a copy the part of the data we need to use, store it somewhere else, so that the raw data remains untouched, safe, and ready to be re-used by ourselves or other members of the team.
* Usually we'd therefore create a directory in the `/scratch` partition of the computer we are in, which usually has a lot more storage but is not necessarily constantly backed up, or may not be backed up at all. This is space where we'd prepare for "things to be lost any time", and it is where we keep all the copied raw data files we are going to use for our projects, as well as the intermediate analyses files that may be large, numerous, but ultimately not so important.
* When we are done with analyses and have *results*, we'd now transfer them to either `/home` or `/data` depending on their size and intended accessibility to others (e.g. some project results need to be accessible to other project members, so doesn't make sense to keep in your own `/home`), so that they will not be lost "any time".
* This is an [old article](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424) on organization of bioinformatics projects that I think still applies today. 
## Conda environment 
Usually for computational biology and bioinformatics projects Miniconda (lightweight) is sufficient. This installer is small and only includes conda and Python. You then install only the packages you need. You can install Miniconda in your `/cluster/home/$username/` directory. 
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
### Installation of Miniconda (Quick Guide)
Use the bash commands above to download Miniconda and install it.
* Review License: Press Enter to view the terms, and then type `yes` to accept.
* Confirm Path: It will ask to install in `~/miniconda3` or `~/anaconda3`. Press `Enter` to accept the default.
* Initialize: At the end, it will ask if you want to `run conda init`. Type `yes`. This adds conda to your shell so it works every time you open the terminal.
* Restart Terminal: Close your terminal window and reopen it (or run source `~/.bashrc`). You should see (`base`) next to your cursor.
## R for plotting 
We will be using R (or python if you so choose - but in case you use R), and need to install R using conda:
```
# Install R into your current conda environment
conda install -c conda-forge r-base r-essentials
# Check if it works
R --version
```
## Bioinformatics/Quantitative genetics tools starter kit 
We can create a directory `~/bin` for all the softwares you will be using in this course. 
Bioinformatics softwares can be downloaded from their individual webpages or, in many instances, be installed directly by conda. Sometimes it is worth it checking their individual webpages for latest versions, which may or may not be updated in conda. See download pages of the following sites:  
* [samtools](https://www.htslib.org/)
* [picardtools](https://broadinstitute.github.io/picard/)
* [Genome analysis toolkit (GATK)](https://gatk.broadinstitute.org/hc/en-us)
You can install all three using conda at the same time: 
```
conda install -c bioconda samtools picard gatk4 plink2 bcftools
```
## File formats 
* [FASTA](https://en.wikipedia.org/wiki/FASTA_format)
* [Sequence alignment formats SAM/BAM](https://samtools.github.io/hts-specs/SAMv1.pdf)
* [Variant call format VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf)
