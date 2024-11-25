# MarchLab-Metatranscriptomics
An updated metatranscriptomics pipeline for the Marchetti Lab at UNC Chapel Hill.
### Getting Started
I will be demonstrating the metatranscriptomics pipeline with sample data I have provided in this folder: /proj/marchlab/projects/MetaT_Example
These are samples for three stations off the California Current System that represent the progression of an upwelling plume from fresh to aged. If trying out this pipeline for the first time, I encourage you to copy the raw reads into your own directory and use it as practice! But if you are working with your own samples, please make sure to adjust directories/paths/emails accordingly.
## Trimming & FastQC
### Trim Galore
Navigate to your work directory. 
```bash
cd /work/users/s/p/speciale
```
Within your work directory, run the command ```nano trimgalore.sh ```. Copy and paste the following code, edit accordingly:
```bash
#!/bin/bash 
# SBATCH --nodes=1 
# SBATCH --time=00-12:00:00 
# SBATCH --mem=100G 
#SBATCH --cpus-per-task=4 
#SBATCH --mail-type=BEGIN,END,FAIL 
#SBATCH --mail-user=speciale@unc.edu 
#SBATCH -J trimgalore
#SBATCH -o trimgalore.%A.%a.out 
#SBATCH -e trimgalore.%A.%a.err

# load necessary modules for trim_galore, auto loads python, cutadapt
module load trim_galore 
module load pigz 

echo 'BEGIN' 
date 
hostname 
trim_galore -j 4 --stringency 1 --illumina --paired $1 $2 -o $3 
echo 'END'
date
```
Run ```nano trimgalore-all.sh ```. Copy and paste the following code, edit accordingly:
```bash
#!/bin/bash

# Set in directory to where raw reads are
indir=/proj/marchlab/projects/MetaT_Example/Raw_Reads
outdir=/proj/marchlab/projects/MetaT_Example/Trimmed_Reads
# Create the out directory if it does not already exist
echo "Checking if ${outdir} exists ..."
if [ ! -d ${outdir} ]
then
	echo "Create directory ... ${outdir}"
	mkdir -p ${outdir}
else
	echo " ... exists"
fi
# Cut path file at R1 for input list
input=`ls ${indir}/*R1*`
# Run trimgalore.sh for each sample, it will create a job for each one
for file in $input
do
	sbatch trimgalore.sh $file ${file::-14}2_001.fastq.gz ${outdir}
done
```
Run ```sbatch trimgalore-all.sh``` to submit the jobs. 
### FastQC
Navigate to your work directory. 
```bash
cd /work/users/s/p/speciale
```
Within your work directory, run the command ```nano fastqc.sh ```. Copy and paste the following code, edit accordingly:
```bash
#!/bin/bash
#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --time=02-0:00:00
#SBATCH --mem=60G
#SBATCH --ntasks=42
#SBATCH -J fastqc
#SBATCH -o fastqc.%A.out
#SBATCH -e fastqc.%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=speciale@unc.edu

module load fastqc
# Set directories for your raw reads and trimmed reads
raw_reads=/proj/marchlab/projects/MetaT_Example/Raw_Reads
trim_reads=/proj/marchlab/projects/MetaT_Example/Trimmed_Reads
# Set an output directory
outdir=/proj/marchlab/projects/MetaT_Example/FastQC
# Create output directory if it does not already exist
echo "Checking if ${outdir} exists ..."
if [ ! -d ${outdir} ]
then
    echo "Create directory ... ${outdir}"
    mkdir -p ${outdir}
else
    echo " ... exists"
fi
# Run fastqc on raw reads and trimmed reads
fastqc -t 42 $raw_reads/* -o ${outdir}
fastqc -t 42 $trim_reads/*.gz -o ${outdir}
# Run multiqc to get overall report of quality control results
cd ${outdir}
module load multiqc
multiqc . 
```
Run ```sbatch fastqc.sh``` to submit the job. This will output a multiqc.html file, download it onto your computer, then open it so it shows up in your browser. It will show many important quality control variables for your samples. Make sure to check that everything looks good before continuing on with the rest of the pipeline.
## Assembly
Assembly refers to the process of reconstructing complete sequences (also known as contigs) from the reads provided by the high-throughput sequencing company. Because we are studying a microbial community with many organisms, we must conduct a de novo assembly, which is the assembly of sequences from scratch without a reference genome. To do this, we first use the tool rna-SPAdes to create a de novo assembly for each sample. We will then combine all these assemblies together into one large clustered assembly (also called a coassembly, mega-assembly, etc.). We use the tool CD-HIT to do this, as it also contains software to look for contigs that are redundant between individual assemblies and cluster them together, thus simplifying the analysis.
### rnaSPAdes
Navigate to your work directory. 
```bash
cd /work/users/s/p/speciale
```
Within your work directory, run the command ```nano rna_spades.sh ```. Copy and paste the following code, edit accordingly. At this point, I have not made a code where you can submit a job to conduct all the individual sample assemblies at once. So for now, you must run this code to make an assembly for each individual sample you have:
```bash
#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 02-00:00:00
#SBATCH --mem=250g
#SBATCH --ntasks=16
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=speciale@unc.edu
#SBATCH -J rna_spades
#SBATCH -o rna_spades.%A.out
#SBATCH -e rna_spades.%A.err

module load spades
# Set output directory
outdir=/proj/marchlab/projects/MetaT_Example/Spades
# Create output directory if it does not already exist
echo "Checking if ${outdir} exists ..."
if [ ! -d ${outdir} ]
then
    echo "Create directory ... ${outdir}"
    mkdir -p ${outdir}
else
    echo " ... exists"
fi
# Run spades for each sample, you want to use the trimmed reads now. Where it says "1-1A" is where you want to input with your own sample name (such as T0-Ctrl, Q18, etc.)
rnades.py \
 --pe1-1 /proj/marchlab/projects/MetaT_Example/Trimmed_Reads/1-1A_R1_001_val_1.fq.gz \
 --pe1-2 /proj/marchlab/projects/MetaT_Example/Trimmed_Reads/1-1A_R2_001_val_2.fq.gz \
 -o $outdir
```

Run ```sbatch rna_spades.sh``` to submit each job. 

### CD-HIT

In order to run CD-HIT, we need to remame and put all of the individual fasta files produced by rnaSPAdes into one directory. I run this code directly in my command line since it doesn't take long.

```bash
# Navigate to your spades directory
cd /proj/marchlab/projects/MetaT_Example/Spades/

# Rename each transcripts.fasta file to something unique, such as the subdirectory name
for dir in */; do
    # Check if the file transcripts.fasta exists in the subdirectory
    if [ -f "${dir}transcripts.fasta" ]; then
        # Extract the subdirectory name (remove trailing slash)
        subdir_name=$(basename "$dir")
        # Rename the file
        cp "${dir}transcripts.fasta" "${dir}${subdir_name}_transcripts.fasta"
    fi
done

# Create a directory for Transcripts directory if it doesn't exist, this will hold all of our rnaSPAdes transcripts file
mkdir -p Transcripts

# Copy each rnaSPAdes transcript file into the Transcripts directory
for dir in */; do
    # Check if the file with the new name exists in the subdirectory
    if [ -f "${dir}${dir%/}_transcripts.fasta" ]; then
        # Copy the file to the Transcripts directory
        cp "${dir}${dir%/}_transcripts.fasta" "Transcripts/"
    fi
done
```

Navigate back to your work directory. 
```bash
cd /work/users/s/p/speciale
```
Within your work directory, run the command ```nano cdhit.sh ```. Copy and paste the following code, edit accordingly. 

```bash
#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 3-00:00:00
#SBATCH -J cdhit_1
#SBATCH -o cdhit_1.out
#SBATCH -e cdhit_1.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=speciale@unc.edu
#SBATCH --mem=40G
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1

indir=/proj/marchlab/projects/MetaT_Example/Spades/Transcripts
outdir=/proj/marchlab/projects/MetaT_Example/CDHit
# ifn is the combined (concatenated) assemblies
input1=`ls ${indir}/*.fasta`
ifn1="${outdir}/mega_assembly.fasta"
ofn1="${outdir}/clustered_assembly.fasta"

echo "Checking if ${outdir} exists ..."
if [ ! -d ${outdir} ]
then
    echo "Create directory ... ${outdir}"
    mkdir -p ${outdir}
else
    echo " ... exists"
fi

echo "Checking if ${ifn1} exists ..."
if [ ! -f "${ifn1}" ]
then
    echo "Create combined assembly ... ${ifn1}"
    cat ${input1} > ${ifn1}
else
    echo " ... exists"
fi


# load default cdit
module load cdhit

echo "${ifn1}"
echo "${ofn1}"

cd-hit-est \
 -i "${ifn1}" \
 -o "${ofn1}" \
 -c .98 -n 10 -d 100 \
 -T ${SLURM_CPUS_PER_TASK} \
 -M 40000 \

# --------------------- 
# sacct -j $SLURM_JOB_ID --format='JobID,user,elapsed, cputime, totalCPU,MaxRSS,MaxVMSize,ncpus,NTasks,ExitCode'

scontrol show job $SLURM_JOB_ID
```
Run ```cdhit.sh``` to submit the jobs. 

## Annotation
### TransDecoder
In order to annotate, we will use TransDecoder to identify the best candidate open reading frame (ORF) for each contig.
First, you need to install the latest version of TransDecoder into the directory of your choice. This can be done in your command line by grabbing the latest version of TransDecoder from GitHub and extracting it:
```
# Navigating to my project directory
cd /proj/marchlab/projects/MetaT_Example
# Downloading the zip file for TransDecoder 5.7.1
wget https://github.com/TransDecoder/TransDecoder/archive/refs/tags/TransDecoder-v5.7.1.tar.gz
--2024-11-07 09:31:46--  https://github.com/TransDecoder/TransDecoder/archive/refs/tags/TransDecoder-v5.7.1.tar.gz
# Extract that zip file, which will create a new folder in MetaT_Example called TransDecoder-TransDecoder-v5.7.1
tar -xzf /proj/marchlab/projects/MetaT_Example/TransDecoder-v5.7.1.tar.gz
```
### MarFERReT
### eggNOG-mapper
## Alignment
### Salmon 
### Tximport
## Analyzing Results
### Merging mapped read counts with annotations
### Taxonomy of mapped reads
### DESeq2
