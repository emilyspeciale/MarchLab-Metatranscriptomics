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
This will output a multiqc.html file, download it onto your computer, then open it so it shows up in your browser. It will show many important quality control variables for your samples. Make sure to check that everything looks good before continuing on with the rest of the pipeline.
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
rnaspades.py \
 --pe1-1 /proj/marchlab/projects/MetaT_Example/Trimmed_Reads/1-1A_R1_001_val_1.fq.gz \
 --pe1-2 /proj/marchlab/projects/MetaT_Example/Trimmed_Reads/1-1A_R2_001_val_2.fq.gz \
 -o $outdir
```
### CD-HIT

## Annotation
### Shortern Contig Headers
### TransDecoder
### MarFERReT
### eggNOG-mapper
## Alignment
### Salmon 
### Tximport
## Analyzing Results
### Merging mapped read counts with annotations
### Taxonomy of mapped reads
### DESeq2
