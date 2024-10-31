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
## Assembly
### rnaSPAdes
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
