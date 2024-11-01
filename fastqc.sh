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
