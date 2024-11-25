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
# Set output directory, make a different out directory for each sample
outdir=/proj/marchlab/projects/MetaT_Example/Spades/3-1C_Spades
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
 --pe1-1 /proj/marchlab/projects/MetaT_Example/Trimmed_Reads/3-1C_R1_001_val_1.fq.gz \
 --pe1-2 /proj/marchlab/projects/MetaT_Example/Trimmed_Reads/3-1C_R2_001_val_2.fq.gz \
 -o $outdir
