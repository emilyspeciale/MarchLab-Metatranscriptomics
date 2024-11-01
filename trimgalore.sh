#!/bin/bash 
#SBATCH --nodes=1 
#SBATCH --time=00-12:00:00 
#SBATCH --mem=100G 
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
