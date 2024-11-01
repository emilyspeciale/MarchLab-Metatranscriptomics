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
