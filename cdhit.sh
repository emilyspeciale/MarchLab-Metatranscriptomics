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

