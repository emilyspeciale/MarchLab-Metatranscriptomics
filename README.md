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
Now copy your clustered assembly into this new TransDecoder v5.7.1 folder using the command ```cp /proj/marchlab/projects/MetaT_Example/CDHit/clustered_assembly.fasta /proj/marchlab/projects/MetaT_Example/TransDecoder-TransDecoder-v5.7.1/```. 

Navigate to this folder via ```cd /proj/marchlab/projects/MetaT_Example/TransDecoder-TransDecoder-v5.7.1```, then run ```nano transdecoder_longorfs.sh```. Copy and paste the following script, edit accordingly:

```bash
#!/bin/bash

#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --time=0-48:00:00
#SBATCH --mem=400G
#SBATCH --ntasks=12
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=speciale@unc.edu
#SBATCH -J transdecoder_1
#SBATCH -o transdecoder_1.%j.out
#SBATCH -e transdecoder_1.%j.err


./TransDecoder.LongOrfs -t clustered_assembly.fasta
```
Then run ```sbatch transdecoder_longorfs.sh```. 

Once ```transdecoder_longorfs.sh``` completes, in the same directory, create a new file, ```nano transdecoder_predict.sh```. Copy and paste the following script, edit accordingly:

```bash
#!/bin/bash

#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --time=0-48:00:00
#SBATCH --mem=400G
#SBATCH --ntasks=12
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=speciale@unc.edu
#SBATCH -J transdecoder_predict_ex
#SBATCH -o transdecoder_predict_ex.%j.out
#SBATCH -e transdecoder_predict_ex.%j.err


./TransDecoder.Predict -t clustered_assembly.fasta
```
Then run ```sbatch transdecoder_predict.sh```.

### MarFERReT
MarFERReT is an open-source reference library for marine eukaryote genes that was developed in 2023. MarFERReT includes taxonomic annotations from NCBI/PR2 and functional annotations from Pfam. I like to use MarFERReT as my primary annotator because 1) it has the most updated and expansive list of genomes for marine eukaryotes, 2) Pfam tends to produce the highest amount of annotations compared to other functional databases, and 3) it is very user friendly. However, there are some drawbacks to MarFERRet, being that 1) it does not include many reference genomes for bacteria/viruses and 2) Pfam annotations are typically for broader protein domains rather than specific genes, and Pfam does not allow for direct mapping to pathways like KEGG does. Thus, I like to use MarFERReT in conjunction with EUKulele for taxonomic annotation and eggNOG for functional annotation. See below for more details!

MarFERReT annotation takes a few steps. I have already installed MarFERReT into the data folder. First, navigate to your work directory using ```cd /work/users/s/p/speciale```. Then, write ```nano marferret.sh```. 

```bash
#!/bin/bash

#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --time=0-48:00:00
#SBATCH --mem=200G
#SBATCH --ntasks=12
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=speciale@unc.edu
#SBATCH -J mft_ex1
#SBATCH -o mft_ex1.%j.out
#SBATCH -e mft_ex1.%j.err

module load diamond

indir=/proj/marchlab/projects/MetaT_Example/TransDecoder-TransDecoder-v5.7.1

outdir=/proj/marchlab/projects/MetaT_Example/Marferret

echo "Checking if ${outdir} exists ..."
if [ ! -d ${outdir} ]
then
    echo "Create directory ... ${outdir}"
    mkdir -p ${outdir}
else
    echo " ... exists"
fi


samples='clustered_assembly.fasta.transdecoder.pep'

for s in `echo $samples`; do

diamond blastp -d /proj/marchlab/data/MarFERReT_v1/MarFERReT.v1.1.1.dmnd \
        -q $indir/${s} \
        -o $outdir/${s}marferret.m8 \
        -p 12 -e 0.000001 -k 1

done
```


### eggNOG-mapper

First, run  ```mkdir -p /proj/marchlab/projects/MetaT_Example/Eggnog/Chunks``` to set up the appropriate directories for eggNOG. 

eggNOG can have trouble running if your grand assembly is too large. If this is the case, it is best to split your grand assembly into a few separate files, run each of them through eggNOG, then combine them again. To do this, you can use this python script, which I run by calling ```module load python/3.9.6```, then executing this code in my command line:

```python
import os

def split_pep(input_file, output_dir, contigs_per_chunk):
    output_prefix = "split_part_"
    file_number = 0
    contig_count = 0
    output_file = None

    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    def open_new_file():
        nonlocal output_file, file_number, contig_count
        if output_file:
            output_file.close()
        file_number += 1
        contig_count = 0
        output_filename = os.path.join(output_dir, f"{output_prefix}{file_number:03d}.pep")
        output_file = open(output_filename, 'w')
        print(f"Creating new file: {output_filename}")

    with open(input_file, 'r') as infile:
        for line in infile:
            if line.startswith('>'):
                if contig_count >= contigs_per_chunk:
                    open_new_file()
                contig_count += 1

            if output_file is None:
                open_new_file()

            output_file.write(line)

    # Close the last output file
    if output_file:
        output_file.close()

# Example usage
input_file = '/proj/marchlab/projects/MetaT_Example/clustered_assembly.fasta.transdecoder.pep'
output_dir = '/proj/marchlab/projects/MetaT_Example/Eggnog/Chunks'

contigs_per_chunk = 1000000
split_pep(input_file, output_dir, contigs_per_chunk)


```
## Alignment
Alignment is the quantification of reads to their associated contigs. This step is not dependent on annotation, as it only requires the trimmed reads and the clustered assembly. We use salmon to conduct our alignment, and the first step of this process is to create an assembly index. 

### Salmon 

First, write ```nano align_part1.sh```, then copy and paste the following code. Edit accordingly (you do not need to make the directory for Alignment ahead of time, the code does it for you):

```bash
#!/bin/bash
#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --time=5-00:00:00
#SBATCH --mem=400G
#SBATCH --ntasks=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=speciale@unc.edu
#SBATCH -J align_pt1
#SBATCH -o align_pt1.%A.out
#SBATCH -e align_pt1.%A.err

module load salmon

salmon index -i /proj/marchlab/projects/MetaT_Example/Alignment/AssemblyIndex \
--transcripts /proj/marchlab/projects/MetaT_Example/CDHit/clustered_assembly.fasta -k 31

```

Run ```sbatch align_part1.sh```. Wait for this code to finish running before proceeding to part 2 (it might take a while). 

Once part 1 finishes, write ```nano align_part2.sh```, then copy and paste the following code. Edit accordingly, specifically in the samples variable, make sure to include every sample you want Salmon to quantify.

```bash
#!/bin/bash

indir=/proj/marchlab/projects/MetaT_Example/Trimmed_Reads
outdir=/proj/marchlab/projects/MetaT_Example/Alignment/Salmon_Quant

mkdir -p /proj/marchlab/projects/MetaT_Example/Alignment/Salmon_Quant

samples='1-1A 1-1B 2-1A 2-1B 2-1C 3-1A 3-1B 3-1C'

for s in $samples; do
    echo ${s}
    R1=`ls -l $indir | grep -o ${s}_R1_001_val_1.fq.gz`
    R2=`ls -l $indir | grep -o ${s}_R2_001_val_2.fq.gz`
    echo ${R1}
    echo ${R2}
    jobfile="salmonquant${s}.sh"
    echo $jobfile
    cat <<EOF > $jobfile
#!/bin/bash
#SBATCH -N 1
#SBATCH -t 05-00:00:00
#SBATCH --mem=250g
#SBATCH -n 16
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=speciale@unc.edu
#SBATCH -J alignpart2
#SBATCH -o alignpart2.%A.out
#SBATCH -e alignpart2.%A.err


module add salmon
echo 'BEGIN'
date
hostname
salmon quant -l A -i  /proj/marchlab/projects/MetaT_Example/Alignment/AssemblyIndex \
        -1 $indir/${R1} \\
        -2 $indir/${R2} \\
        -p 5 --validateMappings \\
        -o /proj/marchlab/projects/MetaT_Example/Alignment/Salmon_Quant/${s}_quant

echo 'END'

date

EOF

    sbatch $jobfile

done
```

Run ```sbatch align_part2.sh```.

### Tximport
Tximport is an R package that can be used to merge all the results from each sample's alignment into an R object. First, write ```nano tximport.r```, then copy and paste the following code. Edit accordingly.
```
BiocManager::install('tximport')
library(tximport)
library(DESeq2)
library(tidyverse)
library(dplyr)
library(stringr)


samples<-list.files(path="/proj/marchlab/projects/MetaT_Example/Alignment/Salmon_Quant", full.names=T)

files<-file.path(samples,"quant.sf")

names(files)<-str_replace(samples, "/proj/marchlab/projects/MetaT_Example/Alignment/Salmon_Quant", "")%>%str_replace(".salmon","")

txi_obj <- tximport(files, type = "salmon", txOut = TRUE)

txi_obj.rds <- saveRDS(txi_obj, file = "/proj/marchlab/projects/MetaT_Example/Alignment/txi_obj.rds")
```

Now create a job file, ```nano tximport.sh```, to run the tximport R code. Edit accordingly.

```bash
#!/bin/bash
#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --time=0-2:00:00
#SBATCH --mem=300G
#SBATCH --ntasks=3
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=speciale@unc.edu
#SBATCH -J txi
#SBATCH -o txi.%A.out
#SBATCH -e txi.%A.err

module load r/4.2.2
cd /work/users/s/p/speciale
outdir=/proj/marchlab/projects/MetaT_Example/Alignment

echo "Checking if ${outdir} exists ..."
if [ ! -d ${outdir} ]
then
    echo "Create directory ... ${outdir}"
    mkdir -p ${outdir}
else
    echo " ... exists"
fi


Rscript tximport.r

```

Run ```sbatch tximport.sh```, and the final product will be an R object that contains your counts matrix. 

## Analyzing Results


