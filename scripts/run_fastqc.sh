#!/bin/bash

#PBS -l select=1:ncpus=4:mem=64gb
#PBS -l walltime=4:00:00
#PBS -M forrest.koch@unsw.edu.au
#PBS -m ae
#PBS -o logs/fastqc/
#PBS -j oe
#PBS -J 1-50

ncpu=8
module load fastqc

#GENOME_DIR=$HOME/vf-scratch/forrest/HumanGenome/genome/STAR_sparse
#GENOME=$HOME/vf-scratch/forrest/HumanGenome/genome/GRCh38.primary_assembly.genome.fa.gz

# This should be submitted from the SCIMETAR folder ...
cd $PBS_O_WORKDIR

# Get the job identifier.
# to be run ...
i=$PBS_ARRAY_INDEX

# Grab the specific folder we are after
sample=$(find data/ -type f -name '*.fastq.gz' | tail -n+$i | head -n1)


mkdir -p qc/fastqc
mkdir -p qc/fastqc-nano

fastqc -t $ncpu -o qc/fastqc/ $sample
fastqc -t $ncpu -o qc/fastqc-nano/ --nano  $sample
