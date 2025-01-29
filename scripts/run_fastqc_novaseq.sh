#!/bin/bash

#PBS -l select=1:ncpus=2:mem=32gb
#PBS -l walltime=1:00:00
#PBS -M forrest.koch@unsw.edu.au
#PBS -m ae
#PBS -o logs/fastqc-novaseq/
#PBS -j oe
#PBS -J 1-1614

ncpu=2
module load fastqc

# This should be submitted from the SCIMETAR folder ...
cd $PBS_O_WORKDIR

# Get the job identifier.
# to be run ...
i=$PBS_ARRAY_INDEX

# Grab the specific folder we are after
sample=$(find data/novaseq -type f -name '*.fastq.gz' | tail -n+$i | head -n1)


mkdir -p qc/fastqc-novaseq

fastqc -t $ncpu -o qc/fastqc-novaseq/ $sample
