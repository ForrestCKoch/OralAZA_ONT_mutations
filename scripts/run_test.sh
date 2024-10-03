#!/bin/bash

#PBS -l select=1:ncpus=12:mem=48gb
#PBS -l walltime=0:30:00
#PBS -M forrest.koch@unsw.edu.au
#PBS -m ae
#PBS -o logs/mutations
#PBS -j oe

NCPU=12
module load python/3.10

# This should be submitted from the SCIMETAR folder ...
cd $PBS_O_WORKDIR

sam="data/nanoranger_output/PGXXXF240300_pass_barcode07/PGXXXF240300_pass_barcode07_genome_tagged.bam"
mutations="data/mutations-of-interest-reduced.csv"
output="test.csv"
python3.10 scripts/get_mutations.py -s $sam -m $mutations -o $output -c $NCPU
