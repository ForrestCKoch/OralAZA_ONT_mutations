#!/bin/bash

#PBS -l select=1:ncpus=16:mem=64gb
#PBS -l walltime=2:00:00
#PBS -M forrest.koch@unsw.edu.au
#PBS -m ae
#PBS -o logs/mutations-v3/
#PBS -j oe
#PBS -J 1-4


mkdir -p logs/mutations-v3/

NCPU=16
module load python/3.10

# This should be submitted from the SCIMETAR folder ...
cd $PBS_O_WORKDIR

i=$PBS_ARRAY_INDEX

sample=$(ls data/nanoranger_output_v2 | tail -n+$i | head -n1)

sam="data/nanoranger_output_v2/$sample/${sample}_genome_tagged.bam"
mutations="data/mutations-of-interest.csv"
outdir="results/mutations-called-v3/"
mkdir -p $outdir

output="$outdir/mutations-called_${sample}.csv"
python3.10 scripts/get_mutations.py -s $sam -m $mutations -o $output -c $NCPU
