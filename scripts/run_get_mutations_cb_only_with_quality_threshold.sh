#!/bin/bash

#PBS -l select=1:ncpus=4:mem=96gb
#PBS -l walltime=2:00:00
#PBS -M forrest.koch@unsw.edu.au
#PBS -m ae
#PBS -o logs/mutations-cb-only-with-quality-threshold-v2/
#PBS -j oe
#PBS -J 1-24
# PBS -J 9-21:4%1

NCPU=4
module load python/3.10

# This should be submitted from the SCIMETAR folder ...
cd $PBS_O_WORKDIR

mkdir -p logs/mutations-cb-only-with-quality-threshold-v2/

args=$(cat args/pool_thresh_2.txt | tail -n+$PBS_ARRAY_INDEX | head -n1)

i=$(echo $args|cut -d' ' -f2)
threshold=$(echo $args|cut -d' ' -f1)

sample=$(ls data/nanoranger_output_v2 | tail -n+$i | head -n1)

sam="data/nanoranger_output_v2/$sample/${sample}_genome_tagged.bam"
mutations="data/mutations-of-interest.csv"
outdir="results/mutations-called-cb-only-with-quality-threshold_v2/"
mkdir -p $outdir

output="$outdir/mutations-called_${sample}_threshold-$threshold.csv"
python3.10 scripts/get_mutations_cb_only_with_quality_threshold_v2.py -s $sam -m $mutations -o $output -c $NCPU -q $threshold
