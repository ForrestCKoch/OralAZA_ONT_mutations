#!/bin/bash

#PBS -l select=1:ncpus=8:mem=60gb
#PBS -l walltime=2:00:00
#PBS -M forrest.koch@unsw.edu.au
#PBS -m ae
#PBS -o logs/nanoranger-jobs
#PBS -j oe
#PBS -J 1-50%1

NCPU=8
module load python/3.11
module load minimap2/2.26
module load star/2.7.9a
module load samtools/1.15.1
module load seqkit/2.5.1

#GENOME_DIR=$HOME/vf-scratch/forrest/HumanGenome/genome/STAR_sparse
GENOME_REF=$HOME/vf-scratch/forrest/HumanGenome/genome/GRCh38.primary_assembly.genome.fa.gz

# This should be submitted from the SCIMETAR folder ...
cd $PBS_O_WORKDIR

# Get the job identifier.
# to be run ...
i=$PBS_ARRAY_INDEX
#i=42
#i=43
#i=47

FASTQ=$(find data/base_called -type f -name '*.fastq.gz' | tail -n+$i | head -n1)
SAMPLE_NAME="$(echo $FASTQ | rev | cut -d'/' -f1 | rev | cut -d'.' -f1)"
OUT_DIR="data/nanoranger_output/$SAMPLE_NAME"
TRANS_REF="data/transcript_reference.fa.gz"
BARCODES="data/barcodes.txt.gz"

mkdir -p $OUT_DIR
mkdir -p logs/nanoranger-pipeline
mkdir -p logs/nanoranger-jobs


python3.11 $HOME/local/builds/nanoranger/pipeline.py \
    --cores $NCPU \
    --infile $FASTQ \
    --outdir $OUT_DIR \
    --expname $SAMPLE_NAME \
    --mode 5p10XGEX \
    --trns_ref $TRANS_REF \
    --genome_ref $GENOME_REF \
    --barcodes $BARCODES >> logs/nanoranger-pipeline/${SAMPLE_NAME}.log 2>&1

