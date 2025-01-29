#!/bin/bash

#PBS -l select=1:ncpus=1:mem=20gb
#PBS -l walltime=0:30:00
#PBS -M forrest.koch@unsw.edu.au
#PBS -m ae
#PBS -o logs/split-patient
#PBS -j oe
#PBS -J 1-4

NCPU=1
module load python/3.11
module load samtools/1.15.1


# This should be submitted from the SCIMETAR folder ...
cd $PBS_O_WORKDIR

# Get the job identifier.
# to be run ...
i=$PBS_ARRAY_INDEX

FASTQ="data/base_called/PGXXXF240300/PGXXXF240300_pass_barcode0$(($i + 4)).fastq.gz"

SAMPLE_NAME="$(echo $FASTQ | rev | cut -d'/' -f1 | rev | cut -d'.' -f1)"

OUT_DIR="data/nanoranger_output_v2/$SAMPLE_NAME/split-by-patient"

mkdir -p $OUT_DIR
mkdir -p logs/split-patient

BARCODES="$OUT_DIR/barcodes.txt"


#zcat data/HSC_obs.csv.gz | grep "HSPC_pool$i" | tail -n+2| cut -d',' -f1,12 | sed 's/-1,/,/' | pigz -9 > $BARCODES
#zcat data/HSC_obs.csv.gz | grep "HSPC_pool$i" | tail -n+2| cut -d',' -f1,12 > $BARCODES

zcat data/HSC_obs.csv.gz | grep "HSPC_pool$i" | tail -n+2| cut -d',' -f1,12 | sed 's/-1,/,/' | tr ',' '\t' > $BARCODES

sinto filterbarcodes -b "$OUT_DIR/../${SAMPLE_NAME}_genome_tagged.bam" -c $BARCODES --barcodetag "CB" --outdir $OUT_DIR

