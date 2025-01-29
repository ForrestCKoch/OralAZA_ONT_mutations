#!/bin/bash
#PBS -l select=1:ncpus=1:mem=40gb
#PBS -l walltime=1:00:00
#PBS -M forrest.koch@unsw.edu.au
#PBS -m ae
#PBS -o logs/aligner
#PBS -j oe
#PBS -J 1-807
# PBS -J 1-8

# Same as v1, but edited not to exclude umi reads, we just trim them out .. 

set -o pipefail

ncpu=1
module load star/2.7.9a
module load samtools
module load fastp
module load subread

GENOME_DIR=$HOME/vf-scratch/forrest/HumanGenome/genome/STAR_sparse
GENOME=$HOME/vf-scratch/forrest/HumanGenome/genome/GRCh38.primary_assembly.genome.fa.gz
GTF=$HOME/vf-scratch/forrest/HumanGenome/genome/gencode.v44.annotation.gtf.gz

# Get the job identifier.
# to be run ...
i=$PBS_ARRAY_INDEX
#i=779

# This should be submitted from the SCIMETAR folder ...
cd $PBS_O_WORKDIR

# Grab the specific folder we are after
#source_dir=$(cat samples.txt | tail -n+$i | head -n1)
source_dir="data/novaseq/ZOU14201"
sample_ID=$(ls $source_dir | grep '.fastq.gz$' | cut -d'_' -f1-6| sort | uniq |tail -n+$i | head -n1)

out_dir=aligned_data/STAR-aligned_UMI-trimmed/$sample_ID/
mkdir -p $out_dir

# Get the read strands
R1=${source_dir}/${sample_ID}_R1_001.fastq.gz
R2=${source_dir}/${sample_ID}_R2_001.fastq.gz

export TMP_DIR=$out_dir/tmp

if [ ! -f $out_dir/trimmed_*.gz ]; then
    fastp -i $R1 -I $R2 -o $out_dir/trimmed_R1.fastq.gz -O $out_dir/trimmed_R2.fastq.gz \
          -g --poly_g_min_len 8 -x --poly_x_min_len 8 --failed_out $out_dir/fastp_failed.fastq.gz \
          -V -z 9 -w $ncpu -j $out_dir/fastp.json -h $out_dir/fastp.html | tee $out_dir/fastp.log


fi

if [ ! -f $out_dir/umi-extracted_*.gz ]; then
    umi_tools extract -I $out_dir/trimmed_R1.fastq.gz --read2-in $out_dir/trimmed_R2.fastq.gz \
        --stdout $out_dir/umi-extracted_R1.fastq.gz --read2-out $out_dir/umi-extracted_R2.fastq.gz --extract-method=regex \
        --filtered-out=$out_dir/filtered_R1.fastq.gz --filtered-out2=$out_dir/filtered_R2.fastq.gz \
        --umi-separator="_UMI_" \
        --bc-pattern="^(?P<discard_1>.*?ATTGCGCAATG)(?P<umi_1>.{8})." 2>&1 | tee $out_dir/umi_tools.log
    zcat $out_dir/umi-extracted_R1.fastq.gz $out_dir/filtered_R1.fastq.gz | pigz -9 > $out_dir/umi-trimmed_R1.fastq.gz
    zcat $out_dir/umi-extracted_R2.fastq.gz $out_dir/filtered_R2.fastq.gz | pigz -9 > $out_dir/umi-trimmed_R2.fastq.gz
fi



if [ ! -f $out_dir/${sample_ID}_Aligned.out.bam ]; then
    STAR --runThreadN $ncpu --quantMode GeneCounts --genomeDir $GENOME_DIR \
         --readFilesCommand pigz -dc --outTmpDir $out_dir/STAR_tmp \
         --readFilesIn $out_dir/umi-trimmed_R1.fastq.gz $out_dir/umi-trimmed_R2.fastq.gz \
         --outFileNamePrefix $out_dir/${sample_ID}_ --outSAMtype BAM Unsorted \
         --outSAMunmapped Within --outSAMattributes All | tee $out_dir/star.log
fi


if [ ! -f $out_dir/${sample_ID}_assigned_sorted.bam ]; then
	featureCounts -a $GTF -o $out_dir/gene_assigned \
		-R BAM $out_dir/${sample_ID}_Aligned.out.bam -p -B -C \
        -T $ncpu | tee $out_dir/featureCounts.log

    samtools sort $out_dir/${sample_ID}_Aligned.out.bam.featureCounts.bam \
			-o $out_dir/${sample_ID}_assigned_sorted.bam
	samtools index $out_dir/${sample_ID}_assigned_sorted.bam
fi
