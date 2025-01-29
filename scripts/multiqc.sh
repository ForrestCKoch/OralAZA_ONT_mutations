#!/bin/bash

#PBS -l select=1:ncpus=4:mem=64gb
#PBS -l walltime=4:00:00
#PBS -M forrest.koch@unsw.edu.au
#PBS -m ae
#PBS -o logs
#PBS -j oe

ncpu=4

cd $PBS_O_WORKDIR

module load multiqc

multiqc -dd 3 -o multiqc/ . 
