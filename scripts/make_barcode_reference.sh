#!/bin/sh

# Call me from the root redirectory of project

# Extracts just the transcripts of interest from the reference file.

REF="data/HSC_obs.csv.gz"
OUT="data/barcodes.txt.gz"

POI=(HSPC_pool1 HSPC_pool2 HSPC_pool3_repeat HSPC_pool4)
POI_REGEX=",$(echo ${POI[@]}|sed 's# #,\\|,#g'),"

zcat $REF | grep $POI_REGEX | cut -d',' -f1 | cut -d'-' -f1 | pigz -9 > $OUT

