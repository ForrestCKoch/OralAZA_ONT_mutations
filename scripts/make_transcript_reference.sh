#!/bin/sh

# Call me from the root redirectory of project

# Extracts just the transcripts of interest from the reference file.

REF="/srv/scratch/vafaeelab/forrest/HumanGenome/genome/gencode.v44.transcripts.fa.gz"
OUT="data/transcript_reference.fa.gz"

TOI=(ASXL1 ATM RUNX1 SRSF2 TET2 TP53)
TOI_REGEX="|$(echo ${TOI[@]}|sed 's# #|\\||#g')|"

#IFS='\n'
#
#for line in $(zcat $REF | tr -d '\n' | tr '>' '\n' | grep $TOI_REGEX | sed 's/^/>/g'); do
#        echo "$(echo $line|cut -d'|' -f1-8)|"
#        echo $line|cut -d'|' -f9
#done | 

zcat $REF | tr -d '\n' | tr '>' '\n' | grep $TOI_REGEX | sed 's/^/>/g' \
    | while read -r line ; do
        echo "$(echo $line|cut -d'|' -f1-8)|"
        echo $line|cut -d'|' -f9
done | pigz -9 > $OUT

