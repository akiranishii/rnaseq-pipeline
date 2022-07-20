#!/bin/bash

module load Bioinformatics
module load star/2.7.5a

mkdir alignments_star
for f in clean/*; do
  fs="${f##*/}"
  mkdir alignments_star/"$fs"
  STAR --runThreadN 8 \
    --genomeDir mouse_star_index \
    --sjdbGTFfile mouse_reference/Mus_musculus.GRCm39.107.gtf \
    --readFilesIn clean/"$fs"/*1.fq.gz clean/"$fs"/*2.fq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix alignments_star/"$fs"/ \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --twopassMode Basic
done
