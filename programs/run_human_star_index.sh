#!/bin/bash

module load Bioinformatics
module load star/2.7.5a

mkdir human_star_index

STAR --runMode genomeGenerate \
  --runThreadN 8 \
  --genomeDir human_star_index \
  --genomeFastaFiles human_reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --sjdbGTFfile human_reference/Homo_sapiens.GRCh38.105.gtf
