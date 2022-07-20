#!/bin/bash

module load Bioinformatics
module load star/2.7.5a

STAR --runMode genomeGenerate \
  --runThreadN 8 \
  --genomeDir mouse_star_index \
  --genomeFastaFiles mouse_reference/Mus_musculus.GRCm39.dna.primary_assembly.fa \
  --sjdbGTFfile mouse_reference/Mus_musculus.GRCm39.107.gtf
