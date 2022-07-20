#!/bin/bash

module load Bioinformatics
module load fastqc/0.11.8

mkdir fastqc

echo "Files for fastqc: " >> track_files.txt 

for f in clean/*; do
  fs="${f##*/}"
  echo "$fs" >> track_files.txt
  mkdir fastqc/"$fs"
  zcat clean/"$fs"/*.fq.gz | fastqc -o fastqc/"$fs" -t 8 stdin:"$fs"
done

