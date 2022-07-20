#!/bin/bash

module load python/3.9.7
pip install multiqc
multiqc fastqc/ -o fastqc/
