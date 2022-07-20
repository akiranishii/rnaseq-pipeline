#!/bin/bash

module load python/3.9.7
pip install multiqc
multiqc alignments_star/ -o alignments_star/
