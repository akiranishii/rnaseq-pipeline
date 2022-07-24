#!/bin/bash
#SBATCH --job-name run_human_rna_seq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200gb
#SBATCH --time=14-00:00:00
#SBATCH --account macdouga99
#SBATCH --partition=largemem
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user anishii@umich.edu

cp -r /nfs/turbo/umms-macdouga/reference/mouse_reference .
cp -r /nfs/turbo/umms-macdouga/reference/human_reference .
. ./programs/run_fastqc.sh
. ./programs/run_multiqc.sh
. ./programs/run_human_star_index.sh
. ./programs/run_human_star.sh
. ./programs/run_multiqc_post_alignment.sh
