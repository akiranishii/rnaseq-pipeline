#!/bin/bash
#SBATCH --job-name run_mouse_rna_seq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200gb
#SBATCH --time=14-00:00:00
#SBATCH --account macdouga99
#SBATCH --partition=largemem
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user log.rnaseq@gmail.com

cp -r /nfs/turbo/umms-macdouga/reference/mouse_reference .
cp -r /nfs/turbo/umms-macdouga/reference/human_reference .
. ./programs/run_fastqc.sh
. ./programs/run_multiqc.sh
. ./programs/run_mouse_star_index.sh
. ./programs/run_mouse_star.sh
. ./programs/run_multiqc_post_alignment.sh
. ./programs/run_mouse_deseq2.sh
