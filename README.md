# MacLab End-to-End RNA-Seq Pipeline Walkthrough

## Introduction

The RNA sequencing (RNA-Seq) pipeline processes sequencing data to estimate transcript abundance and identify differentially expressed transcripts across samples. There are various tools that can be used to achieve similar results, though this program will  make use of tools optimized for use in the MacDougald laboratory at the University of Michigan, Ann Arbor. 

## Extraction of Raw Sequences

The RNA-Seq pipeline supports the input of raw next-generation sequencing (NGS) data in the FASTQ format. The data can be directly obtained from RNA-seq experiments submitted to Illumina or the Beijing Genome Institute or it may be collected from public repositories. Many published NGS studies in the **Gene Expression Omnibus (GEO)** provide direct links to the raw sequence data stored at the Sequence Read Archive (SRA). The sequence data from SRA normally requires decompression and, sometimes, proper splitting to generate the right FASTQ files.

Most sequencing data in the MacDougald lab are derived from Illumina and the Beijing Genome Institute (BGI). Before you look at the raw sequences, you should ask yourself the following questions:

1. What is the **animal model** used (e.g. human, mouse)?
2. Are the reads **paired-end or strand-specific**? If they are strand-specific, is it forward or reversely stranded?

## RNA-seq pipeline outline

The RNA-seq pipeline conducts three analyses in parallel. The first analysis uses STAR (v2.7.5a) for both the alignment and quantification procedures. The second analysis uses STAR (v2.7.5a) for alignment, then FeatureCounts for quantification. The final procedures uses Salmon (v0.11.30) for pseudoalignment. The result from the three analyses are compared after differential expression analysis is conducted. A more detailed outline of the RNA-seq pipeline can be viewed in the image below. 

<img width="823" alt="rnaseq" src="https://user-images.githubusercontent.com/62619033/179835069-b15add0e-0777-460c-a817-dff6d835d4be.png">

## How to use the RNA-seq pipeline
The use of the RNA-seq pipeline requires four steps: (1) Obtaining all necessary files from GitHub, (2) Preparing input files for the RNA-seq pipeline, (3) Uploading all files onto the Great Lakes Server, and (4) Running the RNA-seq pipeline on the Great Lakes Server. All steps are outlined below. 

### STEP 1: Access the Great Lakes Server (https://arc.umich.edu/greatlakes/user-guide/)

1. Install the Cisco VPN client (https://its.umich.edu/enterprise/wifi-networks/vpn) and login using your level 1 password. You may need to fill out a form on the Great Lakes Server website to create an account if you do not have one already.
2. Search “http://ssh uniqname@greatlakes.arc-ts.umich.edu" in the browser.
3. It will navigate you to a log in page. Use your level 1 account to access it.
4. Click files -> home directory.
5. Click "GO TO" on the top bars, and type the path "/nfs/turbo/umms-macdouga/" into it. Click ok.
6. Create a new directory for yourself. This directory will be referred to as [new directory] in later steps. 
7. Upload all files from this repository into [new directory].
8. Create a new directory called "clean" inside [new directory].
9. Upload all of your raw RNA-seq files (.fastq.gz) into the directory labeled "clean". Make sure each sample has its own directory labeled with sample names of your choosing.

### STEP 2: Run RNA-seq analysis through the Great Lakes Server

1. Type the following into your terminal: `ssh [replace with uniqname]@greatlakes.arc-ts.umich.edu`
2. Use your level 1 account to login.
3. Navigate to your folder of interest by using the "cd" command: `cd /nfs/turbo/umms-macdouga/[new directory]` 
4. Run the command: `sbatch run_human_rna_seq` or `sbatch run_mouse_rna_seq` (depending on animal model)
5. To cancel a run, type the command: `scancel [insert Slurm Job ID]`

### Sections currently being added to the guide

* Limitations to the pipeline
* Sample of descriptions to be used in papers
* Past publications that used these tools
* Images for visualization of directions
* File prep and github usage steps
* Last updated
* Detailed explanation of steps
* Citations and additional resources

## Helpful links and resources

* Fastq file format: https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html
* Great Lakes Server: https://arc.umich.edu/greatlakes/user-guide/
* FASTQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
* FASTQC interpretation: https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/qc_fastqc_assessment.html
* MultiQC: https://multiqc.info/
* Cutadapt: https://cutadapt.readthedocs.io/en/stable/#cutadapt
* STAR: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
* FeatureCounts: https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html
* Salmon: https://salmon.readthedocs.io/en/latest/salmon.html
* DESeq2: [https://bioconductor.org/packages/devel/bioc/manuals/DESeq2/man/DESeq2.pdf](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
* Gene set enrichment analysis (GSEA): https://www.gsea-msigdb.org/gsea/index.jsp

## Contact information

Please contact me at anishii@umich.edu if you experience any issues with using the pipeline. 



