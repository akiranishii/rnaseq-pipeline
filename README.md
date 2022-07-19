# MacLab End-to-End RNA-Seq Pipeline Walkthrough

## Introduction

The RNA sequencing (RNA-Seq) pipeline processes sequencing data to estimate transcript abundance and identify differentially expressed transcripts across samples. There are various tools that can be used to achieve similar results, though this program will outline the methods commonly used in the MacDougald laboratory. 

## Extraction of Raw Sequences

The RNA-Seq pipeline supports the input of raw next-generation sequencing (NGS) data in the FASTQ format. The data can be directly obtained from RNA-seq experiments submitted to Illumina or the Beijing Genome Institute or it may be collected from public repositories. Many published NGS studies in the **Gene Expression Omnibus (GEO)** provide direct links to the raw sequence data stored at the Sequence Read Archive (SRA). The sequence data from SRA normally requires decompression and, sometimes, proper splitting to generate the right FASTQ files.

Most sequencing data in the MacDougald lab are derived from Illumina and the Beijing Genome Institute (BGI). Before you look at the raw sequences, you should ask yourself the following questions:

1. What is the **animal model** used (e.g. human, mouse)?
2. Are the reads **paired-end or strand-specific**? If they are strand-specific, is it forward or reversely stranded?

## RNA-seq pipeline outline

The RNA-seq pipeline conducts three analyses in parallel. The first analysis uses STAR (v2.7.5a) for both the alignment and quantification procedures. The second analysis uses STAR (v2.7.5a) for alignment, then FeatureCounts for quantification. The final procedures uses Salmon (v0.11.30) for pseudoalignment. The result from the three analyses are compared after differential expression analysis is conducted. A more detailed outline of the RNA-seq pipeline can be viewed in the image below. 

<img width="823" alt="rnaseq" src="https://user-images.githubusercontent.com/62619033/179835069-b15add0e-0777-460c-a817-dff6d835d4be.png">

## How to use the RNA-seq pipeline
The RNA-seq pipelien 

### STEP 1: Access the Great Lakes Server (https://arc.umich.edu/greatlakes/user-guide/)

1. Install the Cisco VPN client (https://its.umich.edu/enterprise/wifi-networks/vpn) and login using your level 1 password. You may need to fill out a form on the Great Lakes Server website to create an account if you do not have one already.
2. Search “http://ssh uniqname@greatlakes.arc-ts.umich.edu" in the browser.
3. It will navigate you to a log in page. Use your level 1 account to access it.
4. Click files -> home directory.
5. Click "GO TO" on the top bars, and type the path "/nfs/turbo/umms-macdouga/" into it. Click ok.
6. Create a directory for yourself. Inside your directory, create a directory called "clean".
7. Upload all of your raw RNA-seq files (.fastq) into the directory labeled "clean". Make sure each sample has its own directory labeled with the sample name of your choosing.

### STEP 2: Access the Great Lakes Server to the terminal

1. Type the following into your terminal: ssh [replace with uniqname]@greatlakes.arc-ts.umich.edu
2. Use your level 1 account to login.
3. 

