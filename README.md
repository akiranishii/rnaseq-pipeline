# MacLab End-to-End RNA-Seq Pipeline (Last updated: 07/23/2022)

## Introduction

The RNA sequencing (RNA-Seq) pipeline processes sequencing data to estimate transcript abundance and identify differentially expressed transcripts across samples. There are various tools that can be used to achieve similar results, though this program will  make use of tools optimized for use in the MacDougald laboratory at the University of Michigan, Ann Arbor. 

## Extraction of Raw Sequences

The RNA-Seq pipeline supports the input of raw next-generation sequencing (NGS) data in the FASTQ format. The data can be directly obtained from RNA-seq experiments submitted to **Illumina** or the **Beijing Genome Institute** or it may be collected from public repositories. Many published NGS studies in the **Gene Expression Omnibus (GEO)** provide direct links to the raw sequence data stored at the Sequence Read Archive (SRA). The sequence data from SRA normally requires decompression and, sometimes, proper splitting to generate the right FASTQ files.

Most sequencing data in the MacDougald lab are derived from Illumina and the Beijing Genome Institute (BGI). Before you look at the raw sequences, you should ask yourself the following questions:

1. What is the **animal model** used (e.g. human, mouse)?
2. Are the reads **paired-end or strand-specific**? If they are strand-specific, is it forward or reversely stranded?

## RNA-seq pipeline outline

The RNA-seq pipeline conducts three analyses in parallel. The first analysis uses STAR (v2.7.5a) for both the alignment and quantification procedures. The second analysis uses STAR (v2.7.5a) for alignment, then FeatureCounts for quantification. The final procedures uses Salmon (v0.11.30) for pseudoalignment. The result from the three analyses are compared after differential expression analysis is conducted. A more detailed outline of the RNA-seq pipeline can be viewed in the image below. 

<p align="center" width="100%">
<img width="823" alt="rnaseq" src="https://user-images.githubusercontent.com/62619033/179835069-b15add0e-0777-460c-a817-dff6d835d4be.png">
</p>

## Limitations of the RNA-seq pipeline

The MacLab RNA-seq pipeline was made with the most common use cases of RNA-seq in the MacDougald lab in mind. Therefore, the pipeline will not be able to conduct analyses in the following situations:

* The animal model is not human or mouse.
* The researcher is interested in identifying novel trascripts, rather than just identifying differentially expressed genes. 
* The researcher requires customizable plots.
* Adapters and low quality reads were not removed during the sequencing procedure (The UM Bioinformatics Core and BGI automatically conducts these procedures).
* The researcher is interested in controlling for batch effects or other variables (e.g. sex, age).

*If any of the above procedures/analyses are required, please contact anishii@umich.edu with necessary product adjustments.*

## How to use the RNA-seq pipeline
The use of the RNA-seq pipeline requires four steps: (1) Obtaining all necessary files from GitHub, (2) Preparing input files for the RNA-seq pipeline, (3) Uploading all files onto the Great Lakes Server, and (4) Running the RNA-seq pipeline on the Great Lakes Server. All steps are outlined below. 

### STEP 1: Prepare input files and download the raw sequencing data

1. In excel, create a file that contains the following columns: (1) Sample, which contains the sample names, and (2) TreatmentX (replace X with a number between 1 and 10), which contains the treatment categorization of each group. There can be multiple treatment columns, but please make sure they are labeled with different numbers. Once made, please save this file as `samples.csv`. An example of this file can be found below: 

<p align="center" width="100%">
<img width="195" alt="Screen Shot 2022-07-24 at 9 14 41 AM" src="https://user-images.githubusercontent.com/62619033/180648658-6a71cecc-52fc-428d-8c31-ba8e45896e5b.png">
</p>

2. Download the raw sequencing data (.fastq.gz format). Place all of your raw RNA-seq files (.fastq.gz) into the directory labeled "clean". Make sure each sample is placed in its own directory with a unique sample name that matches the sample name on the samples.csv file created in step 1 (this is the default file organization when sequencing data is received from BGI). *PLEASE MAKE SURE YOUR FILES ARE BACKED UP IN AN EXTERNAL HARD DRIVE OR OTHER DATA STORAGE LOCATIONS BEFORE PROCEEDING TO FURTHER STEPS.* 

<p align="center" width="100%">
<img width="1106" alt="Screen Shot 2022-07-24 at 10 19 36 AM" src="https://user-images.githubusercontent.com/62619033/180651317-7780b7c1-0fac-4552-ac02-c02dae5fbd64.png">
</p>

### STEP 1.5 (OPTIONAL): How to download files from AWS on the Great Lakes Server

1. Navigate to the directory you are interested in downloading the files in using the `cd` command
2. Load the AWS command line interface module by using the command: `module load aws-cli/2`
3. Response to prompts (BGI should provide you with this information) after using the command: `aws configure` 
4. Download files of interest into the current directory by using the command: `aws s3 sync s3://[REPLACE WITH BUCKET NAME] .`

*NOTE: If you use this option, please make sure the files and folders are labeled and organized correctly, as described in STEP 1.*

### STEP 2: Access the Great Lakes Server via the browser

1. Install the Cisco VPN client (https://its.umich.edu/enterprise/wifi-networks/vpn) and login using your level 1 password.
2. Search “http://ssh uniqname@greatlakes.arc-ts.umich.edu" in the browser.
3. It will navigate you to a log in page. Use your level 1 account to access it. 
4. Click files -> home directory.
5. Click "GO TO" on the navigation bar, and type the path "/nfs/turbo/umms-macdouga/" into it. Click ok.
6. Create a new directory for yourself. This directory will be referred to as [NEW DIRECTORY NAME] in later steps.
7. Upload the directory labeled "clean" and the file labeled "samples.csv" from STEP 1 into [NEW DIRECTORY NAME].

*NOTE: If you have difficulty with obtaining permissions for any of the steps above, you may need to fill out a form on the Great Lakes Server website to create an account (https://arc.umich.edu/greatlakes/user-guide/).* 

### STEP 3: Access the Great Lakes Server via the terminal

1. Type the following into your terminal: `ssh [REPLACE WITH UNIQNAME]@greatlakes.arc-ts.umich.edu`
2. Use your level 1 account to login.
3. Navigate to new directory by typing the command: `cd /nfs/turbo/umms-macdouga/[NEW DIRECTORY NAME]`
4. Download the necessary reference genome files by using one of the following commands (depending on the animal model of interest):  
`cp -r ../reference/human_reference .`  
`cp -r ../reference/mouse_reference .`  
5. Load the git module onto the Great Lakes Server by using the command: `module load git`
6. Obtain the necessary programs from the RNA-seq pipeline using the command: `git clone https://github.com/akiranishii/rnaseq-pipeline.git`
7. Run the following commands in sequence to reorganize the files for the pipeline:  
`mv -v rnaseq-pipeline/* .`  
`rm -r rnaseq-pipeline`  
8. Run the analysis by using the command: `sbatch run_human_rna_seq.sh` or `sbatch run_mouse_rna_seq.sh` (make sure you choose the correct animal model)
9. To cancel a run, type the command: `scancel [insert Slurm Job ID]`

*NOTE: By default, email notifications about when the run STARTED, ENDED, and FAILED will be sent to rnaseq.log@gmail.com (password: A.kphqXPzEux-c2CYnoxtt_k). This can be changed in the run_human_rna_seq.sh and run_mouse_rna_seq.sh documents, if desired.*

## Output files

## Post-Pipeline File Cleanup

The Great Lakes Server is NOT meant for long-term storage of data files. Therefore, please transfer all files to an external hard drive or other data storage location. Once all files are safely stored at an external location, the files can be removed by using the following commands in sequence:

`cd /nfs/turbo/umms-macdouga/`  
`rm -r [NEW DIRECTORY NAME]`

## Past Publications that use the tools in the RNA-seq pipeline

* https://pubmed.ncbi.nlm.nih.gov/34774798/
* https://pubmed.ncbi.nlm.nih.gov/33979328/
* https://pubmed.ncbi.nlm.nih.gov/32919095/

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



