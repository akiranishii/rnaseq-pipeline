# MacLab End-to-End RNA-Seq Pipeline (Last updated: 07/23/2022)

## Introduction

The RNA sequencing (RNA-Seq) pipeline processes sequencing data to estimate transcript abundance and identify differentially expressed transcripts across samples. There are various tools that can be used to achieve similar results, though this program will  make use of tools optimized for use in the MacDougald laboratory at the University of Michigan, Ann Arbor. 

## Extraction of Raw Sequences

The RNA-Seq pipeline supports the input of raw next-generation sequencing (NGS) data in the FASTQ format. The data can be directly obtained from RNA-seq experiments submitted to Illumina or the Beijing Genome Institute or it may be collected from public repositories. Many published NGS studies in the **Gene Expression Omnibus (GEO)** provide direct links to the raw sequence data stored at the Sequence Read Archive (SRA). The sequence data from SRA normally requires decompression and, sometimes, proper splitting to generate the right FASTQ files.

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

If any of the above procedures/analyses are required, please contact anishii@umich.edu with necessary product adjustments. 

## How to use the RNA-seq pipeline
The use of the RNA-seq pipeline requires four steps: (1) Obtaining all necessary files from GitHub, (2) Preparing input files for the RNA-seq pipeline, (3) Uploading all files onto the Great Lakes Server, and (4) Running the RNA-seq pipeline on the Great Lakes Server. All steps are outlined below. 

### STEP 1: Prepare input files and download the raw sequencing data

1. In excel, create a file that contains the following columns: (1) Sample, which contains the sample names, and (2) TreatmentX (replace X with a number between 1 and 10), which contains the treatment categorization of each group. There can be multiple treatment columns, but please make sure they are labeled with different numbers. Once made, please save this file as `samples.csv`. An example of this file can be found below: 

<p align="center" width="100%">
<img width="195" alt="Screen Shot 2022-07-24 at 9 14 41 AM" src="https://user-images.githubusercontent.com/62619033/180648658-6a71cecc-52fc-428d-8c31-ba8e45896e5b.png">
</p>

2. Download the raw sequencing data (.fastq.gz format). Place all of your raw RNA-seq files (.fastq.gz) into the directory labeled "clean". Make sure each sample is placed in its own directory with a unique sample name that matches the sample name on the samples.csv file created in step 1 (this is the default file organization when sequencing data is received from BGI). PLEASE MAKE SURE YOUR FILES ARE BACKED UP IN AN EXTERNAL HARD DRIVE OR OTHER DATA STORAGE LOCATIONS BEFORE PROCEEDING TO FURTHER STEPS. 

<p align="center" width="100%">
<img width="1106" alt="Screen Shot 2022-07-24 at 10 19 36 AM" src="https://user-images.githubusercontent.com/62619033/180651317-7780b7c1-0fac-4552-ac02-c02dae5fbd64.png">
</p>

### STEP 2: Access the Great Lakes Server

1. Install the Cisco VPN client (https://its.umich.edu/enterprise/wifi-networks/vpn) and login using your level 1 password. You may need to fill out a form on the Great Lakes Server website to create an account if you do not have one already (https://arc.umich.edu/greatlakes/user-guide/).
2. Search “http://ssh uniqname@greatlakes.arc-ts.umich.edu" in the browser.
3. It will navigate you to a log in page. Use your level 1 account to access it.
4. Click files -> home directory.
5. Click "GO TO" on the navigation bar, and type the path "/nfs/turbo/umms-macdouga/" into it. Click ok.
6. Create a new directory for yourself. This directory will be referred to as [NEW DIRECTORY] in later steps.
7. Upload the directory labeled "clean" from STEP 1 inside [NEW DIRECTORY].

### STEP 3: Run RNA-seq analysis through the Great Lakes Server

1. Type the following into your terminal: `ssh [REPLACE WITH UNIQNAME]@greatlakes.arc-ts.umich.edu`
2. Use your level 1 account to login.
3. Navigate to new directory by typing the command: `cd /nfs/turbo/umms-macdouga/[NEW DIRECTORY]`
4. Load the git module onto the Great Lakes Server by using the command: `module load git`
5. Obtain the necessary programs from the RNA-seq pipeline using the command: `git clone https://github.com/akiranishii/rnaseq-pipeline.git`
6. Run the following commands in sequence to reorganize the files, as necessary for the pipeline:  
`mv -v rnaseq-pipeline/* .`  
`rm -r rnaseq-pipeline`  
7. Run the analysis by using the command: `sbatch run_human_rna_seq` or `sbatch run_mouse_rna_seq` (ensure you choose the correct animal model)
8. To cancel a run, type the command: `scancel [insert Slurm Job ID]`

### Sections currently being added to the guide

* Sample of descriptions to be used in papers
* Past publications that used these tools
* Images for visualization of directions
* File prep and github usage steps
* Detailed explanation of steps
* Description of output files/how to access

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



