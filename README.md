# MacLab End-to-End RNA-Seq Pipeline (Last updated: 08/22/2022)

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
<img width="956" alt="rnaseq" src="https://user-images.githubusercontent.com/62619033/186000990-233c89aa-5029-4a54-8942-2fc633405dbb.png">
</p>

## Limitations of the RNA-seq pipeline

The MacLab RNA-seq pipeline was made with the most common use cases of RNA-seq in the MacDougald lab in mind. Therefore, the pipeline will not be able to conduct analyses in the following situations:

* The animal model is not human or mouse.
* The researcher is interested in identifying novel trascripts, rather than just identifying differentially expressed genes. 
* The researcher requires customizable plots.
* Adapters and low quality reads were not removed during the sequencing procedure (The UM Bioinformatics Core and BGI automatically conducts these procedures).
* The researcher is interested in controlling for batch effects or other variables (e.g. sex, age).

*If any of the above procedures/analyses are required, please contact Akira Nishii (anishii@umich.edu), Rebecca Schill (rschill@umich.edu), or Jessica Maung (jmaung@umich.edu) with an explanation of necessary adjustments.*

## How to use the RNA-seq pipeline

The use of the RNA-seq pipeline requires four steps: (1) preparing input files and downloading the raw sequencing data, (2) accessing the Great Lakes Server via the browser, (3) accessing the Great Lakes Server via the terminal, and (4) running the RNA-seq Pipeline. All steps are outlined below. 

### STEP 1: Prepare input files and download the raw sequencing data

1. In excel, create a file that contains the following columns: (1) Sample, which contains the sample names, and (2) TreatmentX (replace X with a number between 1 and 10), which contains the treatment categorization of each group. There can be multiple treatment columns, but please make sure they are labeled with different numbers. Once made, please save this file as `samples.csv`. An example of this file can be found below: 

<p align="center" width="100%">
<img width="195" alt="Screen Shot 2022-07-24 at 9 14 41 AM" src="https://user-images.githubusercontent.com/62619033/180648658-6a71cecc-52fc-428d-8c31-ba8e45896e5b.png">
</p>

2. In excel, create another file that contains the following columns: (1) Group, which contains one of the Treatments from step 1, (2)  Control, which contains the control group for the comparison you are interested in conducting for the Treatment specified in the Group column, and (3) Treatment, which contains the treatment group for the comparison you are interested in conducting for the Treatment specified in the Group column. Once made, please save this file as `comparisons.csv`. PLEASE MAKE SURE THE `samples.csv` AND `comparisons.csv` FILES ARE CONSISTENT IN LABELING OF TREATMENT NAMES AND GROUPS. An example of this file can be found below: 

<p align="center" width="100%">
<img width="349" alt="Screen Shot 2022-07-24 at 5 45 33 PM" src="https://user-images.githubusercontent.com/62619033/180667045-22370a55-00a7-4427-9c1e-dc92eee0ffe5.png">
</p>

3. Download the raw sequencing data (.fastq.gz format). Place all of your raw RNA-seq files (.fastq.gz) into the directory labeled `clean`. Make sure each sample is placed in its own directory with a unique sample name that matches the sample name on the samples.csv file created in step 1 (this is the default file organization when sequencing data is received from BGI). *PLEASE MAKE SURE YOUR FILES ARE BACKED UP IN AN EXTERNAL HARD DRIVE OR OTHER DATA STORAGE LOCATIONS BEFORE PROCEEDING TO FURTHER STEPS.* 

<p align="center" width="100%">
<img width="1106" alt="Screen Shot 2022-07-24 at 10 19 36 AM" src="https://user-images.githubusercontent.com/62619033/180667095-f8f4f34c-04d9-4bf2-a779-c6183742e318.png">
</p>

4. **Optional (for complex comparisons).** For comparisons that require subsetting of data by more than two or more treatment groups, provide a third file saved as `exclusions.csv`. In this file, please provide the following columns: (1) A group column, which is identical to the group column in the `comparisons.csv` file and (2) seperate columns for treatment groups you would like to subset the data by. Once made, please provide information on which variables you would like to exclude from your data before making your comparisons. For example, for the second comparison in the `comparisons.csv` file (Treatment1 B vs C), if you are actually only interested in comparing B vs C among samples who are classified as 'ab' in Treatment 2, you would need to exclude the samples that are labeled 'up'. Therefore, the second row (for the second comparison) and the column labeled "Treatment2" of the `exclusions.csv` file will contain the word 'up'. This translates to "Compare B vs C in Treatment1 after excluding 'up' in the data set". An alternative translation would be, "Compare B vs C for 'ab' only." An example of this file can be found below: 

<p align="center" width="100%">
<img width="717" alt="optional" src="https://user-images.githubusercontent.com/62619033/186008525-b497bbab-0593-4245-8fd2-5dbf49b27541.png">
</p>

### STEP 1.5 (OPTIONAL): How to download files from AWS on the Great Lakes Server

1. Navigate to the directory you are interested in downloading the files in using the `cd` command
2. Load the AWS command line interface module by using the command: `module load aws-cli/2`
3. Response to prompts (BGI should provide you with this information) after using the command: `aws configure` 
4. Download files of interest into the current directory by using the command: `aws s3 sync s3://[REPLACE WITH BUCKET NAME] .`

*NOTE: If you use this option, please make sure the files and folders are labeled and organized correctly, as described in STEP 1.*

### STEP 2: Access the Great Lakes Server via the browser

1. Install the Cisco VPN client (https://its.umich.edu/enterprise/wifi-networks/vpn) and login using your level 1 password. Use the "UMVPN - Only U-M Traffic" option. 
2. Search “http://ssh uniqname@greatlakes.arc-ts.umich.edu" in the browser.
3. It will navigate you to a log in page. Use your level 1 account to access it. 
4. Click files -> home directory.
5. Click "GO TO" on the navigation bar, and type the path "/nfs/turbo/umms-macdouga/" into it. Click ok.
6. Create a new directory for yourself. This directory will be referred to as `[NEW DIRECTORY NAME]` in later steps.
7. Upload the directory labeled `clean` and the files labeled `samples.csv` and `comparisons.csv` from STEP 1 into `[NEW DIRECTORY NAME]`.

*NOTE: If you have difficulty with obtaining permissions for any of the steps above, you may need to fill out a form on the Great Lakes Server website to create an account (https://arc.umich.edu/greatlakes/user-guide/).* 

### STEP 3: Access the Great Lakes Server via the terminal

1. Type the following into your terminal: `ssh [REPLACE WITH UNIQNAME]@greatlakes.arc-ts.umich.edu`
2. Use your level 1 account to login.
3. Navigate to new directory by typing the command: `cd /nfs/turbo/umms-macdouga/[NEW DIRECTORY NAME]`
4. Load the git module onto the Great Lakes Server by using the command: `module load git`
5. Obtain the necessary programs from the RNA-seq pipeline using the command: `git clone https://github.com/akiranishii/rnaseq-pipeline.git`
6. Run the following commands in sequence to reorganize the files for the pipeline:  
`mv -v rnaseq-pipeline/* .`  
`rm -r rnaseq-pipeline`  

### STEP 4: Run the RNA-seq Pipeline

1. Navigate to the correct directory by using the command: `cd /nfs/turbo/umms-macdouga/[NEW DIRECTORY NAME]`
2. Run the analysis by using one of the following commands (make sure you choose the correct animal model):  
`sbatch run_human_rna_seq.sh` *or* `sbatch run_mouse_rna_seq.sh`
3. You may now exit the terminal. The analysis will take approximately **7 hours** for around 10 samples under current settings. To cancel a run at any time, type the command: `scancel [INSERT SLURM JOB ID HERE]`

*NOTE: By default, email notifications about when the run STARTED, ENDED, and FAILED will be sent to rnaseq.log@gmail.com (password: A.kphqXPzEux-c2CYnoxtt_k). This can be changed in the run_human_rna_seq.sh and run_mouse_rna_seq.sh documents, if desired.*

## Output files

All resulting files from the RNA-seq pipeline will be located in a folder called "Output." The "Output" folder will contain three subfolders called "Figures_STAR," "Figures_Salmon," and "Cross_Validation." Each subfolder will contain additional subfolders for each comparison specified for the analysis. "Figures_STAR" will contain all figures generated using the regular alignment procedure, whereas the "Figures_Salmon" will contain all figures generated using the pseudoalignment procedure. The "Cross Validation" file will contain all figures generated for the purpose of cross-checking the results from regular alignment and pseudoalignment.  

There are three types of figures in the "Figures_STAR" and "Figures_Salmon" subfolders. (1) Quality control files, (2) Differential expression files, and (3) Pathway analysis files. This section will outline the location and interpretation of output files that will result from running the RNA-seq pipeline. 

**1. Quality control files are located in a subfolder called "QualityControl". The folder contains the following figures for each comparison:**

* Dispersion plot  

<p align="center" width="100%">
<img width="662" alt="Screen Shot 2022-08-22 at 4 58 23 PM" src="https://user-images.githubusercontent.com/62619033/186018514-f0baa967-ee3f-4a6a-9dc6-5cae13952955.png">
 </p>
 
* Heatmap

<p align="center" width="100%">
<img width="787" alt="Screen Shot 2022-08-22 at 4 59 17 PM" src="https://user-images.githubusercontent.com/62619033/186018550-ebda172a-7899-474f-941d-b4148447785f.png">
 </p>

* Library size

<p align="center" width="100%">
<img width="549" alt="Screen Shot 2022-08-22 at 5 04 55 PM" src="https://user-images.githubusercontent.com/62619033/186019127-f31bbc2f-c709-4e3c-89e7-a80302e8edd0.png">
 </p>
 
* PCA plot

<p align="center" width="100%">
<img width="602" alt="Screen Shot 2022-08-22 at 4 59 44 PM" src="https://user-images.githubusercontent.com/62619033/186018638-261808b6-f7e0-4130-849f-eae6e2c83054.png">
 </p>
 
* MA plot

<p align="center" width="100%">
<img width="434" alt="Screen Shot 2022-08-22 at 5 05 06 PM" src="https://user-images.githubusercontent.com/62619033/186019170-4b4a89fe-b6b8-41cf-b32f-b3e4244e3c15.png">
 </p>

**2. Differential expression files are located in a subfolder called "OutputTables." The folder contains the following data for each comparison:**

* Differential expression results in a .csv file format

<p align="center" width="100%">
<img width="939" alt="Screen Shot 2022-08-22 at 5 02 00 PM" src="https://user-images.githubusercontent.com/62619033/186020334-faf1dddf-6175-485f-bc60-b275faaa1f69.png">
</p>
 
* Depth normalized counts table

 <p align="center" width="100%">
<img width="1067" alt="Screen Shot 2022-08-22 at 5 00 48 PM" src="https://user-images.githubusercontent.com/62619033/186020234-82583278-1ef7-4884-b991-194b54eac77a.png">
</p>

* rlog normalized counts table
  
<p align="center" width="100%">
<img width="862" alt="Screen Shot 2022-08-22 at 5 01 27 PM" src="https://user-images.githubusercontent.com/62619033/186020270-2329a918-06a3-4c17-bfa9-b95c209ddacf.png">
</p>

* raw counts table
 
<p align="center" width="100%">
<img width="820" alt="Screen Shot 2022-08-22 at 5 01 07 PM" src="https://user-images.githubusercontent.com/62619033/186020253-577b3af2-17d0-4ad9-978e-e329e603ca9e.png">
</p>

* MA plot HTML (use this version of the MA plot to see which points correspond to which genes)

<p align="center" width="100%">
<img width="1432" alt="Screen Shot 2022-08-22 at 5 15 54 PM" src="https://user-images.githubusercontent.com/62619033/186020543-b43af33b-a910-4e5e-9cd9-ddfb4ad818f3.png">
</p>


**3. The Pathway analysis files will be located in a subfolder called "PathwayAnalysis." The folder contains the following data for each comparison:**

* `.rnk` file used as input for Gene Set Enrichment Analysis (GSEA)
* GSEA output folders (analysis conducted for 'Hallmark gene set,' 'Curated gene set,' and 'Regulatory target gene set'). More information on GSEA and available gene sets can be accessed here: http://www.gsea-msigdb.org/gsea/msigdb/index.jsp

## Recommended: Long-term File Storage on the Data Den

The Data Den is an archive optimized for data that is not regularly accessed for extended periods of time (weeks to years). Once you are done running the RNA-seq pipeline, it is recommended that you store all raw data files on the Data Den. Several steps must be taken in order to transfer files from the Great Lakes Server to the Data Den.

1. Ask your PI to email arc-support@umich.edu to update the Data Den access list by providing your name and uniqname.
2. Create an user account on the Flux Transfer Service using this link: https://arc.umich.edu/login-request
3. Nagivate to: https://www.globus.org/
4. Login using your umich email
5. Select middle panel in the top right corner of the screen. This should open up two windows on the right and left side of the screen. See image below. 

<img width="229" alt="Screen Shot 2022-08-22 at 7 30 58 PM" src="https://user-images.githubusercontent.com/62619033/186037443-efe45c18-e1fb-41a0-8f93-32a026a23753.png">

6. On the right side of the screen, type `umich#greatlakes` in the "Collection" box. Type `/nfs/turbo/umms-macdouga/` in the "Path" box.
7. On the left side of the screen, type `umich#flux` in the "Collection" box. Type `/nfs/dataden/umms-macdouga/` in the "Path" box.
8. Select the folder you are interested in moving from the Great Lakes Server to the Data Den on the left side of the screen by checking the square box next to the folder. 
9. Click on the blue "start" button to transfer the folder from the Great Lakers Server to the Data Den.

<img width="1111" alt="Screen Shot 2022-08-22 at 7 11 07 PM" src="https://user-images.githubusercontent.com/62619033/186038023-4ce686d1-226c-4550-ac6f-ac7714106af1.png">


## Post-Pipeline File Cleanup

The Great Lakes Server is NOT meant for long-term storage of data files. Therefore, please transfer all files to an external hard drive or other data storage location (such as the UM Data Den). Once all files are safely stored at an external location, the files can be removed by using the following commands in sequence:

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



