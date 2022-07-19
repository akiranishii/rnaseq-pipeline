# MacLab End-to-End RNA-Seq Pipeline Walkthrough

## Introduction

The RNA sequencing (RNA-Seq) pipeline processes sequencing data to estimate transcript abundance and identify differentially expressed transcripts across samples. There are various tools that can be used to achieve similar results, though this program will outline the methods commonly used in the MacDougald laboratory. 

## Extraction of Raw Sequences

The RNA-Seq pipeline supports the input of raw next-generation sequencing (NGS) data in the FASTQ format. The data can be directly obtained from RNA-seq experiments submitted to Illumina or the Beijing Genome Institute or it may be collected from public repositories. Many published NGS studies in the **Gene Expression Omnibus (GEO)** provide direct links to the raw sequence data stored at the Sequence Read Archive (SRA). The sequence data from SRA normally requires decompression and, sometimes, proper splitting to generate the right FASTQ files.

Most sequencing data in the MacDougald lab are derived from Illumina and the Beijing Genome Institute (BGI). Before you look at the raw sequences, you should ask yourself the following questions:

1. What is the **animal model** used (e.g. human, mouse)?
2. Are the reads **paired-end or strand-specific**? If they are strand-specific, is it forward or reversely stranded?

## RNA-seq pipeline outline

The RNA-seq pipeline conducts three analyses in parallel and 
