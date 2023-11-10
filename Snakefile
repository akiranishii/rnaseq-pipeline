# Define global parameters
GENOME_DIR = "path/to/genome_dir"
REFERENCE_DIR = "path/to/reference"
ALIGNMENTS_DIR = "path/to/alignments"

# Define rule for running FastQC
rule fastqc:
    input:
        "clean/{sample}/*.fq.gz"
    output:
        "fastqc/{sample}"
    shell:
        "module load Bioinformatics; "
        "module load fastqc/0.11.8; "
        "mkdir -p {output}; "
        "zcat {input} | fastqc -o {output} -t 8 stdin:{wildcards.sample}"

# Define rule for running MultiQC before alignment
rule multiqc_pre:
    input:
        expand("fastqc/{sample}", sample=SAMPLES)
    output:
        "fastqc/multiqc_report.html"
    shell:
        "module load python/3.9.7; "
        "pip install multiqc; "
        "multiqc fastqc/ -o fastqc/"

# Define rule for STAR index generation
rule star_index:
    input:
        fasta="{reference_dir}/{reference}.dna.primary_assembly.fa",
        gtf="{reference_dir}/{reference}.105.gtf"
    output:
        directory("{genome_dir}/{reference}_star_index")
    shell:
        "module load Bioinformatics; "
        "module load star/2.7.5a; "
        "STAR --runMode genomeGenerate "
        "--runThreadN 8 "
        "--genomeDir {output} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf}"

# Define rule for STAR alignment
rule star_alignment:
    input:
        index=rules.star_index.output,
        fq1="clean/{sample}/*1.fq.gz",
        fq2="clean/{sample}/*2.fq.gz",
        gtf="{reference_dir}/{reference}.105.gtf"
    output:
        "alignments_star/{sample}/"
    shell:
        "module load Bioinformatics; "
        "module load star/2.7.5a; "
        "mkdir -p {output}; "
        "STAR --runThreadN 8 "
        "--genomeDir {input.index} "
        "--sjdbGTFfile {input.gtf} "
        "--readFilesIn {input.fq1} {input.fq2} "
        "--readFilesCommand zcat "
        "--outFileNamePrefix {output} "
        "--outSAMtype BAM SortedByCoordinate "
        "--quantMode GeneCounts "
        "--twopassMode Basic"

# Define rule for MultiQC after alignment
rule multiqc_post:
    input:
        expand("alignments_star/{sample}/", sample=SAMPLES)
    output:
        "alignments_star/multiqc_report.html"
    shell:
        "module load python/3.9.7; "
        "pip install multiqc; "
        "multiqc alignments_star/ -o alignments_star/"

# Define rules for DESeq2 analysis
rule deseq2:
    input:
        alignment=rules.star_alignment.output
    output:
        "deseq2/{sample}_results.txt"
    shell:
        "module load gcc/8.2.0; "
        "module load R; "
        "Rscript run_{reference}_deseq2.R {input.alignment}"
