# Install libraries
install.packages(c("kableExtra", "tidyverse", "ggplot2", "ggrepel", "Cairo","pheatmap", "RColorBrewer"), repos="https://repo.miserver.it.umich.edu/cran")

# Install bioconductor packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos="https://repo.miserver.it.umich.edu/cran")

BiocManager::install("DESeq2")
BiocManager::install("biomaRt")

#Load libraries
suppressPackageStartupMessages(library('knitr', character.only=TRUE))
suppressPackageStartupMessages(library('kableExtra', character.only=TRUE))
suppressPackageStartupMessages(library('tidyverse', character.only=TRUE))
suppressPackageStartupMessages(library('ggrepel', character.only=TRUE))
suppressPackageStartupMessages(library('Cairo', character.only=TRUE))
options(bitmapType = "cairo")

suppressPackageStartupMessages(library('DESeq2', character.only=TRUE))
suppressPackageStartupMessages(library('pheatmap', character.only=TRUE))
suppressPackageStartupMessages(library('RColorBrewer', character.only=TRUE))
suppressPackageStartupMessages(library('biomaRt', character.only=TRUE))
source("deseq2_functions.R")

# Data extraction
data <- read_csv("samples.csv", show_col_types = FALSE)
comparisons <- read_csv("comparisons.csv", show_col_types = FALSE)

#Run deseq2
system("mkdir ./Figures")
for (i in 1:nrow(comparisons)) {
  #Create necessary variables
  group <- comparisons[[i,1]]
  control <- comparisons[[i,2]]
  treatment <- comparisons[[i,3]]
  group_sym <- sym(group)
  Comparison <- paste0(group, "_", treatment, "_vs_", control)
  
  #Paths
  plotPath = paste0("./Figures/", Comparison, "/QualityControl/")
  plotPath2 = paste0("./Figures/", Comparison, "/OutputTables/")
  
  #Generate necessary folders
  generate_folders()
  
  #Subset dataset as necessary
  metadata = eval(expr(subset(data, !!sym(group) == control | !!sym(group) == treatment)))
  kable(metadata, caption = "Samples") %>%
    kable_styling(latex_options = "hold_position")
  
  #Get data from count files
  counts <- get_data()
  
  #Plot library size
  plot_library_size(counts)
  
  #Generate matrix to feed into DESEQ2 program
  cts <- matrix_convert()
  
  #Generate meta data dataframe for DESEQ2 program
  coldata <- coldata_gen()
  
  #Check before proceeding
  print(all(rownames(coldata) == colnames(cts)))
  
  #run deseq2
  res <- run_deseq2_program()
  dds <- res$dds_
  rld <- res$rld_
  
  #Get deseq2 results
  res <- results(dds, name=Comparison)
  summary(res)
  
  #MA plot
  plot_ma()
  
  # PCA plot
  plot_pca()
  
  #Output tables
  generate_output_tables()
  
  #Heatmaps
  generate_heatmaps()
  
  #Convert to mouse symbol
  res_table_merge <- convert_human_symbol()
  
  #Output file
  diff_output()
}