#Function to make necessary directories to store output files
generate_folders <- function() {
  system(paste0("mkdir ./Figures/", Comparison))
  system(paste0("mkdir ./Figures/", Comparison, "/QualityControl"))
  system(paste0("mkdir ./Figures/", Comparison, "/OutputTables"))
}

#Function to get data from count files
read_counts <- function(sample, column = "X4") {
  file <- paste0("alignments_star/", sample, "/ReadsPerGene.out.tab")
  df <- read_tsv(file, skip = 4, col_names = FALSE, show_col_types = FALSE)
  df <- dplyr::select(df, X1, !!sym(column))
  df <- set_names(df, "gene_id", sample)
  return(df)
}

get_data <- function() {
  counts <- map_dfc(metadata$Sample, read_counts)
  counts <- dplyr::select(counts, "gene_id...1", metadata$Sample)
  counts <- dplyr::rename(counts, "gene_id" = `gene_id...1`)
  return(counts)
}

plot_library_size <- function(counts) {
  size <- colSums(dplyr::select(counts, -gene_id))
  size <- as.data.frame(size) %>%
    rownames_to_column("Sample") %>%
    mutate(size = size / 10^6) %>%
    left_join(metadata, by = "Sample") %>%
    mutate(Sample = factor(Sample, levels = metadata$Sample))
  
  # plot
  eval(expr(p <- ggplot(size, aes(x = Sample, y = size, fill = !!group_sym)) +
              geom_hline(aes(yintercept = 30), linetype = "dashed") +
              geom_col(color = "black") +
              geom_text(aes(label = round(size, 0), y = size + 3)) +
              scale_y_continuous(expand = c(0, 0),
                                 limits = c(0, max(size$size) + 5)) +
              theme_void() +
              theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                               color = "black"),
                    plot.title = element_text(hjust = 0.5)) +
              ggtitle("Library size (M)")))
  ggsave(paste0('Figures/', Comparison,'/QualityControl/', 'library_size_', group, '_', control, '_', treatment, '.png'))
}

# function to convert to matrix
matrix_convert <- function() {
  cts <- as.matrix(counts[, c(2:ncol(counts))]) 
  row.names(cts) <- counts$gene_id 
  str(cts)
  return(cts)
}

coldata_gen <- function() {
  
  eval(expr(coldata <- metadata %>% dplyr::select(!!group_sym)))
  a <- expression(coldata$group <- as.factor(coldata$group))
  eval(do.call("substitute", list(a[[1]], list(group = group_sym))))
  rownames(coldata) <- colnames(cts) ##this is a shortcut
  return(coldata)
  
}

#run deseq2
run_deseq2_program <- function() {
  dds <- eval(expr(DESeqDataSetFromMatrix(countData = cts,
                                          colData = coldata,
                                          design = ~!!sym(group))))
  
  a <- expression(dds$group <- relevel(dds$group, ref = control))
  eval(do.call("substitute", list(a[[1]], list(group = group_sym))))
  
  colData(dds) 
  #rlog normalization for subsequent plots
  keep <- rowSums(counts(dds)) >= 1
  dds <- dds[keep, ]
  
  rld <- rlog(dds, blind = FALSE)
  
  # Apply model
  dds <- DESeq(dds)
  resultsNames(dds)
  
  return(list("dds_" = dds, "rld_" = rld))
}

plot_ma <- function() {
  png(paste0(plotPath, Comparison, "_MA_plot.png"))
  plotMA(res, ylim = c(-2,2))
  dev.off()
}


plot_pca <- function() {
  
  pdf(file = paste0(plotPath, 'PCAplot_All_', Comparison, '.pdf'), onefile = TRUE)
  #PCA plot for Rlog-Normalized counts for all samples
  a <- expression(CombinatoricGroup <- factor(coldata$group))
  eval(do.call("substitute", list(a[[1]], list(group = group_sym))))
  SampleName <- factor(row.names(coldata))
  p.all <- plotPCA(rld, intgroup = c(group))
  gp <- ggplot(p.all$data, aes(x = PC1, y = PC2, color = SampleName, shape = CombinatoricGroup)) + xlab(p.all$labels[2]) + ylab(p.all$labels[1]) + scale_shape_manual(values=1:nlevels(CombinatoricGroup), name = 'Combinatoric Group') + geom_point(size=7) + ggtitle(label = as.character('All samples Rlog-Normalized')) + theme(plot.title = element_text(hjust = 0.5)) + guides(colour=guide_legend(nrow=12, title = 'Sample'), legend.key = element_rect(size = 1), legend.key.size = unit(0, 'cm')) + theme_classic(base_size = 10) + theme(legend.margin=margin(t = 0, unit='mm'))
  plot(gp)
  
  #PCA for Depth-Normalized counts for all samples
  DNC <- SummarizedExperiment(log2(counts(dds, normalized = TRUE)), colData=colData(dds))
  p.DNC <- plotPCA(DESeqTransform(DNC), intgroup = group)
  gpDNC <- ggplot(p.DNC$data, aes(x = PC1, y = PC2, color = SampleName, shape = CombinatoricGroup)) + xlab(p.DNC$labels[2]) + ylab(p.DNC$labels[1]) + scale_shape_manual(values=1:nlevels(CombinatoricGroup), name = 'Combinatoric Group') + geom_point(size=2) + ggtitle(label = as.character('All samples Depth-Normalized')) + theme(plot.title = element_text(hjust = 0.5)) + guides(colour=guide_legend(nrow=12, title = 'Sample'), legend.key = element_rect(size = 1), legend.key.size = unit(0, 'cm')) + theme_classic(base_size = 10) + theme(legend.margin=margin(t = 0, unit='mm'))
  plot(gpDNC)
  
  #PCA for Raw counts for all samples
  RC <- SummarizedExperiment(log2(counts(dds, normalized = FALSE)), colData=colData(dds))
  p.RC <- plotPCA(DESeqTransform(RC), intgroup = group)
  gpRC <- ggplot(p.RC$data, aes(x = PC1, y = PC2, color = SampleName, shape = CombinatoricGroup)) + xlab(p.RC$labels[2]) + ylab(p.RC$labels[1]) + scale_shape_manual(values=1:nlevels(CombinatoricGroup), name = 'Combinatoric Group') + geom_point(size=2) + ggtitle(label = as.character('All samples Raw')) + theme(plot.title = element_text(hjust = 0.5)) + guides(colour=guide_legend(nrow=12, title = 'Sample'), legend.key = element_rect(size = 1), legend.key.size = unit(0, 'cm')) + theme_classic(base_size = 10) + theme(legend.margin=margin(t = 0, unit='mm'))
  plot(gpRC)
  
  dev.off()
}

generate_output_tables <- function() {
  countsdds <- counts(dds, normalized = FALSE)
  write.csv(countsdds, file=paste0(plotPath2, "DESeq2_raw_counts.csv"))
  
  countsdds <- counts(dds, normalized = TRUE)
  write.csv(countsdds, file=paste0(plotPath2,"DESeq2_depthNormalized_counts.csv"))
  
  rld2 <- assay(rld) # includes model if 'blind=FALSE'
  write.csv(rld2, file=paste0(plotPath2,"DESeq2_rlogNormalized_counts.csv"))
}

generate_heatmaps <- function() {
  #plot dispersions
  pdf(file = paste0(plotPath,'Dispersion_', Comparison, '.pdf'), onefile = FALSE)
  disp <- plotDispEsts(dds)
  dev.off()
  
  #heatmap of normalized data, sample distibution matrix
  sampleDists <- dist(t(assay(rld))) #rlog normalized data
  sampleDistMatrix <- as.matrix(sampleDists)
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette(rev(brewer.pal(9, 'Blues')))(255)
  pdf(file = paste0(plotPath,'Heatmap_Dispersions_', Comparison, '.pdf'), onefile = FALSE)
  pheatmap(sampleDistMatrix, 
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  dev.off()
  
  #heatmap with top 500 variant or expressed genes, rlog normalized data
  colors <- colorRampPalette(brewer.pal(9, 'Blues'))(255)
  select <- order(rowVars(assay(rld)), decreasing=TRUE)[1:500]
  df <- data.frame(Group = colData(rld)[,c(group)], row.names = rownames(colData(dds)))
  pdf(file = paste0(plotPath,'Heatmap_TopVar_', Comparison, '.pdf'), onefile = FALSE, width=10, height=20)
  pheatmap(assay(rld)[select,], scale="row",  cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=df, fontsize = 7, las = 2, fontsize_row = 7, color = colors, main = '500 Top Variably Expressed Genes Heatmap')
  dev.off()
  
  select <- order(rowMeans(assay(rld)), decreasing=TRUE)[1:500]
  df <- data.frame(Group = colData(rld)[,c(group)], row.names = rownames(colData(dds)))
  pdf(file = paste0(plotPath,'Heatmap_TopExp_', Comparison, '.pdf'), onefile = FALSE, width=10, height=20)
  p <- pheatmap(assay(rld)[select,], scale="row",  cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=df, fontsize = 7, las = 2, fontsize_row = 7, color = colors, main = '500 Top Expressed Genes Heatmap')
  dev.off()
  
}

convert_human_symbol <- function() {
  
  res_table <- cbind(genes=row.names(res), res[ ,c(1:6)])
  res_table <- as.data.frame(res_table)
  
  biomartCacheClear()
  genes <- res_table$genes
  
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                            "hgnc_symbol", "description"),values=genes,mart= mart)
  res_table_merge <- merge(res_table,G_list,by.x="genes",by.y="ensembl_gene_id")
  
}

convert_mouse_symbol <- function() {
  
  res_table <- cbind(genes=row.names(res), res[ ,c(1:6)])
  res_table <- as.data.frame(res_table)
  
  biomartCacheClear()
  genes <- res_table$genes
  
  mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                            "mgi_symbol", "description"),values=genes,mart= mart)
  res_table_merge <- merge(res_table,G_list,by.x="genes",by.y="ensembl_gene_id")
  return(res_table_merge)
}

diff_output <- function() {
  res_table_merge <- arrange(res_table_merge, padj)
  write.csv(res_table_merge, 
            row.names = FALSE,
            na = "",
            file= paste0(plotPath2, Comparison, "_DESeq2.csv"))
}



