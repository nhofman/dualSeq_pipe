set.seed(123)
# Parse arguments
args <- commandArgs(F)
file.dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
counttable_file <- args[match('--counttable', args) + 1]
condition_file <- args[match('--conditions', args) + 1]
comparisons_file <- args[match('--comparisons', args) + 1]
feature_counts_log_file <- args[match('--featcounts-log', args) + 1]
output_folder <- args[match('--output', args) + 1]
threads <- args[match('--threads', args) + 1]
color_file <- args[match('--color', args) + 1]
rRNA_file <- args[match('--rRNA', args) + 1]

# required packages
for(package in c("BiocParallel", "DESeq2", "pheatmap", "ggplot2", "reshape2", "gplots", "tidyr", "gtools", "dplyr", "openxlsx", "org.Hs.eg.db")) {
  if(!(package %in% rownames(installed.packages()))) {
    stop(paste('Package "', package, '" not installed', sep=""))
  } else {
    print(paste("Import:", package))
    library(package, character.only=TRUE)
  }
}

# Run on multiple threads
if("BiocParallel" %in% rownames(installed.packages())) {
  register(MulticoreParam(threads))
}

# Annotate genes
annotate <- function(genes, keytype){
  genes.ann <- AnnotationDbi::select(org.Hs.eg.db, genes,c("UNIPROT","GENENAME","PATH"), keytype)
  genes.ann <- aggregate(genes.ann, by = list(genes.ann$SYMBOL), FUN = function(x) paste(unique(x), collapse = ";"))[,-1]
  return(genes.ann)
}

# Bar charts showing the assignment of alignments to genes (featureCounts statistics)
create_feature_counts_statistics <- function(featureCountsLog) {
  d <- read.table(featureCountsLog, header = T, row.names = 1)
  colnames(d) <- gsub("mapping\\..*\\.(.*)\\.bam", "\\1", colnames(d))
  d <- d[,mixedorder(colnames(d))]
  
  dpct <- t(t(d) / colSums(d))
  
  dm <- melt(t(d))
  dpctm <- melt(t(dpct))
  
  colnames(dm) <- c("Sample", "Group", "Reads")
  dm$Group <- factor(dm$Group, levels = rev(levels(dm$Group)[order(levels(dm$Group))]))
  dm <- separate(dm, "Sample", c("Virus","Time"), "_", F, extra = "merge")
  
  assignment.absolute <- ggplot(dm[dm$Reads > 0,], aes(x = factor(Time, levels = mixedsort(unique(Time))), y = Reads)) +
    geom_bar(aes(fill = Group), stat = "identity", group = 1) +
    facet_wrap(~ Virus, scales = "free_x") + xlab("Time_Replicate") +
    theme(axis.title.x = element_text(size = 20, face = "bold", margin= margin(t=10,r=0,b=0,l=0)), 
          axis.title.y = element_text(size = 20, face = "bold", margin = margin(t=0,r=10,b=0,l=0)),
          axis.text = element_text(size = 16), axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text = element_text(size = 20, face = "bold"),
          legend.text = element_text(size = 16), legend.title = element_text(size = 20))
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  colnames(dpctm) <- c("Sample", "Group", "Reads")
  dpctm$Group = factor(dpctm$Group, levels = rev(levels(dpctm$Group)[order(levels(dpctm$Group))]))
  dpctm <- separate(dpctm, "Sample", c("Virus","Time"), "_", F, extra = "merge")
  assignment.relative <- ggplot(dpctm[dpctm$Reads > 0,], aes(x = factor(Time, levels = mixedsort(unique(Time))), y = Reads)) +
    geom_bar(aes(fill = Group), stat = "identity", group = 1) +
    facet_wrap(~ Virus, scales = "free_x") + xlab("Time_Replicate") +
    theme(axis.title.x = element_text(size = 20, face = "bold", margin= margin(t=10,r=0,b=0,l=0)), 
          axis.title.y = element_text(size = 20, face = "bold", margin = margin(t=0,r=10,b=0,l=0)),
          axis.text = element_text(size = 16), axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text = element_text(size = 20, face = "bold"),
          legend.text = element_text(size = 16), legend.title = element_text(size = 20))
  
  return(list("absolute"=assignment.absolute, "relative"=assignment.relative))
}


# Import count table (featureCounts)
countdata.raw <- read.csv(counttable_file, header = TRUE, row.names = 1, sep = "\t", comment.char = "#")
countdata <- as.matrix(countdata.raw[, c(6 : length(countdata.raw))])
colnames(countdata) <- as.vector(sapply(colnames(countdata), function(x) gsub("mapping\\..*\\.(.*)\\.bam", "\\1", x)))
countdata <- countdata[,mixedorder(colnames(countdata))]
# Import rRNA genes and filter them from countdata
if(rRNA_file != "NULL"){
  rRNA <- read.table(rRNA_file, header = FALSE, stringsAsFactors = FALSE)
  rRNA <- rRNA$V1
  countdata <- countdata[!rownames(countdata) %in% rRNA,] #countdata[-grep("RNA[0-9-]+S", rownames(countdata)),]
}

# Import condition file
conditiontable <- read.csv(condition_file, header = FALSE, row.names = 1, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
colnames(conditiontable) <- c('condition')
conditiontable <- separate(conditiontable, condition, c("Treatment", "Time"), "_", FALSE)
conditiontable <- conditiontable[mixedorder(rownames(conditiontable)),]
conditiontable <- as.data.frame(row.names=colnames(countdata), lapply(conditiontable, as.factor))
condition <- as.factor(conditiontable[, 1])

# Import group comparison file
comparisons.df <- read.csv(comparisons_file, header = FALSE, sep = "\t", comment.char = "#", colClasses = c("character"))
# Create DESeqDataSet
deseqDataset <- DESeqDataSetFromMatrix(countData = countdata, colData = conditiontable, design = ~ condition)

# Write normalized count table
deseqDataset <- estimateSizeFactors(deseqDataset)
countdata.normalized <- counts(deseqDataset, normalized = TRUE)
write.table(countdata.normalized, file = paste(output_folder, "/counts_normalized.tsv", sep = ""), sep = "\t", row.names = TRUE, col.names = NA)

# Import color
if(color_file != "NULL"){
  color.df <- read.table(color_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  color <- color.df[,2]
  names(color) <- color.df[,1]
  
  # Set color for column names of countdata
  colnames.col <- c()
  for(x in names(color)){
    i <- grep(x, colnames(countdata.normalized))
    colnames.col[i] <- color[x]
  }
}else{
  colnames.col <- "lightgrey"
}

# Plot raw countdata and normalized countdata 
pdf(paste(output_folder, "/counts.pdf", sep = ""), width = 30, height = 15)
par(mfrow=c(1,2), mar=c(12, 5.5, 5, 2) + 0.1, lwd=2, cex.lab=1.75)
boxplot(countdata, outline = FALSE, las = 2, ylab = "Raw reads", xlab = "", col = colnames.col)
boxplot(countdata.normalized, outline = FALSE, las = 2, ylab = "DESeq2 normalized reads", xlab = "", col = colnames.col)
dev.off()

png(paste(output_folder, "/counts.png", sep = ""), width = 30, height = 15, units = "in", res = 300)
par(mfrow=c(1,2), mar=c(8, 5.5, 5, 2) + 0.1, lwd=2, cex.lab=1.75)
boxplot(countdata, outline = FALSE, las = 2, ylab = "Raw reads", xlab = "", col = colnames.col)
boxplot(countdata.normalized, outline = FALSE, las = 2, ylab = "DESeq2 normalized reads", xlab = "", col = colnames.col)
dev.off()

# Calculate differential expression analysis steps
deseq.results <- DESeq(object = deseqDataset, parallel = FALSE)

# Calculate variance-stabilized read counts
deseq.results.vst <- vst(deseq.results, blind = FALSE) # or vst()

#deseq.results[ rowSums(counts(deseq.results)) == 0, ] <- 1 # replace rows that have no reads with pseudocount

# Plot PCA
shape <- if(length(unique(conditiontable$Time)) <= 6){ scales::shape_pal()(length(unique(conditiontable$Time))) }else{ c(1:length(unique(conditiontable$Time)))}
names(shape) <- unique(conditiontable$Time)
pca <- plotPCA(deseq.results.vst, intgroup = c("Treatment", "Time"), returnData = TRUE)
plot_PCA <- ggplot(pca, aes(PC1, PC2, color = Treatment, shape = factor(Time, levels = mixedsort(as.character(unique(conditiontable$Time)))))) + 
  geom_point(size=3) + labs(color = "Infection", shape = "Time") + scale_shape_manual(values=shape) + guides(color=guide_legend(override.aes=list(fill=NA))) +
  theme(axis.title = element_text(size = 20, face = "bold"), axis.text = element_text(size = 16, face = "bold"),  
        legend.text = element_text(size = 16), legend.title = element_text(size = 20, face = "bold"), legend.key=element_blank())
if(exists("color")){
  plot_PCA <- plot_PCA + scale_colour_manual(values=color)
}
ggsave("PCA.pdf", plot = plot_PCA, device = "pdf", path = output_folder, width = 10, height = 6)
ggsave("PCA.png", plot = plot_PCA, device = "png", path = output_folder, width = 10, height = 6)

for(time in unique(conditiontable$Time)) {
  pca_time <- plotPCA(deseq.results.vst[,grep(time, colnames(deseq.results.vst))], intgroup = c("Treatment", "Time"), returnData = TRUE)
  plot_PCA <- ggplot(pca_time, aes(PC1, PC2, color = Treatment, shape = Time)) + geom_point(size=5) + 
    labs(color = "Infection", shape = "Time") + ggtitle(main) + 
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), axis.title = element_text(size = 20, face = "bold"), axis.text = element_text(size = 16),  
          legend.text = element_text(size = 16), legend.title = element_text(size = 20, face = "bold"), legend.key=element_blank()) +
    scale_shape_manual(values=shape) + guides(color = guide_legend(order = 2), shape = F) 
  if(exists("color")){
    plot_PCA <- plot_PCA + scale_colour_manual(values=color)
  }
  ggsave(paste0("PCA_", time, ".pdf"), plot = plot_PCA, device = "pdf", path = output_folder)
  ggsave(paste0("PCA_", time, ".png"), plot = plot_PCA, device = "png", path = output_folder)
}

# Heatmap showing correlations between samples
if("pheatmap" %in% rownames(installed.packages())) {
  sample_cor <- cor(assay(deseq.results.vst), method = 'pearson', use = 'pairwise.complete.obs')
  if(exists("color")){
    annColor <- list(Treatment = color)
  }else(
    annColor <- NA
  )
  plot.heat <- pheatmap(sample_cor, annotation_col = conditiontable[,-1], annotation_row = conditiontable[,-1], fontsize=8, annotation_colors = annColor)
  pdf(paste(output_folder, 'correlation_heatmap.pdf', sep = ""), width = 15, height = 15, onefile = FALSE)
  print(plot.heat)
  dev.off()
  png(paste(output_folder, 'correlation_heatmap.png', sep = ""), width = 1300, height = 1300)
  print(plot.heat)
  dev.off()

  for(time in unique(conditiontable$Time)) {
    sample_cor <- cor(assay(deseq.results.vst[,grep(time, colnames(deseq.results.vst))]), method = 'pearson', use = 'pairwise.complete.obs')
    plot.heat <- pheatmap(sample_cor, annotation_col = conditiontable[grep(time, rownames(conditiontable)), -1], annotation_row = conditiontable[grep(time, rownames(conditiontable)), -1], 
                          fontsize=12, annotation_colors = annColor, fontfamily = "")
    pdf(paste(output_folder, 'correlation_heatmap_', time, '.pdf', sep = ""), width = 10, height = 10, onefile = FALSE)
    print(plot.heat)
    dev.off()
    png(paste(output_folder, 'correlation_heatmap_', time, '.png', sep = ""), width = 600, height = 600)
    print(plot.heat)
    dev.off()
  }
}

# Bar charts showing the assignment of alignments to genes (featureCounts statistics)
if(("ggplot2" %in% rownames(installed.packages())) && ("reshape2" %in% rownames(installed.packages()))) {
  plot.list <- create_feature_counts_statistics(feature_counts_log_file)
  invisible(sapply(list("pdf", "png"), function(p){
    sapply(names(plot.list), function(x){
      ggsave(plot=plot.list[[x]],filename=paste0("counts_assignment_", x, ".", p), device=p, path=output_folder, width=18, height=12)
    })}))
}

# Create all DESeq2 comparisons from comparison table
if(!dir.exists(paste(output_folder, "deseq2_comparisons_shrunken", sep = ""))) {
  dir.create(paste(output_folder, "deseq2_comparisons_shrunken", sep = ""))
}
if(!dir.exists(paste(output_folder, "plots", sep = ""))) {
  dir.create(paste(output_folder, "plots", sep = ""))
}

res.list.raw <- list()
res.list <- list()

for(n in 1:nrow(comparisons.df)) {
  #res.list <- mclapply(1:nrow(comparisons.df), function(n){
  print(comparisons.df[n,])
  group_1 <- comparisons.df[n,1] # e.g. infected
  group_2 <- comparisons.df[n,2] # e.g. uninfected
  # calculate log2 fold changes for condition group_1 vs. group_2
  res <- results(deseq.results, contrast = c("condition", group_1, group_2), parallel = FALSE)
  res.list.raw[[paste(group_1, group_2, sep="_vs_")]] <- res 
  write.table(as.data.frame(res), file = paste(output_folder, "deseq2_comparisons_shrunken/deseq2_results_", group_1, "_vs_", group_2, "_unshrunken.tsv", sep = ""), row.names = TRUE, col.names = NA, sep = "\t")
  write.xlsx(as.data.frame(res), file = paste(output_folder, "deseq2_comparisons_shrunken/deseq2_results_", group_1, "_vs_", group_2, "_unshrunken.xlsx", sep = ""), row.names = TRUE)
  
  resLFC <- lfcShrink(deseq.results, contrast = c("condition", group_1, group_2), type = "ashr", res = res)
  
  pdf(paste(output_folder, "plots/MAplot_", group_1, "_vs_", group_2, "_shrunk.pdf", sep = ""))
  plotMA(resLFC, ylim = c(-5, 5), main = paste(group_1, group_2, sep = "_vs_"))
  dev.off()
  
  png(paste(output_folder, "plots/MAplot_", group_1, "_vs_", group_2, "_shrunk.png", sep = ""))
  plotMA(resLFC, ylim = c(-5, 5), main = paste(group_1, group_2, sep = "_vs_"))
  dev.off()
  
  res <- as.data.frame(resLFC)
  res <- cbind(res, `-log10(padj)` = -log10(res$padj))
  
  plot_volcano <- ggplot(res, aes(log2FoldChange, -log10(padj))) + geom_point() + 
    theme(axis.title = element_text(size=18, face = "bold")) +
    geom_hline(yintercept = -log10(0.05), color = "red")
  ggsave(paste("Volcano_", group_1, "_vs_", group_2, "_shrunk.pdf", sep = ""), plot = plot_volcano, device = "pdf", path = paste(output_folder, "plots/", sep = ""))
  ggsave(paste("Volcano_", group_1, "_vs_", group_2, "_shrunk.png", sep = ""), plot = plot_volcano, device = "png", path = paste(output_folder, "plots/", sep = ""))
  
  #res <- cbind(res, countdata.normalized[,grep(paste(as.character(comparisons.df[n,]), collapse="|"),colnames(countdata.normalized))])
  res <- cbind(res, countdata.normalized[,grep(paste(as.character(group_1),sub("_",".*_",group_2), sep="|"),colnames(countdata.normalized))])
  res <- cbind(SYMBOL = rownames(res), res)
  res <- merge(annotate(as.character(res$SYMBOL), "SYMBOL"), res, by = "SYMBOL", all = T)
  colnames(res)[11:(ncol(res))] <- paste0("normalized_", colnames(res)[11:(ncol(res))])
  res_filter <- res[rowSums(res[,grep("normalized",colnames(res))])>0,] # remove genes with no read counts
  res_filter <- res_filter[order(res_filter$log2FoldChange, decreasing = TRUE),] # sort for LFC
  write.table(res_filter, file = paste(output_folder, "deseq2_comparisons_shrunken/deseq2_results_", group_1, "_vs_", group_2, ".tsv", sep = ""), row.names = TRUE, col.names = NA, sep = "\t")
  write.xlsx(res_filter, paste(output_folder, "deseq2_comparisons_shrunken/deseq2_results_", group_1, "_vs_", group_2, ".xlsx", sep = ""), overwrite = T)
  
  res.list[[paste(group_1, group_2, sep="_vs_")]] <- res 
  #return(res)
}#, mc.cores = threads)

#names(res.list.shrunk) <- comparisons.df[,2]

lfc.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL"),lapply(names(res.list), function(x){x.df <- data.frame(res.list[[x]][,c("SYMBOL","log2FoldChange")]); colnames(x.df) <- c("SYMBOL",x); return(x.df)}))
rownames(lfc.df) <- lfc.df$SYMBOL
lfc.df <- lfc.df[,-1, drop = FALSE]
lfc.df <- lfc.df[,mixedorder(colnames(lfc.df))]

padj.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL"),lapply(names(res.list), function(x){x.df <- data.frame(res.list[[x]][,c("SYMBOL","padj")]); colnames(x.df) <- c("SYMBOL",x); return(x.df)}))
rownames(padj.df) <- padj.df$SYMBOL
padj.df <- padj.df[,-1, drop = FALSE]
padj.df <- padj.df[,mixedorder(colnames(padj.df))]
write.xlsx(list(log2FoldChange=lfc.df, padj=padj.df), file = paste(output_folder,"deseq2_comparisons_shrunken/expression_data_all.xlsx", sep = ""), row.names = T, overwrite = T)
write.table(lfc.df, file = paste(output_folder,"deseq2_comparisons_shrunken/LFC_data_all.xlsx", sep = ""), row.names = T, sep = "\t")
write.table(padj.df, file = paste(output_folder,"deseq2_comparisons_shrunken/padj_data_all.xlsx", sep = ""), row.names = T, sep = "\t")

lfc.df[is.na(lfc.df)] <- 0

# Count differentially expressed genes
padj_cut <- 0.05
count.genes <- data.frame()
for(sample in names(res.list)) {
  res <- res.list[[sample]]
  #expr_genes <- res[which(res$log2FoldChange!=0 & !is.na(res$padj)),]
  count_padj <- length(which(res$padj<padj_cut))
  for(LFC.cut in c(0,1,1.5,2,3)){
    count_fc_up <- length(which(res$padj<padj_cut & res$log2FoldChange> LFC.cut & apply(res[,grep("normalized", colnames(res))],1,max) >= 10))
    count_fc_down <- length(which(res$padj<padj_cut & res$log2FoldChange< -LFC.cut & apply(res[,grep("normalized", colnames(res))],1,max) >= 10))
    count.genes <- rbind(count.genes,c(sample,count_fc_up,"up",LFC.cut), stringsAsFactors = F)
    count.genes <- rbind(count.genes,c(sample,-count_fc_down,"down",LFC.cut), stringsAsFactors = F)
  }
  
}
colnames(count.genes) <- c("Sample", "Count", "Direction", "LFC_cutoff")
count.genes <- as.data.frame(count.genes)
write.table(count.genes, file = paste(output_folder, "deseq2_comparisons_shrunken/gene_count.csv", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE)
write.xlsx(count.genes, file = paste(output_folder, "deseq2_comparisons_shrunken/gene_count.xlsx", sep = ""), row.names = FALSE, col.names = TRUE)

count.genes <- separate(separate(count.genes, "Sample", c("Virus","Mock"), "_[^_]*_vs_", F), "Mock", c("Control","Time"), "_")
#count.genes <- separate(count.genes, "Sample", c("Virus","Control"), "_[^_]*_vs_", F)
count.genes$Count <- as.numeric(count.genes$Count) 
count.genes <- count.genes[mixedorder(count.genes$Time),]
for(LFC.cut in unique(count.genes$LFC_cutoff)){
  p <- ggplot(count.genes[count.genes$LFC_cutoff==LFC.cut,], aes(y=Count, x=factor(paste(Control,Time,sep="_"), levels = unique(paste(Control,Time,sep="_"))), group=Direction, fill=Direction)) + 
    geom_bar(stat = "identity", position = "stack", color = "black") + facet_wrap(~Virus, scales = "free_x") + xlab("Time") + 
    scale_y_continuous(breaks = pretty(count.genes$Count[count.genes$LFC_cutoff==LFC.cut], n=10), labels = abs(pretty(count.genes$Count[count.genes$LFC_cutoff==LFC.cut], n=10))) + 
    geom_hline(yintercept = 0) + scale_fill_manual(values=c(up="red", down="green"))
  ggsave(paste0("DEG_count_LFC",LFC.cut,".pdf"), p, "svg", output_folder)
  ggsave(paste0("DEG_count_LFC",LFC.cut,".png"), p, "png", output_folder, width = 14, height = 7)
}

save.image(paste(output_folder, "/deseq2.RData", sep = ""))
