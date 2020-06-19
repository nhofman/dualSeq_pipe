set.seed(123)
# Parse arguments
args <- commandArgs(TRUE)
counttable_file <- args[match('--counttable', args) + 1]
condition_file <- args[match('--conditions', args) + 1]
comparisons_file <- args[match('--comparisons', args) + 1]
feature_counts_log_file <- args[match('--featcounts-log', args) + 1]
output_folder <- args[match('--output', args) + 1]
threads <- args[match('--threads', args) + 1]

# Required packages
for (package in c("DESeq2")) {
  if (!(package %in% rownames(installed.packages()))) {
    library("crayon")
    stop(paste('Package "', package, '" not installed', sep=""))
  } else {
    print(paste("Import:", package))
    library(package, character.only=TRUE)
  }
}

# "Optional" packages
for (package in c("BiocParallel", "pheatmap", "ggplot2", "reshape2", "gplots", "tidyr", "gtools", "clusterProfiler", "ReactomePA",
                  "dplyr", "UniProt.ws", "Glimma", "openxlsx", "pathview", "STRINGdb")) {
  if (!(package %in% rownames(installed.packages()))) {
    stop(paste('Package "', package, '" not installed', sep=""))
  } else {
    print(paste("Import:", package))
    library(package, character.only=TRUE)
  }
}

# Heatmap showing similarities between samples (needs a count table and conditions )
create_correlation_matrix <- function(countdata, conditiontable) {
  countdata.normalized.processed <- as.matrix(countdata)
  countdata.normalized.processed <- countdata.normalized.processed[rowSums(countdata.normalized.processed) >= 10,]
  countdata.normalized.processed <- log2(countdata.normalized.processed + 1)
  sample_cor <- cor(countdata.normalized.processed, method = 'pearson', use = 'pairwise.complete.obs')
  
  return(pheatmap(sample_cor, annotation_col = conditiontable, annotation_row = conditiontable))
}

# Bar charts showing the assignment of allignments to genes (featureCounts statistics)
create_feature_counts_statistics <- function(featureCountsLog) {
  d <- read.table(featureCountsLog, header = T, row.names = 1)
  colnames(d) <- gsub("mapping\\..*\\.(.*)\\.bam", "\\1", colnames(d))
  d <- d[,mixedorder(colnames(d))]
  #d <- d[,grep("CoV229E", colnames(d))]
  
  dpct <- t(t(d) / colSums(d))
  
  dm <- melt(t(d))
  dpctm <- melt(t(dpct))
  
  colnames(dm) <- c("Sample", "Group", "Reads")
  dm$Group <- factor(dm$Group, levels = rev(levels(dm$Group)[order(levels(dm$Group))]))
  
  assignment.absolute <- ggplot(dm[dm$Reads > 0,], aes(x = Sample, y = Reads)) +
    geom_bar(aes(fill = Group), stat = "identity", group = 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  colnames(dpctm) <- c("Sample", "Group", "Reads")
  dpctm$Group = factor(dpctm$Group, levels = rev(levels(dpctm$Group)[order(levels(dpctm$Group))]))
  assignment.relative <- ggplot(dpctm[dpctm$Reads > 0,], aes(x = Sample, y = Reads)) +
    geom_bar(aes(fill = Group), stat = "identity", group = 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+
  #theme(axis.title = element_text(size=18, face = "bold"), axis.text = element_text(size=15, face = "bold"), legend.text = element_text(size=12, face="bold"))
  
  return(list(assignment.absolute, assignment.relative))
}

# Heatmap showing log fold change of one condition versus all others (TODO legend and description)
plotHeatmap2 <- function(x, name = "no_name_set.pdf", row_subset = NA, distMethod = "euclidean", clusterMethod = "complete", clrn = NA){
  if (! is.na(clrn)) {
    # set the custom distance and clustering functions
    hclustfunc <- function(x) hclust(x, method = clusterMethod)
    distfunc <- function(x) dist(x, method = distMethod)
    # perform clustering on rows and columns
    cl.row <- hclustfunc(distfunc(x))
    gr.row <- cutree(cl.row, clrn)
    require(RColorBrewer)
    if (clrn < 3)clrn <- 3
    col1 <- brewer.pal(clrn, "Set1")
  }else {
    gr.row <- NA
  }
  
  pdf(name, width = 25, height = 25)
  nCol <- 40
  mycol2 <- colorpanel(n = nCol, low = "green", mid = "black", high = "red")
  mx <- max(abs(x))
  pairs.breaks <- seq(- mx, mx, by = (2 * mx / nCol))
  if (is.na(row_subset[1])) {
    xx <- x
  }else {
    xx <- x[row_subset,]
  }
  if (is.na(clrn)) {
    heatmap.2(xx,
              distfun = function(x) dist(x, method = distMethod),
              hclustfun = function(x) hclust(x, method = clusterMethod),
              col = mycol2,
              breaks = pairs.breaks,
              margins = c(20, 10),
              trace = "none")
  }else {
    heatmap.2(xx,
              distfun = function(x) dist(x, method = distMethod),
              hclustfun = function(x) hclust(x, method = clusterMethod),
              col = mycol2,
              breaks = pairs.breaks,
              margins = c(20, 10),
              RowSideColors = col1[gr.row],
              trace = "none")
  }
  dev.off()
  return(gr.row)
}

rotate_vector <- function(vec, n=1L){
  x <- seq(1, length(vec))
  while (n > 0) {
    x <- c(x[2 : length(x)], x[1])
    n <- n - 1
  }
  vec[x]
}

# Get Uniprot ID and Proteinname
get_uniprot <- function(gene, human.uniprot){
  gene_bitr <- bitr(gene, fromType="SYMBOL", toType="UNIPROT", "org.Hs.eg.db")
  View(gene_bitr)
  gene_uni <- select(human.uniprot, gene_bitr$UNIPROT, c("PROTEIN-NAMES","REVIEWED"), "UNIPROTKB")
  #gene_bitr <- merge(data.frame(SYMBOL=gene),gene_bitr, all.x = T)
  gene.df <- merge(gene_uni[gene_uni$REVIEWED=="reviewed",c(1,2)], gene_bitr, by.x = "UNIPROTKB", by.y = "UNIPROT")
  #gene.df <- rbind(gene.df, data.frame("UNIPROTKB"=NA, "PROTEIN-NAMES"=NA, "SYMBOL"=setdiff(gene, gene_bitr$SYMBOL), check.names = F))
  return(unique(gene.df))
}

# Over-representation analysis for up- and down-regulated genes
calc_ora <- function(gene, name = "", filename, subdir, GO = T, KEGG = T, REACTOME = F, ont = "BP", pvalue = 0.05){
  if(!dir.exists(paste0(output_folder, subdir))){
    dir.create(paste0(output_folder, subdir))
  }   
  gene_bitr <- bitr(gene, fromType="SYMBOL", toType="ENTREZID", "org.Hs.eg.db")
  if(GO){
    for(o in ont){
      ora_go <- try(enrichGO(gene_bitr$ENTREZID, keyType = "ENTREZID", OrgDb = "org.Hs.eg.db", ont = o, readable = T, pvalueCutoff = pvalue))
      if(nrow(data.frame(ora_go)) > 0){
        plot_go <- dotplot(ora_go, showCategory = 20) + ggtitle(name)
        ggsave(paste(filename,"_GO_",o,"_dotplot.svg", sep = ""), device = "svg", plot = plot_go, path = paste(output_folder, subdir, sep = ""), width = 18, height = 15)
        ora_go <- data.frame(ora_go)
        ora_go$SYMBOL <- sapply(ora_go$geneID, function(x){ x.split <- unlist(strsplit(x,"/")); return(paste(gene_bitr$SYMBOL[gene_bitr$ENTREZID %in% x.split], collapse = "/"))})
        write.table(ora_go, file = paste(output_folder, subdir, "/", filename, "_",o,".csv", sep = ""), sep = "\t", row.names = FALSE)
      }
    }
  }
  if(KEGG){
    ora_kegg <- try(enrichKEGG(gene_bitr$ENTREZID, organism = "hsa", use_internal_data = FALSE, pvalueCutoff = pvalue))
    if(nrow(data.frame(ora_kegg)) > 0){
      plot_kegg <- dotplot(ora_kegg, showCategory = 20) + ggtitle(name)
      ggsave(paste(filename,"_KEGG_dotplot.svg", sep = ""), device = "svg", plot = plot_kegg, path = paste(output_folder, subdir, sep = ""), width = 18, height = 15)
      ora_kegg <- data.frame(ora_kegg)
      ora_kegg$SYMBOL <- sapply(ora_kegg$geneID, function(x){ x.split <- unlist(strsplit(x,"/")); return(paste(gene_bitr$SYMBOL[gene_bitr$ENTREZID %in% x.split], collapse = "/"))})
      write.table(ora_kegg, file = paste(output_folder, subdir, "/", filename, "_KEGG.csv", sep = ""), sep = "\t", row.names = FALSE)
    }
  }
  if(REACTOME){
    ora_reactome <- try(enrichPathway(gene_bitr$ENTREZID, organism = "human", readable = T, pvalueCutoff = pvalue))
    if(nrow(data.frame(ora_reactome)) > 0){
      plot_reactome <- dotplot(ora_reactome, showCategory = 20) + ggtitle(name)
      ggsave(paste(filename,"_REACTOME_dotplot.svg", sep = ""), device = "svg", plot = plot_reactome, path = paste(output_folder, subdir, sep = ""), width = 18, height = 15)
      ora_reactome <- data.frame(ora_reactome)
      ora_reactome$SYMBOL <- sapply(ora_reactome$geneID, function(x){ x.split <- unlist(strsplit(x,"/")); return(paste(gene_bitr$SYMBOL[gene_bitr$ENTREZID %in% x.split], collapse = "/"))})
      write.table(ora_reactome, file = paste(output_folder, subdir, "/", filename, "_REACTOME.csv", sep = ""), sep = "\t", row.names = FALSE)
    }
  }
}

# Gene Set Enrichment Analysis
calc_gsea <- function(res, name, ont = "BP", sort.by = "stat", KEGG = T, GO = T, REACTOME = F, nPerm = 1000, p.cut = 0.05, subfol = "GSEA/"){
  if(!dir.exists(paste0(output_folder, subfol))){
    dir.create(paste0(output_folder, subfol), recursive = T)
  }  
  gsea.list <- list()
  gene_bitr <- bitr(res$SYMBOL, fromType="SYMBOL", toType="ENTREZID", "org.Hs.eg.db")
  res <- merge(res, gene_bitr, by = "SYMBOL")
  res <- res[order(res[,sort.by], decreasing = TRUE),]
  geneset_num <- as.numeric(res[,sort.by])
  names(geneset_num) <- res$ENTREZID
  if(KEGG){
    gsea_kegg <- try(gseKEGG(geneList = geneset_num, organism = "hsa", nPerm = nPerm, pvalueCutoff = p.cut))
    if(nrow(data.frame(gsea_kegg)) > 0){
      plot_kegg <- dotplot(gsea_kegg, showCategory = 20, split = ".sign") + facet_grid(.~.sign) + theme(axis.title = element_text(size=18, face = "bold"), axis.text = element_text(size=15, face = "bold"), strip.text.x = element_text(size = 18, face = "bold"))
      ggsave(paste(name,"_KEGG_dotplot.svg", sep = ""), device = "svg", plot = plot_kegg, path = paste(output_folder, subfol, sep = ""), width = 18, height = 15)
      gsea_kegg <- data.frame(gsea_kegg)
      gsea_kegg$core_enrichment_SYMBOL <- sapply(gsea_kegg$core_enrichment, function(x){ x.split <- unlist(strsplit(x,"/")); return(paste(bitr(x.split, "ENTREZID", "SYMBOL", "org.Hs.eg.db")[["SYMBOL"]], collapse = "/"))})
      write.table(gsea_kegg, file = paste(output_folder, subfol, "/", name, "_KEGG.csv", sep = ""), sep = ",", row.names = FALSE)
      gsea.list[["KEGG"]] <- data.frame(gsea_kegg)
    }
  }
  if(GO){
    for(o in ont){
      gsea_go <- try(gseGO(geneset_num, ont = o, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", nPerm = nPerm, pvalueCutoff = p.cut))
      if(nrow(data.frame(gsea_go)) > 0){
        plot_go <- dotplot(gsea_go, showCategory = 20, split = ".sign") + facet_grid(.~.sign) + theme(axis.title = element_text(size=18, face = "bold"), axis.text = element_text(size=15, face = "bold"), strip.text.x = element_text(size = 18, face = "bold"))
        ggsave(paste(name, "_GO_", o, "_dotplot.svg", sep = ""), device = "svg", plot = plot_go, path = paste(output_folder, subfol, sep = ""), width = 18, height = 15)
        gsea_go <- data.frame(gsea_go)
        gsea_go$core_enrichment_SYMBOL <- sapply(gsea_go$core_enrichment, function(x){ x.split <- unlist(strsplit(x,"/")); return(paste(bitr(x.split, "ENTREZID", "SYMBOL", "org.Hs.eg.db")[["SYMBOL"]], collapse = "/"))})
        write.table(gsea_go, file = paste(output_folder, subfol, "/", name, "_GO_", o, ".csv", sep = ""), sep = ",", row.names = FALSE)
        gsea.list[[o]] <- data.frame(gsea_go)
      }
    }
  }
  if(REACTOME){
    gsea_reactome <- try(gsePathway(geneset_num, nPerm = nPerm, pvalueCutoff = p.cut))
    if(nrow(data.frame(gsea_reactome)) > 0){
      plot_reactome <- dotplot(gsea_reactome, showCategory = 20, split = ".sign") + facet_grid(.~.sign) + theme(axis.title = element_text(size=18, face = "bold"), axis.text = element_text(size=15, face = "bold"), strip.text.x = element_text(size = 18, face = "bold"))
      ggsave(paste(name,"_REACTOME_dotplot.svg", sep = ""), device = "svg", plot = plot_reactome, path = paste(output_folder, subfol, sep = ""), width = 18, height = 15)
      gsea_reactome <- data.frame(gsea_reactome)
      gsea_reactome$core_enrichment_SYMBOL <- sapply(gsea_reactome$core_enrichment, function(x){ x.split <- unlist(strsplit(x,"/")); return(paste(bitr(x.split, "ENTREZID", "SYMBOL", "org.Hs.eg.db")[["SYMBOL"]], collapse = "/"))})
      write.table(gsea_reactome, file = paste(output_folder, subfol, "/", name, "_REACTOME.csv", sep = ""), sep = ",", row.names = FALSE)
      gsea.list[["REACTOME"]] <- data.frame(gsea_reactome)
    }
  }
  return(gsea.list)
}

loadKEGG <- function(geneID){
  query_names <- paste0("hsa:", geneID)
  query <- try(keggGet(query_names))
  query.unlist <- unlist(query, recursive=F)
}

getPATHWAY <- function(genes, kegg){
  pathway.list <- lapply(genes, function(x) try(kegg[[x]]$PATHWAY))
  pathway.list <- lapply(1:length(pathway.list), function(x){if(class(pathway.list[[x]])=="try-error"){pathway.list[[x]]<-NULL};return(pathway.list[[x]])})
  #rownames(pathway.list) <- genes
  return(pathway.list)
}

# Run on multiple threads
if ("BiocParallel" %in% rownames(installed.packages())) {
  register(MulticoreParam(threads))
}

# Import count table (featureCounts)
countdata.raw <- read.csv(counttable_file, header = TRUE, row.names = 1, sep = "\t", comment.char = "#")
countdata <- as.matrix(countdata.raw[, c(6 : length(countdata.raw))])
colnames(countdata) <- as.vector(sapply(colnames(countdata), function(x) gsub("mapping\\..*\\.(.*)\\.bam", "\\1", x)))
countdata <- countdata[,mixedorder(colnames(countdata))]
print(head(countdata))

# Import condition file
conditiontable <- read.csv(condition_file, header = FALSE, row.names = 1, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
colnames(conditiontable) <- c('condition')
conditiontable <- separate(conditiontable, condition, c("treatment", "time"), "_", FALSE)
conditiontable <- conditiontable[mixedorder(rownames(conditiontable)),]
conditiontable <- as.data.frame(row.names=colnames(countdata), lapply(conditiontable, as.factor))
#colnames(conditiontable) <- c('condition')
#conditiontable$treatment <- strsplit(conditiontable$condition, "_", fixed = TRUE)[0]
#conditiontable$time <- strsplit(conditiontable$condition, "_", fixed = TRUE)[1]
condition <- as.factor(conditiontable[, 1])
#conditiontable$virus <- sub("_.*","",conditiontable$condition) #remain everything before '_' -> virus name
#print(class(condition))
print(conditiontable)

# Import group comparison file
comparisons.df <- read.csv(comparisons_file, header = FALSE, sep = "\t", comment.char = "#", colClasses = c("character"))
#virus_names <- unique(sub("_.*", "", comparisons.df[,2]))
print(comparisons.df)

deseqDataset <- DESeqDataSetFromMatrix(countData = countdata, colData = conditiontable, design = ~ condition)

# Write normalized count table
deseqDataset <- estimateSizeFactors(deseqDataset)
countdata.normalized <- counts(deseqDataset, normalized = TRUE)
write.table(countdata.normalized, file = paste(output_folder, "counts_normalized.txt", sep = ""), sep = "\t", row.names = TRUE, col.names = NA)
for(virus in unique(conditiontable$treatment)){
  write.table(c(comparisons.df[grep(virus, comparisons.df[,2]),1], comparisons.df[grep(virus, comparisons.df[,2]),2]), paste("/nfs/sfb1021/SFB1021_Virus/conditions_", virus, ".tsv", sep = ""), row.names = F, col.names = F)
  countdata.normalized.virus <- countdata.normalized[,grep(paste(comparisons.df[grep(virus, comparisons.df[,2]),1], comparisons.df[grep(virus, comparisons.df[,2]),2], collapse = "|", sep = "|"), colnames(countdata.normalized))]
  countdata.normalized.virus <- as.data.frame(countdata.normalized.virus)
  countdata.normalized.virus <- cbind(Gene=rownames(countdata.normalized.virus), countdata.normalized.virus)
  write.table(countdata.normalized.virus, file = paste(output_folder, "counts_normalized_", virus, ".tsv", sep = ""), sep = "\t", row.names = F)
}
# plot raw countdata and normalized countdata 
pdf(paste(output_folder, "/counts.pdf", sep = ""), width = 25, height=15)
par(mfrow=c(1,2), mar=c(20, 4, 5, 2) + 0.1, lwd=2)
boxplot(countdata, outline = FALSE, las = 2, ylab = "raw reads", xlab = "")
boxplot(countdata.normalized, outline = FALSE, las = 2, ylab = "DESeq2 normalized reads", xlab = "")
dev.off()

deseq.results <- DESeq(object = deseqDataset, parallel = TRUE)

deseq.results.vst <- vst(deseq.results, blind = FALSE) # or vst()

#deseq.results[ rowSums(counts(deseq.results)) == 0, ] <- 1 # replace rows that have no reads with pseudocount

# Plot PCA
pca <- plotPCA(deseq.results.vst, intgroup = c("treatment", "time"), returnData = TRUE) 
#shape <- c(1:length(unique(pca$time)))
plot_PCA <- ggplot(pca, aes(PC1, PC2, color = treatment, shape = time)) + geom_point(size=3)
ggsave("PCA.svg", plot = plot_PCA, device = "svg", path = output_folder)

for (time in unique(conditiontable$time)) {
  pca_time <- plotPCA(deseq.results.vst[,grep(time, colnames(deseq.results.vst))], intgroup = c("treatment"), returnData = TRUE)
  plot_PCA <- ggplot(pca_time, aes(PC1, PC2, color = treatment)) + geom_point(size=3)
  ggsave(paste("PCA_", time, ".svg"), plot = plot_PCA, device = "svg", path = output_folder)
  print(time)
}

# Heatmap showing similarities between samples (needs a count table and conditions )
if ("pheatmap" %in% rownames(installed.packages())) {
  sample_cor <- cor(assay(deseq.results.vst), method = 'pearson', use = 'pairwise.complete.obs')
  pdf(paste(output_folder, 'correlation_heatmap.pdf', sep = ""), width = 15, height = 15, onefile = FALSE)
  #print(create_correlation_matrix(countdata.normalized, conditiontable))
  pheatmap(sample_cor, annotation_col = conditiontable[,-1], annotation_row = conditiontable[,-1], fontsize=5)
  dev.off()
}

for (time in unique(conditiontable$time)) {
  sample_cor <- cor(assay(deseq.results.vst[,grep(time, colnames(deseq.results.vst))]), method = 'pearson', use = 'pairwise.complete.obs')
  svg(paste(output_folder, 'correlation_heatmap_', time, '.svg', sep = ""), width = 10, height = 10, onefile = FALSE)
  pheatmap(sample_cor, annotation_col = conditiontable[grep(time, rownames(conditiontable)),], annotation_row = conditiontable[grep(time, rownames(conditiontable)),], fontsize=8)
  dev.off()
}


# Bar charts showing the assignment of allignments to genes (featureCounts statistics)
if (("ggplot2" %in% rownames(installed.packages())) && ("reshape2" %in% rownames(installed.packages()))) {
  pdf(paste(output_folder, 'counts_assignment.pdf', sep = ""), width = 20, height = 10)
  invisible(lapply(create_feature_counts_statistics(feature_counts_log_file), print))
  dev.off()
}

# Create all DESeq2 comparisons from comparison table
if (!file.exists(paste(output_folder, "deseq2_comparisons_shrunken", sep = ""))) {
  dir.create(paste(output_folder, "deseq2_comparisons_shrunken", sep = ""))
}
if (!file.exists(paste(output_folder, "plots", sep = ""))) {
  dir.create(paste(output_folder, "plots", sep = ""))
}


load(paste0(output_folder, "hsa_kegg.RData"))
human.uniprot <- UniProt.ws(taxId=9606)
genes.uniprot <- get_uniprot(rownames(countdata.normalized), human.uniprot)
save(human.uniprot, genes.uniprot, file = paste(output_folder,"uniprot.RData", sep = ""))
load(paste(output_folder,"uniprot.RData", sep = ""))


res.list <- list()
res.list.filter <- list()
lfc <- NULL
#lfc_filtered_reg <- list()
#lfc_filtered_not <- list()
lfc_names <- NULL
#for (n in 1:1) {
res.list.BPL <- mclapply(1:nrow(comparisons.df), function(n){
  print(comparisons.df[n,])
  res <- results(deseq.results, contrast = c("condition", comparisons.df[n,2], comparisons.df[n,1]), parallel = FALSE)
  res.list[[comparisons.df[n,2]]] <- res 
  resLFC <- lfcShrink(deseq.results, contrast = c("condition", comparisons.df[n,2], comparisons.df[n,1]), type = "ashr", res = res)
  #pdf(paste(output_folder, "plots/MA-plot/MAplot_", comparisons.df[n,2], "_Vs_", comparisons.df[n,1], ".pdf", sep = ""))
  #plotMA(res, ylim = c(-5, 5), main = paste(comparisons.df[n,2], comparisons.df[n,1], sep = "_Vs_"))
  #dev.off()
  pdf(paste(output_folder, "plots/MAplot_", comparisons.df[n,2], "_Vs_", comparisons.df[n,1], "_shrunk.pdf", sep = ""))
  plotMA(resLFC, ylim = c(-5, 5), main = paste(comparisons.df[n,2], comparisons.df[n,1], sep = "_Vs_"))
  dev.off()
  
  res <- as.data.frame(resLFC)
  res <- cbind(res, `-log10(padj)` = -log10(res$padj))
  
  #plot_volcano <- ggplot(res, aes(log2FoldChange, -log10(padj))) + geom_point() + 
  #  theme(axis.title = element_text(size=18, face = "bold")) +
  #  geom_hline(yintercept = -log10(0.05), color = "red")
  #ggsave(paste("Volcano_", comparisons.df[n,2], "_Vs_", comparisons.df[n,1], "_shrunk.pdf", sep = ""), plot = plot_volcano, device = "pdf", path = paste(output_folder, "plots/Volcano/", sep = ""))
  
  
  write.table(res, file = paste(output_folder, "deseq2_comparisons_shrunken/deseq2_results_", comparisons.df[n,2], "_Vs_", comparisons.df[n,1], "_forHTML.tsv", sep = ""), row.names = TRUE, col.names = NA, sep = "\t")
  #res <- cbind(res, countdata.normalized[,grep(paste(as.character(comparisons.df[n,]), collapse="|"),colnames(countdata.normalized))])
  res <- cbind(res, countdata.normalized[,grep(paste(as.character(comparisons.df[n,2]),sub("_",".*_",comparisons.df[n,1]), sep="|"),colnames(countdata.normalized))])
  res <- cbind(SYMBOL = rownames(res), res)
  genes.uniprot.filter <- genes.uniprot[genes.uniprot$SYMBOL %in% res$SYMBOL,] %>% group_by(SYMBOL) %>% summarise_at(c("UNIPROTKB","PROTEIN-NAMES"),paste, collapse = ";")
  res <- merge(genes.uniprot.filter, res, by = "SYMBOL", all.y = T)
  for(x in as.character(res$SYMBOL)){try(res[res$SYMBOL == x,"PATHWAY"] <- paste(query.unlist[[x]]$PATHWAY, collapse = ";"))} # add associated pathways for each gene
  colnames(res)[10:(ncol(res)-1)] <- paste0("normalized_", colnames(res)[10:(ncol(res)-1)])
  res_filter <- res[rowSums(res[,grep("normalized",colnames(res))])>0,] # remove genes with no read counts
  res_filter <- res_filter[order(res_filter$log2FoldChange, decreasing = TRUE),] # sort for LFC
  write.csv(res_filter, file = paste(output_folder, "deseq2_comparisons_shrunken/deseq2_results_", comparisons.df[n,2], "_Vs_", comparisons.df[n,1], ".csv", sep = ""), row.names = TRUE, col.names = NA)
  write.xlsx(data.frame(lapply(res_filter, function(x) if(is.numeric(x)){formatC(x, digits = 4, format = "g")}else{x})), file = paste(output_folder, "deseq2_comparisons_shrunken/deseq2_results_", comparisons.df[n,2], "_Vs_", comparisons.df[n,1], ".xlsx", sep = ""))
  #write.csv(data.frame(lapply(res_filter, function(x) if(is.numeric(x)){formatC(x, digits = 4, format = "g")}else{x})), file = paste(output_folder, "deseq2_comparisons_shrunken/deseq2_results_", comparisons.df[n,2], "_Vs_", comparisons.df[n,1], ".csv", sep = ""), row.names = TRUE, col.names = NA)
  #write.xlsx(res_filter, paste(output_folder, "deseq2_comparisons_shrunken/deseq2_results_", comparisons.df[n,2], "_Vs_", comparisons.df[n,1], ".xlsx", sep = ""))
  
  res[is.na(res$`-log10(padj)`),"-log10(padj)"] <- 0
  res[is.na(res$padj),"padj"] <- 1
  glXYPlot(x = res$log2FoldChange, y = res$`-log10(padj)`, counts = res[,c(10:(ncol(res)-1))], groups = sub("_[[:digit:]]$","",colnames(res[,c(10:(ncol(res)-1))])), 
           status = ifelse(res$log2FoldChange>1 & res$padj<0.05, 1, ifelse(res$log2FoldChange<(-1) & res$padj<0.05, -1, 0)), 
           xlab = "log2FoldChange", ylab = "-log10(padj)", anno = res, side.main = "SYMBOL",
           display.columns = c("SYMBOL", "UNIPROTKB", "PROTEIN.NAMES", "PATHWAY", "padj"), 
           html = paste("Volcano_", comparisons.df[n,2], "_Vs_", comparisons.df[n,1], sep = ""),
           folder = "glimma_plots", path = output_folder, launch = F)
  
  return(res)
}, mc.cores = threads)
names(res.list) <- comparisons.df[,2]

res.list.filter <- sapply(res.list, function(x){x[which(apply(x[,grep("normalized", colnames(x))],1,max) >= 10 & x$padj < 0.05 & abs(x$log2FoldChange) > 1),]})

lfc.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL"),lapply(names(res.list), function(x){x.df <- data.frame(res.list[[x]][,c("SYMBOL","log2FoldChange")]); colnames(x.df) <- c("SYMBOL",paste0(x,".log2FoldChange")); return(x.df)}))
rownames(lfc.df) <- lfc.df$SYMBOL
lfc.df <- lfc.df[,-1]
lfc.df <- lfc.df[,mixedorder(colnames(lfc.df))]

padj.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL"),lapply(names(res.list), function(x){x.df <- data.frame(res.list[[x]][,c("SYMBOL","padj")]); colnames(x.df) <- c("SYMBOL",paste0(x,".padj")); return(x.df)}))
rownames(padj.df) <- padj.df$SYMBOL
padj.df <- padj.df[,-1]
padj.df <- padj.df[,mixedorder(colnames(padj.df))]
write.xlsx(list(log2FoldChange=data.frame(lapply(lapply(lfc.df, formatC, digits = 4, format = "g"), as.numeric), row.names = row.names(lfc.df)), padj=data.frame(lapply(lapply(padj.df, formatC, digits = 4, format = "g"), as.numeric), row.names = row.names(lfc.df))), file = paste(output_folder,"deseq2_comparisons_shrunken/expression_data_all.xlsx", sep = ""), row.names = T)

lfc.df[is.na(lfc.df)] <- 0
save.image(paste(output_folder, "/deseq2.RData", sep = ""))

# Count differentially expressed genes
count.genes <- data.frame()
padj_cut <- 0.05
for (sample in names(res.list)) {
  res <- res.list[[sample]]
  #expr_genes <- res[which(res$log2FoldChange!=0 & !is.na(res$padj)),]
  count_padj <- length(which(res$padj<padj_cut))
  count_fc_3_up <- length(which(res$padj<padj_cut & res$log2FoldChange> 3))
  count_fc_3_down <- length(which(res$padj<padj_cut & res$log2FoldChange< -3)) 
  count_fc_2_up <- length(which(res$padj<padj_cut & res$log2FoldChange> 2))
  count_fc_2_down <- length(which(res$padj<padj_cut & res$log2FoldChange< -2))
  count_fc_1_up <- length(which(res$padj<padj_cut & res$log2FoldChange> 1))
  count_fc_1_down <- length(which(res$padj<padj_cut & res$log2FoldChange< -1))
  count.genes <- rbind(count.genes,c(sample,list(count_padj,count_fc_1_up,count_fc_1_down,count_fc_2_up,count_fc_2_down,count_fc_3_up,count_fc_3_down)), stringsAsFactors = F)
}
colnames(count.genes) <- c("Sample", paste("padj<",padj_cut,sep = ""), "LFC>1", "LFC<-1", "LFC>2", "LFC<-2", "LFC>3", "LFC<-3")
count.genes <- as.data.frame(count.genes)
count.genes <- count.genes[mixedorder(count.genes$Sample),]
write.table(count.genes, file = paste(output_folder, "deseq2_comparisons_shrunken/gene_count.csv", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE)


# create boxplot over LFC (for every virus separatly) -> show course of infection 
par(mfrow = c(3,4))
for (virus in unique(conditiontable$treatment)) {
  if (!grepl("Mock", virus)){
    print(virus)
    #svg(paste(output_folder, "Boxplot_", virus, "_Vs_Control_counts", ".svg", sep = ""))#, width = 1500, height = 1000)
    #par(mar=c(10,4,4,2))
    #boxplot(countdata.normalized[,grep(virus, colnames(countdata.normalized))], las = 2, outline = FALSE)
    #dev.off()
    lfc <- do.call(cbind,lapply(res.list[grep(virus, names(res.list))],function(res){res$log2FoldChange}))
    lfc <- lfc[,mixedorder(colnames(lfc))]
    svg(paste(output_folder, "Boxplot_", virus, "_Vs_Control", ".svg", sep = ""))#, width = 1500, height = 1000)
    par(mar=c(10,4,4,2))
    boxplot(lfc, las = 2)
    dev.off()
    #svg(paste(output_folder, "Boxplot_", virus, "_Vs_Control_LFC>3", ".svg", sep = ""))#, width = 1500, height = 1000)
    #par(mar=c(10,4,4,2))
    #boxplot(lfc_filtered_reg[grep(virus,names(lfc_filtered_reg))], las = 2)
    #dev.off()
    #png(paste(output_folder, "Boxplot_", virus, "_Vs_Control_LFC<1", ".png", sep = ""), width = 1500, height = 1000)
    #boxplot(lfc_filtered[grep(virus,names(lfc_filtered))])
    #dev.off()
  }
}


# over-representation/GSEA analysis
#mclapply(1:nrow(comparisons.df), 
#    function(n){
gsea.list <- list()
for(n in names(res.list)){
  print(n)
  res <- res.list[[n]]
  res <- res[rowSums(res[,grep("normalized",colnames(res))])>0,] # remove genes with no read counts
  #res <- cbind(SYMBOL = rownames(res), res)
  gsea.list[[n]] <- calc_gsea(res, n, sort.by = "log2FoldChange", REACTOME = T, ont = c("CC", "MF", "BP"), nPerm = 10000,
                              p.cut = 0.05, subfol = "GSEA")
}
#	calc_ora(res$SYMBOL[res$padj<0.05 & res$log2FoldChange>1], "", comparisons.df[n,2], "up/ORA/")
#	calc_ora(res$SYMBOL[res$padj<0.05 & res$log2FoldChange<(-1)], "", comparisons.df[n,2], "down/ORA)
#    }#,
#    mc.cores = threads
#)



# dir.create(paste(output_folder, "summary", sep = ""))
# sapply(1 : nrow(comparisons.df), function(control_i) {
#     control = comparisons.df[control_i,2]
#     vs_condition <- comparisons.df[control_i,1]
#     print(paste("Control:\t",control,"\nvs_condition:\t",vs_condition,sep=""))
# 
#    log2foldChange <- Reduce(
#        function(d1, d2){
#            merge(d1, d2, by = "gene_id", all = TRUE)
#        },
#        lapply(
#            lapply(vs_condition, function(x) c(control, x)),
#            function(y){
#                df <- as.data.frame(results(deseq.results, addMLE = FALSE, contrast = c("condition", y[1], y[2])));
#                colnames(df) <- paste(y[2], colnames(df), sep = ".");
#                df$gene_id <- rownames(df);
#                df
#            }
#        )
#    )
#    print("log2foldChange:")
#    print(head(log2foldChange))
# 
#    if (length(vs_condition) > 2) {
#        log2foldChangeOnly <- log2foldChange[, c(1,
#                                                 seq(from = 3, to = length(log2foldChange), by = 6),
#                                                 seq(from = 7, to = length(log2foldChange), by = 6),
#                                                 seq(from = 6, to = length(log2foldChange), by = 6))
#                                            ]
#        log2foldChangeOnly$minFC <- apply(log2foldChangeOnly[, 2:(length(vs_condition) + 1)], 1, min, na.rm = T)
#        log2foldChangeOnly$maxFC <- apply(log2foldChangeOnly[, 2:(length(vs_condition) + 1)], 1, max, na.rm = T)
#        fc2 <- log2foldChangeOnly[log2foldChangeOnly$minFC <= -2 | log2foldChangeOnly$maxFC >= 2,]
#        fc2$minFC <- NULL
#        fc2$maxFC <- NULL
#    } else {
#        log2foldChangeOnly <- log2foldChange[, c(7, 2, 5, 6)]
#        fc2 <- log2foldChangeOnly[log2foldChangeOnly[, 2] <= - 2 | log2foldChangeOnly[, 2] >= 2,]
#        fc2 <- fc2[complete.cases(fc2[, 2]),]
#    }
# 
#    rownames(fc2) <- fc2$gene_id
#    fc2 <- merge(fc2, counts(deseqDataset, normalized = TRUE), by = 0)
#    fc2$Row.names <- NULL
#    rownames(fc2) <- fc2$gene_id
#    print("fc2:")
#    print(head(fc2))
#    write.table(x = fc2[, c(1,
#                            rotate_vector(2:(1 + (length(vs_condition) * 3)), length(vs_condition)),
#                            (length(vs_condition) * 3 + 2):length(fc2))],
#                paste(output_folder, 'summary/', control, ".tsv", sep = ""),
#                row.names = F, col.names = T, sep = "\t")
#    if (length(vs_condition) > 2) {
#        m <- as.matrix(fc2[, 2 : (length(vs_condition) + 1)])
#        colnames(m) <- sapply(colnames(m), function(x) gsub("\\..{1,}$", '', x))
#        plotHeatmap2(m, name = paste(output_folder, 'summary/', control, ".pdf", sep = ""))
#    }
# })