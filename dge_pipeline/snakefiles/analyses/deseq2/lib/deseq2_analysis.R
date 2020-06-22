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
                  "dplyr", "UniProt.ws", "htmlwidgets", "Glimma", "openxlsx", "pathview")) {
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

# Get Uniprot ID and Proteinname
get_uniprot <- function(gene, human.uniprot){
  gene_bitr <- bitr(gene, fromType="SYMBOL", toType="UNIPROT", "org.Hs.eg.db")
  gene_uni <- select(human.uniprot, gene_bitr$UNIPROT, c("PROTEIN-NAMES","REVIEWED"), "UNIPROTKB")
  #gene_bitr <- merge(data.frame(SYMBOL=gene),gene_bitr, all.x = T)
  gene.df <- merge(gene_uni[gene_uni$REVIEWED=="reviewed",c(1,2)], gene_bitr, by.x = "UNIPROTKB", by.y = "UNIPROT")
  #gene.df <- rbind(gene.df, data.frame("UNIPROTKB"=NA, "PROTEIN-NAMES"=NA, "SYMBOL"=setdiff(gene, gene_bitr$SYMBOL), check.names = F))
  return(unique(gene.df))
}

# Over-representation analysis for a set of genes
calc_ora <- function(gene, main = "", filename, subdir, GO = T, KEGG = T, REACTOME = F, ont = "BP", p.cut = 0.05){
    gene_bitr <- bitr(gene, fromType="SYMBOL", toType="ENTREZID", "org.Hs.eg.db")
    ora.list <- list()
    if(GO){
      for(o in ont){
        ora_go <- try(enrichGO(gene_bitr$ENTREZID, keyType = "ENTREZID", OrgDb = "org.Hs.eg.db", ont = o, readable = T, pvalueCutoff = p.cut))
        if(nrow(data.frame(ora_go)) > 0){
        	plot_go <- dotplot(ora_go, showCategory = 20) + ggtitle(main)
        	ggsave(paste(filename,"_",o,"_dotplot.svg", sep = ""), device = "svg", plot = plot_go, path = paste(output_folder, subdir, sep = ""), width = 18, height = 15)
        	write.table(data.frame(ora_go), file = paste(output_folder, subdir, "/", filename, "_",o,".csv", sep = ""), sep = "\t", row.names = FALSE)
        }
        ora.list[[o]] <- ora_go
      }
    }
    if(KEGG){
      ora_kegg <- try(enrichKEGG(gene_bitr$ENTREZID, organism = "hsa", use_internal_data = FALSE, pvalueCutoff = p.cut))
      if(nrow(data.frame(ora_kegg)) > 0){
        ora_kegg <- setReadable(ora_kegg, org.Hs.eg.db, keyType = "ENTREZID")
  	    plot_kegg <- dotplot(ora_kegg, showCategory = 20) + ggtitle(main)
  	    ggsave(paste(filename,"_KEGG_dotplot.svg", sep = ""), device = "svg", plot = plot_kegg, path = paste(output_folder, subdir, sep = ""), width = 18, height = 15)
  	    write.table(data.frame(ora_kegg), file = paste(output_folder, subdir, "/", filename, "_KEGG.csv", sep = ""), sep = "\t", row.names = FALSE)
      }
      ora.list[["KEGG"]] <- ora_kegg
    }
    if(REACTOME){
      ora_reactome <- try(enrichPathway(gene_bitr$ENTREZID, organism = "human", readable = T, pvalueCutoff = p.cut))
      if(nrow(data.frame(ora_reactome)) > 0){
        plot_reactome <- dotplot(ora_reactome, showCategory = 20) + ggtitle(main)
        ggsave(paste(filename,"_REACTOME_dotplot.svg", sep = ""), device = "svg", plot = plot_reactome, path = paste(output_folder, subdir, sep = ""), width = 18, height = 15)
        write.table(data.frame(ora_reactome), file = paste(output_folder, subdir, "/", filename, "_REACTOME.csv", sep = ""), sep = "\t", row.names = FALSE)
      }
      ora.list[["REACTOME"]] <- ora_reactome
    }
    return(ora.list)
}

# Gene Set Enrichment Analysis
calc_gsea <- function(res, name, ont = "BP", sort.by = "stat", KEGG = T, GO = T, REACTOME = F, nPerm = 1000, p.cut = 0.05, subfol = "GSEA/"){
    if(!dir.exists(paste0(output_folder, subfol))){
      dir.create(paste0(output_folder, subfol))
    }
    gsea.list <- list()
    gene_bitr <- bitr(res$SYMBOL, fromType="SYMBOL", toType="ENTREZID", "org.Hs.eg.db")
    res <- merge(res, gene_bitr, by = "SYMBOL")
    res <- res[order(res[,sort.by], decreasing = TRUE),]
    geneset_num <- as.numeric(res[,sort.by])
    names(geneset_num) <- res$ENTREZID
    if(KEGG){
      gsea_kegg <- try(gseKEGG(geneList = geneset_num, organism = "hsa", nPerm = nPerm, pvalueCutoff = p.cut, seed = T))
      if(nrow(data.frame(gsea_kegg)) > 0){
      	plot_kegg <- dotplot(gsea_kegg, showCategory = 20, split = ".sign") + facet_grid(.~.sign) + theme(axis.title = element_text(size=18, face = "bold"), axis.text = element_text(size=15, face = "bold"), strip.text.x = element_text(size = 18, face = "bold"))
      	ggsave(paste(name,"_KEGG_dotplot.svg", sep = ""), device = "svg", plot = plot_kegg, path = paste(output_folder, subfol, sep = ""), width = 18, height = 15)
      	gsea_kegg.df <- data.frame(gsea_kegg)
      	gsea_kegg.df$core_enrichment_SYMBOL <- sapply(gsea_kegg.df$core_enrichment, function(x){ x.split <- unlist(strsplit(x,"/")); return(paste(bitr(x.split, "ENTREZID", "SYMBOL", "org.Hs.eg.db")[["SYMBOL"]], collapse = "/"))})
      	write.table(gsea_kegg.df, file = paste(output_folder, subfol, "/", name, "_KEGG.csv", sep = ""), sep = ",", row.names = FALSE)
      	gsea.list[["KEGG"]] <- gsea_kegg
      }
    }
    if(GO){
      for(o in ont){
        gsea_go <- try(gseGO(geneset_num, ont = o, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", nPerm = nPerm, pvalueCutoff = p.cut, seed = T))
        if(nrow(data.frame(gsea_go)) > 0){
          plot_go <- dotplot(gsea_go, showCategory = 20, split = ".sign") + facet_grid(.~.sign) + theme(axis.title = element_text(size=18, face = "bold"), axis.text = element_text(size=15, face = "bold"), strip.text.x = element_text(size = 18, face = "bold"))
          ggsave(paste(name, "_GO_", o, "_dotplot.svg", sep = ""), device = "svg", plot = plot_go, path = paste(output_folder, subfol, sep = ""), width = 18, height = 15)
          gsea_go.df <- data.frame(gsea_go)
          gsea_go.df$core_enrichment_SYMBOL <- sapply(gsea_go.df$core_enrichment, function(x){ x.split <- unlist(strsplit(x,"/")); return(paste(bitr(x.split, "ENTREZID", "SYMBOL", "org.Hs.eg.db")[["SYMBOL"]], collapse = "/"))})
          write.table(gsea_go.df, file = paste(output_folder, subfol, "/", name, "_GO_", o, ".csv", sep = ""), sep = ",", row.names = FALSE)
          gsea.list[[o]] <- gsea_go
          # Problem: Bei gleichem value von by werden alle Pathways mit diesem Wert ausgegeben!
          # Loesung: by = pvalue?
          gsea.list[[paste0(o,"_simplify")]] <- simplify(gsea_go, cutoff = 0.5, by = "pvalue")
        }
      }
    }
    if(REACTOME){
      gsea_reactome <- try(gsePathway(geneset_num, nPerm = nPerm, pvalueCutoff = p.cut, seed = T))
      if(nrow(data.frame(gsea_reactome)) > 0){
        plot_reactome <- dotplot(gsea_reactome, showCategory = 20, split = ".sign") + facet_grid(.~.sign) + theme(axis.title = element_text(size=18, face = "bold"), axis.text = element_text(size=15, face = "bold"), strip.text.x = element_text(size = 18, face = "bold"))
        ggsave(paste(name,"_REACTOME_dotplot.svg", sep = ""), device = "svg", plot = plot_reactome, path = paste(output_folder, subfol, sep = ""), width = 18, height = 15)
        gsea_reactome.df <- data.frame(gsea_reactome)
        gsea_reactome.df$core_enrichment_SYMBOL <- sapply(gsea_reactome.df$core_enrichment, function(x){ x.split <- unlist(strsplit(x,"/")); return(paste(bitr(x.split, "ENTREZID", "SYMBOL", "org.Hs.eg.db")[["SYMBOL"]], collapse = "/"))})
        write.table(gsea_reactome.df, file = paste(output_folder, subfol, "/", name, "_REACTOME.csv", sep = ""), sep = ",", row.names = FALSE)
        gsea.list[["REACTOME"]] <- gsea_reactome
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
condition <- as.factor(conditiontable[, 1])

# Import group comparison file
comparisons.df <- read.csv(comparisons_file, header = FALSE, sep = "\t", comment.char = "#", colClasses = c("character"))

# Create DESeqDataSet
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

# Plot raw countdata and normalized countdata 
pdf(paste(output_folder, "/counts.pdf", sep = ""), width = 25, height=15)
par(mfrow=c(1,2), mar=c(20, 4, 5, 2) + 0.1, lwd=2)
boxplot(countdata, outline = FALSE, las = 2, ylab = "raw reads", xlab = "")
boxplot(countdata.normalized, outline = FALSE, las = 2, ylab = "DESeq2 normalized reads", xlab = "")
dev.off()

# Calculate differential expression analysis steps
deseq.results <- DESeq(object = deseqDataset, parallel = TRUE)

# Calculate variance-stabilized read counts
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

# Heatmap showing correlations between samples
if ("pheatmap" %in% rownames(installed.packages())) {
    sample_cor <- cor(assay(deseq.results.vst), method = 'pearson', use = 'pairwise.complete.obs')
    pdf(paste(output_folder, 'correlation_heatmap.pdf', sep = ""), width = 15, height = 15, onefile = FALSE)
    pheatmap(sample_cor, annotation_col = conditiontable[,-1], annotation_row = conditiontable[,-1], fontsize=5)
    dev.off()
}

for (time in unique(conditiontable$time)) {
    sample_cor <- cor(assay(deseq.results.vst[,grep(time, colnames(deseq.results.vst))]), method = 'pearson', use = 'pairwise.complete.obs')
    svg(paste(output_folder, 'correlation_heatmap_', time, '.svg', sep = ""), width = 10, height = 10, onefile = FALSE)
    pheatmap(sample_cor, annotation_col = conditiontable[grep(time, rownames(conditiontable)),], annotation_row = conditiontable[grep(time, rownames(conditiontable)),], fontsize=8)
    dev.off()
}


# Bar charts showing the assignment of a  lignments to genes (featureCounts statistics)
if (("ggplot2" %in% rownames(installed.packages())) && ("reshape2" %in% rownames(installed.packages()))) {
    pdf(paste(output_folder, 'counts_assignment.pdf', sep = ""), width = 20, height = 10)
    invisible(lapply(create_feature_counts_statistics(feature_counts_log_file), print))
    dev.off()
}

# Proteinname and Uniprot ID for each gene
human.uniprot <- UniProt.ws(taxId=9606)
genes.uniprot <- get_uniprot(rownames(countdata.normalized), human.uniprot)
save(human.uniprot, genes.uniprot, file = paste(output_folder,"uniprot.RData", sep = ""))
load(paste(output_folder,"uniprot.RData", sep = ""))

# Gene associations with KEGG pathways
load(paste0(output_folder, "hsa_kegg.RData"))
# Gene associations with GO terms
genes.GO <- select(org.Hs.eg.db, rownames(countdata.normalized), "GO", "SYMBOL")
genes.GO$Description <- select(GO.db,genes.GO$GO, "TERM")$TERM
genes.GO.group <- data.frame(genes.GO %>% group_by(SYMBOL,ONTOLOGY) %>% summarise_at("Description", paste, collapse = ";"))
genes.GO.df <- Reduce(rbind,mclapply(unique(genes.GO.group$SYMBOL), function(x){
  y <- genes.GO.group[genes.GO.group$SYMBOL==x,]
  if(!is.na(y[,"ONTOLOGY"])[1]){
    data.frame(SYMBOL=x,
               BP=ifelse(nrow(y[y$ONTOLOGY=="BP",])>0,y[y$ONTOLOGY=="BP","Description"],NA),
               MF=ifelse(nrow(y[y$ONTOLOGY=="MF",])>0,y[y$ONTOLOGY=="MF","Description"],NA),
               CC=ifelse(nrow(y[y$ONTOLOGY=="CC",])>0,y[y$ONTOLOGY=="CC","Description"],NA))
  }else{
    data.frame(SYMBOL=x,BP=NA,MF=NA,CC=NA)
  }
}, mc.cores = 4))
save(genes.GO, file = "/home/nina/Documents/Virus_project/genes2GO.RData")

# Create all DESeq2 comparisons from comparison table
if (!dir.exists(paste(output_folder, "deseq2_comparisons_shrunken", sep = ""))) {
    dir.create(paste(output_folder, "deseq2_comparisons_shrunken", sep = ""))
}
if (!dir.exists(paste(output_folder, "plots", sep = ""))) {
  dir.create(paste(output_folder, "plots", sep = ""))
}

res.list <- list()
res.list.filter <- list()
#for (n in 1:1) {
res.list.shrunk <- mclapply(1:nrow(comparisons.df), function(n){
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
    res <- merge(res, genes.GO.df, by = "SYMBOL") # add associated GO terms
    res_filter <- res[rowSums(res[,grep("normalized",colnames(res))])>0,] # remove genes with no read counts
    res_filter <- res_filter[order(res_filter$log2FoldChange, decreasing = TRUE),] # sort for LFC
    write.csv(res_filter, file = paste(output_folder, "deseq2_comparisons_shrunken/deseq2_results_", comparisons.df[n,2], "_Vs_", comparisons.df[n,1], ".csv", sep = ""), row.names = TRUE, col.names = NA)
    write.xlsx(res_filter, paste(output_folder, "deseq2_comparisons_shrunken/deseq2_results_", comparisons.df[n,2], "_Vs_", comparisons.df[n,1], ".xlsx", sep = ""))
    
    res[is.na(res$`-log10(padj)`),"-log10(padj)"] <- 0
    res[is.na(res$padj),"padj"] <- 1
    glXYPlot(x = res$log2FoldChange, y = res$`-log10(padj)`, counts = res[,c(10:(ncol(res_filter)-1))], groups = sub("_[[:digit:]]$","",colnames(res[,c(10:(ncol(res_filter)-1))])), 
             status = ifelse(res$log2FoldChange>1 & res$padj<0.05, 1, ifelse(res$log2FoldChange<(-1) & res$padj<0.05, -1, 0)), 
             xlab = "log2FoldChange", ylab = "-log10(padj)", anno = res, side.main = "SYMBOL",
             display.columns = c("SYMBOL", "UNIPROTKB", "PROTEIN.NAMES", "PATHWAY", "padj"), 
             html = paste("Volcano_", comparisons.df[n,2], "_Vs_", comparisons.df[n,1], sep = ""),
             folder = "glimma_plots", path = output_folder, launch = F)
    
    return(res)
}, mc.cores = threads)
names(res.list.shrunk) <- comparisons.df[,2]

res.list.shrunk.filter <- sapply(res.list.shrunk, function(x){x[which(apply(x[,grep("normalized", colnames(x))],1,max) >= 10 & abs(x$log2FoldChange) > 1),]})

lfc.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL"),lapply(names(res.list), function(x){x.df <- data.frame(res.list[[x]][,c("SYMBOL","log2FoldChange")]); colnames(x.df) <- c("SYMBOL",x); return(x.df)}))
rownames(lfc.df) <- lfc.df$SYMBOL
lfc.df <- lfc.df[,-1]
lfc.df <- lfc.df[,mixedorder(colnames(lfc.df))]

padj.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL"),lapply(names(res.list), function(x){x.df <- data.frame(res.list[[x]][,c("SYMBOL","padj")]); colnames(x.df) <- c("SYMBOL",x); return(x.df)}))
rownames(padj.df) <- padj.df$SYMBOL
padj.df <- padj.df[,-1]
padj.df <- padj.df[,mixedorder(colnames(padj.df))]
write.xlsx(list(log2FoldChange=lfc.df, padj=padj.df), file = paste(output_folder,"deseq2_comparisons_shrunken/expression_data_all.xlsx", sep = ""), row.names = T)

lfc.df[is.na(lfc.df)] <- 0
save.image(paste(output_folder, "/deseq2.RData", sep = ""))

# Count differentially expressed genes
padj_cut <- 0.05
count.genes <- data.frame()
for (sample in names(res.list)) {
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
colnames(count.genes) <- c("Sample", "Count","Direction","LFC_cutoff")
count.genes <- as.data.frame(count.genes)
write.table(count.genes, file = paste(output_folder, "deseq2_comparisons_shrunken/gene_count.csv", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE)

count.genes <- separate(count.genes, "Sample", c("Virus","Time"), "_", F)
count.genes$Count <- as.numeric(count.genes$Count) 
for(LFC.cut in unique(count.genes$LFC_cutoff)){
  p <- ggplot(count.genes[count.genes$LFC_cutoff==LFC.cut,], aes(y=Count, x=factor(Time, levels = c("3h","6h","12h","24h","48h","BPL")), group=Direction, fill=Direction)) + 
    geom_bar(stat = "identity", position = "stack", color = "black") + facet_wrap(~Virus) + xlab("Time") + 
    scale_y_continuous(breaks = pretty(count.genes$Count[count.genes$LFC_cutoff==LFC.cut], n=10), labels = abs(pretty(count.genes$Count[count.genes$LFC_cutoff==LFC.cut], n=10))) + 
    geom_hline(yintercept = 0) + scale_fill_manual(values=c(up="red", down="green"))
  ggsave(paste0("DEG_count_LFC",LFC.cut,".svg"), p, "svg", output_folder)
  ggsave(paste0("DEG_count_LFC",LFC.cut,".png"), p, "png", output_folder, width = 10, height = 7)
}

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


# GSEA analysis
gsea.list <- list()
for(n in names(res.list)){
  print(n)
	res <- res.list[[n]]
	res <- res[rowSums(res[,grep("normalized",colnames(res))])>0,] # remove genes with no read counts
	gsea.list[[n]] <- calc_gsea(res, n, sort.by = "log2FoldChange", REACTOME = T, ont = c("CC", "MF", "BP"), nPerm = 10000,
  	p.cut = 0.05, subfol = "GSEA")
	gc()
}
calc_ora(gene = res$SYMBOL[res$padj<padj.cut & res$log2FoldChange>LFC.cut], main = "", filename = paste0(comparisons.df[n,2],"up"), subdir = "ORA/",)
