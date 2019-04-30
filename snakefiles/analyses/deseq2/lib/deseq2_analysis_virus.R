args <- commandArgs(TRUE)
counttable_file <- args[match('--counttable', args) + 1]
condition_file <- args[match('--conditions', args) + 1]
comparisons_file <- args[match('--comparisons', args) + 1]
counttable_host_file <- args[match('--counttable_host', args) + 1]
feature_counts_log_file <- args[match('--featcounts', args) + 1]
output_folder <- args[match('--output', args) + 1]
color <- args[match('--color', args) + 1]
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
for (package in c("BiocParallel", "pheatmap", "ggplot2", "reshape2", "gplots", "tidyr", "gtools")) {
    if (!(package %in% rownames(installed.packages()))) {
        stop(paste('Package "', package, '" not installed', sep=""))
    } else {
        print(paste("Import:", package))
        library(package, character.only=TRUE)
    }
}



create_correlation_matrix <- function(countdata, conditiontable) {
    countdata.normalized.processed <- as.matrix(countdata)
    countdata.normalized.processed <- countdata.normalized.processed[rowSums(countdata.normalized.processed) >= 10,]
    countdata.normalized.processed <- log2(countdata.normalized.processed + 1)
    sample_cor <- cor(countdata.normalized.processed, method = 'pearson', use = 'pairwise.complete.obs')

    return(pheatmap(sample_cor, annotation_col = conditiontable, annotation_row = conditiontable))
}


create_feature_counts_statistics <- function(featureCountsLog) {
    d <- read.table(featureCountsLog, header = T, row.names = 1)
    colnames(d) <- gsub(".bam", "", colnames(d))

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
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    return(list(assignment.absolute, assignment.relative))
}

plotHeatmap2 <- function(x, name = "no_name_set.pdf", row_subset = NA, distMethod = "euclidean", clusterMethod = "complete", clrn = NA, ...){
    require('gplots')
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

countdata.host <- read.csv(counttable_host_file, header = TRUE, row.names = 1, sep = "\t", comment.char = "#")
countdata.host <- as.matrix(countdata.host[, c(6 : length(countdata.host))])
colnames(countdata.host) <- as.vector(sapply(colnames(countdata.host), function(x) gsub("mapping\\..*\\.(.*)\\.bam", "\\1", x)))
countdata.host <- countdata.host[,mixedorder(colnames(countdata.host))]

# Import condition file
conditiontable <- read.csv(condition_file, header = FALSE, row.names = 1, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
colnames(conditiontable) <- c('condition')
conditiontable <- separate(conditiontable, condition, c("treatment", "time"), "_", FALSE)
conditiontable <- conditiontable[mixedorder(rownames(conditiontable)),]
condition <- as.factor(conditiontable[, 1])
print(class(condition))
print(condition)

# Import group comparison file
comparisons.df <- read.csv(comparisons_file, header = FALSE, sep = "\t", comment.char = "#", colClasses = c("character"))
#virus_names <- unique(sub("_.*", "", comparisons.df[,2]))
comparisons.df <- comparisons.df[grep(conditiontable$treatment, comparisons.df[,1]),]
print(comparisons.df)

#deseqDataset <- DESeqDataSetFromMatrix(countData = countdata, colData = conditiontable, design = ~ condition)

# Write normalized count table
#deseqDataset <- estimateSizeFactors(deseqDataset)
#countdata.normalized <- counts(deseqDataset, normalized = TRUE)
#write.table(countdata.normalized, file = paste(output_folder, "counts_normalized.txt", sep = ""), sep = "\t", row.names = TRUE, col.names = NA)

#deseq.results <- DESeq(object = deseqDataset, parallel = FALSE)

#deseq.results.vst <- rlog(deseq.results, blind = FALSE) # or vst()

#deseq.results <- deseq.results[ rowSums(counts(deseq.results)) > 0, ] # remove rows that have no reads

# Plot PCA
#pca <- plotPCA(deseq.results.vst, intgroup = c("treatment", "time"), returnData = TRUE) 
#plot_PCA <- ggplot(pca, aes(PC1, PC2, color = treatment, shape = time)) + geom_point(size=3)
#ggsave("PCA.png", plot = plot_PCA, device = "png", path = output_folder)

#pdf(paste(output_folder, 'correlation_heatmap.pdf', sep = ""), width = 15, height = 15, onefile = FALSE)
#print(create_correlation_matrix(countdata.normalized, conditiontable))
#dev.off()

pdf(paste(output_folder, 'counts_assignment.pdf', sep = ""))
invisible(lapply(create_feature_counts_statistics(feature_counts_log_file), print))
dev.off()

rpkm.virus <- NULL
countdata.mean <- NULL
for(virus in unique(conditiontable$condition)) {
    print(virus)
    countdata.sub <- countdata[, grep(virus, colnames(countdata))]
    print(countdata.sub)
    rpkm <- NULL
    for(i in 1:ncol(countdata.sub)){
        rpkm <- cbind(rpkm,countdata.sub[,i]/((countdata.raw[, "Length"]/1000)*((sum(countdata.host[,grep(colnames(countdata.sub)[i], colnames(countdata.host))])+sum(countdata.sub[,i]))/1000000)))
    }
    #rpkm <- t(t(countdata.sub) / ((countdata.raw[, "Length"]/1000)*((colSums(countdata.host[,grep(virus, colnames(countdata.host))])+colSums(countdata.sub[,grep(virus, colnames(countdata.sub))]))/1000000)))
    #print(countdata.raw[, "Length"]/1000)
    #print(colSums(countdata.host[,grep(virus, colnames(countdata.host))]))
    #print(colSums(countdata.sub[,grep(virus, colnames(countdata.sub))]))
    print(rpkm)
    countdata.mean <- cbind(countdata.mean, rowMeans(rpkm)) 
}
colnames(countdata.mean) <- unique(conditiontable$condition)
#rownames(rpk.virus) <- rownames(countdata)
svg(paste(output_folder, 'boxplot_mean_rpkm.svg', sep=""))
par(mar=c(10,4,4,2))
boxplot(countdata.mean, log = "", col = color, las = 2, ylim = c(0,300000))
dev.off()

#svg(paste(output_folder, 'boxplot_mean_norm.svg', sep=""))
#par(mar=c(10,4,4,2))
#boxplot(countdata.mean/(countdata.raw[, "Length"]/1000), col = color, las = 2, ylim = c(0,8000000))
#dev.off()

# Create all DESeq2 comparisons from comparison table
#dir.create(paste(output_folder, "deseq2_comparisons", sep = ""))

#lfc <- NULL
#lfc_filtered_reg <- list()
#lfc_filtered_not <- list()
#lfc_names <- NULL
#for (n in 1:nrow(comparisons.df)) {
#    res <- results(deseq.results, contrast = c("condition", comparisons.df[n,2], comparisons.df[n,1]))
#    res <- as.data.frame(res)
#    lfc <- cbind(lfc, res$log2FoldChange)
    #lfc_filtered_reg[n] <- subset(res, abs(log2FoldChange)>=1 & padj<0.05, "log2FoldChange")
    #lfc_filtered_reg[n] <- subset(res, abs(log2FoldChange)<1 & padj<0.05, "log2FoldChange")
#    lfc_names <- c(lfc_names,paste(comparisons.df[n,2], comparisons.df[n,1], sep = "_Vs_"))
    #png(paste(output_folder, "deseq2_comparisons/MAplot_", comparisons.df[n,2], "_Vs_", comparisons.df[n,1], ".png", sep = ""), width = 1500, height = 1000)
    #plotMA(res, ylim = c(-5, 5), main = paste(comparisons.df[n,2], comparisons.df[n,1], sep = "_Vs_"))
    #dev.off()
#    write.table(res, file = paste(output_folder, "deseq2_comparisons/deseq2_results_", comparisons.df[n,2], "_Vs_", comparisons.df[n,1], ".csv", sep = ""), sep = "\t", row.names = TRUE, col.names = NA)
#}
#colnames(lfc) <- lfc_names
#names(lfc_filtered_reg) <- lfc_names
#print(head(lfc))

save.image(paste(output_folder, "/deseq2.RData", sep = ""))

# create boxplot over LFC (for every virus separatly) -> show course of infection 
#for (virus in unique(conditiontable$treatment)) {
#    if (!grepl("Mock", virus)){
#	print(virus)
#	png(paste(output_folder, "Boxplot_", virus, "_Vs_Control", ".png", sep = ""), width = 1500, height = 1000)
#	boxplot(countdata[,grep(virus, colnames(countdata))])
#	dev.off()
#	png(paste(output_folder, "Boxplot_", virus, "_Vs_Control.log", ".png", sep = ""), width = 1500, height = 1000)
#	boxplot(log(countdata[,grep(virus, colnames(countdata))]))
#	dev.off()
	#png(paste(output_folder, "Boxplot_", virus, "_Vs_Control_LFC>1", ".png", sep = ""), width = 1500, height = 1000)
	#boxplot(lfc_filtered_reg[grep(virus,names(lfc_filtered_reg))])
	#dev.off()
	#png(paste(output_folder, "Boxplot_", virus, "_Vs_Control_LFC<1", ".png", sep = ""), width = 1500, height = 1000)
	#boxplot(lfc_filtered[grep(virus,names(lfc_filtered))])
	#dev.off()
#    }
#}

