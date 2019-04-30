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

# Import condition file
conditiontable <- read.csv(condition_file, header = FALSE, row.names = 1, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
colnames(conditiontable) <- c('condition')
conditiontable <- separate(conditiontable, condition, c("treatment", "time"), "_", FALSE)
conditiontable <- conditiontable[mixedorder(rownames(conditiontable)),]
#colnames(conditiontable) <- c('condition')
#conditiontable$treatment <- strsplit(conditiontable$condition, "_", fixed = TRUE)[0]
#conditiontable$time <- strsplit(conditiontable$condition, "_", fixed = TRUE)[1]
condition <- as.factor(conditiontable[, 1])
#conditiontable$virus <- sub("_.*","",conditiontable$condition) #remain everything before '_' -> virus name
print(class(condition))
print(condition)

# Import group comparison file
comparisons.df <- read.csv(comparisons_file, header = FALSE, sep = "\t", comment.char = "#", colClasses = c("character"))
#virus_names <- unique(sub("_.*", "", comparisons.df[,2]))
print(comparisons.df)

deseqDataset <- DESeqDataSetFromMatrix(countData = countdata, colData = conditiontable, design = ~ condition)

# Write normalized count table
deseqDataset <- estimateSizeFactors(deseqDataset)
countdata.normalized <- counts(deseqDataset, normalized = TRUE)
write.table(countdata.normalized, file = paste(output_folder, "counts_normalized.txt", sep = ""), sep = "\t", row.names = TRUE, col.names = NA)

pdf(paste(output_folder, "/counts.pdf", sep = ""), width = 25, height=15)
par(mfrow=c(1,2), mar=c(20, 4, 5, 2) + 0.1, lwd=2)
boxplot(countdata, outline = FALSE, las = 2, ylab = "raw reads", xlab = "")
boxplot(countdata.normalized, outline = FALSE, las = 2, ylab = "DESeq2 normalized reads", xlab = "")
dev.off()

deseq.results <- DESeq(object = deseqDataset, parallel = FALSE)

deseq.results.vst <- vst(deseq.results, blind = FALSE) # or vst()

deseq.results <- deseq.results[ rowSums(counts(deseq.results)) > 0, ] # remove rows that have no reads

# Plot PCA
pca <- plotPCA(deseq.results.vst, intgroup = c("treatment", "time"), returnData = TRUE) 
plot_PCA <- ggplot(pca, aes(PC1, PC2, color = treatment, shape = time)) + geom_point(size=3)
ggsave("PCA_vst.png", plot = plot_PCA, device = "png", path = output_folder)

pdf(paste(output_folder, 'correlation_heatmap.pdf', sep = ""), width = 15, height = 15, onefile = FALSE)
print(create_correlation_matrix(countdata.normalized, conditiontable))
dev.off()

pdf(paste(output_folder, 'counts_assignment.pdf', sep = ""), width = 15, height = 15)
invisible(lapply(create_feature_counts_statistics(feature_counts_log_file), print))
dev.off()

# Create all DESeq2 comparisons from comparison table
dir.create(paste(output_folder, "deseq2_comparisons", sep = ""))

res.list <- list()
lfc <- NULL
lfc_filtered_reg <- list()
#lfc_filtered_not <- list()
lfc_names <- NULL
for (n in 1:nrow(comparisons.df)) {
    res <- results(deseq.results, contrast = c("condition", comparisons.df[n,2], comparisons.df[n,1]))
    res <- as.data.frame(res)
    res.list[[comparisons.df[n,2]]] <- res
    lfc <- cbind(lfc, res$log2FoldChange)
    lfc_filtered_reg[n] <- subset(res, abs(log2FoldChange)>=4 & padj<0.05, "log2FoldChange")
    #lfc_filtered_reg[n] <- subset(res, abs(log2FoldChange)<1 & padj<0.05, "log2FoldChange")
    lfc_names <- c(lfc_names,paste(comparisons.df[n,2], comparisons.df[n,1], sep = "_Vs_"))
    #png(paste(output_folder, "deseq2_comparisons/MAplot_", comparisons.df[n,2], "_Vs_", comparisons.df[n,1], ".png", sep = ""), width = 1500, height = 1000)
    #plotMA(res, ylim = c(-5, 5), main = paste(comparisons.df[n,2], comparisons.df[n,1], sep = "_Vs_"))
    #dev.off()
    write.table(res, file = paste(output_folder, "deseq2_comparisons/deseq2_results_", comparisons.df[n,2], "_Vs_", comparisons.df[n,1], ".csv", sep = ""), sep = "\t", row.names = TRUE, col.names = NA)
}
colnames(lfc) <- lfc_names
names(lfc_filtered_reg) <- lfc_names
print(head(lfc_filtered_reg))

save.image(paste(output_folder, "/deseq2.RData", sep = ""))

# create boxplot over LFC (for every virus separatly) -> show course of infection 
for (virus in unique(conditiontable$treatment)) {
    if (!grepl("Mock", virus)){
	print(virus)
	png(paste(output_folder, "Boxplot_", virus, "_Vs_Control", ".png", sep = ""), width = 1500, height = 1000)
	boxplot(lfc[,grep(virus, colnames(lfc))])
	dev.off()
	png(paste(output_folder, "Boxplot_", virus, "_Vs_Control_LFC>4", ".png", sep = ""), width = 1500, height = 1000)
	boxplot(lfc_filtered_reg[grep(virus,names(lfc_filtered_reg))])
	dev.off()
	#png(paste(output_folder, "Boxplot_", virus, "_Vs_Control_LFC<1", ".png", sep = ""), width = 1500, height = 1000)
	#boxplot(lfc_filtered[grep(virus,names(lfc_filtered))])
	#dev.off()
    }
}

human.uniprot <- UniProt.ws(taxId=9606)
intersections_up <- list()
intersections_down <- list()
for(time in unique(conditiontable$time)[-5]){
    print(time)
    for(cutoff in c(0,1,2,3,4)){
	print(cutoff)
	res.list.sub.up <- lapply(lapply(res.list[grep(time,names(res.list))], subset, padj<0.05 & log2FoldChange>cutoff), rownames)
	res.list.sub.down <- lapply(lapply(res.list[grep(time,names(res.list))], subset, padj<0.05 & log2FoldChange<(-cutoff)), rownames)
	#print(lapply(res.list.sub.down, nrow))
	#print(names(res.list.sub.down))
	overlap <- venn.diagram(res.list.sub.up, 
			filename = paste0("dge_analyses2/",output_folder,"overlap_",time,"_LFC>",cutoff,".svg"),
			imagetype = "svg", fill = c("red","blue","green","orange","purple"), margin = 0.1, width = 15, height = 15,
			cex = 2, cat.cex = 2, lwd = 1)
	overlap <- venn.diagram(res.list.sub.down, 
				filename = paste0("dge_analyses2/",output_folder,"overlap_",time,"_LFC<-",cutoff,".svg"),
				imagetype = "svg", fill = c("red","blue","green","orange","purple"), margin = 0.1, width = 15, height = 15,
				cex = 2, cat.cex = 2, lwd = 1)
	overlap_up <- attr(venn(res.list.sub.up, show.plot = FALSE),"intersections")
	overlap_down <- attr(venn(res.list.sub.down, show.plot = FALSE),"intersections")
	intersections_up[[paste0(time,"_",cutoff)]] <- overlap_up[paste(names(res.list.sub.up),collapse=":")]
	intersections_down[[paste0(time,"_",cutoff)]] <- overlap_down[paste(names(res.list.sub.up),collapse=":")]
	#intersections_up <- sapply(intersections_up, "[", seq(max(sapply(intersections_up, length))))
	#print(intersections)
	#write.csv(intersections, file = paste0("/nfs/sfb1021/SFB1021_Virus/dge_analyses2/analyses/host/deseq2/overlap_|LFC|>",cutoff,".csv"))
    }
}
bitr_24h_3 <- bitr(intersections_up$`24h_3`[[1]], fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
res <- select(human.uniprot, bitr_24h_3$ENTREZID, c("PROTEIN-NAMES", "GO", "REVIEWED"), "ENTREZ_GENE")
res_merge <- merge(x=bitr_24h_3, y=res, by.y="ENTREZ_GENE", by.x="ENTREZID")
write.table(res_merge[res_merge$REVIEWED=="reviewed",-c(1,5)],file="/nfs/sfb1021/SFB1021_Virus/dge_analyses2/analyses/host/deseq2/24h_3_up.txt", row.names=F, col.names=T)

sapply(names(intersections_up),function(name){write(unlist(intersections_up[[name]]), file=paste("/nfs/sfb1021/SFB1021_Virus/dge_analyses2/analyses/host/deseq2/overlap/",name,"_up.csv", sep=""), sep="\t")})
sapply(names(intersections_down),function(name){write(unlist(intersections_down[[name]]), file=paste("/nfs/sfb1021/SFB1021_Virus/dge_analyses2/analyses/host/deseq2/overlap/",name,"_down.csv", sep=""), sep="\t")})

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
