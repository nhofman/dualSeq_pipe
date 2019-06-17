args <- commandArgs(TRUE)
#counttable_file <- args[match('--counttable', args) + 1]
condition_file <- args[match('--conditions', args) + 1]
counttable_host_file <- args[match('--counttable_host', args) + 1]
counttable_sum_host_file <- args[match('--counttable_sum_host', args) + 1]
#feature_counts_log_file <- args[match('--featcounts', args) + 1]
output_folder <- args[match('--output', args) + 1]
#color <- args[match('--color', args) + 1]
threads <- args[match('--threads', args) + 1]
gff <- args[match('--gff', args) + 1]
bam <- args[match('--bam', args) + 1]
mapping_stat_host_file <- args[match('--stat_host', args) + 1]

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
for (package in c("BiocParallel", "pheatmap", "ggplot2", "reshape2", "gplots", "tidyr", "gtools", "Rsamtools")) {
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

save.image(paste(output_folder, "/deseq2.RData", sep = ""))
# Import count table (featureCounts)
#countdata.raw <- read.csv(counttable_file, header = TRUE, row.names = 1, sep = "\t", comment.char = "#")
#countdata <- as.matrix(countdata.raw[, c(6 : length(countdata.raw))])
#colnames(countdata) <- as.vector(sapply(colnames(countdata), function(x) gsub("mapping\\..*\\.(.*)\\.bam", "\\1", x)))
#countdata <- countdata[,mixedorder(colnames(countdata))]
#print(head(countdata))

# Import gff file
#gff_files <- lapply(strsplit(gff, ",")[[1]],read.table,comment.char="#", header=F, sep="\t", col.names=c("seqname","source","feature","start","end","score","strand","frame","attribute"))
#gff <- read.table(gff_file, comment.char="#", header=F, sep="\t")
#colnames(gff) <- c("seqname","source","feature","start","end","score","strand","frame","attribute")
#print(gff_files) #[grep("^region",gff$feature),c(1:8)])

# Import bam files
print(bam)
bam_files <- unlist(strsplit(bam,","))
print(bam_files)
bam_list <- BamFileList(bam_files)
bam_count <- countBam(bam_list, param=ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE, isSecondaryAlignment=FALSE)))
rownames(bam_count) <- sapply(rownames(bam_count),function(x){substr(x,1,nchar(x)-4)})
print(bam_count)
header <- scanBamHeader(bam_files, what="targets")
header_LN <- sapply(header,function(x){sum(x$targets)})

# Import host countdata
countdata.host <- read.csv(counttable_host_file, header = TRUE, row.names = 1, sep = "\t", comment.char = "#")
#countdata.host <- as.matrix(countdata.host[, c(6 : length(countdata.host))])
colnames(countdata.host) <- as.vector(sapply(colnames(countdata.host), function(x) gsub("mapping\\..*\\.(.*)\\.bam", "\\1", x)))
#countdata.host <- countdata.host[,mixedorder(colnames(countdata.host))]
countdata.summary.host <- read.csv(counttable_sum_host_file, header = TRUE, row.names = 1, sep = "\t")
colnames(countdata.summary.host) <- as.vector(sapply(colnames(countdata.summary.host), function(x) gsub("mapping\\..*\\.(.*)\\.bam", "\\1", x)))
mapping_stat_host <- read.table(mapping_stat_host_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
print(mapping_stat_host)

# Import condition file
conditiontable <- read.csv(condition_file, header = FALSE, row.names = 1, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
colnames(conditiontable) <- c('condition')
conditiontable <- separate(conditiontable, condition, c("treatment", "time"), "_", FALSE)
conditiontable <- conditiontable[mixedorder(rownames(conditiontable)),]
condition <- as.factor(conditiontable[, 1])
print(class(condition))
print(condition)


#pdf(paste(output_folder, 'correlation_heatmap.pdf', sep = ""), width = 15, height = 15, onefile = FALSE)
#print(create_correlation_matrix(countdata.normalized, conditiontable))
#dev.off()

#pdf(paste(output_folder, 'counts_assignment.pdf', sep = ""))
#invisible(lapply(create_feature_counts_statistics(feature_counts_log_file), print))
#dev.off()

#colSums(countdata.host[,c(6:length(countdata.host))]/countdata.host[,'Length'])

rpkm.virus <- NULL
countdata.mean <- NULL
mean.ratio <- NULL
mean.tpm <- NULL
for(treatment in unique(conditiontable$treatment)) { 
    print(treatment)
    for(time in unique(conditiontable$time[conditiontable$treatment==treatment])) {
	print(time)
	countdata.sub <- bam_count[grep(paste(treatment,time,sep="_"), rownames(bam_count)),]
	ratio <- NULL
	tpm <- NULL
	for(i in 1:nrow(countdata.sub)){
	    #rpkm <- cbind(rpkm,countdata.sub[,i]/((countdata.raw[, "Length"]/1000)*((sum(countdata.host[,grep(colnames(countdata.sub)[i], colnames(countdata.host))])+sum(countdata.sub[,i]))/1000000)))
	    print(countdata.sub[i,"records"])
	    print(mapping_stat_host$uniquely_mapped[grep(paste(treatment, time, i, sep="_"), mapping_stat_host$sample)])
	    ratio <- cbind(ratio,countdata.sub[i,"records"]/mapping_stat_host$uniquely_mapped[grep(paste(treatment, time, i, sep="_"), mapping_stat_host$sample)])
	    print(ratio)#sum(countdata.summary.host[c(1,10,12),grep(rownames(countdata.sub)[i], colnames(countdata.summary.host))]))
	    rpk_virus <- (countdata.sub[i,"records"]*mapping_stat_host$avg_mapped_read_length[mapping_stat_host$sample==paste(treatment, time, i, sep="_")])/(unique(header_LN[grep(paste(treatment,time,sep="_"),names(header_LN))]))
	    scaling <- (sum((countdata.host[,grep(rownames(countdata.sub)[i], colnames(countdata.host))]*mapping_stat_host$avg_mapped_read_length[grep(paste(treatment, time, i, sep="_"), mapping_stat_host$sample)])/countdata.host[,'Length'])+rpk_virus)/1000000
	    tpm <- cbind(tpm,rpk_virus/scaling)
	}
	mean.ratio <- rbind(mean.ratio, c(rowMeans(ratio),time,treatment)) 
	mean.tpm <- rbind(mean.tpm, c(rowMeans(tpm),time,treatment)) 
    }
}
colnames(mean.tpm) <- c("tpm","time","virus")
mean.tpm <- as.data.frame(mean.tpm)
print(mean.tpm)
colnames(mean.ratio) <- c("ratio","time","virus")
mean.ratio <- as.data.frame(mean.ratio)
print(mean.ratio)
plot_ratio <- ggplot(mean.ratio, aes(x=factor(time,level=unique(mixedsort(conditiontable$time))),y=as.numeric(as.character(ratio)),color=virus,group=virus))+geom_line()+labs(x="Time",y="Ratio virus/human")
ggsave("ratio.svg", plot_ratio, "svg", "analyses/virus/deseq2/")
plot_tpm <- ggplot(mean.tpm, aes(x=factor(time,level=unique(mixedsort(conditiontable$time))),y=as.numeric(as.character(tpm)),color=virus,group=virus))+geom_line()+labs(x="Time",y="TPM")
ggsave("tpm.svg", plot_tpm, "svg", "analyses/virus/deseq2/")
#colnames(countdata.mean) <- unique(conditiontable$condition)
#rownames(rpk.virus) <- rownames(countdata)
#svg(paste(output_folder, 'boxplot_mean_rpkm.svg', sep=""))
#par(mar=c(5,7,4,2))
#boxplot(countdata.mean, log = "", col = color, las = 2, ylim = c(0,300000))
#boxplot(countdata.mean, ylim = c(0,300000), las = 1, cex.axis = 2, names = unique(conditiontable$time))
#legend("topleft", unique(conditiontable$treatment), cex = 2)
#dev.off()

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

