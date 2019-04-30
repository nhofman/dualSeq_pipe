library("ggplot2")

# calculate RPKM with total number of reads in each sample (incl. unmapped reads)
calculate_rpkm <- function(feature_counts, feature_counts_log_file){
    counts <- read.table(feature_counts, header = TRUE, row.names = 1)
    fc_log <- read.table(feature_counts_log_file, header = TRUE, row.names = 1)
    
    rpkm.virus <- NULL
    for(i in 6:ncol(counts)){
	rpkm <- counts[, i] / ((counts[, "Length"]/1000) * (sum(counts[, i]) / 1000000))
	rpkm.virus <- cbind(rpkm.virus, rpkm)
	colnames(rpkm.virus)[i] <- colnames(countdata.virus)[i]
    }
    rownames(rpkm.virus) <- counts$Geneid
    rpkm.virus.t <- cbind(as.data.frame(t(rpkm.virus)),sample=substr(colnames(rpkm.virus),4, nchar(colnames(rpkm.virus))-5))
    rpkm.virus.t.mean <- aggregate(rpkm.virus.t[,-c(8)], by = list(rpkm.virus.t$sample), FUN = mean)

    boxplot(t(rpkm.virus.t.mean[,-1]), names = rpkm.virus.t.mean$Group.1, las = 2)

    
}

args <- commandArgs(TRUE)
#counttable_file <- args[match('--counttable', args) + 1]
condition_file <- args[match('--conditions', args) + 1]
feature_counts <- args[match('--featcounts', args) + 1]
feature_counts_log_file <- args[match('--featcounts-log', args) + 1]
feature_counts_log_file_star <-args[match('--featcounts-log-star', args) + 1]
output_folder <- args[match('--output', args) + 1]
threads <- args[match('--threads', args) + 1]

conditiontable <- read.csv(condition_file, header = FALSE, row.names = 1, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
colnames(conditiontable) <- c('condition')
conditiontable <- separate(conditiontable, condition, c("treatment", "time"), "_", FALSE)





