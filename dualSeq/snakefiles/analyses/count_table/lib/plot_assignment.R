library(ggplot2)
library(reshape2)

# Parse arguments
args <- commandArgs(F)
feature_counts_log_file <- args[match('--featurecounts-log', args) + 1]

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

# Bar charts showing the assignment of alignments to genes (featureCounts statistics)
plot.list <- create_feature_counts_statistics(feature_counts_log_file)
invisible(sapply(list("pdf", "png"), function(p){
  sapply(names(plot.list), function(x){
    ggsave(plot=plot.list[[x]],filename=paste0("counts_assignment_", x, ".", p), device=p, path=output_folder, width=18, height=12)
  })}))
