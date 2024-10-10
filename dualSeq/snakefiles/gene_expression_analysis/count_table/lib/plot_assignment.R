library(ggplot2)
library(reshape2)
library(gtools)
library(tidyr)

# Parse arguments
args <- commandArgs(F)
feature_counts_log_file <- args[match('--featurecounts-log', args) + 1]
output_folder <- args[match('--output', args) + 1]

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
    geom_bar(aes(fill = Group), stat = "identity", group = 1, color = "black") +
    facet_wrap(~ Virus, scales = "free_x") + xlab("Time_Replicate") +
    theme(text = element_text(face = "bold"), line = element_line(linewidth = 0.25),
          axis.line.x.top = element_blank(), axis.line.x.bottom = element_line(color = "black", linewidth = 0.25),
          axis.line.y.right = element_blank(), axis.line.y.left = element_line(color = "black", linewidth = 0.25),
          panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "gray58"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
          axis.title.x = element_text(size = 12, margin = margin(t=1,r=0,b=0,l=0)), 
          axis.title.y = element_text(size = 12, margin = margin(t=0,r=1,b=0,l=0)),
          axis.text = element_text(size = 10), axis.text.x = element_text(angle = 90), 
          strip.text = element_text(size = 12), strip.background = element_rect(fill = NA, color = NA), #panel.spacing = unit(2, "lines"),
          legend.text = element_text(size = 10, face = "plain"), legend.title = element_blank())
  
  colnames(dpctm) <- c("Sample", "Group", "Reads")
  dpctm$Group = factor(dpctm$Group, levels = rev(levels(dpctm$Group)[order(levels(dpctm$Group))]))
  dpctm <- separate(dpctm, "Sample", c("Virus","Time"), "_", F, extra = "merge")
  assignment.relative <- ggplot(dpctm[dpctm$Reads > 0,], aes(x = factor(Time, levels = mixedsort(unique(Time))), y = Reads)) +
    geom_bar(aes(fill = Group), stat = "identity", group = 1, color = "black") +
    facet_wrap(~ Virus, scales = "free_x") + xlab("Time_Replicate") +
    theme(text = element_text(face = "bold"), line = element_line(linewidth = 0.25),
          axis.line.x.top = element_blank(), axis.line.x.bottom = element_line(color = "black", linewidth = 0.25),
          axis.line.y.right = element_blank(), axis.line.y.left = element_line(color = "black", linewidth = 0.25),
          panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "gray58"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
          axis.title.x = element_text(size = 12, margin = margin(t=1,r=0,b=0,l=0)), 
          axis.title.y = element_text(size = 12, margin = margin(t=0,r=1,b=0,l=0)),
          axis.text = element_text(size = 10), axis.text.x = element_text(angle = 90), 
          strip.text = element_text(size = 12), strip.background = element_rect(fill = NA, color = NA), #panel.spacing = unit(2, "lines"),
          legend.text = element_text(size = 10, face = "plain"), legend.title = element_blank())
  
  return(list("absolute"=assignment.absolute, "relative"=assignment.relative))
}

# Bar charts showing the assignment of alignments to genes (featureCounts statistics)
plot.list <- create_feature_counts_statistics(feature_counts_log_file)
invisible(sapply(list("pdf", "png"), function(p){
  sapply(names(plot.list), function(x){
    ggsave(plot=plot.list[[x]],filename=paste0("counts_assignment_", x, ".", p), device=p, path=output_folder, width=20, height=12)
  })}))
