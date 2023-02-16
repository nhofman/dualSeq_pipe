stat_host <- "Documents/Virus_project/analyses/host/deseq2/mapping/host_stats/" 
stat_virus <- "Documents/Virus_project/analyses/host/deseq2/mapping/virus_stats/" 

for (package in c("ggplot2", "tidyr", "gtools", "dplyr")) {
  if (!(package %in% rownames(installed.packages()))) {
    stop(paste('Package "', package, '" not installed', sep=""))
  } else {
    print(paste("Import:", package))
    library(package, character.only=TRUE)
  }
}

color_file <- "Documents/Virus_project/analyses/host/deseq2/colors.txt"
color.df <- read.table(color_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
color <- color.df[,2]
names(color) <- color.df[,1]
color <- c("host:uniquely mapped"="gray60", "host:multimapped"="gray80", color)

virus.levels <- c("H1N1","H5N1","RVFV","SFSV","RSV","NiV","EBOV","MARV","LASV")

stats.df <- Reduce(rbind,lapply(c(list.files(stat_host, ".*_stats.txt", full.names = T),list.files(stat_virus, ".*_stats.txt", full.names = T)), function(x){
  tmp <- read.table(x, header = T, sep = "\t")
  tmp$multimapped <- tmp$mapped - tmp$uniquely_mapped - tmp$mapped_to_both
  tmp$mapped <- tmp$mapped - tmp$mapped_to_both
  tmp$mapped_to_both <- sum(tmp$mapped_to_both)
  return(tmp)
}))
write.table(stats.df[,c(1,2,3,4,7,5,6)], paste0(output_folder, "stats/mapping_statistic.tsv"), sep = "\t", row.names = F)
for(x in unique(stats.df$sample)){
  tmp <- sum(stats.df[stats.df$sample==x,"mapped_to_both"])
  stats.df[stats.df$sample==x,"mapped_to_both"] <- tmp
}
stats.print <- merge(stats.df[stats.df$group=="host",-6], stats.df[stats.df$group=="virus",-c(5,6)], by = c("sample","total_reads"), all = T)
colnames(stats.print) <- c("Sample","# reads in sample","Host: mapped","Host: uniquely mapped","Mapped to both","Host: multimapped","Virus: mapped","Virus: uniquely mapped","Virus: multimapped")
stats.print <- stats.print[,c(1,2,3,4,6,7,8,9,5)]
write.table(stats.print, paste0(output_folder, "mapping_statistic_paper.tsv"), sep = "\t", row.names = F)

# transpose data frame for plotting
stats.t <- data.frame("Sample"=character(), "Total"=integer(), "Count"=integer(), "Organism"=character(), "Category"=character())
stats.t <- data.frame(Map(c, stats.t, data.frame(stats.df[,c(1,2,3,6)], "mapped")))
stats.t <- data.frame(Map(c, stats.t, data.frame(stats.df[,c(1,2,4,6)], "uniquely mapped")))
stats.t <- data.frame(Map(c, stats.t, data.frame(stats.df[,c(1,2,7,6)], "multimapped")))
stats.t <- data.frame(Map(c, stats.t, data.frame(stats.df[stats.df$group=="host",c(1,2,5,6)], "both")))

stats.t <- separate(stats.t, "Sample", c("Virus","Time", "Rep"),"_",F,extra = "merge")
stats.t$Count_percent <- stats.t$Count/stats.t$Total
stats.t$Category_2 <- paste(stats.t$Organism, stats.t$Category, sep = ":")
stats.t$Time_Rep <- paste(stats.t$Time, stats.t$Rep, sep=":")
# calculate mean of Count_percent for all replicates
stats.t.mean <- stats.t %>% group_by(Virus,Time,Organism,Category) %>% summarise_at("Count_percent", "mean")
stats.t.mean$Category_2 <- paste(stats.t.mean$Organism, stats.t.mean$Category, sep = ":")
stats.t.mean$Category_new <- NA
for(i in 1:nrow(stats.t.mean)){
  stats.t.mean$Category_new[i] <- ifelse(stats.t.mean$Organism[i]=="virus", sub("virus.*", stats.t.mean$Virus[i], stats.t.mean$Category_2[i]), stats.t.mean$Category_2[i])
}

ncol.facet <- ceiling(sqrt(length(unique(stats.t.mean$Virus[!grepl("Mock",stats.t.mean$Virus)]))))
nrow.facet <- ceiling(length(unique(stats.t.mean$Virus[!grepl("Mock",stats.t.mean$Virus)]))/ncol.facet)
times <- length(unique(stats.t.mean$Time[stats.t.mean$Virus==stats.t.mean$Virus[1]]))

# plot mapping statistic for host + virus samples
p <- ggplot(stats.t.mean[stats.t.mean$Category!="mapped" & !grepl("Mock",stats.t.mean$Virus),], aes(x=factor(Time, levels = unique(mixedsort(Time))), y=Count_percent, group=Category_2, fill=Category_2)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Virus, scales = "free_x", ncol = ncol.facet) +
  xlab("Time") + ylab("Ratio of mapped reads") + labs(fill="Category") +
  scale_fill_manual(labels = c("Mapped to both", "Host: Multimapped", "Host: Uniquely mapped", "Virus: Multimapped", "Virus: Uniquely mapped"), 
                    values = c("gray10", "gray80", "gray40", "lightsteelblue", "steelblue")) + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 15), strip.text = element_text(size = 15, face = "bold", vjust = 0.3),
        legend.text = element_text(size = 15), legend.title = element_text(size = 18), strip.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = NA), panel.grid = element_line(colour = "grey"))
ggsave("mapping_statistic.svg", p, "svg", output_folder, width = 0.6*times*(ncol.facet+1), height = 2.5*(nrow.facet+1))
ggsave("mapping_statistic.pdf", p, "pdf", output_folder, width = 0.6*times*(ncol.facet+1), height = 2.5*(nrow.facet+1))
system(paste0("inkscape -l ", output_folder, "mapping_statistic.svg ", output_folder, "mapping_statistic.pdf"))

p <- ggplot(stats.t.mean[stats.t.mean$Category!="mapped" & !stats.t.mean$Category_2%in%c("virus:multimapped","host:both") & !grepl("Mock",stats.t.mean$Virus),], 
            aes(x=factor(Time, levels = unique(mixedsort(Time))), y=Count_percent, fill=factor(Category_new, levels = c("host:multimapped", "host:uniquely mapped", virus.levels)))) +
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(~ factor(Virus, levels = virus.levels), scales = "free_x", ncol = ncol.facet) + scale_x_discrete(labels=c("3 h", "6 h", "12 h", "24 h", "BPL")) +
  xlab("Time") + ylab("Ratio of mapped reads") + labs(fill="") + scale_fill_manual(values=color, breaks = c("host:uniquely mapped", "host:multimapped"), labels = c("Host: uniquely mapped", "Host: multimapped")) + 
  theme(text = element_text(family = "Helvetica", face = "bold"), axis.title = element_text(size = 35), axis.text = element_text(size = 28), 
        axis.title.x = element_text(margin = margin(7, 0, 0, 0, "mm")), axis.title.y = element_text(margin = margin(0, 7, 0, 0, "mm")), 
        strip.text = element_text(size = 35, face = "bold", vjust = 0.3),
        legend.text = element_text(size = 28), legend.title = element_text(size = 35), #strip.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"), panel.grid = element_blank(), axis.ticks.x = element_blank())
#ggsave("mapping_statistic_modified.svg", p, "svg", output_folder, width = 0.6*times*(ncol.facet+1), height = 2.5*(nrow.facet+1))
ggsave("mapping_statistic_modified.pdf", p, "pdf", output_folder, width = 23, height = 12) #width = 0.6*times*(ncol.facet+1), height = 2.5*(nrow.facet+1))
system(paste0("inkscape -l ", output_folder, "mapping_statistic_modified.svg ", output_folder, "mapping_statistic_modified.pdf"))
