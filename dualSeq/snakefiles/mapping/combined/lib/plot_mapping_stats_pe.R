# Parse arguments
args <- commandArgs(T)
stat_host_file <- args[match('--stat-host', args) + 1]
stat_virus_file <- args[match('--stat-virus', args) + 1]

for (package in c("ggplot2", "tidyr", "gtools", "dplyr")) {
  if (!(package %in% rownames(installed.packages()))) {
    stop(paste('Package "', package, '" not installed', sep=""))
  } else {
    print(paste("Import:", package))
    library(package, character.only=TRUE)
  }
}

# Mapping statistics for host
mapping.stat <- read.table(stat_host_file, header = T, stringsAsFactors = F)
mapping.stat <- separate(mapping.stat, "sample", c("virus","time","rep"), "_", F)
mapping.stat$type <- "Host"
mapping.stat$uniquely_mapped_percent <- as.numeric(mapping.stat$uniquely_mapped)/as.numeric(mapping.stat$total_reads)
# Mapping statistics for virus
mapping.stat.virus <- read.table(stat_virus_file, header = T, stringsAsFactors = F)
mapping.stat.virus <- merge(mapping.stat.virus, mapping.stat[,c(1,5)], by = "sample")
mapping.stat.virus$uniquely_mapped_conc_percent <- as.numeric(mapping.stat.virus$aligned_conc_1_time)/as.numeric(mapping.stat.virus$total_read)
mapping.stat.virus$uniquely_mapped_disconc_percent <- as.numeric(mapping.stat.virus$aligned_disconc_1_time)/as.numeric(mapping.stat.virus$total_read)
mapping.stat.virus$uniquely_mapped_mixed_percent <- as.numeric(mapping.stat.virus$mates_aligned_1_time)/as.numeric(mapping.stat.virus$total_read)
mapping.stat.virus <- separate(mapping.stat.virus, "sample", c("virus","time","rep"), "_", F)
mapping.stat.combi <- mapping.stat.virus[,c("virus","time","rep","uniquely_mapped_conc_percent","total_reads")]
colnames(mapping.stat.combi)[colnames(mapping.stat.combi)=="uniquely_mapped_conc_percent"] <- "uniquely_mapped_percent"
mapping.stat.combi$type <- "Virus: concordant"
tmp <- mapping.stat.virus[,c("virus","time","rep","uniquely_mapped_disconc_percent","total_reads")]
colnames(tmp)[colnames(tmp)=="uniquely_mapped_disconc_percent"] <- "uniquely_mapped_percent"
tmp$type <- "Virus: discordant"
mapping.stat.combi <- rbind(mapping.stat.combi, tmp)
tmp <- mapping.stat.virus[,c("virus","time","rep","uniquely_mapped_mixed_percent","total_reads")]
colnames(tmp)[colnames(tmp)=="uniquely_mapped_mixed_percent"] <- "uniquely_mapped_percent"
tmp$type <- "Virus: mixed"
mapping.stat.combi <- rbind(mapping.stat.combi, tmp)
rm(tmp)
mapping.stat.combi <- rbind(mapping.stat.combi,
                            mapping.stat[,c("virus","time","rep","uniquely_mapped_percent","total_reads","type")])
mapping.stat.combi <- mapping.stat.combi %>% group_by(virus,time,type) %>% summarise_at("uniquely_mapped_percent", "mean")

plot.map <- ggplot(mapping.stat.combi[!grepl("Mock",mapping.stat.combi$virus),], aes(x=factor(time,levels=unique(mixedsort(time))), y=uniquely_mapped_percent, fill=type, group=type)) + 
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ virus, scales = "free_x") + xlab("Time") + ylab("Ratio of mapped reads") + ylim(c(0,1)) +
  scale_fill_manual(values=c("darkorange", "darkblue", "blue", "lightblue3")) + labs(fill="Type") +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12), strip.text = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 15), legend.title = element_text(size = 18, face = "bold"))
ggsave("mapping_stat_dual.svg", plot.map, "svg", "mapping/", width = 12, height = 8)

plot.map <- ggplot(mapping.stat.combi[grepl("Mock",mapping.stat.combi$virus),], aes(x=factor(time,levels=unique(mixedsort(time))), y=uniquely_mapped_percent, fill=organism, group=organism)) + 
  geom_bar(stat = "identity", position = "stack") + 
  facet_wrap(~ virus) + xlab("Time") + ylab("Ratio of mapped reads") + ylim(c(0,1)) +
  scale_fill_manual(values=c("darkorange")) + labs(fill = "Organism") +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12), strip.text = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 15), legend.title = element_text(size = 18, face = "bold"))
ggsave("mapping_stat_mock.svg", plot.map, "svg", "mapping/")
