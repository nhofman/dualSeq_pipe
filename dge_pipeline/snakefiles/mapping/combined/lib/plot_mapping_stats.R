# Parse arguments
args <- commandArgs(T)
stat_host <- args[match('--stat-host', args) + 1]
stat_virus <- args[match('--stat-virus', args) + 1]
output_folder <- args[match('--output', args) + 1] 

for (package in c("ggplot2", "tidyr", "gtools", "dplyr")) {
  if (!(package %in% rownames(installed.packages()))) {
    stop(paste('Package "', package, '" not installed', sep=""))
  } else {
    print(paste("Import:", package))
    library(package, character.only=TRUE)
  }
}

# read mapping statistics for every file
stats.df <- Reduce(rbind,lapply(c(list.files(stat_host, ".*_stats.txt", full.names = T),list.files(stat_virus, ".*_stats.txt", full.names = T)), function(x){
  tmp <- read.table(x, header = T, sep = "\t")
  tmp$multimapped <- tmp$mapped - tmp$uniquely_mapped - tmp$mapped_to_both
  tmp$mapped <- tmp$mapped - tmp$mapped_to_both
  tmp$mapped_to_both <- sum(tmp$mapped_to_both)
  return(tmp)
}))
write.table(stats.df[,c(1,2,3,4,7,5,6)], paste0(output_folder, "stats/mapping_statistic.tsv"), sep = "\t", row.names = F)

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

ncol.facet <- ceiling(sqrt(length(unique(stats.t.mean$Virus[!grepl("Mock",stats.t.mean$Virus)]))))
nrow.facet <- ceiling(length(unique(stats.t.mean$Virus[!grepl("Mock",stats.t.mean$Virus)]))/ncol.facet)
times <- length(unique(stats.t.mean$Time[stats.t.mean$Virus==stats.t.mean$Virus[1]]))

# plot mapping statistic for host + virus samples
p <- ggplot(stats.t.mean[stats.t.mean$Category!="mapped" & !grepl("Mock",stats.t.mean$Virus),], aes(x=factor(Time, levels = unique(mixedsort(Time))), y=Count_percent, group=Category_2, fill=Category_2)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Virus, scales = "free_x", ncol = ncol.facet) +
  xlab("Time") + ylab("Ratio of mapped reads") + labs(fill="Category") +
  scale_fill_manual(labels = c("Mapped to Both", "Host: Multimapped", "Host: Uniquely mapped", "Virus: Multimapped", "Virus: Uniquely mapped"), 
                    values = c("red", "peachpuff2", "darkorange", "lightblue", "darkblue")) + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 15), strip.text = element_text(size = 15, face = "bold", vjust = 0.3),
        legend.text = element_text(size = 15), legend.title = element_text(size = 18), strip.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = NA), panel.grid = element_line(colour = "grey"))
ggsave("mapping_statistic.svg", p, "svg", paste0(output_folder,"stats/"), width = 0.6*times*(ncol.facet+1), height = 2.5*(nrow.facet+1))
ggsave("mapping_statistic.png", p, "png", paste0(output_folder,"stats/"), width = 0.6*times*(ncol.facet+1), height = 2.5*(nrow.facet+1))

ncol.facet <- ceiling(sqrt(length(unique(stats.t.mean$Virus[grepl("Mock",stats.t.mean$Virus)]))))
nrow.facet <- ceiling(length(unique(stats.t.mean$Virus[grepl("Mock",stats.t.mean$Virus)]))/ncol.facet)

# plot mapping statistic for control samples
p <- ggplot(stats.t.mean[!stats.t.mean$Category%in%c("mapped","both") & grepl("Mock",stats.t.mean$Virus),], aes(x=factor(Time, levels = unique(mixedsort(Time))), y=Count_percent, group=Category_2, fill=Category_2)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Virus, scales = "free_x") +
  xlab("Time") + ylab("Ratio of mapped reads") + labs(fill="Category") +
  scale_fill_manual(labels = c("Host: Multimapped", "Host: Uniquely mapped"), 
                    values = c("peachpuff2", "darkorange")) + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 15), strip.text = element_text(size = 15, face = "bold", vjust = 0.3),
        legend.text = element_text(size = 15), legend.title = element_text(size = 18), strip.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = NA), panel.grid = element_line(colour = "grey"))
ggsave("mapping_statistic_mock.svg", p, "svg", paste0(output_folder,"stats/"), width = 0.6*times*(ncol.facet+1), height = 2.5*(nrow.facet+1))
ggsave("mapping_statistic_mock.png", p, "png", paste0(output_folder,"stats/"), width = 0.6*times*(ncol.facet+1), height = 2.5*(nrow.facet+1))
