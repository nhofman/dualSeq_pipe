# Parse arguments
args <- commandArgs(T)
stat_host <- args[match('--stat-host', args) + 1]
stat_virus <- args[match('--stat-virus', args) + 1]
output_folder <- args[match('--output', args) + 1] 

for (package in c("ggplot2", "tidyr", "gtools", "dplyr", "openxlsx")) {
  if (!(package %in% rownames(installed.packages()))) {
    stop(paste('Package "', package, '" not installed', sep=""))
  } else {
    print(paste("Import:", package))
    library(package, character.only=TRUE)
  }
}

# read mapping statistics for every file
stats.df <- Reduce(rbind,lapply(c(list.files(stat_host, ".*_stats.csv", full.names = T),list.files(stat_virus, ".*_stats.csv", full.names = T)), function(x){
  tmp <- read.csv(x, header = T)
  tmp$multimapped <- tmp$mapped - tmp$uniquely_mapped
  return(tmp)
}))
write.csv(stats.df[,c(1,2,3,4,6,5)], paste0(output_folder, "stats/mapping_statistic.csv"), row.names = F)
write.xlsx(stats.df[,c(1,2,3,4,6,5)], paste0(output_folder, "stats/mapping_statistic.xlsx"))

# Replace virus total_reads with host total_reads, since the host count represents the total number of reads in the sample (Bowtie2 used the unmapped reads from STAR)
for(name in unique(stats.df$sample)){
  if(nrow(stats.df[stats.df$sample==name & stats.df$group=="virus", ])>0){
    stats.df[stats.df$sample==name & stats.df$group=="virus", "total_reads"] <- stats.df[stats.df$sample==name & stats.df$group=="host", "total_reads"]
  }
}

# transpose data frame for plotting
stats.t <- data.frame("Sample"=character(), "Total"=integer(), "Count"=integer(), "Organism"=character(), "Category"=character())
stats.t <- data.frame(Map(c, stats.t, data.frame(stats.df[,c(1,2,3,5)], "mapped")))
stats.t <- data.frame(Map(c, stats.t, data.frame(stats.df[,c(1,2,4,5)], "uniquely mapped")))
stats.t <- data.frame(Map(c, stats.t, data.frame(stats.df[,c(1,2,6,5)], "multimapped")))

stats.t <- separate(stats.t, "Sample", c("Virus","Time", "Rep"),"_",F,extra = "merge")
stats.t$Count_percent <- stats.t$Count/stats.t$Total
stats.t$Category_2 <- paste(stats.t$Organism, stats.t$Category, sep = ":")
stats.t$Time_Rep <- paste(stats.t$Time, stats.t$Rep, sep=":")
# calculate mean of Count_percent for all replicates
stats.t.mean <- stats.t %>% group_by(Virus,Time,Organism,Category) %>% summarise_at("Count_percent", "mean")
stats.t.mean$Category_2 <- paste(stats.t.mean$Organism, stats.t.mean$Category, sep = ":")

ncol.facet <- ceiling(sqrt(length(unique(stats.t.mean$Virus))))
ncol.facet <- ifelse(ncol.facet==1, 2, ncol.facet)
nrow.facet <- ceiling(length(unique(stats.t.mean$Virus))/ncol.facet)
nrow.facet <- ifelse(nrow.facet==1, 2, nrow.facet)
times <- length(unique(stats.t.mean$Time[stats.t.mean$Virus==stats.t.mean$Virus[1]]))
times <- ifelse(times<3, 3, times)

# plot mapping statistic for host + virus samples
p <- ggplot(stats.t.mean[stats.t.mean$Category!="mapped",], aes(x=factor(Time, levels = unique(mixedsort(Time))), y=Count_percent, group=Category_2, fill=Category_2)) +
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(~ Virus, scales = "free_x", ncol = ncol.facet) +
  xlab("Time") + ylab("Ratio of mapped reads") + labs(fill="Category") +
  scale_fill_manual(labels = c("Host: Multimapped", "Host: Uniquely mapped (proper pair)", "Virus: Multimapped", "Virus: Uniquely mapped (concordantly)"), 
                    values = c("peachpuff2", "darkorange", "lightblue", "darkblue")) + 
  theme(text = element_text(face = "bold"), line = element_line(linewidth = 0.25),
        axis.line.x.top = element_blank(), axis.line.x.bottom = element_line(color = "black", linewidth = 0.25),
        axis.line.y.right = element_blank(), axis.line.y.left = element_line(color = "black", linewidth = 0.25),
        panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "gray58"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        axis.title.x = element_text(size = 12, margin = margin(t=8,r=0,b=0,l=0)), 
        axis.title.y = element_text(size = 12, margin = margin(t=0,r=8,b=0,l=0)),
        axis.text = element_text(size = 10, face = "plain"), axis.text.x = element_text(angle = 0), 
        strip.text = element_text(size = 12), strip.background = element_rect(fill = NA, color = NA), #panel.spacing = unit(2, "lines"),
        legend.text = element_text(size = 10, face = "plain"), legend.title = element_blank())
ggsave("mapping_statistic.pdf", p, "pdf", paste0(output_folder,"stats/"), width = 0.6*times*(ncol.facet+1), height = 2.5*(nrow.facet+1))
ggsave("mapping_statistic.png", p, "png", paste0(output_folder,"stats/"), width = 0.6*times*(ncol.facet+1), height = 2.5*(nrow.facet+1))

