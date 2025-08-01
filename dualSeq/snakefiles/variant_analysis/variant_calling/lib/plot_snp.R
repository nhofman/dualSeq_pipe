library(vcfR)
library(ggplot2)
library(tidyr)
library(gtools)
library(Gviz)
library(lattice)

options(ucscChromosomeNames=FALSE)
#vcf.files <- list.files("Documents/Virus_project/variant_calling/", pattern = ".*h.*.vcf", full.names = T)
args <- commandArgs()
vcf <- args[match("--vcf", args) + 1]
gff_file <- args[match("--gff", args) + 1]
outdir <- args[match("--out", args) + 1]
# read vcf
vcf.list <- sapply(vcf.files, function(f){
  vcf <- read.vcfR(f)
  vcf.df <- cbind((vcf@fix),data.frame(extract_info_tidy(vcf)))
  #vcf.filter.df <- vcf.df[vcf.df$AF>0.01,]
  vcf.split <- vcf.df %>% group_split(CHROM)
  names(vcf.split) <- sapply(vcf.split, function(x) x[[1,"CHROM"]])
  return(vcf.split)
})
names(vcf.list) <-  gsub("\\..*","",basename(names(vcf.list)))

# create SNP dataframe
df.SNP <- Reduce(rbind, sapply(names(vcf.list), function(n){
  print(n)
  tmp.df <- Reduce(rbind, sapply(names(vcf.list[[n]]), function(x){
    tmp <- data.frame("Sample"=character(), "Chrom"=character(), "AF"=numeric(), "POS"=character(), "DP"=character(), "SB"=character(), "INDEL"=character())
    if(x %in% names(vcf.list[[n]])){
      tmp <- cbind("Sample"=n, vcf.list[[n]][[x]])
      #tmp <- data.frame("Sample"=n, "Chrom"=x, "AF"=vcf.list[[n]][[x]]$AF, "POS"=vcf.list[[n]][[x]]$POS, "DP"=vcf.list[[n]][[x]]$DP, 
      #                  "SNP"=paste(vcf.list[[n]][[x]]$CHROM, vcf.list[[n]][[x]]$POS, vcf.list[[n]][[x]]$REF, vcf.list[[n]][[x]]$ALT, sep = "_"),
      #                  "SB"=vcf.list[[n]][[x]]$SB, "INDEL"=vcf.list[[n]][[x]]$INDEL)
      tmp$INDEL <- ifelse(tmp$INDEL, ifelse(nchar(tmp$REF) > nchar(tmp$ALT), "Deletion", "Insertion"),  "SNP")
      #tmp.no <- rbind(tmp.no, data.frame("Sample"=n, "Chrom"=x, "AF"=vcf.list.nofilter[[n]][[x]]$DP, "Dedup"="no", "Type"="DP"))
    }
    return(tmp)
  }, simplify = F))
  return(tmp.df)
}, simplify = F))
df.SNP <- separate(df.SNP, "Sample", c("Virus", "Time"), "_", F, extra = "merge")
df.SNP$Time <- factor(df.SNP$Time, levels = gtools::mixedsort(unique(df.SNP$Time)))
#df.SNP$INDEL <- ifelse(df.SNP$INDEL, "INDEL", "SNP")
df.SNP$mutation <- paste(df.SNP$CHROM, df.SNP$POS, df.SNP$REF, df.SNP$ALT, sep = "_")

plot_theme <- theme(text = element_text(face = "plain"), line = element_line(linewidth = 0.25),
                    axis.line.x.top = element_blank(), axis.line.x.bottom = element_line(color = "black", linewidth = 0.25),
                    axis.line.y.right = element_blank(), axis.line.y.left = element_line(color = "black", linewidth = 0.25),
                    panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "gray58"), 
                    panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
                    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t=8,r=0,b=0,l=0)), 
                    axis.title.y = element_text(size = 12, face = "bold", margin = margin(t=0,r=8,b=0,l=0)),
                    axis.text = element_text(size = 10), axis.text.x = element_text(angle = 0), 
                    strip.text = element_text(size = 12, face = "bold"), strip.background = element_rect(fill = NA, color = NA), #panel.spacing = unit(2, "lines"),
                    legend.text = element_text(size = 10, face = "plain"), legend.title = element_blank())

plot_count <- ggplot(df.SNP, aes(Time, group = INDEL, fill = INDEL)) + geom_bar(position = "stack", color = "black") + 
  facet_wrap(~Virus) + xlab("Time_Replicate") + ylab("Number of SNP/INDEL") +
  scale_fill_manual(values = c("SNP" = "darkblue", "Insertion" = "grey70", "Deletion" = "lightblue")) + plot_theme
ggplot2::ggsave("Count_SNP.pdf", plot_count, "pdf", outdir, width = 16, height = 8)

plot_AF <- ggplot(df.SNP, aes(color=Time, x=AF)) + geom_freqpoly(binwidth=0.00025) + facet_wrap(~Virus) + xlim(c(0,1)) + 
  xlab("Frequency") + ylab("Count") + scale_color_brewer(palette = "Paired") + plot_theme
ggplot2::ggsave("Frequency_density.pdf", plot_AF, "pdf", out.dir, width = 15, height = 8)

plot_AF <- ggplot(df.SNP, aes(color=Time, x=AF)) + geom_freqpoly(binwidth=0.00025) + facet_wrap(~Virus) + xlim(c(0,0.01)) + 
  xlab("Frequency") + ylab("Count") + scale_color_brewer(palette = "Paired") + plot_theme
ggplot2::ggsave("Frequency_density_zoom.pdf", plot_AF, "pdf", out.dir, width = 15, height = 8)