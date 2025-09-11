library(vcfR)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gtools)
library(gggenomes)
library(egg)

#out.dir <- "Documents/Virus_project/variant_calling_new/plots/"
#vcf.dir <- "~/Documents/Virus_project/variant_calling_new/vcf"
#gff.file <- paste0("Documents/Virus_project/gff/", virus, ".gff")

# Parse arguments
args <- commandArgs(F)
out.dir <- args[match('--output', args) + 1] 
vcf.dir <- args[match('--vcf', args) + 1]
sample2gff.file <- args[match('--gff', args) + 1]
  
vcf.files <- list.files(vcf.dir, ".*[1|2].vcf", full.names = T)
vcf.list <- sapply(vcf.files, function(f){
  vcf <- read.vcfR(f)
  vcf.df <- cbind((vcf@fix),data.frame(extract_info_tidy(vcf)))
  #vcf.filter.df <- vcf.df[vcf.df$AF>0.01,]
  vcf.split <- vcf.df %>% group_split(CHROM)
  names(vcf.split) <- sapply(vcf.split, function(x) x[[1,"CHROM"]])
  return(vcf.split)
})
names(vcf.list) <-  gsub("\\..*","",basename(names(vcf.list)))

#vcf <- readVcf("Documents/Virus_project/variant_calling_new/vcf/EBOV_24h_1.vcf")

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

# Assign sample to gff
sample2gff <- read.table(sample2gff.file, header = F, sep = ",")
colnames(sample2gff) <- c("Sample", "GFF")
sample2gff$Sample <- sub("_*", "", sample2gff$Sample)
#sample2gff <- separate(sample2gff, "Sample", c("Virus", "Time"), "_", F, extra = "merge")
sample2gff <- unique(sample2gff)

# plot theme for all plots
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

# Plot overview of allele frequencies and number of variants 
plot_AF <- ggplot(df.SNP, aes(color=Time, x=AF)) + geom_freqpoly(binwidth=0.00025) + facet_wrap(~Virus) + xlim(c(0,1)) + 
  xlab("Frequency") + ylab("Count") + scale_color_brewer(palette = "Paired") + plot_theme
ggplot2::ggsave("Frequency_density.pdf", plot_AF, "pdf", out.dir, width = 15, height = 8)
plot_AF <- ggplot(df.SNP, aes(color=Time, x=AF)) + geom_freqpoly(binwidth=0.00025) + facet_wrap(~Virus) + xlim(c(0,0.01)) + 
  xlab("Frequency") + ylab("Count") + scale_color_brewer(palette = "Paired") + plot_theme
ggplot2::ggsave("Frequency_density_zoom.pdf", plot_AF, "pdf", out.dir, width = 15, height = 8)
plot_count <- ggplot(df.SNP, aes(Time, group = INDEL, fill = INDEL)) + geom_bar(position = "stack", color = "black") + 
  facet_wrap(~Virus) + xlab("Time_Replicate") + ylab("Number of SNP/INDEL") +
  scale_fill_manual(values = c("SNP" = "darkblue", "Insertion" = "grey70", "Deletion" = "lightblue")) + plot_theme
ggplot2::ggsave("Count_SNP.pdf", plot_count, "pdf", out.dir, width = 16, height = 8)

# Positional plots of virus-specific variants 
for(virus in unique(df.SNP$Virus)){
  # gggenomes package
  #gbk <- read_gbk(paste0("Documents/Virus_project/gff/feature_viewer/genbank_mod/", virus, ".gb"))
  gff.df <- read_gff3(sample2gff[virus, "GFF"])
  gff.df$length <- gff.df$end - gff.df$start + 1 
  colnames(gff.df)[1] <- "chrom"
  gff.df$seq_id <- "ann"
  #gff.df <- gff.df[rep(seq_len(nrow(gff.df)), 5),]
  #gff.df$seq_id <- rep(unique(df.SNP.filter$seq_id), each = nrow(gff.df)/5)
  
  df.SNP.filter <- df.SNP[df.SNP$Virus==virus,] # & df.SNP$AF>0.01
  df.SNP.filter <- separate(df.SNP.filter, "Time", c("Time_point", "Rep"), "_", F)
  df.SNP.filter$Rep <- paste("Rep", df.SNP.filter$Rep)
  df.SNP.filter$POS <- as.numeric(df.SNP.filter$POS)
  colnames(df.SNP.filter)[4] <- "seq_id"
  df.SNP.filter <- df.SNP.filter[mixedorder(df.SNP.filter$seq_id),]
  colnames(df.SNP.filter)[19] <- "type"
  df.SNP.filter$start <- df.SNP.filter$POS
  df.SNP.filter$end <- df.SNP.filter$POS + nchar(df.SNP.filter$ALT) - 1 
  df.SNP.filter$length <- nchar(df.SNP.filter$ALT)
  df.SNP.filter$CHROM <- sub("\\..*", "", df.SNP.filter$CHROM)
  for(chrom in unique(df.SNP.filter$CHROM)){
    gff.genes <- gff.df[gff.df$type %in% c("gene") & gff.df$chrom == chrom,]
    gff.genes$seq_id <- gff.genes$chrom
    i <- sapply(gff.genes$introns, function(x){0 %in% x}, simplify = T)
    if(TRUE %in% i){
      for(j in which(i)){
        gff.genes[[j,"introns"]] <- list(NULL)
      }
    }
    gff.seq <- gff.df[gff.df$type == "region" & gff.df$start == 1 & gff.df$chrom == chrom,]
    gff.seq$seq_id <- gff.seq$chrom
    
    gg.snp.manual <- ggplot(df.SNP.filter[df.SNP.filter$CHROM==chrom,], aes(x = POS, y = AF)) + 
      geom_point(aes(alpha = 1, color = Rep, shape = type), size = 1, alpha = 0.7) + 
      facet_grid(rows = vars(factor(seq_id, levels = mixedsort(unique(df.SNP.filter$seq_id)))), 
                 cols = vars(CHROM), scales = "free_x") + #scale_shape_discrete(na.translate = FALSE) +
      ylim(c(0,1)) + scale_x_continuous(expand = c(0.01,0.01), limits = c(0, max(gff.seq$end))) + 
      scale_color_manual(values = c("darkorchid4", "gold2"), na.translate = FALSE) +
      scale_shape_manual(values = c(15, 17, 19)) + ylab("Frequency") +
      theme(text = element_text(face = "plain"), line = element_line(linewidth = 0.25),
            axis.line.x.top = element_blank(), axis.line.x.bottom = element_line(color = "black", linewidth = 0.25),
            axis.line.y.right = element_blank(), axis.line.y.left = element_line(color = "black", linewidth = 0.25),
            panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "gray95"), 
            panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
            axis.title.x = element_blank(), #element_text(size = 12, face = "bold", margin = margin(t=8,r=0,b=0,l=0)), 
            axis.title.y = element_text(size = 12, face = "bold"), #, margin = margin(t=0,r=8,b=0,l=0)),
            axis.text = element_text(size = 10), axis.text.x = element_blank(), 
            strip.text = element_text(size = 12, face = "bold"), strip.background = element_rect(fill = NA, color = NA), 
            panel.spacing = unit(2, "lines"),
            legend.text = element_text(size = 10, face = "plain"), legend.title = element_blank())
    
    gg.genome <- gggenomes(genes = gff.genes, seqs = gff.seq) + geom_gene(data = genes(.gene_types = c("gene")), aes(fill = name), position = "pile") + 
      geom_gene_text(data = genes(.gene_types = c("gene")), aes(label = name), size = 3, vjust = -0.75, hjust = 0.5, angle = 0, position = "pile", nudge_y = 0) + 
      geom_seq(arrow = 1) + scale_x_bp(accuracy = NULL, expand = c(0.01,0.01), limits = c(0, max(gff.seq$end))) + 
      theme(legend.position = "none",
            text = element_text(face = "plain"), line = element_line(linewidth = 0.25),
            axis.title.y = element_blank(), 
            axis.text.x = element_text(size = 10))
    SNP_filename <- paste0(out.dir, virus, "_", chrom, ".pdf")
    pdf(SNP_filename, width = 8, height = 10, bg = "white")
    ggarrange(gg.snp.manual, gg.genome, heights = c(10,1), widths = c(1,1), newpage = F)
    dev.off()
    system(paste0("inkscape -lo ", SNP_filename, ".svg ", SNP_filename))
  }
}
