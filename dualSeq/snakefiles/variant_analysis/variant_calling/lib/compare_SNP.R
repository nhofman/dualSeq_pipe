library(vcfR)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gtools)
library(gggenomes)
library(egg)

# Parse arguments
args <- commandArgs(F)
out.dir <- args[match('--output', args) + 1] 
vcf.dir <- args[match('--vcf', args) + 1]
sample2gff.file <- args[match('--gff', args) + 1]

# Read vcf files
vcf.files <- list.files(vcf.dir, ".*.vcf", full.names = T)
print(vcf.files)
vcf.list <- sapply(vcf.files, function(f){
  print(f)
  vcf <- read.vcfR(f)
  vcf.df <- cbind((vcf@fix),data.frame(extract_info_tidy(vcf)))
  print(vcf.df)
  #vcf.filter.df <- vcf.df[vcf.df$AF>0.01,]
  vcf.split <- vcf.df %>% group_split(CHROM)
  names(vcf.split) <- sapply(vcf.split, function(x) x[[1,"CHROM"]])
  print(vcf.split)
  return(vcf.split)
}, simplify = F)
print(names(vcf.list))
names(vcf.list) <-  gsub("\\..*","",basename(names(vcf.list)))

# Create data frame of variants
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
colnames(sample2gff) <- c("Sample", "GFF", "FASTA")
sample2gff$Sample <- sub("_*", "", sample2gff$Sample)
#sample2gff <- separate(sample2gff, "Sample", c("Virus", "Time"), "_", F, extra = "merge")
sample2gff <- unique(sample2gff)

# Plot theme for all plots
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
save.image(paste0(out.dir, "/variant.RData"))

# Function to assign different tracks to overlapping features
assign_tracks <- function(feat.data, track.y, base){
  tracks <- c()
  ends <- c()
  for (i in 1:nrow(feat.data)) {
    pos.start <- feat.data[[i, "start"]]
    pos.end <- feat.data[[i, "end"]]
    placed <- FALSE
    for(e in seq_along(ends)){
      if(pos.start < ends[e]){
        tracks[i] <- tracks[i-1]+track.y
        placed <- TRUE
        break
      }
    }
    if(!placed){
      tracks[i] <- base
      ends <- c(ends, pos.end)
    }
  }
  tracks
}

# Plot positions of virus-specific variants 
for(virus in unique(df.SNP$Virus)){
  # gggenomes package
  #gbk <- read_gbk(paste0("Documents/Virus_project/gff/feature_viewer/genbank_mod/", virus, ".gb"))
  gff.df <- read_gff3(sample2gff[grep(virus, sample2gff$Sample), "GFF"][1])
  gff.df$length <- gff.df$end - gff.df$start + 1 
  colnames(gff.df)[1] <- "chrom"
  gff.df$seq_id <- gff.df$chrom
  fa.seqs <- read_seqs(sample2gff[grep(virus, sample2gff$Sample), "FASTA"][1])
  
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
  #df.SNP.filter$CHROM <- sub("\\..*", "", df.SNP.filter$CHROM)
  for(chrom in unique(df.SNP.filter$CHROM)){
    gff.chrom <- gff.df[gff.df$chrom == chrom,]
    gff.chrom$seq_id <- gff.chrom$chrom
    gff.genes <- gff.chrom[gff.chrom$type == "gene",]
    gff.cds <- gff.chrom[gff.chrom$type == "CDS",]
    gene.track <- 1
    if(nrow(gff.genes)>0){
      gff.genes$track <- assign_tracks(gff.genes, 0.3, 1)
      gene.track <- max(gff.genes$track)
      if(!"name" %in% colnames(gff.genes)){
        gff.genes$name <- NA
      }
    }
    if(nrow(gff.cds)>0){
      gff.cds$track <- assign_tracks(gff.cds, 0.3, gene.track+0.3)
      if(!"name" %in% colnames(gff.cds)){
        gff.cds$name <- NA
      }
    }
    gff.seq <- fa.seqs[fa.seqs$seq_id == chrom,]
    
    gg.snp.manual <- ggplot(df.SNP.filter[df.SNP.filter$CHROM==chrom,], aes(x = POS, y = AF)) + 
      geom_point(aes(alpha = 1, color = Rep, shape = type), size = 1, alpha = 0.7) + 
      facet_grid(rows = vars(factor(seq_id, levels = mixedsort(unique(df.SNP.filter$seq_id)))), 
                 cols = vars(CHROM), scales = "free_x") + #scale_shape_discrete(na.translate = FALSE) +
      ylim(c(0,1)) + scale_x_continuous(expand = c(0.01,0.01), limits = c(0, max(gff.seq$length))) + 
      #scale_color_manual(values = c("darkorchid4", "gold2"), na.translate = FALSE) +
      scale_color_manual(values = viridis::viridis_pal()(length(unique(df.SNP.filter$Rep))), na.translate = FALSE) +
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
    
    gg.genome <- gggenomes(genes = rbind(gff.genes, gff.cds), seqs = gff.seq) + geom_seq(arrow = 1) +
      geom_gene(data = genes(.gene_types = "gene"), aes(y = track, fill = type)) + 
      geom_gene_text(data = genes(.gene_types = "gene"), aes(y = track, label = name), size = 2.5, vjust = -1, hjust = 0.5, angle = 0, nudge_y = 0, check_overlap = F) + 
      geom_gene(data = genes(.gene_types = "CDS"), aes(y = track, fill = type)) + 
      geom_gene_text(data = genes(.gene_types = "CDS"), aes(y = track, label = name), size = 2.5, vjust = -1, hjust = 0.5, angle = 0, nudge_y = 0, check_overlap = F) + 
      scale_x_bp(accuracy = NULL, expand = c(0.01,0.01), limits = c(0, max(gff.seq$length))) + labs(fill = "Type") +
      theme(text = element_text(face = "plain"), line = element_line(linewidth = 0.25),
            axis.title.y = element_blank(), 
            axis.text.x = element_text(size = 10))
    SNP_filename <- paste0(out.dir, "/", virus, "_", chrom, ".pdf")
    pdf(SNP_filename, width = 10, height = 12, bg = "white")
    ggarrange(gg.snp.manual, gg.genome, heights = c(10,3), widths = c(1,1), newpage = F)
    dev.off()
    #system(paste0("inkscape -lo ", SNP_filename, ".svg ", SNP_filename))
  }
}
