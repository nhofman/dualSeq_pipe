library(vcfR)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gtools)
library(gggenomes)
library(egg)

vcf.files <- list.files("~/Documents/Virus_project/variant_calling_new/vcf", ".*[1|2].vcf", full.names = T)
vcf.list <- sapply(vcf.files, function(f){
  vcf <- read.vcfR(f)
  vcf.df <- cbind((vcf@fix),data.frame(extract_info_tidy(vcf)))
  #vcf.filter.df <- vcf.df[vcf.df$AF>0.01,]
  vcf.split <- vcf.df %>% group_split(CHROM)
  names(vcf.split) <- sapply(vcf.split, function(x) x[[1,"CHROM"]])
  return(vcf.split)
})
names(vcf.list) <-  gsub("\\..*","",basename(names(vcf.list)))

gff_file <- paste0("Documents/Virus_project/gff/", virus, ".gff")
gff <- read.table(gff_file, sep = "\t", comment.char = "#", stringsAsFactors = F, quote = "\"")
colnames(gff) <- c("chromosome","source","feature","start","end","score","strand","phase","attributes")
gff <- gff[-1,]
i <- 1
for(l in lapply(as.character(gff$attributes), FUN = strsplit, split=";", fixed=TRUE)){
  for(att in l[[1]]){
    if(grepl("Name", att)){
      gff$id[i] <- strsplit(att, "=", fixed=TRUE)[[1]][2]
    }
  }
  i <- i+1
}
gff$length <- gff$end - gff$start + 1
gff$seq_id <- "BPL"

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

out.dir <- "Documents/Virus_project/variant_calling_new/plots/"

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

df.SNP.total <- df.SNP %>% group_by(Sample) %>% summarise("Count" = n())
df.SNP.total <- merge(df.SNP.total, counts_noSplice.total)
df.SNP.total$perc <- df.SNP.total$Count/df.SNP.total$count_noSplice
df.SNP.total <- separate(df.SNP.total, "Sample", c("Virus", "Time_Rep"), "_", F, extra = "merge")
plot_SNP <- ggplot(df.SNP.total, aes(x=factor(Time_Rep, levels = mixedsort(unique(Time_Rep))), y=perc)) + geom_col() + 
  facet_wrap(~Virus) + xlab("Time_Rep") + ylab("# of SNP/# of viral Reads") + plot_theme
ggsave("Count_SNP_perc.pdf", plot_SNP, "pdf", "Documents/Virus_project/variant_calling_new/plots/", width = 14, height = 8)

ggplot(df.SNP, aes(x=Time, y=AF, group=mutation)) + 
  geom_point(show.legend = F) + geom_line() + ylim(c(0,0.01)) +
  facet_wrap(~Virus, scales = "free_x") + plot_theme
ggplot(df.SNP, aes(x=Time, y=DP)) + geom_boxplot() + facet_wrap(~Virus)
ggplot(df.SNP, aes(x=Time, y=AF)) + geom_violin() + facet_wrap(~Virus) + ylim(c(0,0.01))
plot_AF <- ggplot(df.SNP, aes(color=Time, x=AF)) + geom_freqpoly(binwidth=0.00025) + facet_wrap(~Virus) + xlim(c(0,0.01)) + 
  xlab("Frequency") + ylab("Count") + scale_color_brewer(palette = "Paired") + plot_theme
ggplot2::ggsave("Frequency_density_zoom.pdf", plot_AF, "pdf", out.dir, width = 15, height = 8)
plot_count <- ggplot(df.SNP, aes(Time, group = INDEL, fill = INDEL)) + geom_bar(position = "stack", color = "black") + 
  facet_wrap(~Virus) + xlab("Time_Replicate") + ylab("Number of SNP/INDEL") +
  scale_fill_manual(values = c("SNP" = "darkblue", "Insertion" = "grey70", "Deletion" = "lightblue")) + plot_theme
ggplot2::ggsave("Count_SNP.pdf", plot_count, "pdf", out.dir, width = 16, height = 8)

sapply(unique(df.SNP$Virus), function(x){
  #tmp <- as.data.frame(df.SNP[df.SNP$Virus==x,] %>% dplyr::select(Time, SNP) %>% dplyr::count(SNP, Time) %>% 
  #  pivot_wider(names_from = Time, values_from = n, values_fill = list(n = 0)))
  print(x)
  tmp <- as.data.frame(df.SNP[df.SNP$Virus==x,] %>% dplyr::select(Time, SNP, AF) %>% 
                         pivot_wider(names_from = Time, values_from = AF, values_fill = list(AF = NA)) %>% 
                         mutate(AvgFreq=rowMeans(pick(c(-1)), na.rm=T)) %>% 
                         mutate(across(2:last_col(1), ~ifelse(is.na(.), 0, 1))))
  #row.names(tmp) <- tmp$SNP
  #tmp <- tmp[,-1]
  tmp <- tmp[rowSums(tmp[,2:(ncol(tmp)-1)])>1,]
  n_intersect <- nrow(count(tmp[,2:(ncol(tmp)-1)], pick(everything())))
  pdf(paste0(out.dir, x, "_UpSet_new_BPL.pdf"), width = n_intersect/4, height = 8, family = "Helvetica")
  print(upset(tmp, mixedsort(colnames(tmp)[-c(1, length(colnames(tmp)))]), base_annotations = list(
    'intersection_size'=intersection_size(text = list(size = 3), mapping = aes(width = 1))),
    annotations = list("AvgFreq"=upset_annotate('AvgFreq', geom_boxplot(na.rm = TRUE))), 
    width_ratio = 0.1, sort_sets=FALSE), newpage=F)
  dev.off()
})

for(virus in virus.levels){
  # gggenomes package
  gbk <- read_gbk(paste0("Documents/Virus_project/gff/feature_viewer/genbank_mod/", virus, ".gb"))
  #gff.ggg <- read_gff3(paste0("Documents/Virus_project/gff/", virus, ".gff"))
  gbk$length <- gbk$end - gbk$start + 1 
  colnames(gbk)[1] <- "chrom"
  gbk$seq_id <- "ann"
  #gbk <- gbk[rep(seq_len(nrow(gbk)), 5),]
  #gbk$seq_id <- rep(unique(df.SNP.filter$seq_id), each = nrow(gbk)/5)
  
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
    gff.genes <- gbk[gbk$type %in% c("gene") & gbk$chrom == chrom,]
    gff.genes$seq_id <- gff.genes$chrom
    gff.genes$type <- "CDS" 
    i <- sapply(gff.genes$introns, function(x){0 %in% x}, simplify = T)
    if(TRUE %in% i){
      for(j in which(i)){
        gff.genes[[j,"introns"]] <- list(NULL)
      }
    }
    gff.seq <- gbk[gbk$type == "region" & gbk$start == 1 & gbk$chrom == chrom,]
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
    
    gg.genome <- gggenomes(genes = gff.genes, seqs = gff.seq) + geom_gene(data = genes(.gene_types = c("CDS")), aes(fill = name), position = "pile") + 
      geom_gene_text(data = genes(.gene_types = c("CDS")), aes(label = name), size = 3, vjust = -0.75, hjust = 0.5, angle = 0, position = "pile", nudge_y = 0) + 
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

gg.gene <- gggenomes(genes = gbk[gbk$type %in% c("gene", "CDS", "region"),]) + 
  #geom_seq() + 
  geom_gene(data = genes(.gene_types = c("gene","CDS", "region")), position = position_strandpile(grouped = TRUE, gap = 1)) + 
  scale_x_bp(expand = c(0.05,0.05), accuracy = NULL) +  #, limits = c(0, max(gbk$end)))  + 
  facet_grid(~chrom, scales = "free_x") +
  geom_gene_label(aes(label = gene), size = 3, angle = , nudge_y = 0, hjust = 0.5, vjust = 0) +
  theme(strip.background = element_blank(), strip.text = element_blank(), panel.spacing = unit(2, "lines"),
        text = element_text(face = "plain"), line = element_line(linewidth = 0.25),
        axis.title.y = element_blank(), axis.line.y.left = element_line(color = "black", linewidth = 0.25),
        axis.text.x = element_text(size = 10))

pdf(paste(out.dir, "SNP_Pos_", virus, "_filter.pdf"), width = 8*length(unique(gbk$chrom)), height = 10)
egg::ggarrange(gg.snp.manual, gg.gene, heights = c(10,2), widths = c(1,1), newpage = T)
dev.off()
