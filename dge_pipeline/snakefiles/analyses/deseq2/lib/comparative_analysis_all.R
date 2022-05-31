set.seed(123)
library(tidyr)
library(gtools)
library(openxlsx)
library(DESeq2)
library(reshape2)
library(VennDiagram)

#source("/homes/nhofmann/Virus_project/virus_pipeline/dge_pipeline/snakefiles/analyses/deseq2/lib/plot_heatmap.R")
#source("/homes/nhofmann/Virus_project/virus_pipeline/dge_pipeline/snakefiles/analyses/deseq2/lib/enrichment.R")
source("plot_heatmap.R")
source("enrichment.R")
source("STRINGdb.R")
source("/homes/nhofmann/Virus_project/virus_pipeline/dge_pipeline/snakefiles/analyses/deseq2/lib/venn.R")
source("/homes/nhofmann/Virus_project/virus_pipeline/dge_pipeline/snakefiles/analyses/deseq2/lib/upset.R")

# common subset of diff. expressed genes
# - separate up- and down-regulated genes
# - exclude MARV, LASV
# Heatmap
# Upset
# ORA
# STRING network
# “How many genes are regulated (up or down) during how may time points?”

output_folder <- "/vol/sfb1021/SFB1021_Virus/dge_analyses_antisense_new/analyses/host/deseq2/customized/"
if(!dir.exists(output_folder)){
  dir.create(output_folder)
}

virus.levels <- c("H1N1","H5N1","RVFV","SFSV","RSV","NiV","EBOV","MARV","LASV")

lfc.df.all <- Reduce(function(x,y)merge(x,y,by="SYMBOL"),lapply(names(res.list)[grep(".*h", names(res.list))], function(x){x.df <- data.frame(res.list[[x]]$SYMBOL,res.list[[x]]$log2FoldChange); colnames(x.df) <- c("SYMBOL",x); return(x.df)}))
rownames(lfc.df.all) <- lfc.df.all$SYMBOL
lfc.df.all <- lfc.df.all[,-1]
lfc.df.all <- lfc.df.all[,mixedorder(colnames(lfc.df.all))]
lfc.df.all <- lfc.df.all[,unlist(sapply(virus.levels, function(v){grep(v,colnames(lfc.df.all))}, simplify = F, USE.NAMES = F))]
#colnames(lfc.df.all) <- sub("3h","_3h",colnames(lfc.df.all))
#colnames(lfc.df.all) <- sub("6h","_6h",colnames(lfc.df.all))
#lfc.df <- lfc.df[,grep(".*h.*Mock", colnames(lfc.df))]
lfc.df <- lfc.df.all[, -grep("vs_.*BPL$", colnames(lfc.df.all))]
lfc.df.mod <- Reduce(function(x,y)merge(x,y,by="SYMBOL"),lapply(names(res.list)[-grep(".*BPL", names(res.list))], function(x){
  y <- res.list[[x]]
  y$log2FoldChange <- ifelse(y$padj<0.05 & apply(y[,grep("normalized", colnames(y))],1,max) >= 10, y$log2FoldChange, NA)
  x.df <- data.frame(y[,c("SYMBOL","log2FoldChange")])
  colnames(x.df) <- c("SYMBOL",x)
  return(x.df)
}))
rownames(lfc.df.mod) <- lfc.df.mod$SYMBOL
lfc.df.mod <- lfc.df.mod[, -1]

res.list <- res.list[mixedorder(names(res.list))]

# plot number of differentially expressed genes - sorted by customized order
LFC.cut <- 1
count.genes <- count.genes[-grep("BPL",count.genes$Time),]
p <- ggplot(count.genes[count.genes$LFC_cutoff==LFC.cut,], aes(y=Count, x=factor(Time, levels = c("3h","6h","12h","24h")), group=Direction, fill=factor(Direction, labels = c("Down","Up")))) + 
  geom_bar(stat = "identity", position = "stack", width = 0.8) + facet_wrap(~factor(Virus, levels = virus.levels), scales = "free_x") + 
  scale_x_discrete(labels=c("3h"="3","6h"="6","12h"="12","24h"="24")) + 
  scale_y_continuous(breaks = pretty(count.genes$Count[count.genes$LFC_cutoff==LFC.cut], n=5), labels = abs(pretty(count.genes$Count[count.genes$LFC_cutoff==LFC.cut], n=5))) + 
  geom_hline(yintercept = 0) + scale_fill_manual(values=c(Up="red", Down="blue"), guide = guide_legend(reverse=T)) + 
  xlab("Time after infection") + ylab("Number of genes") +
  theme(axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold", margin = margin(t=10,r=0,b=0,l=0)), 
        axis.title.y = element_text(family = "Helvetica", size = 28, face = "bold", margin = margin(t=0,r=10,b=0,l=0)),
        axis.text = element_text(family = "Helvetica", size = 30), axis.text.x = element_text(family = "Helvetica", angle = 0), 
        strip.text = element_text(family = "Helvetica", size = 30, face = "bold"),
        legend.text = element_text(family = "Helvetica", size = 28), legend.title = element_blank())
ggsave(paste0("DEG_count_LFC",LFC.cut,".pdf"), p, "pdf", output_folder, width = 15, height = 10)
system(paste0("inkscape -l ", output_folder, "DEG_count_LFC", LFC.cut,".svg ", output_folder, "DEG_count_LFC", LFC.cut, ".pdf"))

# Plot PCA
shape <- if(length(unique(conditiontable$time)) <= 6){ scales::shape_pal()(length(unique(conditiontable$time))) }else{ c(1:length(unique(conditiontable$time)))}
names(shape) <- unique(conditiontable$time)
#color.df <- read.table("/vol/sfb1021/SFB1021_Virus/colors.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
color.df <- read.table("Documents/Virus_project/analyses/host/deseq2_antisense/new/colors.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
color <- color.df[,2]
names(color) <- color.df[,1]
pca <- plotPCA(deseq.results.vst, intgroup = c("treatment", "time"), returnData = TRUE)
plot_PCA <- ggplot(pca, aes(PC1, PC2, color = factor(treatment, levels = c(virus.levels, "Mock")), shape = factor(time, levels = mixedsort(as.character(unique(conditiontable$time)))))) + 
  geom_point(size=3) + labs(color = "Infection", shape = "Time") + scale_shape_manual(values=shape) + guides(color=guide_legend(override.aes=list(fill=NA))) +
  theme(axis.title = element_text(family = "Helvetica", size = 20, face = "bold"), axis.text = element_text(family = "Helvetica", size = 16),  
        legend.text = element_text(family = "Helvetica", size = 16), legend.title = element_text(family = "Helvetica", size = 20, face = "bold"), legend.key=element_blank()) + 
  guides(color = guide_legend(order = 2), shape = guide_legend(order = 1))
if(exists("color")){
  plot_PCA <- plot_PCA + scale_colour_manual(values=color)
}
ggsave("PCA.pdf", plot = plot_PCA, device = "pdf", path = output_folder, width = 10, height = 8)
system(paste0("inkscape -l ", output_folder, "PCA.svg ", output_folder, "PCA.pdf"))

deseq.results.vst.mod <- deseq.results.vst[,!grepl("Mock", colnames(deseq.results.vst))]
for (time in unique(conditiontable$time[!conditiontable$treatment%in%c("Mock")])) {
  pca_time <- plotPCA(deseq.results.vst.mod[,grep(paste(virus.levels[!virus.levels %in% c("LASV")], time, collapse = "|", sep="_"), colnames(deseq.results.vst.mod))], intgroup = c("treatment"), returnData = TRUE)
  if(time!="BPL"){
    main <- paste(time, "p.i.")
  }else{
    main <- time
  }
  plot_PCA <- ggplot(pca_time, aes(PC1, PC2, color = factor(treatment, levels = virus.levels), shape = time)) + geom_point(size=5) + 
    labs(color = "Infection", shape = "Time") + ggtitle(main) + 
    theme(plot.title = element_text(family = "Helvetica", size = 20, face = "bold", hjust = 0.5), axis.title = element_text(family = "Helvetica", size = 20, face = "bold"), axis.text = element_text(family = "Helvetica", size = 16),  
          legend.text = element_text(family = "Helvetica", size = 16), legend.title = element_text(family = "Helvetica", size = 20, face = "bold"), legend.key=element_blank()) +
    scale_shape_manual(values=shape) + guides(color = guide_legend(order = 2), shape = "none") 
  if(exists("color")){
    plot_PCA <- plot_PCA + scale_colour_manual(values=color[!grepl("LASV|Mock",names(color))])
  }
  ggsave(paste0("PCA_", time, ".pdf"), plot = plot_PCA, device = "pdf", path = output_folder)
  system(paste0("inkscape -l ", output_folder, "PCA_", time, ".svg ", output_folder, "PCA_", time, ".pdf"))
  
}

pca_time <- plotPCA(deseq.results.vst[,grep("LASV", colnames(deseq.results.vst))], intgroup = c("treatment","time"), returnData = TRUE)
plot_PCA <- ggplot(pca_time, aes(PC1, PC2, color = factor(treatment, levels = virus.levels), shape = factor(time, levels = mixedsort(as.character(unique(conditiontable$time)))))) + geom_point(size=3) + labs(color = "Infection", shape = "Time") +
  theme(axis.title = element_text(family = "Helvetica", size = 20, face = "bold"), axis.text = element_text(family = "Helvetica", size = 16),  
        legend.text = element_text(family = "Helvetica", size = 16), legend.title = element_text(family = "Helvetica", size = 20, face = "bold"), legend.key=element_blank()) +
  scale_shape_manual(values=shape) + guides(color = guide_legend(order = 2), shape = guide_legend(order = 1))
if(exists("color")){
  plot_PCA <- plot_PCA + scale_colour_manual(values=color)
}

ggsave(paste0("PCA_LASV.pdf"), plot = plot_PCA, device = "pdf", path = output_folder)
system(paste0("inkscape -l ", output_folder, "PCA_LASV.svg ", output_folder, "PCA_LASV.pdf"))
svg(paste0(output_folder,"PCA_LASV.svg"), family = "Helvetica", width = 15, height = 10)
plot(plot_PCA)
dev.off()

# Violin plot
normalized.stack <- melt(countdata.normalized) 
normalized.stack$Var2 <- as.character(normalized.stack$Var2)
normalized.stack <- separate(normalized.stack, "Var2", c("Virus", "Time"), "_", remove=F, extra="merge")
normalized.stack$Virus <- sub("Mock.*","Mock",normalized.stack$Virus)
normalized.stack$Time[grep("MockMR",normalized.stack$Var2)] <- sub("_1","_3",normalized.stack$Time[grep("MockMR",normalized.stack$Var2)])
normalized.stack$Time[grep("MockMR",normalized.stack$Var2)] <- sub("_2","_4",normalized.stack$Time[grep("MockMR",normalized.stack$Var2)])
normalized.stack <- separate(normalized.stack, Time, c("hours","rep"), "_", remove = F)
plot.violin <- ggplot(normalized.stack, aes(factor(Time, levels = mixedsort(unique(Time))), log2(value), fill=Virus)) + geom_violin() + 
  facet_wrap(~factor(Virus, levels = c(virus.levels, "MockGI", "MockMR")), scales = "free_x") + #geom_boxplot(width=0.1) #stat_summary(fun=mean, geom="point", size=1) +
  labs(x="Time", y="log2(normalized counts)") + scale_fill_manual(values = color) +
  theme(axis.title = element_text(family = "Helvetica", size = 20, face = "bold"), axis.text = element_text(family = "Helvetica", size = 16, face = "bold"),
        axis.text.x=element_text(family = "Helvetica", angle = 90, hjust = 1.25, vjust = 0.5), 
        strip.text = element_text(family = "Helvetica", size = 30, face = "bold"),
        legend.position = "None", plot.margin = margin(t = 0.5, r = 1.1, b = 0.5, l = 0.5, "cm"))
ggsave("Violin_plot_counts_alt.pdf", plot.violin, "pdf", output_folder, width = 20, height = 9)
system(paste0("inkscape -l ", output_folder, "Violin_plot_counts_alt.svg ", output_folder, "Violin_plot_counts_alt.pdf"))
svg(paste0(output_folder,"Violin_plot_counts.svg"), family = "Helvetica", width = 15, height = 10)
plot(plot.violin)
dev.off()

plot.box <- ggplot(normalized.stack, aes(factor(Time, levels = mixedsort(unique(Time))), log2(value), fill=Virus)) + geom_boxplot() + 
  facet_wrap(~factor(Virus, levels = c(virus.levels, "MockGI", "MockMR")), scales = "free_x") +
  labs(x="Time", y="log2(normalized counts)") + scale_fill_manual(values = color) +
  theme(axis.title = element_text(family = "Helvetica", size = 20, face = "bold"), axis.text = element_text(family = "Helvetica", size = 16, face = "bold"), 
        axis.text.x=element_text(family = "Helvetica", angle=90, hjust = 1.25, vjust=0.5), 
        strip.text = element_text(family = "Helvetica", size = 30, face = "bold"),
        legend.position = "None", plot.margin = margin(t = 0.5, r = 0.8, b = 0.5, l = 0.5, "cm"))
ggsave("Boxplot_counts_alt.pdf", plot.box, "pdf", output_folder, width = 18, height = 9)
system(paste0("inkscape -l ", output_folder, "Boxplot_counts_alt.svg ", output_folder, "Boxplot_counts_alt.pdf"))
svg(paste0(output_folder,"Boxplot_counts.svg"), family = "Helvetica", width = 15, height = 10)
plot(plot.box)
dev.off()


LFC.cut <- 1
res.list.filter <- sapply(names(res.list)[-grep(".*BPL", names(res.list))], function(n){
  x <- res.list[[n]]
  return(x[which(apply(x[,grep("normalized", colnames(x))],1,max) >= 10 & x$padj < 0.05 & abs(x$log2FoldChange) > LFC.cut),])
})
names(res.list.filter) <- sub("_"," ",sub("_vs.*", "", names(res.list.filter)))
 

# “How many genes are regulated (up or down) during how may time points?”
# UpSet plot 
virus.dge.binary <- list2binary(lapply(res.list.filter, "[[", "SYMBOL"))#, paste0(out.dir,"/all_genes_binary.csv"))
for(virus in virus.levels){
  print(virus)
  virus.sub <- virus.dge.binary[,grep(virus,colnames(virus.dge.binary))]
  virus.sub <- virus.sub[rowSums(virus.sub)>0,]
  #svglite::svglite(paste0(output_folder, "UpSet_",virus,"_Mock.svg"), width = 14, height = 8)
  pdf(paste0(output_folder, "UpSet_",virus,"_Mock.pdf"), width = 15, height = 8, family = "Helvetica")
  #svg(paste0(output_folder, "UpSet_", virus, "_Mock.svg"), width = 14, height = 8, family = "Helvetica")
  print(upset(virus.sub, order.by = "freq", nsets = 4, nintersects = NA, text.scale = c(3,3,2,2,3,3), line.size = 1, point.size = 4,
              set_size.scale_max = round(max(colSums(virus.sub))+1300, -3)+100, set_size.show = T,
              sets = colnames(virus.sub[,colSums(virus.sub)>0])[mixedorder(colnames(virus.sub[,colSums(virus.sub)>0]), decreasing = T)], keep.order = TRUE), newpage=F)
  dev.off()
  system(paste0("inkscape -l ", output_folder, "UpSet_", virus, "_Mock.svg ", output_folder, "UpSet_", virus, "_Mock.pdf"))
}

# Common genes between viruses at any time point 
virus.levels <- c("H1N1","H5N1","RVFV","SFSV","RSV","NiV","EBOV")
virus.dge <- sapply(virus.levels,function(virus){unique(unlist(sapply(res.list.filter[grep(virus, names(res.list.filter))], function(x){return(as.character(x$SYMBOL))})))}, USE.NAMES = T)
out.dir <- paste0(output_folder,"/common_pattern/")
if(!dir.exists(out.dir)){
  dir.create(out.dir)
}
genes.common <- Reduce(intersect, virus.dge)
# UpSet plot of gene sets
pdf(paste0(out.dir, "/UpSet_minus_LASV_MARV.pdf"), width = 20, height = 8, family = "Helvetica")
#svg(paste0(out.dir, "/UpSet_minus_LASV_MARV_svg.svg"), width = 20, height = 8, family = "")
print(upset(fromList(virus.dge), nsets = 7, nintersects = 40, order.by = "freq", text.scale = c(2,2,1.5,1.5,2,2),
            point.size = 3, line.size = 1, number.angles = 0, set_size.show = T, set_size.scale_max = round(max(sapply(virus.dge, length))+500, -3)+100,
            queries = list(list(query = intersects, params = list(virus.levels), active = T, color = "red"))))
dev.off()
system(paste0("inkscape -l ", out.dir, "UpSet_minus_LASV_MARV.svg ", out.dir, "UpSet_minus_LASV_MARV.pdf"))

# plot heatmap of common gene set
virus.heat <- plotHeatmap(lfc.df.all[, -grep("MARV|LASV", colnames(lfc.df.all))], filename = paste0(out.dir,"/Heatmap_common_genes_LFC",LFC.cut,"_all.pdf"), 
                          row_subset = genes.common, 
                          colClust = F, clusterMethod = "ward.D2", legend.limit = 1, clrn = 1,
                          fontsize_row = 3.5, fontsize_col = 3.5, height = 7, border_col = NA)
virus.heat <- plotHeatmap(lfc.df[, -grep("MARV|LASV", colnames(lfc.df))], filename = paste0(out.dir,"/Heatmap_common_genes_LFC",LFC.cut,".pdf"), 
                          row_subset = genes.common, 
                          colClust = F, clusterMethod = "ward.D2", legend.cut = 1, clrn = 1,
                          fontsize_row = 3.5, fontsize_col = 3.5, height = 7, border_col = NA, family = "Helvetica")

max.lfc <- max(lfc.df[genes.common, -grep("HCV|MARV|LASV", colnames(lfc.df))])
min.lfc <- min(lfc.df[genes.common, -grep("HCV|MARV|LASV", colnames(lfc.df))])
genes.common.sum <- sign(lfc.df[genes.common, -grep("MARV|LASV|H1N1|H5N1", colnames(lfc.df))])
genes.common.list <- list("up" = rownames(genes.common.sum[rowSums(genes.common.sum)>0,]), "down"= rownames(genes.common.sum[rowSums(genes.common.sum)<0,]))

for(reg in names(genes.common.list)){
  virus.heat.cl <- plotHeatmap(lfc.df[, -grep("HCV|MARV|LASV", colnames(lfc.df))],
                               filename = paste0(out.dir,"/Heatmap_common_genes_LFC",LFC.cut,"_",reg,".pdf"), 
                               row_subset = genes.common.list[[reg]], colNames = F, cellwidth = 6.5, cellheight = 3.5,
                               colClust = F, clusterMethod = "ward.D2", clrn = 1, legend.limit.up = max.lfc, legend.limit.down = min.lfc, 
                               fontsize_row = 3.5, fontsize_col = 6, border_col = NA, column_title_rot = 0, family = "Helvetica")
  #ora <- calc_ora(genes.common.list[[reg]], filename = paste0("ORA_common_",reg), out.dir = paste0(out.dir, "/ORA"), GO = T, REACTOME = F, ont = c("CC","BP","MF"), 
   #               p.cut = 0.05, label.size = 30, legend.size = 25, legend.title.size = 20, width = 15, imagetype = "pdf")
}
write.xlsx(res.list[[1]][res.list[[1]]$SYMBOL %in% genes.common, c("SYMBOL","UNIPROT","GENENAME", "PATH")], paste0(out.dir,"/common_genes.xlsx"))
write.xlsx(genes.common, paste0(out.dir,"/common_genes.xlsx"))
virus.dge.BPL <- Reduce(function(x,y)merge(x,y,by="SYMBOL"),lapply(virus.levels[!virus.levels %in% c("HCV","LASV","MARV")],function(virus){
  tmp <- res.list.filter[[grep(paste0(virus,".*h.*BPL"), names(res.list.filter))]]
  tmp.df <- data.frame(genes.common, ifelse(genes.common %in% tmp$SYMBOL, "yes", "-"))
  colnames(tmp.df) <- c("SYMBOL", virus)
  return(tmp.df)
}))
write.table(virus.dge.BPL, paste0(out.dir,"/Common_genes_vs_inaktive.tsv"), sep = "\t", row.names = F)
rownames(virus.dge.BPL) <- virus.dge.BPL$SYMBOL
virus.dge.BPL <- virus.dge.BPL[,-1]
pheatmap::pheatmap(virus.dge.BPL)

# over-representation analysis
ora <- calc_ora(genes.common, filename = "ORA_common", out.dir = paste0(out.dir, "/ORA/"), GO = T, REACTOME = T, ont = c("CC","BP","MF"), 
                p.cut = 0.05, label.size = 30, legend.size = 25, legend.title.size = 20, imagetype = "pdf", width = 22, height = 17, family = "Helvetica")

# network analysis using STRING
string_ppi(string_db, gene.df = data.frame("SYMBOL"=genes.common), filename = "common_genes", out.dir = paste0(out.dir, "/STRING"), link = F, 
           cluster = T, required_score = 400)


# Transcription factor binding site analysis
#Convert gtf to GRanges object
gtf.human <- "Documents/Virus_project/genome/genes.gtf" #"/vol/sfb1021/SFB1021_Virus/genomes/hg38/genes.gtf"
human.GRanges <- getTranscriptsfromGFF(gtf.human)
human.GRanges <- keepStandardChromosomes(human.GRanges, pruning.mode = "coarse")
genes.bg <- rownames(countdata[rowSums(countdata)>0,])
genes.bg <- genes.bg[!genes.bg %in% genes.common]

# TFBS analysis with MEME
path2meme <- "/home/nina/meme/meme-5.4.1/bin/"
motifDB <- "/home/nina/Documents/Virus_project/motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant.meme"
motifDB <- "/home/nina/Documents/Virus_project/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme"
fasta.file <- "Documents/Virus_project/genome/genome.fa" #"/vol/sfb1021/SFB1021_Virus/genomes/hg38/genome.fa"
out.dir <- paste0(out.dir, "/TFBS/")
promotor_prim <- "promotor.common.fa"
promotor_back <- "promotor.bg.fa"
writeGRanges2Fasta(human.GRanges, out.dir, promotor_prim, promotor_back, fasta.file, genes.common, genes.bg, upstream = 1000, downstream = 0)
meme.res <- meme(path2meme, outdir = out.dir, paste0(out.dir,promotor_prim), paste0(out.dir,promotor_back), nmotif = 10, alph = "dna", objfun = "de")
tomtom(path2meme, motifDB, motif_file = meme.res, outdir = out.dir)
ame(path2meme, motifDB, paste0(out.dir,promotor_prim), paste0(out.dir,promotor_back), out.dir)

# Compare Mock and BPL
LFC.cut <- 1
res.list.filter.24h <- sapply(names(res.list)[grep("24h", names(res.list))], function(n){
  x <- res.list[[n]]
  return(x[which(apply(x[,grep("normalized", colnames(x))],1,max) >= 10 & x$padj < 0.05 & abs(x$log2FoldChange) > LFC.cut),])
})
names(res.list.filter.24h) <- sub("_BPL",":BPL_24h",sub("24h_vs_", "vs_", names(res.list.filter.24h)))
sapply(virus.levels, function(v){
  venn.plot <- venn.diagram(lapply(res.list.filter.24h[grep(v,names(res.list.filter.24h))],"[[", "SYMBOL"), filename = NULL, fontfamily = "Helvetica", cat.fontfamily = "Helvetica",
                            cex = 2, cat.cex = 1.5, margin = 0.1, cat.dist = c(0.125,0.125), ext.text = T, ext.percent = c(0.2,0.2,0.2), disable.logging = T, euler.d = T, scaled = T)
  ggsave(paste0("Venn_24h_",v,".pdf"), venn.plot, "pdf", output_folder)
  system(paste0("inkscape -l ", output_folder, "Venn_24h_",v,".svg ", output_folder, "Venn_24h_",v,".pdf"))
  #svg(paste0(output_folder,"Venn_24h_",v,".svg"), family = "Helvetica", width = 15, height = 10)
  #plot(venn.plot)
  #dev.off()
})

### Compare 2 groups of viruses ###

# Compare two sets and find set specific and common values
compare_geneset <- function(set1.list, set2.list, set1.name, set2.name){
  # get genes common in set1 or set2
  set1 <- Reduce(intersect, set1.list)
  set2 <- Reduce(intersect, set2.list)
  # get all genes that are in at least one sample of set1 or set2 
  set1.all <- unique(unlist(set1.list))
  set2.all <- unique(unlist(set2.list))
  # calculate genes that are specific for set1 or set2
  set1.spec <- set1[!set1 %in% set2.all]
  set2.spec <- set2[!set2 %in% set1.all]
  spec.list <- list(set1.spec, set2.spec, intersect(set1,set2))
  names(spec.list) <- c(set1.name, set2.name, paste(set1.name, set2.name, sep = "_"))
  return(spec.list)
}

# comparison of extremly pathogenic vs low pathogenic viruses
# extremly pathogenic: EBOV, NIV
# low pathogenic: 229E, SFSV
out.dir <- paste0(output_folder, "/weak_vs_high")
high <- virus.levels[grep("EBOV|NiV", virus.levels)]
low <- virus.levels[grep("HCoV-229E|SFSV", virus.levels)]
lfc.sub <- lfc.df[,grep(paste(high, low, sep = "|", collapse = "|"), colnames(lfc.df))]
lfc.sub <- lfc.sub[,unlist(sapply(c(high,low), function(v){grep(v,colnames(lfc.sub))}, simplify = F, USE.NAMES = F))]

venn.intersect <- draw_venn(virus.dge[grep(paste(high, low, sep = "|", collapse = "|"), names(virus.dge))], 
                            out_name = paste(out.dir, paste(high, low, sep = "_", collapse = "_"), sep = "/"), 
                            imagetype = "pdf", fill = color[names(virus.dge[grep(paste(high, low, sep = "|", collapse = "|"), names(virus.dge))])], margin = 0.1)
# venn.intersect <- venn.intersect[grep(paste(paste(low,collapse = ":"),paste(high,collapse = ":"),names(venn.intersect)[grep(".+:.+:.+:.+",names(venn.intersect))],sep = "$|^"), names(venn.intersect))]
list.intersect <- compare_geneset(set1.list = virus.dge[grep(paste0(high,collapse = "|"),names(virus.dge))], set1.name = "high",
                                  set2.list = virus.dge[grep(paste0(low,collapse = "|"),names(virus.dge))], set2.name = "low")
#lapply(names(venn.intersect), function(x){
sapply(names(list.intersect), function(x){ 
  virus.heat <- plotHeatmap(lfc.sub, filename = paste0(out.dir,"/Heatmap_",x,"_LFC",LFC.cut,".pdf"), 
                            row_subset = list.intersect[[x]], colNames = T,
                            colClust = F, clusterMethod = "ward.D2", legend.cut = 1, clrn = 1,
                            fontsize_row = 3.5, fontsize_col = 14, height = 7, border_col = NA, rowNames = F)
})
list.intersect.ENTREZ <- lapply(venn.intersect, function(x)bitr(x, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID)
calc_compareCluster(dataset = list.intersect, filename = "Compare_high_low_pathogenic_all", out.dir = out.dir, GO = F, 
                    ont = c("CC","BP","MF"), label.size = 20, w = 30, REACTOME = F, keytype = "SYMBOL", legend.size = 15)
calc_ora(c(list.intersect.ENTREZ$`CoV229E:SFSV`), filename = "ORA_low_specific", out.dir = out.dir, GO = F, ont = "BP", 
         p.cut = 0.05, label.size = 20, keytype = "ENTREZID")

# negative sense vs positive sense 
# (+)sense: CoV229E, MERS, HCV
# (-)sense: Infuenza, EBOV, MARV, NIV, RSV, (RVFV, SFSV)
out.dir <- paste0(output_folder, "positive_vs_negative")
positive <- virus.levels[grep("HCoV-229E|MERS-CoV", virus.levels)]
negative <- virus.levels[grep("H1N1|H5N1|EBOV|NiV|RSV|RVFV|SFSV", virus.levels)]
lfc.sub <- lfc.df[,grep(paste(positive, negative, sep = "|", collapse = "|"), colnames(lfc.df))]
lfc.sub <- lfc.sub[,unlist(sapply(c(positive,negative), function(v){grep(v,colnames(lfc.sub))}, simplify = F, USE.NAMES = F))]
annCol <- data.frame("Sense"=ifelse(grepl(paste0(positive,collapse = "|"),colnames(lfc.sub)),"positive-sense","negative-sense"), row.names = colnames(lfc.sub))
my_color <- list(Sense=c("positive-sense"="dark blue", "negative-sense"="light blue"))

svglite::svglite(paste(out.dir, "UpSet.svg", sep = "/"), width = 14, height = 8)
print(upset(fromList(virus.dge[grep(paste(c(positive, negative), collapse = "|"), names(virus.dge))]), keep.order = T,
            intersections = list(list(positive), list(negative)), text.scale = c(1.75,1.75,1.5,1.5,1.5,1.75)))
dev.off()

list.intersect <- compare_geneset(set1.list = virus.dge[grep(paste0(positive,collapse = "|"),names(virus.dge))], set1.name = "positive",
                                  set2.list = virus.dge[grep(paste0(negative,collapse = "|"),names(virus.dge))], set2.name = "negative")
sapply(names(list.intersect), function(x){ 
  virus.heat <- plotHeatmap(lfc.sub, filename = paste0(out.dir,"/Heatmap_",x,"_LFC",LFC.cut,".pdf"), 
                            row_subset = list.intersect[[x]], colNames = F,
                            colClust = F, clusterMethod = "ward.D2", legend.limit = 1, clrn = 1,
                            fontsize_row = 4.5, fontsize_col = 6, border_col = NA)
})
list.intersect.ENTREZ <- lapply(list.intersect, function(x)bitr(x, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID)
cc <- compareCluster(list.intersect.ENTREZ, fun = "enrichKEGG")
dotplot(cc)
