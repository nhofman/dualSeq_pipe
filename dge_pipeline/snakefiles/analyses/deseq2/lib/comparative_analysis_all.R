set.seed(123)
library(tidyr)
library(gtools)
library(openxlsx)
library(DESeq2)
library(reshape2)
library(gplots)
library(ggplot2)
source("plot_heatmap.R")
source("enrichment.R")
source("STRINGdb.R")
source("venn.R")
source("tfbs.R")

# common subset of diff. expressed genes
# - separate up- and down-regulated genes
# - exclude HCV, MARV, LASV
# Heatmap
# Upset
# ORA
# STRING network
# Transcriptionfactor analyses
# “How many genes are regulated (up or down) during how may time points?”

virus.levels <- c("H1N1","H5N1","RVFV","SFSV","RSV","NiV","EBOV","MARV","LASV")
virus.levels <- c("H1N1"="H1N1","H5N1"="H5N1","MERS"="MERS-CoV","CoV229E"="HCoV-229E","RVFV"="RVFV","SFSV"="SFSV","RSV"="RSV","NIV"="NiV","EBOV"="EBOV","MARV"="MARV","LASV"="LASV","HCV"="HCV")

lfc.df.all <- Reduce(function(x,y)merge(x,y,by="SYMBOL"),lapply(names(res.list), function(x){
  x.df <- res.list[[x]][,c("SYMBOL","log2FoldChange")]; colnames(x.df) <- c("SYMBOL",x); return(x.df)}))
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
count.genes <- separate(count.genes, "Sample", c("Virus","Time","Vs"), "_", F, extra = "merge")
count.genes$Vs <- sub("vs_","",count.genes$Vs)
count.genes$x <- ifelse(grepl("Mock", count.genes$Vs), count.genes$Time, "active")
count.genes <- count.genes[-grep("BPL",count.genes$Vs),]
count.genes.mod$Time <- ifelse(count.genes.mod$Control=="inactive", "BPL", count.genes.mod$Time)
p <- ggplot(count.genes[count.genes$LFC_cutoff==LFC.cut,], aes(y=Count, x=factor(x, levels = c("3h","6h","12h","24h","48h")), group=Direction, fill=factor(Direction, labels = c("Down","Up")))) + 
  geom_bar(stat = "identity", position = "stack") + facet_wrap(~factor(Virus, levels = virus.levels), scales = "free_x") + 
  scale_x_discrete(labels=c("3h"="3","6h"="6","12h"="12","24h"="24","48h"="48")) + 
  scale_y_continuous(breaks = pretty(count.genes$Count[count.genes$LFC_cutoff==LFC.cut], n=5), labels = abs(pretty(count.genes$Count[count.genes$LFC_cutoff==LFC.cut], n=5))) + 
  geom_hline(yintercept = 0) + scale_fill_manual(values=c(Up="red", Down="blue"), guide = guide_legend(reverse=F)) + 
  xlab("Time after infection") + ylab("Number of genes") +
  theme(axis.title.x = element_text(size = 28, face = "bold", margin = margin(t=10,r=0,b=0,l=0)), axis.title.y = element_text(size = 28, face = "bold", margin = margin(t=0,r=10,b=0,l=0)),
        axis.text = element_text(size = 30), axis.text.x = element_text(angle = 0), 
        strip.text = element_text(size = 30, face = "bold"),
        legend.text = element_text(size = 28), legend.title = element_blank(), 
        panel.background = element_rect(fill = NA), panel.grid = element_line(colour = "grey"))
ggsave(paste0("DEG_count_LFC",LFC.cut,"_wobg.pdf"), p, "pdf", output_folder, width = 12, height = 10)

# Plot PCA
shape <- if(length(unique(conditiontable$time)) <= 6){ scales::shape_pal()(length(unique(conditiontable$time))) }else{ c(1:length(unique(conditiontable$time)))}
names(shape) <- unique(conditiontable$time)
pca <- plotPCA(deseq.results.vst, intgroup = c("treatment", "time"), returnData = TRUE)
plot_PCA <- ggplot(pca, aes(PC1, PC2, color = factor(treatment, levels = c(names(virus.levels), "Mock"), labels = c(virus.levels, "Mock")), shape = factor(time, levels = mixedsort(as.character(unique(conditiontable$time)))))) + 
  geom_point(size=3) + labs(color = "Infection", shape = "Time") + scale_shape_manual(values=shape) + guides(color=guide_legend(override.aes=list(fill=NA))) +
  theme(axis.title = element_text(size = 20, face = "bold"), axis.text = element_text(size = 16),  
        legend.text = element_text(size = 16), legend.title = element_text(size = 20, face = "bold"), legend.key=element_blank()) + 
  guides(color = guide_legend(order = 2), shape = guide_legend(order = 1))
if(exists("color")){
  plot_PCA <- plot_PCA + scale_colour_manual(values=color)
}
ggsave("PCA.svg", plot = plot_PCA, device = "svg", path = output_folder, width = 10, height = 8)

for (time in unique(conditiontable$time[!conditiontable$treatment%in%c("HCV","LASV","Mock")])) {
  pca_time <- plotPCA(deseq.results.vst[,grep(paste(names(virus.levels)[!names(virus.levels) %in% c("HCV","LASV")], time, collapse = "|", sep="_"), colnames(deseq.results.vst))], intgroup = c("treatment"), returnData = TRUE)
  if(time!="BPL"){
    main <- paste(time, "p.i.")
  }else{
    main <- time
  }
  plot_PCA <- ggplot(pca_time, aes(PC1, PC2, color = factor(treatment, levels = names(virus.levels), labels = virus.levels), shape = time)) + geom_point(size=5) + 
  labs(color = "Infection", shape = "Time") + ggtitle(main) + 
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), axis.title = element_text(size = 20, face = "bold"), axis.text = element_text(size = 16),  
          legend.text = element_text(size = 16), legend.title = element_text(size = 20, face = "bold"), legend.key=element_blank()) +
  scale_shape_manual(values=shape) + guides(color = guide_legend(order = 2), shape = F) 
  if(exists("color")){
    plot_PCA <- plot_PCA + scale_colour_manual(values=color)
  }
  ggsave(paste0("PCA_", time, ".pdf"), plot = plot_PCA, device = "pdf", path = output_folder)
}

pca_time <- plotPCA(deseq.results.vst[,grep("HCV|LASV", colnames(deseq.results.vst))], intgroup = c("treatment","time"), returnData = TRUE)
plot_PCA <- ggplot(pca_time, aes(PC1, PC2, color = factor(treatment, levels = names(virus.levels), labels = virus.levels), shape = factor(time, levels = mixedsort(as.character(unique(conditiontable$time)))))) + geom_point(size=3) + labs(color = "Infection", shape = "Time") +
  theme(axis.title = element_text(size = 20, face = "bold"), axis.text = element_text(size = 16),  
        legend.text = element_text(size = 16), legend.title = element_text(size = 20, face = "bold"), legend.key=element_blank()) +
  scale_shape_manual(values=shape) + guides(color = guide_legend(order = 2), shape = guide_legend(order = 1))
if(exists("color")){
  plot_PCA <- plot_PCA + scale_colour_manual(values=color)
}
ggsave(paste0("PCA_HCV_LASV.pdf"), plot = plot_PCA, device = "pdf", path = output_folder)

# Violin plot
#normalized.stack <- stack(countdata.normalized[, grep("NiV", colnames(countdata.normalized))])
normalized.stack <- melt(countdata.normalized) 
normalized.stack$Var2 <- as.character(normalized.stack$Var2)
normalized.stack <- separate(normalized.stack, "Var2", c("Virus", "Time"), "_", remove=F, extra="merge")
plot.violin <- ggplot(normalized.stack, aes(factor(Time, levels = mixedsort(unique(Time))), log2(value), fill=Virus)) + geom_violin() + 
  facet_wrap(~factor(Virus, levels = c(virus.levels, "MockGI", "MockMR"))) + #geom_boxplot(width=0.1) #stat_summary(fun=mean, geom="point", size=1) +
  labs(x="Time", y="log2(normalized counts)") + scale_fill_manual(values = color) +
  theme(axis.title = element_text(size = 20, face = "bold"), axis.text = element_text(size = 16, face = "bold"), axis.text.x=element_text(angle=-30, hjust = 0.2, vjust=0.4), 
        strip.text = element_text(size = 30, face = "bold"),
        legend.position = "None", plot.margin = margin(t = 0.5, r = 1.1, b = 0.5, l = 0.5, "cm"))
ggsave("Violin_plot_counts.svg", plot.violin, "svg", output_folder, width = 16, height = 8)
ggsave("Violin_plot_counts.pdf", plot.violin, "pdf", output_folder, width = 18, height = 9)

plot.box <- ggplot(normalized.stack, aes(factor(Time, levels = mixedsort(unique(Time))), log2(value), fill=Virus)) + geom_boxplot() + 
  facet_wrap(~factor(Virus, levels = c(virus.levels, "MockGI", "MockMR"))) +
  labs(x="Time", y="log2(normalized counts)") + scale_fill_manual(values = color) +
  theme(axis.title = element_text(size = 20, face = "bold"), axis.text = element_text(size = 16, face = "bold"), axis.text.x=element_text(angle=-30, hjust = 0.2, vjust=0.4), 
        strip.text = element_text(size = 30, face = "bold"),
        legend.position = "None", plot.margin = margin(t = 0.5, r = 0.8, b = 0.5, l = 0.5, "cm"))
ggsave("Boxplot_counts.svg", plot.box, "svg", output_folder, width = 18, height = 9)
ggsave("Boxplot_counts.pdf", plot.box, "pdf", output_folder, width = 18, height = 9)
ggsave("Boxplot_counts.png", plot.box, "png", output_folder, width = 16, height = 8)
 
LFC.cut <- 1
res.list.filter <- sapply(names(res.list)[-grep(".*BPL", names(res.list))], function(n){
  x <- res.list[[n]]
  return(x[which(apply(x[,grep("normalized", colnames(x))],1,max) >= 10 & x$padj < 0.05 & abs(x$log2FoldChange) > LFC.cut),])
})
names(res.list.filter) <- sub("_"," ",sub("_vs.*", "", names(res.list.filter)))

# “How many genes are regulated (up or down) during how may time points?”
# UpSet plot or bar chart
virus.dge.binary <- list2binary(lapply(res.list.filter, "[[", "SYMBOL"))#, paste0(out.dir,"/all_genes_binary.csv"))
for(virus in virus.levels){
  print(virus)
  virus.sub <- virus.dge.binary[,grep(virus,colnames(virus.dge.binary))]
  virus.sub <- virus.sub[rowSums(virus.sub)>0,]
  #svglite::svglite(paste0(output_folder, "UpSet_",virus,"_Mock.svg"), width = 14, height = 8)
  pdf(paste0(output_folder, "UpSet_",virus,"_Mock.pdf"), width = 14, height = 8)
  print(upset(virus.sub, order.by = "freq", nsets = 4, nintersects = NA, text.scale = c(3,3,2,2,3,3), set_size.scale_max = round(max(colSums(virus.sub))+500, -3)+100,
              sets = colnames(virus.sub[,colSums(virus.sub)>0])[mixedorder(colnames(virus.sub[,colSums(virus.sub)>0]), decreasing = T)], keep.order = TRUE), newpage=F)
  dev.off()
}

# Common genes between viruses at any time point 
out.dir <- paste0(output_folder,"/common_pattern")
dir.create(out.dir)
virus.levels <- c("H1N1","H5N1","RVFV","SFSV","RSV","NiV","EBOV")
virus.dge <- sapply(unname(virus.levels),function(virus){unique(unlist(sapply(res.list.filter[grep(virus, names(res.list.filter))], function(x){return(as.character(x$SYMBOL))})))}, USE.NAMES = T)
genes.common <- Reduce(intersect, virus.dge)
# Venn diagram or better UpSet plot of gene sets
virus.dge.binary <- list2binary(virus.dge)#, paste0(out.dir,"/all_genes_binary.csv"))
svglite::svglite(paste0(out.dir, "/UpSet_minus_LASV_HCV_MARV.svg"), width = 20, height = 8)
pdf(paste0(out.dir, "/UpSet_minus_LASV_HCV_MARV.pdf"), width = 20, height = 8)
print(upset(virus.dge.binary, nsets = 9, nintersects = 40, order.by = "freq", text.scale = c(2,2,1.5,2,2,2), number.angles = 0,
            queries = list(list(query = intersects, params = as.list(virus.levels), active = T, color = "red"))))
dev.off()
upset_json(file = "UpSet_minus_LASV_HCV_MARV", out.dir = output_folder, name = "DGE Virus", start = 1, end = 12)

# plot heatmap of common gene set
virus.heat <- plotHeatmap(lfc.df.all[, -grep("HCV|MARV|LASV", colnames(lfc.df.all))], filename = paste0(out.dir,"/Heatmap_common_genes_LFC",LFC.cut,"_all.pdf"), 
                          row_subset = genes.common, 
                          colClust = F, clusterMethod = "ward.D2", legend.limit = 1, clrn = 1,
                          fontsize_row = 3.5, fontsize_col = 3.5, height = 7, border_col = NA)
virus.heat <- plotHeatmap(lfc.df[, -grep("HCV|MARV|LASV", colnames(lfc.df))], filename = paste0(out.dir,"/Heatmap_common_genes_LFC",LFC.cut,".pdf"), 
                          row_subset = genes.common, 
                          colClust = F, clusterMethod = "ward.D2", legend.limit = 1, clrn = 1,
                          fontsize_row = 3.5, fontsize_col = 3.5, height = 7, border_col = NA)

virus.heat <- plotHeatmap(lfc.df[, -grep("HCV|MARV|LASV", colnames(lfc.df))], 
                          row_subset = genes.common, plot.fig = F,
                          colClust = F, clusterMethod = "ward.D2", legend.cut = 1, clrn = 2,
                          fontsize_row = 3.5, fontsize_col = 3.5, height = 7, border_col = NA)
max.lfc <- max(lfc.df[genes.common, -grep("HCV|MARV|LASV", colnames(lfc.df))])
min.lfc <- min(lfc.df[genes.common, -grep("HCV|MARV|LASV", colnames(lfc.df))])
for(i in unique(virus.heat$row_cluster)){
  genes.cl <- names(virus.heat$row_cluster[virus.heat$row_cluster==i])
  if(sum(lfc.df[genes.cl,] > 0) > sum(lfc.df[genes.cl,] < 0)){
    reg <- "up"
  }else{
    reg <- "down"
  }
  virus.heat.cl <- plotHeatmap(lfc.df[, -grep("HCV|MARV|LASV", colnames(lfc.df))],
                            filename = paste0(out.dir,"/Heatmap_common_genes_LFC",LFC.cut,"_",reg,".pdf"), 
                            row_subset = genes.cl, colNames = F, cellwidth = 9.5, cellheight = 5.2,
                            colClust = F, clusterMethod = "ward.D2", clrn = 1, legend.limit.up = max.lfc, legend.limit.down = min.lfc, 
                            fontsize_row = 5.5, fontsize_col = 6, border_col = NA, column_title_rot = 0)
  ora <- calc_ora(genes.cl, filename = paste0("ORA_common_",reg), out.dir = paste0(out.dir, "/ORA"), GO = F, REACTOME = F, ont = c("CC","BP","MF"), 
                  p.cut = 0.05, label.size = 30, legend.size = 25, legend.title.size = 20, width = 15, imagetype = "pdf")
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
ora <- calc_ora(genes.common, filename = "ORA_common", out.dir = paste0(out.dir, "/ORA"), GO = T, REACTOME = F, ont = c("CC","BP","MF"), 
                p.cut = 0.05, label.size = 30, legend.size = 25, legend.title.size = 20)

# network analysis using STRING
string_ppi(string_db, gene.df = data.frame("SYMBOL"=genes.common), filename = "common_genes", out.dir = paste0(out.dir, "/STRING"), link = F, 
           cluster = T, required_score = 400)

# Transcription factor binding site analysis
#Convert gtf to GRanges object
gtf.human <- "/vol/sfb1021/SFB1021_Virus/genomes/hg38/genes.gtf"
human.GRanges <- getTranscriptsfromGFF(gtf.human)
human.GRanges <- keepStandardChromosomes(human.GRanges, pruning.mode = "coarse")

# TFBS analysis with MEME
path2meme <- "/home/nina/meme/bin/"
motifDB <- "/home/nina/Documents/Virus_project/motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant.meme"
motifDB <- "/home/nina/Documents/Virus_project/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme"
fasta.file <- "/vol/sfb1021/SFB1021_Virus/genomes/hg38/genome.fa"
out.dir <- paste0(out.dir, "/TFBS/")
promotor_prim <- "promotor.common.fa"
promotor_back <- "promotor.bg.fa"
writeGRanges2Fasta(human.GRanges, out.dir, promotor_prim, promotor_back, fasta.file, genes.common, upstream = 1000, downstream = 100)
meme.res <- meme(path2meme, outdir = out.dir, paste0(out.dir,promotor_prim), paste0(out.dir,promotor_back), nmotif = 10, alph = "dna", objfun = "de")
tomtom(path2meme, motifDB, motif_file = meme.res, outdir = out.dir)
ame(path2meme, motifDB, paste0(output_folder,"/common_pattern/TFBS/",promotor_prim), paste0(output_folder,"/common_pattern/TFBS/",promotor_back), out.dir)

# Compare Mock and BPL
LFC.cut <- 1
res.list.filter.24h <- sapply(names(res.list)[grep("24h", names(res.list))], function(n){
  x <- res.list[[n]]
  return(x[which(apply(x[,grep("normalized", colnames(x))],1,max) >= 10 & x$padj < 0.05 & abs(x$log2FoldChange) > LFC.cut),])
})
names(res.list.filter.24h) <- sub("_BPL",":BPL_24h",sub("24h_vs_", "vs_", names(res.list.filter.24h)))
intersect.24h <- Reduce(rbind, sapply(virus.levels, function(v){
  virus.list <- sapply(res.list.filter.24h[grep(v,names(res.list.filter.24h))],function(x){
    if(nrow(x)>0){
      x <- x[, "SYMBOL"]
    }
  })
  names(virus.list) <- gsub(v,"Virus",names(virus.list))
  virus.list <- virus.list[order(names(virus.list))]
  if(sum(sapply(virus.list, length)>0)){
    data <- attr(venn(virus.list, show.plot = F),"intersections")
  }else{
    data <- list(character(), character(), character())
    names(data) <- c(names(virus.list)[1], names(virus.list)[2], paste0(names(virus.list)[1], ":",names(virus.list)[2]))
  }
  #names(data) <- gsub("_vs_:*","",gsub(v,"",names(data)))
  data.df <- data.frame("Virus"=v,"Group"=names(data), "Count"=sapply(data, length), row.names = NULL)
  data.df$Sum <- sum(data.df$Count)
  data.df$Percentage <- data.df$Count/data.df$Sum
  return(data.df)
}, simplify = F))
intersect.24h$Group <- factor(intersect.24h$Group, levels = c("BPL_24h", "Mock_24h:BPL_24h", "Mock_24h"))
intersect.24h$Virus <- factor(intersect.24h$Virus, levels = virus.levels)
plot_count <- ggplot(intersect.24h, aes(x=Count, y=Virus, group=Group, fill=Group)) + geom_col() + xlab("# Genes") +
  scale_fill_manual(values = c("Virus_vs_Virus:BPL_24h"="grey40", "Virus_vs_Mock_24h"="steelblue", "Virus_vs_Mock_24h:Virus_vs_Virus:BPL_24h"="lightsalmon")) +
  scale_y_discrete(limits = rev) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 20), axis.text = element_text(size = 30), 
        axis.title.y = element_blank(), axis.title.x = element_text(size = 28, face = "bold"), 
        panel.background = element_rect(fill = NA), panel.grid = element_line(colour = "grey"))
ggsave("24h_Mock_vs_BPL_counts.pdf", plot_count, "pdf", output_folder, width = 14, height = 8)

plot_percent <- ggplot(intersect.24h, aes(x=Percentage, y=Virus, group=Group, fill=Group)) + geom_col() + 
  scale_fill_manual(values = c("Virus_vs_Virus:BPL_24h"="grey40", "Virus_vs_Mock_24h"="steelblue", "Virus_vs_Mock_24h:Virus_vs_Virus:BPL_24h"="lightsalmon")) +
  scale_x_continuous(labels = scales::percent) + scale_y_discrete(limits = rev) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 20), axis.text = element_text(size = 30), 
        axis.title.y = element_blank(), axis.title.x = element_blank(), 
        panel.background = element_rect(fill = NA), panel.grid = element_line(colour = "grey"))
ggsave("24h_Mock_vs_BPL_percentage.pdf", plot_percent, "pdf", output_folder, width = 12, height = 8)

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
