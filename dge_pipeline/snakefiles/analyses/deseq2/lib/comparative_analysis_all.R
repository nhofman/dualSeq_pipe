set.seed(123)
library(tidyr)
library(gtools)
library(openxlsx)
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

virus.levels <- c("H1N1","H5N1","MERS","CoV229E","RVFV","SFSV","RSV","NIV","EBOV","MARV","HCV","LASV")
#lfc.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL"),lapply(names(res.list)[grep(".*h.*Mock", names(res.list))], function(x){x.df <- data.frame(res.list[[x]]$SYMBOL,res.list[[x]]$log2FoldChange); colnames(x.df) <- c("SYMBOL",x); return(x.df)}))
lfc.df.all <- Reduce(function(x,y)merge(x,y,by="SYMBOL"),lapply(names(res.list)[grep(".*h", names(res.list))], function(x){x.df <- data.frame(res.list[[x]]$SYMBOL,res.list[[x]]$log2FoldChange); colnames(x.df) <- c("SYMBOL",x); return(x.df)}))
rownames(lfc.df.all) <- lfc.df.all$SYMBOL
lfc.df.all <- lfc.df.all[,-1]
lfc.df.all <- lfc.df.all[,mixedorder(colnames(lfc.df.all))]
lfc.df.all <- lfc.df.all[,unlist(sapply(virus.levels, function(v){grep(v,colnames(lfc.df.all))}, simplify = F, USE.NAMES = F))]
colnames(lfc.df.all) <- sub("3h","_3h",colnames(lfc.df.all))
colnames(lfc.df.all) <- sub("6h","_6h",colnames(lfc.df.all))
#lfc.df <- lfc.df[,grep(".*h.*Mock", colnames(lfc.df))]
lfc.df <- lfc.df.all[, -grep("vs_.*BPL$", colnames(lfc.df.all))]

res.list <- res.list[mixedorder(names(res.list))]

LFC.cut <- 1
res.list.filter <- sapply(names(res.list)[grep(".*h", names(res.list))], function(n){
  x <- res.list[[n]]
  return(x[which(apply(x[,grep("normalized", colnames(x))],1,max) >= 10 & x$padj < 0.05 & abs(x$log2FoldChange) > LFC.cut),])
})

# Plot PCA
shape <- if(length(unique(conditiontable$time)) <= 6){ scales::shape_pal()(length(unique(conditiontable$time))) }else{ c(1:length(unique(conditiontable$time)))}
names(shape) <- unique(conditiontable$time)
pca <- plotPCA(deseq.results.vst, intgroup = c("treatment", "time"), returnData = TRUE)
plot_PCA <- ggplot(pca, aes(PC1, PC2, color = factor(treatment, levels = c(virus.levels, "Mock")), shape = factor(time, levels = mixedsort(as.character(unique(conditiontable$time)))))) + 
  geom_point(size=3) + labs(color = "Treatment", shape = "Time") + scale_shape_manual(values=shape) + guides(color=guide_legend(override.aes=list(fill=NA))) +
  theme(axis.title = element_text(size = 20, face = "bold"), axis.text = element_text(size = 16),  
        legend.text = element_text(size = 16), legend.title = element_text(size = 20, face = "bold"), legend.key=element_blank()) + 
  guides(color = guide_legend(order = 2), shape = guide_legend(order = 1))
if(exists("color")){
  plot_PCA <- plot_PCA + scale_colour_manual(values=color)
}
#plot_PCA <- ggplot(pca, aes(PC1, PC2, color = treatment, shape = factor(time, levels = mixedsort(levels(pca$time))))) + geom_point(size=3) + 
#  labs(color = "Treatment", shape = "Time")
ggsave("PCA.svg", plot = plot_PCA, device = "svg", path = output_folder, width = 10, height = 8)

for (time in unique(conditiontable$time)) {
  pca_time <- plotPCA(deseq.results.vst[,grep(time, colnames(deseq.results.vst))], intgroup = c("treatment"), returnData = TRUE)
  plot_PCA <- ggplot(pca_time, aes(PC1, PC2, color = factor(treatment, levels = c(virus.levels, "Mock")), shape = time)) + geom_point(size=3) + labs(color = "Treatment", shape = "Time") +
    theme(axis.title = element_text(size = 20, face = "bold"), axis.text = element_text(size = 16),  
          legend.text = element_text(size = 16), legend.title = element_text(size = 20, face = "bold"), legend.key=element_blank()) +
  scale_shape_manual(values=shape) + guides(color = guide_legend(order = 2), shape = guide_legend(order = 1))
  if(exists("color")){
    plot_PCA <- plot_PCA + scale_colour_manual(values=color)
  }
  ggsave(paste0("PCA_", time, ".svg"), plot = plot_PCA, device = "svg", path = output_folder)
}

# plot number of diffeerentially expressed genes - sorted by customized order
count.genes <- separate(count.genes, "Sample", c("Virus","Time","Vs"), "_", F, extra = "merge")
count.genes$Vs <- sub("vs_","",count.genes$Vs)
count.genes$x <- ifelse(grepl("Mock", count.genes$Vs), count.genes$Time, "active")
count.genes.mod <- count.genes[count.genes$Vs!="Mock_BPL",]
count.genes.mod$Time <- ifelse(count.genes.mod$Control=="inactive", "BPL", count.genes.mod$Time)
p <- ggplot(count.genes[count.genes$LFC_cutoff==LFC.cut,], aes(y=Count, x=factor(x, levels = c("3h","6h","12h","24h","active","48h","BPL")), group=Direction, fill=factor(Direction, labels = c("Down","Up")))) + 
  geom_bar(stat = "identity", position = "stack", color = "black") + facet_wrap(~factor(Virus, levels = virus.levels), scales = "free_x") + 
  scale_y_continuous(breaks = pretty(count.genes$Count[count.genes$LFC_cutoff==LFC.cut], n=5), labels = abs(pretty(count.genes$Count[count.genes$LFC_cutoff==LFC.cut], n=5))) + 
  geom_hline(yintercept = 0) + scale_fill_manual(values=c(Up="red", Down="green"), guide = guide_legend(reverse=T)) + 
  xlab("Time after infection") + ylab("Number of genes") +
  theme(axis.title.x = element_text(size = 20, face = "bold", margin = margin(t=10,r=0,b=0,l=0)), axis.title.y = element_text(size = 20, face = "bold", margin = margin(t=0,r=10,b=0,l=0)),
        axis.text = element_text(size = 20), axis.text.x = element_text(angle = 0), 
        strip.text = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 20), legend.title = element_text(size = 24, face = "bold")) + labs(fill="Diff. expression")
ggsave(paste0("DEG_count_LFC",LFC.cut,".svg"), p, "svg", output_folder, width = 20, height = 10)
 

# “How many genes are regulated (up or down) during how may time points?”
# UpSet plot or bar chart
virus.dge.binary <- list2binary(lapply(res.list.filter, "[[", "SYMBOL"))#, paste0(out.dir,"/all_genes_binary.csv"))
for(v in virus.levels){
  print(v)
  virus.sub <- virus.dge.binary[,grep(paste0(v, ".*Mock"),colnames(virus.dge.binary))]
  virus.sub <- virus.sub[rowSums(virus.sub)>0,]
  svglite::svglite(paste0(output_folder, "UpSet_",v,"_Mock.svg"), width = 10, height = 8)
  print(upset(virus.sub, order.by = "freq", nsets = 4, nintersects = NA, text.scale = c(2,2,2,1.75,2,2), 
              sets = colnames(virus.sub[,colSums(virus.sub)>0])[mixedorder(colnames(virus.sub[,colSums(virus.sub)>0]), decreasing = T)], keep.order = TRUE))
  dev.off()
  virus.sub <- virus.dge.binary[,grep(v,colnames(virus.dge.binary))]
  virus.sub <- virus.sub[rowSums(virus.sub)>0,]
  svglite::svglite(paste0(output_folder, "UpSet_",v,".svg"), width = 10, height = 8)
  print(upset(virus.sub, order.by = "freq", nsets = 4, nintersects = NA, text.scale = c(2,2,2,1.75,2,2), 
              sets = colnames(virus.sub[,colSums(virus.sub)>0])[mixedorder(colnames(virus.sub[,colSums(virus.sub)>0]), decreasing = T)], keep.order = TRUE))
  dev.off()
  #p <- ggplot(data.frame("sum"=rowSums(virus.sub)), aes(x=sum)) + geom_bar() + labs(x = "Number of time points", y = "# genes") + ggtitle(v)
  #ggsave(paste0("Time_count_",v,".svg"), p, "svg", output_folder)
}

# Comparison Mock vs inactivatd viruses
if(!dir.exists(paste0(output_folder, "/Comparison_controls"))){
  dir.create(paste0(output_folder, "/Comparison_controls"))
}

for(v in virus.levels){
  print(v)
  res.list.v <- res.list.filter[grep(paste0(v, ".*24h"), names(res.list.filter))]
  names(res.list.v) <- ifelse(grepl("BPL",names(res.list.v)), paste0(sub("_vs_.*","",names(res.list.v)),"_active"), sub("_vs_.*","",names(res.list.v)))
  res.list.v.split <- lapply(res.list.v, function(r){
    #r <- res.list.v[[x]]
    up <- as.character(r[r$log2FoldChange>0, "SYMBOL"])
    down <- as.character(r[r$log2FoldChange<0, "SYMBOL"])
    return(list("up"=up, "down"=down))
  })
  res.list.v.split <- unlist(res.list.v.split, recursive = F)
  venn.intersect <- draw_venn(lapply(res.list.v, "[[", "SYMBOL"), 
                              out_name = paste(output_folder, "Comparison_controls","Venn", paste("Venn", v, sep = "_"), sep = "/"), 
                              imagetype = "svg", fill = color[v], margin = 0.1, cex = 2, cat.cex = 2, ext.text = T, cat.dist = 0.04, cat.pos = c(-150,150))
  calc_compareCluster(venn.intersect, v, paste0(output_folder, "Comparison_controls"), GO = T,  ont = c("BP","MF"), label.size = 20, w = 30, REACTOME = F, keytype = "SYMBOL", legend.size = 15)
  venn.intersect <- draw_venn(res.list.v.split[grep("up", names(res.list.v.split))], 
                              out_name = paste(output_folder, "Comparison_controls", paste("Venn", v, "up", sep = "_"), sep = "/"), 
                              imagetype = "svg", fill = color[v], margin = 0.3, cex = 1, cat.cex = 1)
  venn.intersect <- draw_venn(res.list.v.split[grep("down", names(res.list.v.split))], 
                              out_name = paste(output_folder, "Comparison_controls", paste("Venn", v, "down", sep = "_"), sep = "/"), 
                              imagetype = "svg", fill = color[v], margin = 0.3, cex = 1, cat.cex = 1)
}

# Common genes between viruses at any time point 
virus.dge <- sapply(virus.levels,function(virus){unique(unlist(sapply(res.list.filter[grep(paste0(virus,".*h.*Mock"), names(res.list.filter))], function(x){return(x$SYMBOL)})))}, USE.NAMES = T)
out.dir <- paste0(output_folder,"/common_pattern")
genes.common <- Reduce(intersect, virus.dge[grep("CoV229E|MERS|H1N1|H5N1|RSV|RVFV|EBOV|NIV|SFSV", names(virus.dge))])
# Venn diagram or better UpSet plot of gene sets
virus.dge.binary <- list2binary(virus.dge)#, paste0(out.dir,"/all_genes_binary.csv"))
svglite::svglite(paste0(out.dir, "/UpSet_minus_LASV_HCV_MARV.svg"), width = 20, height = 8)
print(upset(virus.dge.binary, nsets = 9, nintersects = 40, order.by = "freq", text.scale = c(2,2,1.5,1.5,2,2),
            queries = list(list(query = intersects, params = list("CoV229E","MERS","H1N1","H5N1","RVFV","SFSV","RSV","NIV","EBOV"), active = T, color = "red"))))
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
heat.gg <- ggplotify::as.ggplot(virus.heat$plot$gtable)
heat.gg <- heat.gg$theme$legend.text + theme(legend.text = element_text(size = 20))
ggsave("Heatmap_test.svg",heat.gg,"svg","Documents/Virus_project/analyses/host/deseq2_new/common_pattern/")

virus.heat <- plotHeatmap(lfc.df[, -grep("HCV|MARV|LASV", colnames(lfc.df))], filename = paste0(out.dir,"/Heatmap_common_genes_LFC",LFC.cut,"_split.pdf"), 
                          row_subset = genes.common, plot.fig = F,
                          colClust = F, clusterMethod = "ward.D2", legend.limit = 1, clrn = 2,
                          fontsize_row = 3.5, fontsize_col = 3.5, height = 7, border_col = NA)
for(i in unique(virus.heat$row_cluster)){
  genes.cl <- names(virus.heat$row_cluster[virus.heat$row_cluster==i])
  if(sum(lfc.df[genes.cl,] > 0) > sum(lfc.df[genes.cl,] < 0)){
    reg <- "up"
  }else{
    reg <- "down"
  }
  virus.heat.cl <- plotHeatmap(lfc.df[, -grep("HCV|MARV|LASV", colnames(lfc.df))],
                            filename = paste0(out.dir,"/Heatmap_common_genes_LFC",LFC.cut,"_Cluster_",i,"_",reg,".pdf"), 
                            row_subset = genes.cl, 
                            colClust = F, clusterMethod = "ward.D2", legend.limit = 1, clrn = 1,
                            fontsize_row = 3.5, fontsize_col = 3.5, height = 7, border_col = NA)
  ora <- calc_ora(genes.cl, filename = paste0("ORA_common_Cluster_",i,"_",reg), out.dir = paste0(out.dir, "/ORA"), GO = T, REACTOME = T, ont = c("CC","BP","MF"), 
                  p.cut = 0.05, label.size = 18, legend.size = 18, legend.title.size = 20)
}
write.xlsx(res.list[[1]][res.list[[1]]$SYMBOL %in% genes.common, c("SYMBOL","UNIPROTKB","PROTEIN-NAMES")], paste0(out.dir,"/common_genes.xlsx"))
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
ora <- calc_ora(genes.common, filename = "ORA_common", out.dir = paste0(out.dir, "/ORA"), GO = T, REACTOME = T, ont = c("CC","BP","MF"), 
                p.cut = 0.05, label.size = 18, legend.size = 18, legend.title.size = 20)

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
motifDB <- "/home/nina/Documents/Virus_project/motif_databases/JASPAR/JASPAR2020_CORE_vertebrates_redundant.meme"
motifDB <- "/home/nina/Documents/Virus_project/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme"
fasta.file <- "/vol/sfb1021/SFB1021_Virus/genomes/hg38/genome.fa"
out.dir <- paste0(out.dir, "/TFBS/")
promotor_prim <- "promotor.common.fa"
promotor_back <- "promotor.bg.fa"
writeGRanges2Fasta(human.GRanges, out.dir, promotor_prim, promotor_back, fasta.file, genes.common, upstream = 1000, downstream = 100)
meme.res <- meme(path2meme, outdir = out.dir, paste0(out.dir,promotor_prim), paste0(out.dir,promotor_back), nmotif = 10, alph = "dna", objfun = "de")
tomtom(path2meme, motifDB, motif_file = meme.res, outdir = out.dir)
ame(path2meme, motifDB, paste0(output_folder,"/common_pattern/TFBS/",promotor_prim), paste0(output_folder,"/common_pattern/TFBS/",promotor_back), out.dir)

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
high <- virus.levels[grep("EBOV|NIV", virus.levels)]
low <- virus.levels[grep("CoV229E|SFSV", virus.levels)]
lfc.sub <- lfc.df[,grep(paste(high, low, sep = "|", collapse = "|"), colnames(lfc.df))]
lfc.sub <- lfc.sub[,unlist(sapply(c(high,low), function(v){grep(v,colnames(lfc.sub))}, simplify = F, USE.NAMES = F))]

venn.intersect <- draw_venn(virus.dge[grep(paste(high, low, sep = "|", collapse = "|"), names(virus.dge))], 
                            out_name = paste(out.dir, paste(high, low, sep = "_", collapse = "_"), sep = "/"), 
                            imagetype = "svg", fill = color[names(virus.dge[grep(paste(high, low, sep = "|", collapse = "|"), names(virus.dge))])], margin = 0.1)
# venn.intersect <- venn.intersect[grep(paste(paste(low,collapse = ":"),paste(high,collapse = ":"),names(venn.intersect)[grep(".+:.+:.+:.+",names(venn.intersect))],sep = "$|^"), names(venn.intersect))]
list.intersect <- compare_geneset(set1.list = virus.dge[grep(paste0(high,collapse = "|"),names(virus.dge))], set1.name = "high",
                                  set2.list = virus.dge[grep(paste0(low,collapse = "|"),names(virus.dge))], set2.name = "low")
#lapply(names(venn.intersect), function(x){
sapply(names(list.intersect), function(x){ 
  virus.heat <- plotHeatmap(lfc.sub, filename = paste0(out.dir,"/Heatmap_",x,"_LFC",LFC.cut,".pdf"), 
                            row_subset = list.intersect[[x]], 
                            colClust = F, clusterMethod = "ward.D2", legend.limit = 1, clrn = 1,
                            fontsize_row = 3.5, fontsize_col = 6, height = 7, border_col = NA)
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
positive <- virus.levels[grep("CoV229E|MERS", virus.levels)]
negative <- virus.levels[grep("H1N1|H5N1|EBOV|NIV|RSV|RVFV|SFSV", virus.levels)]
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
                            row_subset = list.intersect[[x]], 
                            colClust = F, clusterMethod = "ward.D2", legend.limit = 1, clrn = 1,
                            fontsize_row = 4.5, fontsize_col = 6, height = 7, border_col = NA,
                            annCol = annCol, annotation_colors = my_color)
})
list.intersect.ENTREZ <- lapply(list.intersect, function(x)bitr(x, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID)
cc <- compareCluster(list.intersect.ENTREZ, fun = "enrichKEGG")
dotplot(cc)
