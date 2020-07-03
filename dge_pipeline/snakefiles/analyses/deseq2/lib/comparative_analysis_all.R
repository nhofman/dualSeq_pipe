set.seed(123)
library(tidyr)
library(gtools)
source("plot_heatmap.R")
source("enrichment.R")
source("STRINGdb.R")
source("venn.R")

# common subset of diff. expressed genes
# - separate up- and down-regulated genes
# - exclude HCV, MARV, LASV
# Heatmap
# Upset
# ORA
# STRING network
# Transcriptionfactor analyses
# “How many genes are regulated (up or down) during how may time points?”

virus.levels <- c("H1N1","H5N1","MERS","CoV229E","RVFV","RSV","NIV","SFSV","EBOV","MARV","HCV","LASV")
lfc.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL"),lapply(names(res.list), function(x){x.df <- data.frame(res.list[[x]]$SYMBOL,res.list[[x]]$log2FoldChange); colnames(x.df) <- c("SYMBOL",x); return(x.df)}))
rownames(lfc.df) <- lfc.df$SYMBOL
lfc.df <- lfc.df[,-1]
lfc.df <- lfc.df[,mixedorder(colnames(lfc.df))]
lfc.df <- lfc.df[,unlist(sapply(virus.levels, function(v){grep(v,colnames(lfc.df))}, simplify = F, USE.NAMES = F))]
colnames(lfc.df) <- sub("3h","_3h",colnames(lfc.df))
colnames(lfc.df) <- sub("6h","_6h",colnames(lfc.df))
lfc.df <- lfc.df[,grep(".*h$", colnames(lfc.df))]

res.list <- res.list[mixedorder(names(res.list))]

LFC.cut <- 1
res.list.filter <- sapply(names(res.list), function(n){
  x <- res.list[[n]]
  return(x[which(apply(x[,grep("normalized", colnames(x))],1,max) >= 10 & x$padj < 0.05 & abs(x$log2FoldChange) > LFC.cut),])
})

# plot number of diffeerentially expressed genes - sorted by customized order
p <- ggplot(count.genes[count.genes$LFC_cutoff==LFC.cut,], aes(y=Count, x=factor(Time, levels = mixedsort(unique(count.genes$Time))), group=Direction, fill=Direction)) + 
  geom_bar(stat = "identity", position = "stack", color = "black") + facet_wrap(~factor(Virus, levels = virus.levels), scales = "free_x") + xlab("Time") + 
  scale_y_continuous(breaks = pretty(count.genes$Count[count.genes$LFC_cutoff==LFC.cut], n=10), labels = abs(pretty(count.genes$Count[count.genes$LFC_cutoff==LFC.cut], n=10))) + 
  geom_hline(yintercept = 0) + scale_fill_manual(values=c(up="red", down="green"))
ggsave(paste0("DEG_count_LFC",LFC.cut,".svg"), p, "svg", output_folder)
 

# “How many genes are regulated (up or down) during how may time points?”
# UpSet plot or bar chart
virus.dge.binary <- list2binary(lapply(res.list.filter, "[[", "SYMBOL"))#, paste0(out.dir,"/all_genes_binary.csv"))
for(v in virus.levels){
  virus.sub <- virus.dge.binary[,grep(v,colnames(virus.dge.binary))]
  virus.sub <- virus.sub[rowSums(virus.sub)>0,]
  p <- ggplot(data.frame("sum"=rowSums(virus.sub)), aes(x=sum)) + geom_bar() + labs(x = "Number of time points", y = "# genes") + ggtitle(v)
  ggsave(paste0("Time_count_",v,".svg"), p, "svg", output_folder)
  svglite::svglite(paste0(output_folder, "UpSet_",v,".svg"), width = 10, height = 8)
  print(upset(virus.sub, order.by = "freq", nsets = 4, nintersects = NA))
  dev.off()
}

  # Common genes between viruses at any time point 
virus.dge <- sapply(virus.levels,function(virus){unique(unlist(sapply(res.list.filter[grep(paste0(virus,".*h"), names(res.list.filter))], function(x){return(x$SYMBOL)})))}, USE.NAMES = T)
out.dir <- paste0(output_folder,"/common_pattern")
genes.common <- Reduce(intersect, virus.dge[grep("CoV229E|MERS|H1N1|H5N1|RSV|RVFV|EBOV|NIV|SFSV", names(virus.dge))])
# Venn diagram or better UpSet plot of gene sets
virus.dge.binary <- list2binary(virus.dge)#, paste0(out.dir,"/all_genes_binary.csv"))
svglite::svglite(paste0(output_folder, "UpSet_minus_LASV_HCV_MARV.svg"), width = 15, height = 8)
print(upset(virus.dge.binary, nsets = 9, nintersects = 40, order.by = "freq"))
dev.off()
upset_json(file = "UpSet_minus_LASV_HCV_MARV", out.dir = output_folder, name = "DGE Virus", start = 1, end = 12)

# plot heatmap of common gene set
virus.heat <- plotHeatmap(lfc.df, filename = paste0(out.dir,"/Heatmap_common_genes_LFC",LFC.cut,".pdf"), 
                          row_subset = genes.common, 
                          colClust = F, clusterMethod = "ward.D2", legend.limit = 1, clrn = 1,
                          fontsize_row = 3.5, fontsize_col = 3.5, height = 7, border_col = NA)
write.xlsx(res.list[[1]][res.list[[1]]$SYMBOL %in% genes.common, c("SYMBOL","UNIPROTKB","PROTEIN-NAMES")], paste0(out.dir,"/common_genes.xlsx"))
# separate up- and down-regulated genes
# How to define up and down since we have multiple time points?

# over-representation analysis
ora <- calc_ora(genes.common, filename = "ORA_common_", out.dir = paste0(out.dir, "/ORA"), GO = T, REACTOME = T, ont = c("CC","BP","MF"), p.cut = 0.05, label.size = 20)

# network analysis using STRING
string_ppi(string_db, gene.df = data.frame("SYMBOL"=genes.common), filename = "common_genes", out.dir = paste0(out.dir, "/STRING"), link = F, 
           cluster = T, required_score = 400)

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
subdir <- "weak_vs_high"
high <- virus.levels[grep("EBOV|NIV", virus.levels)]
low <- virus.levels[grep("CoV229E|SFSV", virus.levels)]
lfc.sub <- lfc.df[,grep(paste(high, low, sep = "|", collapse = "|"), colnames(lfc.df))]
lfc.sub <- lfc.sub[,unlist(sapply(c(high,low), function(v){grep(v,colnames(lfc.sub))}, simplify = F, USE.NAMES = F))]
annCol <- data.frame("Pathogenicity"=ifelse(grepl(paste0(high,collapse = "|"),colnames(lfc.sub)),"high pathogenic","low pathogenic"), row.names = colnames(lfc.sub))
my_color <- list(Pathogenicity=c("high pathogenic"="dark blue", "low pathogenic"="light blue"))

venn.intersect <- draw_venn(virus.dge[grep(paste(high, low, sep = "|", collapse = "|"), names(virus.dge))], 
                            out_name = paste(out.dir, subdir, paste(high, low, sep = "_", collapse = "_"), sep = "/"), 
                            imagetype = "svg", fill = c("red", "blue", "green", "orange"), margin = 0.1)
# venn.intersect <- venn.intersect[grep(paste(paste(low,collapse = ":"),paste(high,collapse = ":"),names(venn.intersect)[grep(".+:.+:.+:.+",names(venn.intersect))],sep = "$|^"), names(venn.intersect))]
virus.heat <- plotHeatmap(lfc.sub, filename = paste0(out.dir,"/",subdir,"/Heatmap_",paste(high, low, sep = "_", collapse = "_"),"_LFC",LFC.cut,".pdf"), 
                          row_subset = Reduce(union, virus.dge[grep(paste(high, low, sep = "|", collapse = "|"), names(virus.dge))]), 
                          colClust = F, clusterMethod = "ward.D2", legend.limit = 1, clrn = 1,
                          fontsize_row = 3.5, fontsize_col = 3.5, height = 7, border_col = NA)
list.intersect <- compare_geneset(set1.list = virus.dge[grep(paste0(high,collapse = "|"),names(virus.dge))], set1.name = "high",
                                  set2.list = virus.dge[grep(paste0(low,collapse = "|"),names(virus.dge))], set2.name = "low")
#lapply(names(venn.intersect), function(x){
sapply(names(list.intersect), function(x){ 
  virus.heat <- plotHeatmap(lfc.sub, filename = paste0(out.dir,"/",subdir,"/Heatmap_",x,"_LFC",LFC.cut,".pdf"), 
                            row_subset = list.intersect[[x]], 
                            colClust = F, clusterMethod = "ward.D2", legend.limit = 1, clrn = 1,
                            fontsize_row = 3.5, fontsize_col = 6, height = 7, border_col = NA, 
                            annCol = annCol, annotation_colors = my_color)
})
list.intersect.ENTREZ <- lapply(list.intersect, function(x)bitr(x, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID)
cc <- compareCluster(list.intersect.ENTREZ, fun = "enrichKEGG", pvalueCutoff = 0.1)
cc.plot <- dotplot(cc) 

# negative sense vs positive sense 
# (+)sense: CoV229E, MERS, HCV
# (-)sense: Infuenza, EBOV, MARV, NIV, RSV, (RVFV, SFSV)
subdir <- "positive_vs_negative"
positive <- virus.levels[grep("CoV229E|MERS", virus.levels)]
negative <- virus.levels[grep("H1N1|H5N1|EBOV|NIV|RSV|RVFV|SFSV", virus.levels)]
lfc.sub <- lfc.df[,grep(paste(positive, negative, sep = "|", collapse = "|"), colnames(lfc.df))]
lfc.sub <- lfc.sub[,unlist(sapply(c(positive,negative), function(v){grep(v,colnames(lfc.sub))}, simplify = F, USE.NAMES = F))]
annCol <- data.frame("Sense"=ifelse(grepl(paste0(positive,collapse = "|"),colnames(lfc.sub)),"positive-sense","negative-sense"), row.names = colnames(lfc.sub))
my_color <- list(Sense=c("positive-sense"="dark blue", "negative-sense"="light blue"))

list.intersect <- compare_geneset(set1.list = virus.dge[grep(paste0(positive,collapse = "|"),names(virus.dge))], set1.name = "positive",
                                  set2.list = virus.dge[grep(paste0(negative,collapse = "|"),names(virus.dge))], set2.name = "negative")
sapply(names(list.intersect), function(x){ 
  virus.heat <- plotHeatmap(lfc.sub, filename = paste0(out.dir,"/",subdir,"/Heatmap_",x,"_LFC",LFC.cut,".pdf"), 
                            row_subset = list.intersect[[x]], 
                            colClust = F, clusterMethod = "ward.D2", legend.limit = 1, clrn = 1,
                            fontsize_row = 4.5, fontsize_col = 6, height = 7, border_col = NA,
                            annCol = annCol, annotation_colors = my_color)
})
list.intersect.ENTREZ <- lapply(list.intersect, function(x)bitr(x, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID)
cc <- compareCluster(list.intersect.ENTREZ, fun = "enrichKEGG")
dotplot(cc)
