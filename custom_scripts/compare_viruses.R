set.seed(123)
library(tidyr)
library(gtools)
library(gplots)
library(UpSetR)
library(enrichplot)

#source("plot_heatmap.R")
#source("enrichment.R")

# Convert lists of DEG to binary table
list2binary <- function(data.list, filename){
  data.binary <- Reduce(function(x,y)merge(x,y,by="SYMBOL",all=T),sapply(names(data.list),function(n){
    if(length(data.list[[n]])>0){
      data.df <- data.frame(data.list[[n]],1)
    }else{
      data.df <- data.frame(matrix(ncol=2,nrow=0)) 
    }
    colnames(data.df) <- c("SYMBOL", n)
    return(data.df)
  }, USE.NAMES = T, simplify = F))
  data.binary[is.na(data.binary)] <- 0
  if(!missing(filename)){
    write.csv(data.binary, filename, row.names = F)
  }
  return(data.binary)
}

#load("deseq2.RData") # load data from DESeq2 analysis
virus.levels <- c("H1N1","H5N1","RVFV","SFSV","RSV","NiV","EBOV","MARV","LASV")
res.list <- res.list[mixedorder(names(res.list))]

lfc.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL"),lapply(names(res.list)[grep(".*h_vs_Mock.*", names(res.list))], function(x){x.df <- data.frame(res.list[[x]]$SYMBOL,res.list[[x]]$log2FoldChange); colnames(x.df) <- c("SYMBOL",x); return(x.df)}))
rownames(lfc.df) <- lfc.df$SYMBOL
lfc.df <- lfc.df[,-1]
lfc.df <- lfc.df[,mixedorder(colnames(lfc.df))]
lfc.df <- lfc.df[,unlist(sapply(virus.levels, function(v){grep(v,colnames(lfc.df))}, simplify = F, USE.NAMES = F))]

# Get significantly differentially expressed genes
LFC.cut <- 1
res.list.filter <- sapply(names(res.list)[-grep(".*BPL", names(res.list))], function(n){
  x <- res.list[[n]]
  return(x[which(apply(x[,grep("normalized", colnames(x))],1,max) >= 10 & x$padj < 0.05 & abs(x$log2FoldChange) > LFC.cut),])
})
names(res.list.filter) <- sub("_"," ",sub("_vs.*", "", names(res.list.filter)))

# Common genes between viruses at any time point 
virus.levels <- c("H1N1","H5N1","RVFV","SFSV","RSV","NiV","EBOV")
virus.dge <- sapply(virus.levels,function(virus){unique(unlist(sapply(res.list.filter[grep(virus, names(res.list.filter))], function(x){return(as.character(x$SYMBOL))})))}, USE.NAMES = T)
out.dir <- paste0(output_folder,"/common_pattern/")
if(!dir.exists(out.dir)){
  dir.create(out.dir)
}
genes.common <- Reduce(intersect, virus.dge)

# UpSet plot of gene sets
virus.dge.binary <- list2binary(virus.dge)
pdf(paste0(out.dir, "/UpSet_minus_LASV_MARV.pdf"), width = 20, height = 8, family = "Helvetica")
#svg(paste0(out.dir, "/UpSet_minus_LASV_MARV_svg.svg"), width = 20, height = 8, family = "")
print(upset(fromList(virus.dge), nsets = 7, nintersects = 40, order.by = "freq", text.scale = c(2,2,2,2,2,2),
            point.size = 3, line.size = 1, number.angles = 0, set_size.show = F, set_size.numbers_size = 8,
            set_size.scale_max = round(max(sapply(virus.dge, length))+500, -3)+100,
            queries = list(list(query = intersects, params = list(virus.levels), active = T, color = "red"))), newpage = F)
dev.off()
#system(paste0("inkscape --export-type=svg ", out.dir, "UpSet_minus_LASV_MARV.pdf"))

lfc.df.common <- lfc.df[genes.common, -grep("MARV|LASV", colnames(lfc.df))]

# Create heatmap of common genes across viruses and time points
genes.clust <- hclust(dist(lfc.df.common, method = "euclidean"), method = "complete")
distMethod <- "euclidean"
clusterMethod <- "ward.D2"
virus.heat <- plotHeatmap(lfc.df.common, filename = paste0(out.dir,"/Heatmap_common_genes_",distMethod,"_",clusterMethod,"_split2.pdf"), 
                          row.dend = T, distMethod = distMethod, plot.fig = T,
                          colClust = F, clusterMethod = clusterMethod, legend.cut = 1, clrn = 1, cellheight = 5,
                          fontsize_row = 4, fontsize_col = 7, height = 13, border_col = NA, family = "Helvetica")
dend.1 <- genes.common[row_order(virus.heat)[[1]]]
dend.2 <- genes.common[row_order(virus.heat)[[2]]]
length(dend.1) <- length(dend.2)
write.csv(cbind(dend.1, dend.2), paste0(out.dir, "Heatmap_split.csv"), quote = F, na = "", row.names = F)

# Split heatmap in two parts
genes.clust <- row_order(virus.heat)
virus.heat1 <- plotHeatmap(lfc.df.common, filename = paste0(out.dir,"/Heatmap_common_genes_LFC",LFC.cut,"_part1.pdf"), 
                           row_subset = genes.clust[1:(length(genes.clust)/2)], row.dend = F, rowClust = F, clrn = 1,
                           colClust = F, clusterMethod = "ward.D2", legend.cut = 1, legend.limit.up = max(lfc.df.common), legend.limit.down = min(lfc.df.common), 
                           fontsize_row = 5.5, fontsize_col = 7.5, height = 10, border_col = NA, family = "Helvetica")
virus.heat2 <- plotHeatmap(lfc.df.common, filename = paste0(out.dir,"/Heatmap_common_genes_LFC",LFC.cut,"_part2.pdf"), 
                           row_subset = genes.clust[((length(genes.clust)/2)+1):length(genes.clust)], row.dend = F, rowClust = F, clrn = 1,
                           colClust = F, clusterMethod = "ward.D2", legend.cut = 1, legend.limit.up = max(lfc.df.common), legend.limit.down = min(lfc.df.common), 
                           fontsize_row = 5.5, fontsize_col = 7.5, height = 10, border_col = NA, family = "Helvetica")

# Over-representation analysis of common genes
ora <- calc_ora(genes.common, filename ="ORA_common", GO = T, REACTOME = T, KEGG = T, ont = c("CC","BP","MF"), #legendlimit=c(0.003, 0.0425),
                p.cut = 0.05, label.size = 8, legend.size = 6, title.size = 8, imagetype = "pdf", ink = "--export-type=svg",
                width = 10, height = 7, family = "Helvetica", out.dir = paste0(output_folder, "common_pattern"), 
                label_format = function(x) stringr::str_wrap(x, width = 50), category_top = 20)

# Plot heatplot for ORA result
offspring.tbl_tree <- utils::getFromNamespace("offspring.tbl_tree", "tidytree")
offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")

lfc.common <- lfc.df[genes.common,]
for(db in names(ora)){
  print(db)
  if(db != "CC"){
    #options(enrichplot.colours = c("#327eba","#e06663"))
    #options(enrichplot.colours = c("blue","white","red"))
    categories <- nrow(data.frame(ora[[db]]))
    height <- ifelse(categories > 20, 20, categories)
    fc <- rowMedians(as.matrix(lfc.common[, grep("24", colnames(lfc.common))]))
    fc <- rowMedians(as.matrix(lfc.common))
    names(fc) <- rownames(lfc.common)
    width_heat <- length(unique(unlist(strsplit(data.frame(ora[[db]])[1:20, "geneID"], "/", fixed = T)))) 
    heat_plot <- heatplot(ora[[db]], showCategory=50, foldChange = fc, label_format = function(x) stringr::str_wrap(x, width = 30)) + 
      scale_fill_gradient2(limits = NULL, low = "blue", high = "red", name = "LFC") +
      theme(text = element_text(family = "Helvetica", face = "bold", size = 10), legend.text = element_text(face = "plain"), 
            legend.title = element_text(size = 10), axis.title = element_text(size = 1), axis.text.y = element_text(hjust = 1), 
            axis.title.x = element_text(margin = margin(8, 0, 0, 0, "mm")), axis.line.x.top = element_blank(), 
            axis.line.y.right = element_blank(), axis.line.x.bottom = element_line(color = "black"), 
            axis.line.y.left = element_line(color = "black"), panel.border = element_blank()) 
    ggsave(paste("ORA_common_",db,"_heatplot_median.", "pdf", sep = ""), device = "pdf", path = paste0(output_folder, "common_pattern"), 
           plot = egg::set_panel_size(heat_plot, width = unit(width_heat*0.4, "cm"), height = unit(height-2, "cm")), 
           width = width_heat*0.6, height = height, units = "cm") 
  }
}

