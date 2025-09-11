library(UpSetR)
library(ggplot2)
library(ggpattern)

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

virus.levels <- c("H1N1","H5N1","RVFV","SFSV","RSV","NiV","EBOV","MARV","LASV")
res.list <- res.list[mixedorder(names(res.list))]

# Get significantly differentially expressed genes
LFC.cut <- 1
res.list.filter <- sapply(names(res.list)[-grep(".*BPL", names(res.list))], function(n){
  x <- res.list[[n]]
  return(x[which(apply(x[,grep("normalized", colnames(x))],1,max) >= 10 & x$padj < 0.05 & abs(x$log2FoldChange) > LFC.cut),])
})
names(res.list.filter) <- sub("_"," ",sub("_vs.*", "", names(res.list.filter)))


# “How many genes are regulated (up or down) during how may time points?”
# UpSet plot 
virus.dge.binary <- list2binary(lapply(res.list.filter, "[[", "SYMBOL"))#, paste0(out.dir,"/all_genes_binary.csv"))
par(mfrow=c(3,3))
for(virus in virus.levels){
  print(virus)
  virus.sub <- virus.dge.binary[,grep(virus,colnames(virus.dge.binary))]
  virus.sub <- virus.sub[rowSums(virus.sub)>0,]
  #svglite::svglite(paste0(output_folder, "UpSet_",virus,"_Mock.svg"), width = 14, height = 8)
  pdf(paste0(output_folder, "UpSet_",virus,"_Mock.pdf"), width = 15, height = 8, family = "Helvetica")
  #svg(paste0(output_folder, "UpSet_", virus, "_Mock.svg"), width = 14, height = 8, family = "Arial")
  print(upset(virus.sub, order.by = "freq", nsets = 4, nintersects = NA, text.scale = c(3,3,2,2,3,3), line.size = 1, point.size = 4,
              set_size.scale_max = round(max(colSums(virus.sub))+1300, -3)+100, set_size.show = T,
              sets = colnames(virus.sub[,colSums(virus.sub)>0])[mixedorder(colnames(virus.sub[,colSums(virus.sub)>0]), decreasing = T)], keep.order = TRUE), newpage=F)
  dev.off()
  #system(paste0("inkscape -l ", output_folder, "UpSet_", virus, "_Mock.svg ", output_folder, "UpSet_", virus, "_Mock.pdf"))
  #system(paste0("inkscape --export-type=svg ", output_folder, "UpSet_", virus, "_Mock.pdf"))
}

# BPL
res.list.filter.BPL <- sapply(names(res.list)[grep("BPL_vs_Mock", names(res.list))], function(n){
  x <- res.list[[n]]
  return(x[which(apply(x[,grep("normalized", colnames(x))],1,max) >= 10 & x$padj < 0.05 & abs(x$log2FoldChange) > 0),])
}, simplify = F)

# Compare Mock and BPL

p <- ggplot(count.genes[count.genes$LFC_cutoff==LFC.cut & grepl("24h|BPL", count.genes$Sample),], aes(y=Count, x=factor(Sample, levels = unique(Sample)), group=Direction, fill=factor(Direction, labels = c("Down","Up")))) + 
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.25, width = 0.8) + facet_wrap(~Virus, scales = "free_x") + xlab("Time") + 
  scale_y_continuous(breaks = pretty(count.genes$Count[count.genes$LFC_cutoff==LFC.cut], n=5), labels = abs(pretty(count.genes$Count[count.genes$LFC_cutoff==LFC.cut], n=5))) + 
  geom_hline(yintercept = 0, linewidth = 0.25) + scale_fill_manual(values=c(Up="red", Down="blue"), guide = guide_legend(reverse=T)) +
  xlab("Time after infection") + ylab("Number of genes") +
  theme(text = element_text(face = "plain"), line = element_line(linewidth = 0.25),
        axis.line.x.top = element_blank(), axis.line.x.bottom = element_line(color = "black", linewidth = 0.25),
        axis.line.y.right = element_blank(), axis.line.y.left = element_line(color = "black", linewidth = 0.25),
        panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "gray58"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(t=8,r=0,b=0,l=0)), 
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(t=0,r=8,b=0,l=0)),
        axis.text = element_text(size = 10), axis.text.x = element_text(angle = 20, vjust = 0.5), 
        strip.text = element_text(size = 12, face = "bold"), strip.background = element_rect(fill = NA, color = NA), #panel.spacing = unit(2, "lines"),
        legend.text = element_text(size = 10, face = "plain"), legend.title = element_blank())
ggsave(paste0("DEG_count_MockvsBPL_LFC",LFC.cut,".pdf"), p, "pdf", output_folder, width = 16, height = 7)

LFC.cut <- 1
res.list.filter.24h <- sapply(names(res.list)[grep("24h", names(res.list))], function(n){
  x <- res.list[[n]]
  return(x[which(apply(x[,grep("normalized", colnames(x))],1,max) >= 10 & x$padj < 0.05 & abs(x$log2FoldChange) > LFC.cut),])
})
intersect.list <- sapply(virus.levels, function(x){
  tmp <- res.list.filter.24h[grep(x, names(res.list.filter.24h))]
  if(sum(sapply(tmp, nrow))>0){
    intersect.tmp <- attr(venn(sapply(tmp, "[[", "SYMBOL"), show.plot = T), "intersections")
    venn(list("Mock"=intersect.tmp[[paste0(x, "_24h_vs_Mock_24h")]], "BPL"=res.list.filter.BPL[[grep(x, names(res.list.filter.BPL))]]$SYMBOL))
  #intersect.list <- attr(venn(list("vs_Mock_24h"=res.list.filter.BPL[[paste0(x,"_BPL_vs_Mock_BPL")]]$SYMBOL, "vs_Mock_BPL"=res.list.filter.BPL[[paste0(x,"_BPL_vs_Mock_24h")]]$SYMBOL)), "intersections")
    return(intersect.tmp)
  }
})

#names(res.list.filter.24h) <- sub("_BPL",":BPL_24h",sub("24h_vs_", "vs_", names(res.list.filter.24h)))
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
  names(data) <- gsub("_24h","",names(data))
  data.df <- data.frame("Virus"=v,"Group"=names(data), "Count"=sapply(data, length), row.names = NULL)
  data.df$Sum <- sum(data.df$Count)
  data.df$Percentage <- data.df$Count/data.df$Sum
  return(data.df)
}, simplify = F))
intersect.24h[is.nan(intersect.24h$Percentage), "Percentage"] <- 0
intersect.24h$Group2 <- gsub("_", " ", sub(".*:.*","intersect",intersect.24h$Group))
intersect.24h$Group2 <- sub("Virus BPL","BPL-Virus",intersect.24h$Group2)
intersect.24h$Group2 <- factor(intersect.24h$Group2, levels = c("Virus vs BPL-Virus", "intersect", "Virus vs Mock"))
#intersect.24h$Group <- factor(intersect.24h$Group, levels = c("BPL_24h", "Mock_24h:BPL_24h", "Mock_24h"))
intersect.24h$Virus <- factor(intersect.24h$Virus, levels = virus.levels)
plot_24h <- ggplot(intersect.24h, aes(x=Percentage, y=Virus, group=Group2, fill=Group2)) 
plot_24h <- ggplot(intersect.24h[intersect.24h$Virus!="LASV",], aes(x=Count, y=Virus, group=Group2, fill=Group2)) 
plot_24h <- plot_24h + 
  geom_col_pattern(aes(pattern=Group2), pattern_fill = "grey40", pattern_color = "grey40", pattern_density = 0.5, pattern_spacing = 0.025) + 
  xlab("Ratio of DEG") + scale_pattern_manual(values = c("none", "stripe", "none")) + 
  scale_fill_manual(values = c("Virus vs BPL-Virus"="grey40", "Virus vs Mock"="grey80", "intersect"="grey80"),
                    breaks = c("Virus vs Mock", "Virus vs BPL-Virus"), #labels = c("Virus vs Mock", "Virus vs BPL-Virus"),
                    guide = guide_legend(override.aes = list(pattern = "none"))) +
  guides(pattern = "none") + scale_y_discrete(limits = rev) + scale_x_continuous(labels = scales::percent) +
  theme(text = element_text(family = "Helvetica", face = "bold"),
        legend.title = element_blank(), legend.text = element_text(size = 28, face = "plain"), axis.text = element_text(size = 28), 
        axis.title.x = element_text(size = 35, margin = margin(7, 0, 0, 0, "mm")), 
        axis.title.y = element_text(size = 35, margin = margin(0, 7, 0, 0, "mm")), 
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.background = element_rect(fill = NA)) #, panel.grid = element_blank())
ggsave("24h_Mock_vs_BPL_percent_new.pdf", plot_24h, "pdf", output_folder, width = 23, height = 12)
system(paste0("inkscape -l ", output_folder, "24h_Mock_vs_BPL_percentage_new.svg ", output_folder, "24h_Mock_vs_BPL_percentage_new.pdf"))
system(paste0("inkscape --export-type svg ", output_folder, "24h_Mock_vs_BPL_percentage_new.pdf"))