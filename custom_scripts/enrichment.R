library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(ggplot2)
library(openxlsx)


# Over-representation analysis for a set of genes
calc_ora <- function(geneset, main = "", filename, out.dir = "ORA", ink = "-l", GO = T, KEGG = T, REACTOME = F, ont = "BP", p.cut = 0.05, legendlimit = NULL, label.size = 12, dpi = 300,
                     legend.size = 8, title.size = 8, keytype = "SYMBOL", width = 18, height = 15, category_top = 20, imagetype = "svg", family = family, label_format = 30){
  if(!dir.exists(out.dir)){
    dir.create(out.dir)
  }
  if(keytype != "ENTREZID"){
    gene_bitr <- bitr(geneset, fromType=keytype, toType="ENTREZID", "org.Hs.eg.db")
  }else{
    gene_bitr <- data.frame(ENTREZID=geneset)
  }
  theme_plot <- theme(text = element_text(family = family, face = "bold", size = label.size),
                  legend.text = element_text(size = legend.size, face = "plain"), legend.title = element_text(size = title.size),
                  axis.title = element_text(size = title.size), axis.text.y = element_text(hjust = 1),
                  axis.title.x = element_text(margin = margin(7, 0, 0, 0, "mm")),
                  axis.line.x.top = element_blank(), axis.line.y.right = element_blank(), axis.line.x.bottom = element_line(color = "black"), 
                  axis.line.y.left = element_line(color = "black"), panel.border = element_blank())
  ora.list <- list()
  if(GO){
    for(o in ont){
      ora_go <- try(enrichGO(gene_bitr$ENTREZID, keyType = "ENTREZID", OrgDb = "org.Hs.eg.db", ont = o, readable = T, pvalueCutoff = p.cut))
      if(class(ora_go)!="try-error"){
        categories <- nrow(data.frame(ora_go))
        if(categories > 0){
          height <- ifelse(categories > 20, 20, categories)
          plot_go <- dotplot(ora_go, showCategory = category_top, font.size = label.size, label_format = label_format) + ggtitle(main) +
            scale_fill_gradient(limits = legendlimit, low = "red", high = "blue") + theme_plot
          ggsave(paste(filename,"_",o,"_dotplot.", imagetype, sep = ""), device = imagetype, path = paste(out.dir, sep = ""),
                 plot = egg::set_panel_size(plot_go, width = unit(width, "cm"), height = unit(height, "cm")), width = width+10, height = height+2, units = "cm", dpi = dpi)
          write.table(data.frame(ora_go), file = paste(out.dir, "/", filename, "_", o, ".csv", sep = ""), sep = ",", row.names = FALSE)
        }
        ora.list[[o]] <- ora_go
      }
    }
  }
  if(KEGG){
    ora_kegg <- try(enrichKEGG(gene_bitr$ENTREZID, organism = "hsa", use_internal_data = FALSE, pvalueCutoff = p.cut))
    if(class(ora_kegg)!="try-error"){
      categories <- nrow(data.frame(ora_kegg))
      if(categories > 0){
        height <- ifelse(categories > 20, 20, categories)
        ora_kegg <- setReadable(ora_kegg, org.Hs.eg.db, keyType = "ENTREZID")
        plot_kegg <- dotplot(ora_kegg, showCategory = category_top, font.size = label.size, label_format = label_format) + ggtitle(main) +
          scale_fill_gradient(limits = legendlimit, low = "red", high = "blue") + theme_plot
        ggsave(paste(filename,"_KEGG_dotplot.", imagetype, sep = ""), device = imagetype, path = paste(out.dir, sep = ""),
               plot = egg::set_panel_size(plot_kegg,width = unit(width, "cm"), height = unit(height, "cm")), width = width+10, height = height+2, units = "cm", dpi = dpi)
        write.table(data.frame(ora_kegg), file = paste(out.dir, "/", filename, "_KEGG.csv", sep = ""), sep = ",", row.names = FALSE)
      }
      ora.list[["KEGG"]] <- ora_kegg
    }
  }
  if(REACTOME){
    ora_reactome <- try(enrichPathway(gene_bitr$ENTREZID, organism = "human", readable = T, pvalueCutoff = p.cut))
    if(class(ora_reactome)!="try-error"){
      categories <- nrow(data.frame(ora_reactome))
      if(categories > 0){
        height <- ifelse(categories > 20, 20, categories)
        plot_reactome <- dotplot(ora_reactome, showCategory = category_top, font.size = label.size, label_format = label_format) + ggtitle(main) +
          scale_fill_gradient(limits = legendlimit, low = "red", high = "blue") + theme_plot
        ggsave(paste(filename,"_REACTOME_dotplot.", imagetype, sep = ""), device = imagetype, path = paste(out.dir, sep = ""),
               plot = egg::set_panel_size(plot_reactome, width = unit(width, "cm"), height = unit(height, "cm")), width = width+10, height = height+2, units = "cm", dpi = dpi)
        write.table(data.frame(ora_reactome), file = paste(out.dir, "/", filename, "_REACTOME.csv", sep = ""), sep = ",", row.names = FALSE)
      }
      ora.list[["REACTOME"]] <- ora_reactome
    }
  }
  return(ora.list)
}


# Gene Set Enrichment Analysis
calc_gsea <- function(res, name, ont = "BP", sort.by = "stat", KEGG = T, GO = T, REACTOME = F, nPerm = 1000, p.cut = 0.05, out.dir = "GSEA/"){
  if(!dir.exists(out.dir)){
    dir.create(out.dir)
  }
  gsea.list <- list()
  gene_bitr <- bitr(res$SYMBOL, fromType="SYMBOL", toType="ENTREZID", "org.Hs.eg.db")
  res <- merge(res, gene_bitr, by = "SYMBOL")
  res <- res[order(res[,sort.by], decreasing = TRUE),]
  geneset_num <- as.numeric(res[,sort.by])
  names(geneset_num) <- res$ENTREZID
  if(KEGG){
    gsea_kegg <- try(gseKEGG(geneList = geneset_num, organism = "hsa", pvalueCutoff = p.cut, seed = T))
    if(class(gsea_kegg) != "try-error" & nrow(data.frame(gsea_kegg)) > 0){
      plot_kegg <- dotplot(gsea_kegg, showCategory = 20, split = ".sign", font.size=10, label_format = function(x) stringr::str_wrap(x, width=30)) + 
        facet_grid(.~.sign) + theme(axis.title = element_text(size=18, face = "bold"), axis.text = element_text(size=10, face = "plain"),
                                    strip.text.x = element_text(size = 18, face = "bold"))
      ggsave(paste(name,"_KEGG_dotplot.pdf", sep = ""), device = "pdf", plot = plot_kegg, path = paste(out.dir, sep = ""), width = 18, height = 18)
      ggsave(paste(name,"_KEGG_dotplot.png", sep = ""), device = "png", plot = plot_kegg, path = paste(out.dir, sep = ""), width = 18, height = 18)
      gsea_kegg.df <- data.frame(gsea_kegg)
      gsea_kegg.df$core_enrichment_SYMBOL <- sapply(gsea_kegg.df$core_enrichment, function(x){ x.split <- unlist(strsplit(x,"/")); return(paste(bitr(x.split, "ENTREZID", "SYMBOL", "org.Hs.eg.db")[["SYMBOL"]], collapse = "/"))})
      gsea.list[["KEGG"]] <- gsea_kegg
    }else{
      gsea_kegg.df <- data.frame()
    }
    write.csv(gsea_kegg.df, file = paste(out.dir, "/", name, "_KEGG.csv", sep = ""), row.names = FALSE)
    write.xlsx(gsea_kegg.df, file = paste(out.dir, "/", name, "_KEGG.xlsx", sep = ""), rowNames = FALSE)
  }
  if(GO){
    for(o in ont){
      gsea_go <- try(gseGO(geneset_num, ont = o, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", pvalueCutoff = p.cut, seed = T))
      if(class(gsea_go) != "try-error" & nrow(data.frame(gsea_go)) > 0){
        plot_go <- dotplot(gsea_go, showCategory = 20, split = ".sign", font.size=10, label_format = function(x) stringr::str_wrap(x, width=30)) + 
          facet_grid(.~.sign) + theme(axis.title = element_text(size=18, face = "bold"), axis.text = element_text(size=15), 
                                      strip.text.x = element_text(size = 18, face = "bold"))
        ggsave(paste(name, "_GO_", o, "_dotplot.pdf", sep = ""), device = "pdf", plot = plot_go, path = paste(out.dir, sep = ""), width = 18, height = 18)
        ggsave(paste(name, "_GO_", o, "_dotplot.png", sep = ""), device = "png", plot = plot_go, path = paste(out.dir, sep = ""), width = 18, height = 18)
        gsea_go.df <- data.frame(gsea_go)
        gsea_go.df$core_enrichment_SYMBOL <- sapply(gsea_go.df$core_enrichment, function(x){ x.split <- unlist(strsplit(x,"/")); return(paste(bitr(x.split, "ENTREZID", "SYMBOL", "org.Hs.eg.db")[["SYMBOL"]], collapse = "/"))})
        gsea.list[[o]] <- gsea_go
      }else{
        gsea_go.df <- data.frame()
      }
      write.csv(gsea_go.df, file = paste(out.dir, "/", name, "_GO_", o, ".csv", sep = ""), row.names = FALSE)
      write.xlsx(gsea_go.df, file = paste(out.dir, "/", name, "_GO_", o, ".xlsx", sep = ""), rowNames = FALSE)
    }
  }
  if(REACTOME){
    gsea_reactome <- try(gsePathway(geneset_num, pvalueCutoff = p.cut, seed = T))
    if(class(gsea_reactome) != "try-error" & nrow(data.frame(gsea_reactome)) > 0){
      plot_reactome <- dotplot(gsea_reactome, showCategory = 20, split = ".sign", font.size=10, label_format = function(x) stringr::str_wrap(x, width=30)) + 
        facet_grid(.~.sign) + theme(axis.title = element_text(size=18, face = "bold"), axis.text = element_text(size=15), 
                                    strip.text.x = element_text(size = 18, face = "bold"))
      ggsave(paste(name,"_REACTOME_dotplot.pdf", sep = ""), device = "pdf", plot = plot_reactome, path = paste(out.dir, sep = ""), width = 18, height = 18)
      ggsave(paste(name,"_REACTOME_dotplot.png", sep = ""), device = "png", plot = plot_reactome, path = paste(out.dir, sep = ""), width = 18, height = 18)
      gsea_reactome.df <- data.frame(gsea_reactome)
      gsea_reactome.df$core_enrichment_SYMBOL <- sapply(gsea_reactome.df$core_enrichment, function(x){ x.split <- unlist(strsplit(x,"/")); return(paste(bitr(x.split, "ENTREZID", "SYMBOL", "org.Hs.eg.db")[["SYMBOL"]], collapse = "/"))})
      gsea.list[["REACTOME"]] <- gsea_reactome
    }else{
      gsea_reactome.df <- data.frame()
    }
    write.csv(gsea_reactome.df, file = paste(out.dir, "/", name, "_REACTOME.csv", sep = ""), row.names = FALSE)
    write.xlsx(gsea_reactome.df, file = paste(out.dir, "/", name, "_REACTOME.xlsx", sep = ""), rowNames = FALSE)
  }
  return(gsea.list)
}
