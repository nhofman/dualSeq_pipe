library(clusterProfiler)

# Over-representation analysis for a set of genes
calc_ora <- function(gene, main = "", filename, out.dir = "ORA", GO = T, KEGG = T, REACTOME = F, ont = "BP", p.cut = 0.05){
  if(!dir.exists(out.dir)){
    dir.create(out.dir)
  }
  gene_bitr <- bitr(gene, fromType="SYMBOL", toType="ENTREZID", "org.Hs.eg.db")
  ora.list <- list()
  if(GO){
    for(o in ont){
      ora_go <- try(enrichGO(gene_bitr$ENTREZID, keyType = "ENTREZID", OrgDb = "org.Hs.eg.db", ont = o, readable = T, pvalueCutoff = p.cut))
      if(nrow(data.frame(ora_go)) > 0){
        plot_go <- dotplot(ora_go, showCategory = 20) + ggtitle(main)
        ggsave(paste(filename,"_",o,"_dotplot.svg", sep = ""), device = "svg", plot = plot_go, path = paste(out.dir, sep = ""), width = 18, height = 15)
        write.table(data.frame(ora_go), file = paste(out.dir, "/", filename, "_",o,".csv", sep = ""), sep = "\t", row.names = FALSE)
      }
      ora.list[[o]] <- ora_go
    }
  }
  if(KEGG){
    ora_kegg <- try(enrichKEGG(gene_bitr$ENTREZID, organism = "hsa", use_internal_data = FALSE, pvalueCutoff = p.cut))
    if(nrow(data.frame(ora_kegg)) > 0){
      ora_kegg <- setReadable(ora_kegg, org.Hs.eg.db, keyType = "ENTREZID")
      plot_kegg <- dotplot(ora_kegg, showCategory = 20) + ggtitle(main)
      ggsave(paste(filename,"_KEGG_dotplot.svg", sep = ""), device = "svg", plot = plot_kegg, path = paste(out.dir, sep = ""), width = 18, height = 15)
      write.table(data.frame(ora_kegg), file = paste(out.dir, "/", filename, "_KEGG.csv", sep = ""), sep = "\t", row.names = FALSE)
    }
    ora.list[["KEGG"]] <- ora_kegg
  }
  if(REACTOME){
    ora_reactome <- try(enrichPathway(gene_bitr$ENTREZID, organism = "human", readable = T, pvalueCutoff = p.cut))
    if(nrow(data.frame(ora_reactome)) > 0){
      plot_reactome <- dotplot(ora_reactome, showCategory = 20) + ggtitle(main)
      ggsave(paste(filename,"_REACTOME_dotplot.svg", sep = ""), device = "svg", plot = plot_reactome, path = paste(out.dir, sep = ""), width = 18, height = 15)
      write.table(data.frame(ora_reactome), file = paste(out.dir, "/", filename, "_REACTOME.csv", sep = ""), sep = "\t", row.names = FALSE)
    }
    ora.list[["REACTOME"]] <- ora_reactome
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
    gsea_kegg <- try(gseKEGG(geneList = geneset_num, organism = "hsa", nPerm = nPerm, pvalueCutoff = p.cut, seed = T))
    if(nrow(data.frame(gsea_kegg)) > 0){
      plot_kegg <- dotplot(gsea_kegg, showCategory = 20, split = ".sign") + facet_grid(.~.sign) + theme(axis.title = element_text(size=18, face = "bold"), axis.text = element_text(size=15, face = "bold"), strip.text.x = element_text(size = 18, face = "bold"))
      ggsave(paste(name,"_KEGG_dotplot.svg", sep = ""), device = "svg", plot = plot_kegg, path = paste(out.dir, sep = ""), width = 18, height = 15)
      gsea_kegg.df <- data.frame(gsea_kegg)
      gsea_kegg.df$core_enrichment_SYMBOL <- sapply(gsea_kegg.df$core_enrichment, function(x){ x.split <- unlist(strsplit(x,"/")); return(paste(bitr(x.split, "ENTREZID", "SYMBOL", "org.Hs.eg.db")[["SYMBOL"]], collapse = "/"))})
      write.table(gsea_kegg.df, file = paste(out.dir, "/", name, "_KEGG.csv", sep = ""), sep = ",", row.names = FALSE)
      gsea.list[["KEGG"]] <- gsea_kegg
    }
  }
  if(GO){
    for(o in ont){
      gsea_go <- try(gseGO(geneset_num, ont = o, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", nPerm = nPerm, pvalueCutoff = p.cut, seed = T))
      if(nrow(data.frame(gsea_go)) > 0){
        plot_go <- dotplot(gsea_go, showCategory = 20, split = ".sign") + facet_grid(.~.sign) + theme(axis.title = element_text(size=18, face = "bold"), axis.text = element_text(size=15, face = "bold"), strip.text.x = element_text(size = 18, face = "bold"))
        ggsave(paste(name, "_GO_", o, "_dotplot.svg", sep = ""), device = "svg", plot = plot_go, path = paste(out.dir, sep = ""), width = 18, height = 15)
        gsea_go.df <- data.frame(gsea_go)
        gsea_go.df$core_enrichment_SYMBOL <- sapply(gsea_go.df$core_enrichment, function(x){ x.split <- unlist(strsplit(x,"/")); return(paste(bitr(x.split, "ENTREZID", "SYMBOL", "org.Hs.eg.db")[["SYMBOL"]], collapse = "/"))})
        write.table(gsea_go.df, file = paste(out.dir, "/", name, "_GO_", o, ".csv", sep = ""), sep = ",", row.names = FALSE)
        gsea.list[[o]] <- gsea_go
        # Problem: Bei gleichem value von by werden alle Pathways mit diesem Wert ausgegeben!
        # Loesung: by = pvalue?
        gsea.list[[paste0(o,"_simplify")]] <- simplify(gsea_go, cutoff = 0.5, by = "pvalue")
      }
    }
  }
  if(REACTOME){
    gsea_reactome <- try(gsePathway(geneset_num, nPerm = nPerm, pvalueCutoff = p.cut, seed = T))
    if(nrow(data.frame(gsea_reactome)) > 0){
      plot_reactome <- dotplot(gsea_reactome, showCategory = 20, split = ".sign") + facet_grid(.~.sign) + theme(axis.title = element_text(size=18, face = "bold"), axis.text = element_text(size=15, face = "bold"), strip.text.x = element_text(size = 18, face = "bold"))
      ggsave(paste(name,"_REACTOME_dotplot.svg", sep = ""), device = "svg", plot = plot_reactome, path = paste(out.dir, sep = ""), width = 18, height = 15)
      gsea_reactome.df <- data.frame(gsea_reactome)
      gsea_reactome.df$core_enrichment_SYMBOL <- sapply(gsea_reactome.df$core_enrichment, function(x){ x.split <- unlist(strsplit(x,"/")); return(paste(bitr(x.split, "ENTREZID", "SYMBOL", "org.Hs.eg.db")[["SYMBOL"]], collapse = "/"))})
      write.table(gsea_reactome.df, file = paste(out.dir, "/", name, "_REACTOME.csv", sep = ""), sep = ",", row.names = FALSE)
      gsea.list[["REACTOME"]] <- gsea_reactome
    }
  }
  return(gsea.list)
}