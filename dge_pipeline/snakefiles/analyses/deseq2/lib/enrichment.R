library(clusterProfiler)
#library(ReactomePA)
library(org.Hs.eg.db)
library(ggplot2)
#library(stringr)

# Over-representation analysis for a set of genes
calc_ora <- function(geneset, main = "", filename, out.dir = "ORA", GO = T, KEGG = T, REACTOME = F, ont = "BP", p.cut = 0.05, label.size = 12, 
                     legend.size = 8, legend.title.size = 8, keytype = "SYMBOL", width = 18, height = 15, imagetype = "svg"){
  if(!dir.exists(out.dir)){
    dir.create(out.dir)
  }
  if(keytype != "ENTREZID"){
    gene_bitr <- bitr(geneset, fromType=keytype, toType="ENTREZID", "org.Hs.eg.db")
  }else{
    gene_bitr <- data.frame(ENTREZID=geneset)
  }
  ora.list <- list()
  if(GO){
    for(o in ont){
      ora_go <- try(enrichGO(gene_bitr$ENTREZID, keyType = "ENTREZID", OrgDb = "org.Hs.eg.db", ont = o, readable = T, pvalueCutoff = p.cut))
      if(class(ora_go)!="try-error"){
        if(nrow(data.frame(ora_go)) > 0){
          plot_go <- enrichplot::dotplot(ora_go, showCategory = 20, font.size = label.size) + ggtitle(main) +
            theme(legend.text = element_text(size = legend.size, family = "sans"), legend.title = element_text(size = legend.title.size, face = "bold", family = "sans"),
                  axis.title = element_text(size = legend.title.size, face = "bold", family = "sans")) #+ scale_y_discrete(labels=function(x)str_wrap(x, width = 20))
          ggsave(paste(filename,"_",o,"_dotplot.", imagetype, sep = ""), device = imagetype, plot = plot_go, path = paste(out.dir, sep = ""), width = width, height = height)
          write.table(data.frame(ora_go), file = paste(out.dir, "/", filename, "_",o,".csv", sep = ""), sep = "\t", row.names = FALSE)
        }
        ora.list[[o]] <- ora_go
      }
    }
  }
  if(KEGG){
    ora_kegg <- try(enrichKEGG(gene_bitr$ENTREZID, organism = "hsa", use_internal_data = FALSE, pvalueCutoff = p.cut))
    if(class(ora_kegg)!="try-error"){
      if(nrow(data.frame(ora_kegg)) > 0){
        ora_kegg <- setReadable(ora_kegg, org.Hs.eg.db, keyType = "ENTREZID")
        plot_kegg <- enrichplot::dotplot(ora_kegg, showCategory = 20, font.size = label.size) + ggtitle(main) +
          theme(legend.text = element_text(size = legend.size, family = "sans"), legend.title = element_text(size = legend.title.size, face = "bold", family = "sans"),
                axis.title = element_text(size = label.size, face = "bold", family = "sans")) #+ scale_y_discrete(labels=function(x)str_wrap(x, width = 20))
        ggsave(paste(filename,"_KEGG_dotplot.", imagetype, sep = ""), device = imagetype, plot = plot_kegg, path = paste(out.dir, sep = ""), width = width, height = height)
        write.table(data.frame(ora_kegg), file = paste(out.dir, "/", filename, "_KEGG.csv", sep = ""), sep = "\t", row.names = FALSE)
      }
      ora.list[["KEGG"]] <- ora_kegg
    }
  }
  if(REACTOME){
    ora_reactome <- try(enrichPathway(gene_bitr$ENTREZID, organism = "human", readable = T, pvalueCutoff = p.cut))
    if(class(ora_reactome)!="try-error"){
      if(nrow(data.frame(ora_reactome)) > 0){
        plot_reactome <- enrichplot::dotplot(ora_reactome, showCategory = 20, font.size = label.size) + ggtitle(main) +
          theme(legend.text = element_text(size = legend.size, family = "sans"), legend.title = element_text(size = legend.title.size, face = "bold", family = "sans"),
                axis.title = element_text(size = legend.title.size, face = "bold", family = "sans")) #+ scale_y_discrete(labels=function(x)str_wrap(x, width = 20))
        ggsave(paste(filename,"_REACTOME_dotplot.", imagetype, sep = ""), device = imagetype, plot = plot_reactome, path = paste(out.dir, sep = ""), width = width, height = height)
        write.table(data.frame(ora_reactome), file = paste(out.dir, "/", filename, "_REACTOME.csv", sep = ""), sep = "\t", row.names = FALSE)
      }
      ora.list[["REACTOME"]] <- ora_reactome
    }
  }
  return(ora.list)
}

# Gene Set Enrichment Analysis
calc_gsea <- function(geneset.df, filename, ont = "BP", sort.by = "stat", KEGG = T, GO = T, REACTOME = F, nPerm = 1000, p.cut = 0.05, 
                      out.dir = "GSEA/", keytype = "SYMBOL", key.colname = "SYMBOL"){
  if(!dir.exists(out.dir)){
    dir.create(out.dir)
  }
  gsea.list <- list()
  if(keytype != "ENTREZID"){
    gene_bitr <- bitr(geneset.df[,key.colname], fromType=keytype, toType="ENTREZID", "org.Hs.eg.db")
    geneset.df <- merge(geneset.df, gene_bitr, by.x = key.colname, by.y = keytype)
  }else{
    colnames(geneset.df)[colnames(geneset.df)==key.colname] <- "ENTREZID"
  }
  geneset.df <- geneset.df[order(geneset.df[,sort.by], decreasing = TRUE),]
  geneset_num <- as.numeric(geneset.df[,sort.by])
  names(geneset_num) <- geneset.df$ENTREZID
  
  if(KEGG){
    gsea_kegg <- try(gseKEGG(geneList = geneset_num, organism = "hsa", nPerm = nPerm, pvalueCutoff = p.cut, seed = T))
    if(class(gsea_kegg)!="try-error"){
      if(nrow(data.frame(gsea_kegg)) > 0){
        plot_kegg <- enrichplot::dotplot(gsea_kegg, showCategory = 20, split = ".sign") + facet_grid(.~.sign) + 
          theme(axis.title = element_text(size=18, face = "bold"), axis.text = element_text(size=15, face = "bold"), 
                strip.text.x = element_text(size = 18, face = "bold"))
        ggsave(paste(filename,"_KEGG_dotplot.svg", sep = ""), device = "svg", plot = plot_kegg, path = paste(out.dir, sep = ""), width = 18, height = 15)
        gsea_kegg.df <- data.frame(gsea_kegg)
        gsea_kegg.df$core_enrichment_SYMBOL <- sapply(gsea_kegg.df$core_enrichment, function(x){ x.split <- unlist(strsplit(x,"/")); return(paste(bitr(x.split, "ENTREZID", "SYMBOL", "org.Hs.eg.db")[["SYMBOL"]], collapse = "/"))})
        write.table(gsea_kegg.df, file = paste(out.dir, "/", filename, "_KEGG.csv", sep = ""), sep = ",", row.names = FALSE)
        gsea.list[["KEGG"]] <- gsea_kegg
      }
    }
  }
  if(GO){
    for(o in ont){
      gsea_go <- try(gseGO(geneset_num, ont = o, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", nPerm = nPerm, pvalueCutoff = p.cut, seed = T))
      if(class(gsea_go)!="try-error"){
        if(nrow(data.frame(gsea_go)) > 0){
          plot_go <- enrichplot::dotplot(gsea_go, showCategory = 20, split = ".sign") + facet_grid(.~.sign) + 
            theme(axis.title = element_text(size=18, face = "bold"), axis.text = element_text(size=15, face = "bold"), 
                  strip.text.x = element_text(size = 18, face = "bold"))
          ggsave(paste(filename, "_GO_", o, "_dotplot.svg", sep = ""), device = "svg", plot = plot_go, path = paste(out.dir, sep = ""), width = 18, height = 15)
          gsea_go.df <- data.frame(gsea_go)
          gsea_go.df$core_enrichment_SYMBOL <- sapply(gsea_go.df$core_enrichment, function(x){ x.split <- unlist(strsplit(x,"/")); return(paste(bitr(x.split, "ENTREZID", "SYMBOL", "org.Hs.eg.db")[["SYMBOL"]], collapse = "/"))})
          write.table(gsea_go.df, file = paste(out.dir, "/", filename, "_GO_", o, ".csv", sep = ""), sep = ",", row.names = FALSE)
          gsea.list[[o]] <- gsea_go
          # Problem: Bei gleichem value von by werden alle Pathways mit diesem Wert ausgegeben!
          # Loesung: by = pvalue?
          gsea.list[[paste0(o,"_simplify")]] <- simplify(gsea_go, cutoff = 0.5, by = "pvalue")
        }
      }
    }
  }
  if(REACTOME){
    gsea_reactome <- try(gsePathway(geneset_num, nPerm = nPerm, pvalueCutoff = p.cut, seed = T))
    if(class(gsea_reactome)!="try-error"){
      if(nrow(data.frame(gsea_reactome)) > 0){
        plot_reactome <- enrichplot::dotplot(gsea_reactome, showCategory = 20, split = ".sign") + facet_grid(.~.sign) + 
          theme(axis.title = element_text(size=18, face = "bold"), axis.text = element_text(size=15, face = "bold"), 
                strip.text.x = element_text(size = 18, face = "bold"))
        ggsave(paste(filename,"_REACTOME_dotplot.svg", sep = ""), device = "svg", plot = plot_reactome, path = paste(out.dir, sep = ""), width = 18, height = 15)
        gsea_reactome.df <- data.frame(gsea_reactome)
        gsea_reactome.df$core_enrichment_SYMBOL <- sapply(gsea_reactome.df$core_enrichment, function(x){ x.split <- unlist(strsplit(x,"/")); return(paste(bitr(x.split, "ENTREZID", "SYMBOL", "org.Hs.eg.db")[["SYMBOL"]], collapse = "/"))})
        write.table(gsea_reactome.df, file = paste(out.dir, "/", filename, "_REACTOME.csv", sep = ""), sep = ",", row.names = FALSE)
        gsea.list[["REACTOME"]] <- gsea_reactome
      }
    }
  }
  return(gsea.list)
}

# Comparative overrepresentation analysis
calc_compareCluster <- function(dataset, filename, out.dir = "ORA", GO = T, KEGG = T, REACTOME = F, ont = "BP", p.cut = 0.05, label.size = 12, 
                                legend.size = 8, main = "", w = 18, h = 15, keytype = "SYMBOL"){
  if(!dir.exists(out.dir)){
    dir.create(out.dir)
  }
  if(keytype != "ENTREZID"){
    if(class(dataset)=="list"){
      dataset <- lapply(dataset, function(x)bitr(x, keytype, "ENTREZID", org.Hs.eg.db)$ENTREZID)
    }
  }
  if(GO){
    for(o in ont){
      cc_go <- try(compareCluster(dataset, fun = "enrichGO", keyType = "ENTREZID", OrgDb = "org.Hs.eg.db", ont = o, readable = T, pvalueCutoff = p.cut))
      if(class(cc_go)!="try-error"){
        if(nrow(data.frame(cc_go)) > 0){
          plot_go <- enrichplot::dotplot(cc_go, showCategory = 20, font.size = label.size) + ggtitle(main) +
            theme(legend.text = element_text(size = legend.size), legend.title = element_text(size = legend.size, face = "bold")) #+ scale_y_discrete(labels=function(x)str_wrap(x, width = 20))
          ggsave(paste(filename,"_",o,"_dotplot.svg", sep = ""), device = "svg", plot = plot_go, path = out.dir, width = w, height = h)
          write.table(data.frame(cc_go), file = paste(out.dir, "/", filename, "_",o,".csv", sep = ""), sep = "\t", row.names = FALSE)
        }
      }
    }
  }
  if(KEGG){
    cc_kegg <- try(compareCluster(dataset, "enrichKEGG", organism = "hsa", use_internal_data = FALSE, pvalueCutoff = p.cut))
    if(class(cc_kegg)!="try-error"){
      if(nrow(data.frame(cc_kegg)) > 0){
        plot_kegg <- enrichplot::dotplot(cc_kegg, showCategory = 20, font.size = label.size) + ggtitle(main) +
          theme(legend.text = element_text(size = legend.size), legend.title = element_text(size = legend.size, face = "bold")) #+ scale_y_discrete(labels=function(x)str_wrap(x, width = 20))
        ggsave(paste(filename,"_KEGG_dotplot.svg", sep = ""), device = "svg", plot = plot_kegg, path = out.dir, width = w, height = h)
        write.table(data.frame(cc_kegg), file = paste(out.dir, "/", filename, "_KEGG.csv", sep = ""), sep = "\t", row.names = FALSE)
      }
    }
  }
  if(REACTOME){
    cc_reactome <- try(compareCluster(dataset, "enrichPathway", organism = "human", readable = T, pvalueCutoff = p.cut))
    if(class(cc_reactome)!="try-error"){
      if(nrow(data.frame(cc_reactome)) > 0){
        plot_reactome <- enrichplot::dotplot(cc_reactome, showCategory = 20, font.size = label.size) + ggtitle(main) +
          theme(legend.text = element_text(size = legend.size), legend.title = element_text(size = legend.size, face = "bold")) #+ scale_y_discrete(labels=function(x)str_wrap(x, width = 20))
        ggsave(paste(filename,"_REACTOME_dotplot.svg", sep = ""), device = "svg", plot = plot_reactome, path = out.dir, width = w, height = h)
        write.table(data.frame(cc_reactome), file = paste(out.dir, "/", filename, "_REACTOME.csv", sep = ""), sep = "\t", row.names = FALSE)
      }
    }
  }
}
