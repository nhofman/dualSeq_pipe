set.seed(123)
# Parse arguments
args <- commandArgs(F)
file.dir <- dirname(sub("--file=","",args[grep("--file=",args)]))
rdata <- args[match('--rdata', args) + 1]
threads <- args[match('--threads', args) + 1]

load(rdata)

source(paste0(file.dir, "/enrichment.R"))
source(paste0(file.dir, "/STRINGdb.R"))

gsea.list <- list()
for(n in names(res.list)){
  print(n)
  res <- res.list[[n]]
  res <- res[rowSums(res[,grep("normalized",colnames(res))])>0,] # remove genes with no read counts
  res <- res[order(res$log2FoldChange, decreasing = T),]
  # GSEA analysis
  gsea.list[[n]] <- calc_gsea(res, n, sort.by = "log2FoldChange", REACTOME = T, ont = c("CC", "MF", "BP"),
  p.cut = 0.05, out.dir = paste0(output_folder,"/GSEA"))
  # Protein-protein interaction analysis with STRING
  if(nrow(data.frame("SYMBOL"=res$SYMBOL[res$log2FoldChange>0 & res$padj < 0.05])) > 0){
          try(string_ppi(string_db, gene.df = data.frame("SYMBOL"=res$SYMBOL[res$log2FoldChange>0 & res$padj < 0.05]), filename = paste0(n, "_up"), out.dir = paste0(output_folder, "/STRING")))
  }
  if(nrow(data.frame("SYMBOL"=res$SYMBOL[res$log2FoldChange<0 & res$padj < 0.05])) > 0){
          try(string_ppi(string_db, gene.df = data.frame("SYMBOL"=res$SYMBOL[res$log2FoldChange<0 & res$padj < 0.05]), filename = paste0(n, "_down"), out.dir = paste0(output_folder, "/STRING")))
  }
  res <- res[order(abs(res$log2FoldChange), decreasing = T),]
  if(nrow(res) > 0){
          try(string_ppi(string_db, gene.df = data.frame("SYMBOL"=res$SYMBOL), filename = paste0(n, "_top_absolut"), out.dir = paste0(output_folder, "/STRING")))
  }
  gc()
}

