set.seed(123)
# Parse arguments
args <- commandArgs(F)
file.dir <- dirname(sub("--file=","",args[grep("--file=",args)]))
input.file <- args[match('--datatable', args) + 1]
threads <- args[match('--threads', args) + 1]
db <- args[match('--db', args) + 1]
ont <- args[match('--ont', args) + 1]
id.type <- args[match('--id_type', args) + 1]
output_folder <- args[match('--output_folder', args) + 1]

source(paste0(file.dir, "/enrichment.R"))

db <- strsplit(db, ",")[[1]]
KEGG <- "KEGG" %in% db
GO <- "GO" %in% db
REACTOME <- "REACTOME" %in% db
ont <- strsplit(ont, ",")[[1]]

# sort genes by LFC and call GSEA
n <- gsub("deseq2_results_(.*)\\.csv", "\\1", basename(input.file))
res <- read.csv(input.file, header = T, row.names = 1) 
res <- res[rowSums(res[,grep("normalized", colnames(res))])>0,] # remove genes with no read counts
res <- res[order(res$log2FoldChange, decreasing = T),]
# GSEA analysis
gsea <- calc_gsea(res, n, sort.by = "log2FoldChange", REACTOME = REACTOME, KEGG = KEGG, GO = GO, ont = ont, keytype = id.type,
p.cut = 0.05, out.dir = paste0(output_folder, "/GSEA"))
  
