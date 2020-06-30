library("STRINGdb")

# besser cairo_pdf, funktioniert auf Laptop nicht...
string_ppi <- function(stringdb, gene.df, filename, top = 400, cluster = FALSE, cluster.algorithm = "fastgreedy", link = TRUE, out.dir = "", 
                            min.clust = 5, required_score = 0){
    if(!dir.exists(out.dir)){
    dir.create(out.dir, recursive = T)
  }
  gene.df <- string_db$map(gene.df, "SYMBOL", removeUnmappedRows = T)
  if(length(gene.df$STRING_id) > 0){
    pdf(paste(out.dir, "/PPI_enrichment_", filename, ".pdf", sep = ""))
    if(length(gene.df$STRING_id) > 1000){
      string_db$plot_ppi_enrichment(gene.df$STRING_id[1:1000])
    }else{
      if(length(gene.df$STRING_id) > 20){
        string_db$plot_ppi_enrichment(gene.df$STRING_id)        
      }else{
        string_db$plot_ppi_enrichment(gene.df$STRING_id, sliceWindow = length(gene.df$STRING_id))
      }
    }
    dev.off()
    if(length(gene.df$STRING_id) > top){
      hits <- gene.df$STRING_id[1:top]
    }else{
      hits <- gene.df$STRING_id
    }
    pdf(paste(out.dir, "/", filename, ".pdf", sep = ""))
    string_db$plot_network(hits, add_link = link, required_score = required_score)
    dev.off()
    if(cluster){
      string.clust <- stringdb$get_clusters(hits, algorithm = cluster.algorithm)
      pdf(paste(out.dir, "/", filename, "_Cluster.pdf", sep = ""), onefile = T)
      for(i in 1:length(string.clust)){
        if(length(string.clust[[i]]) > min.clust){
          string_db$plot_network(string.clust[[i]], add_link = link, required_score = required_score)
        }
      }
      dev.off()
    }
    #enrichmentGO <- string_db$get_enrichment(gene.df$STRING_ID, category = "Process", methodMT = "fdr", iea = FALSE)
    #enrichmentKEGG <- string_db$get_enrichment(gene.df$STRING_ID, category = "KEGG", methodMT = "fdr", iea = FALSE)
    #return(list(GO=enrichmentGO, KEGG=enrichmentKEGG, DF=gene.df))
  }
}

if(!dir.exists(paste(output_folder, "STRINGdb", sep = ""))){
  dir.create(paste(output_folder, "STRINGdb", sep = ""))
}
string_db <- STRINGdb$new(version="10", species=9606, input_directory=paste0(output_folder,"STRINGdb/"))


  