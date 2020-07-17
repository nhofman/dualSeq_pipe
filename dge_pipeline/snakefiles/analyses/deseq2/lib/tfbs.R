library("GenomicFeatures")
library("rtracklayer")
library("Biostrings")
library("Rsamtools")

# calculate GRanges object from hg38 annotation
getTranscriptsfromGFF <- function(gffFilePath){
  # create transcript GRanges object
  txdb <- makeTxDbFromGFF(gffFilePath) # get txdb object from gff file
  transcripts.GRanges <- transcripts(txdb, columns = c("tx_id", "tx_name", "gene_id")) # create GRanges object of gff file
  tlen <- transcriptLengths(txdb, with.cds_len = T, with.utr5_len = T, with.utr3_len = T)
  values(transcripts.GRanges) <- cbind(values(transcripts.GRanges), tlen[, -c(1,2,3)])
  # remove duplicated gene entries
  x <- unlist(values(transcripts.GRanges)$gene_id)
  unique_gene_ids_groups <- lapply(levels(as.factor(x)), function(y) which(y==x)) # get list with indexes of all unique gene_ids
  index_max_unique_gene_ids <- lapply(unique_gene_ids_groups, function(x) which(values(transcripts.GRanges)[x, "tx_len"]==max(values(transcripts.GRanges)[x, "tx_len"])))# get max genelength index for all unique gene_ids indexes
  index <- vector("numeric", length(index_max_unique_gene_ids))
  for(i in 1:length(index_max_unique_gene_ids)){
    ind <- index_max_unique_gene_ids[[i]]
    ind <- unique(ind)
    index <- c(index, unique_gene_ids_groups[[i]][ind[1]])
  }
  transcripts.GRanges <- transcripts.GRanges[index, ]
  return(transcripts.GRanges)
}

# Get promotor regions of genes in geneset and write them to fasta-file
writeGRanges2Fasta <- function(GRanges, outdir, promotor_prim, promotor_back, fasta.file, geneset, upstream = 1000, downstream = 100){
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }
  human.promoter <- promoters(human.GRanges, upstream = upstream, downstream = downstream)
  human.promoter.select <- subset(human.promoter, gene_id %in% geneset)
  human.promoter.bg <- subset(human.promoter, !(gene_id %in% geneset))
  
  #fa <- open.FaFile(FaFile(fasta.file))
  
  human.promoter.select.seq <- getSeq(FaFile(fasta.file), param = human.promoter.select)
  names(human.promoter.select.seq) <- paste(names(human.promoter.select.seq), as.character(human.promoter.select$gene_id), sep=":")
  writeXStringSet(human.promoter.select.seq, filepath = paste(outdir,promotor_prim,sep="/"))
  
  human.promoter.bg.seq <- getSeq(FaFile(fasta.file), param = human.promoter.bg)
  names(human.promoter.bg.seq) <- paste(names(human.promoter.bg.seq), as.character(human.promoter.bg$gene_id), sep=":")
  writeXStringSet(human.promoter.bg.seq, filepath = paste(outdir,promotor_back,sep="/"))
  
  #close.FaFile(FaFile(fasta.file))
}


# call meme
meme <- function(path2meme, outdir, promotor_prim, promotor_back, nmotif = 1){
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }
  outdir <- paste0(outdir, "/meme")
  bash_file <- "meme_command.sh"
  #if(file.exists(bash_file)) file.remove(bash_file)
  f <- file(bash_file, open = "w")
  meme_command <- paste(paste0(path2meme, "meme"),
         promotor_prim,
         paste0("-neg ", promotor_back),
         paste0("-oc ", outdir),
         "-objfun de",
         "-dna",
         paste0("-nmotifs ", nmotif), sep = " ")
  cat(paste("echo \"", meme_command, "\"\n"), file=f)
  cat(meme_command, file=f) 
  close(f)
  meme_cmd_out <- system2("bash", args = bash_file, stderr = TRUE, stdout = TRUE)
  cat(paste(meme_cmd_out, collapse = "\n"))
  cat(paste("mv", bash_file, outdir)) 
  return(paste0(outdir, "/meme.html"))
}

# call tomtom
tomtom <- function(path2meme, motifDB, motif_file, outdir){
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }
  outdir <- paste0(outdir, "/tomtom")
  bash_file <- "tomtom_command.sh"
  #if(file.exists(bash_file)) file.remove(bash_file)
  f <- file(bash_file, open = "w")
  tomtom_command <- paste(paste0(path2meme, "tomtom"),
         paste0("-oc ", outdir),
         motif_file,
         motifDB, sep = " ")
  cat(paste("echo \"", tomtom_command, "\"\n"), file=f)
  cat(tomtom_command, file=f) 
  close(f)
  tomtom_cmd_out <- system2("bash", args = bash_file, stderr = TRUE, stdout = TRUE)
  cat(paste(tomtom_cmd_out, collapse = "\n"))
  system(paste("mv", bash_file, outdir)) 
}

# call ame
ame <- function(path2meme, motifDB, promotor_prim, promotor_back, outdir){
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }
  outdir <- paste0(outdir, "/ame")
  bash_file <- "ame_command.sh"
  #if(file.exists(bash_file)) file.remove(bash_file)
  f <- file(bash_file, open = "w")
  ame_command <- paste(paste0(path2meme, "ame"),
         paste0("-oc ", outdir),
         paste0("--control ",promotor_back),
         promotor_prim,
         motifDB, sep = " ")
  cat(paste("echo \"", ame_command, "\"\n"), file=f)
  cat(ame_command, file=f) 
  close(f)
  ame_cmd_out <- system2("bash", args = bash_file, stderr = TRUE, stdout = TRUE)
  cat(paste(ame_cmd_out, collapse = "\n"))
  cat(paste("mv", bash_file, outdir)) 
}
