library(vcfR)
library(ggplot2)
library(tidyr)
library(gtools)
library(Gviz)
library(lattice)

options(ucscChromosomeNames=FALSE)
#vcf.files <- list.files("Documents/Virus_project/variant_calling/", pattern = ".*h.*.vcf", full.names = T)
args <- commandArgs()
vcf <- args[match("--vcf", args) + 1]
gff_file <- args[match("--gff", args) + 1]
#genome_length <- gff[gff[,3]=="region",c(1,5)]
#vcf <- "Documents/Virus_project/variant_calling/HCV_12h_1.vcf,Documents/Virus_project/variant_calling/HCV_12h_2.vcf,Documents/Virus_project/variant_calling/HCV_24h_1.vcf,Documents/Virus_project/variant_calling/HCV_24h_2.vcf"
print(vcf)
vcf.files <- unlist(strsplit(vcf, ","))
print(vcf.files)
vcf.list <- lapply(vcf.files, read.vcfR, verbose = F)
names(vcf.list) <- sub(".vcf", "", basename(vcf.files))

gff <- read.table(gff_file, sep = "\t", comment.char = "#", stringsAsFactors = F, quote = "\"")
colnames(gff) <- c("chromosome","source","feature","start","end","score","strand","phase","attributes")
gff <- gff[-1,]
i <- 1
for(l in lapply(as.character(gff$attributes), FUN = strsplit, split=";", fixed=TRUE)){
  for(att in l[[1]]){
    if(grepl("Name", att)){
      gff$id[i] <- strsplit(att, "=", fixed=TRUE)[[1]][2]
    }
  }
  i <- i+1
}
gff <- gff[gff$feature=="gene",]

#sapply(unique(sub("_\\d+$", "", names(vcf.list))),function(x){
ymaxs <- c()
dt.all <- c()
snp.all <- c()
for(x in mixedsort(unique(sub("_\\d+$", "", names(vcf.list))))){
  #print(x)
  genome_name <- sub("_.*", "", x)
  gtrack <- GenomeAxisTrack(add35 = T, cex = 1)
  #ann <- AnnotationTrack(gff, feature = as.character(gff$type)) #, group = gff$type)
  #ann <- AnnotationTrack(gff[which(gff$feature=="CDS"),], chromosome = "HCV_Jc1", genome = "HCV", name = "Annotation", fill = "darkblue", shape = "arrow",
  #                       id = gff[which(gff$feature=="CDS"),"id"], cex = 0.75)
  ann <- AnnotationTrack(gff,  genome = genome_name, name = "Annotation", fill = "darkblue", shape = "arrow", chromosome = gff$chromosome[1],
                         id = gff[,"id"], cex = 0.75)
  dt.list <- list()
  snp.list <- list()
  sample_names <- names(vcf.list)[grep(x, names(vcf.list))]
  for(name in sample_names){
    bigwig_file <- paste0("Documents/Virus_project/bigWig/", name, ".bigWig")
    bw <- rtracklayer::import(bigwig_file)
    ymaxs <- c(ymaxs, max(bw$score))
    dt <- DataTrack(range = bw, name = paste("Coverage", x, sep = "\n"), genome = genome_name, type = "l", legend = TRUE, fontsize.legend = 14,
                     groups = factor(name, levels = sample_names), cex.axis = 1, cex.title = 1, chromosome = gff$chromosome[1])
    dt.list <- c(dt.list, dt)
    vcf.df <- as.data.frame(vcf.list[[name]]@fix, stringsAsFactors = F)
    vcf.df$AF <- extract.info(vcf.list[[name]], "AF", as.numeric = T)
    colnames(vcf.df)[1:2] <- c("chromosome", "start")
    vcf.df$start <- as.numeric(vcf.df$start)
    vcf.df$end <- vcf.df$start
    vcf.df <- vcf.df[,c("chromosome", "start", "end", "AF")]
    if(nrow(vcf.df) > 0){
      #snp_track <- DataTrack(data = extract.info(vcf.list[[name]], "AF", as.numeric = T), start = as.numeric(as.character(vcf.df$POS)), width = 1, 
      #                       name = "SNP allele frequency", genome = genome_name, type = "p", cex = 1, legend = TRUE, fontsize.legend = 14,
      #                       groups = factor(name, levels = sample_names), cex.axis = 1, cex.title = 1.5, chromosome = gff$chromosome)
      snp_track <- DataTrack(vcf.df, name = paste("SNP allele frequency", x, sep = "\n"), genome = genome_name, groups = factor(name, levels = sample_names), 
                             type = "p", cex = 1, legend = TRUE, fontsize.legend = 14, cex.axis = 1, cex.title = 1)
      snp.list <- c(snp.list, snp_track)
    }
  }
  dt_overlay <- OverlayTrack(trackList = dt.list)
  snp_overlay <- OverlayTrack(trackList = snp.list)
  
  dt.all <- c(dt.all, dt_overlay)
  snp.all <- c(snp.all, snp_overlay)
  #svglite::svglite(paste0("Documents/Virus_project/variant_calling/Gviz_", x, ".svg"), width = 15, height = 10)
  #print(xyplot(1 ~ chromosome | chromosome, data = chroms, panel = function(x) {
  #  plotTracks(list(snp_overlay, dt_overlay, ann, gtrack), featureAnnotation = "id", 
  #             chromosome = x, add = TRUE, showId = F) },
  #  scales = list(draw = FALSE), xlab = NULL, ylab = NULL))
  #dev.off()

} #)
chroms <- data.frame(chromosome = unique(gff$chromosome), stringsAsFactors = F)
snp.all <- c(snp.all, ann, gtrack)
svglite::svglite(paste0("Documents/Virus_project/variant_calling/", genome_name, "_SNP.svg"), width = 15, height = 10)
print(xyplot(1 ~ chromosome | chromosome, data = chroms, panel = function(x) {
  plotTracks(snp.all, featureAnnotation = "id", ylim = c(0,1), 
             chromosome = x, add = TRUE, showId = F) },
  scales = list(draw = FALSE), xlab = NULL, ylab = NULL))
dev.off()

dt.all <- c(dt.all, ann, gtrack)
svglite::svglite(paste0("Documents/Virus_project/variant_calling/", genome_name, "_Coverage.svg"), width = 15, height = 10)
print(xyplot(1 ~ chromosome | chromosome, data = chroms, panel = function(x) {
  plotTracks(dt.all, featureAnnotation = "id", ylim = c(0,max(ymaxs)),
             chromosome = x, add = TRUE, showId = F) },
  scales = list(draw = FALSE), xlab = NULL, ylab = NULL))
dev.off()

    # bigwig_file_name <- sub(".bigWig", "", basename(bigwig_file))
  # bigwig_file2 <- paste0("Documents/Virus_project/variant_calling/", x, "_2.bigWig")
  # bigwig_file_name2 <- sub(".bigWig", "", basename(bigwig_file2))
  # dt2 <- DataTrack(range = bigwig_file2, genome = "HCV", name = "Coverage", chromosome = "HCV_Jc1", type = "l", legend = TRUE, fontsize.legend = 14,
  #                  groups = factor(bigwig_file_name2, levels = c(bigwig_file_name, bigwig_file_name2)), cex.axis = 1, cex.title = 1.5)
  # 
  # vcf.df <- as.data.frame(vcf.list[[bigwig_file_name]]@fix)
  # snp_track2 <- DataTrack(data = extract.info(vcf.list[[bigwig_file_name2]], "AF", as.numeric = T), start = as.numeric(as.character(vcf.df2$POS)), width = 1, strand = "+", 
  #                        genome = "HCV", name = "SNP allele frequency", chromosome = "HCV_Jc1", type = "p", cex = 1, legend = TRUE, fontsize.legend = 14,
  #                        groups = factor(bigwig_file_name2, levels = c(bigwig_file_name, bigwig_file_name2)), cex.axis = 1, cex.title = 1.5)
  # ylims <- extendrange(range(c(values(snp_track), values(snp_track2))))

# compare SNP between replicates and time points

# vcf.AF <- Reduce(rbind,lapply(names(vcf.list), function(x){
#   print(x)
#   snp.df <- Reduce(rbind, lapply(1:nrow(genome_length), function(x){df <- data.frame("CHROM" = genome_length[x,1], "POS" = 1:genome_length[x,2], "AF" = rep(NA, genome_length[x,2]))}))
#   af <- data.frame("CHROM"= getCHROM(vcf.list[[x]]), "AF" = extract.info(vcf.list[[x]], "AF", as.numeric = T), "POS" = getPOS(vcf.list[[x]]), stringsAsFactors = F)
#   if(nrow(af)>0){
#     for(i in 1:nrow(af)){
#       snp.df[snp.df$CHROM==af$CHROM[i] & snp.df$POS==af$POS[i],"AF"] <- af$AF[i]
#     }
#   }
#   snp.df$sample <- x
#   return(snp.df)
# }))
# 
# vcf.AF <- separate(vcf.AF, sample, c("Virus", "Time", "Rep"), "_", F)
# p <- ggplot(vcf.AF, aes(x = POS, y = AF, shape = Rep)) + geom_point() + facet_wrap(vars(factor(Time, levels = mixedsort(unique(vcf.AF$Time))), CHROM), nrow = 4)
# ggsave(paste0(unique(vcf.AF$Virus), ".svg"), p,"svg", "/home/nina/Documents/Virus_project/variant_calling/", width = 16, height = 8)
# #plot(genome.AF, type = "p", ylab = "SNP frequency", xlab = "Position", main = "CoV229E_24h_1")
# 
# #dna <- ape::read.dna("Documents/Virus_project/igv.js/test_igv.js/CoV229E.fasta", format = "fasta")
# #gff <- read.table("Documents/Virus_project/igv.js/test_igv.js/CoV229E.gff", sep = "\t", comment.char = "#")
# #chrom <- create.chromR(vcf, name = "CoV229E", seq = dna, ann = gff)
# #plot(chrom)
# #chromoqc(chrom)
# #chrom.proc  <- proc.chromR(chrom, win.size = 100)
# #chromoqc(chrom.proc)