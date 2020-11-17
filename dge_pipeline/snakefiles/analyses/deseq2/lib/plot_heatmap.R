library(pheatmap)

# Heatmap 
plotHeatmap <- function(x, filename = "no_name_set.pdf", row_subset = NA, distMethod = "euclidean", clusterMethod = "complete", clrn = 1, clcn = 1,
                        rowClust = T, colClust = T, fontsize_row = 0.8, fontsize_col = 10, annCol = NA, annRow = NA, border_col = "grey60", plot.fig = T,
                        legend.limit = 1, filter_col = NA, annotation_colors = NA, height = NA, width = NA, break_step = 0.1, display_numbers = F){
  if(!dir.exists(dirname(filename))){
    dir.create(dirname(filename), recursive = T)
  }
  if(is.na(row_subset[1])){
    xx <- x
  }else{
    xx <- x[row_subset,, drop = F]
  }
  #xx[xx==0] <- 0.000001
  xx <- na.omit(xx)
  if (!is.na(clrn) & !is.na(clcn)) {
    # set the custom distance and clustering functions
    hclustfunc <- function(x) hclust(x, method = clusterMethod)
    distfunc <- function(x) dist(x, method = distMethod)
    # perform clustering on rows and columns
    cl.row <- hclustfunc(distfunc(xx))
    cl.col <- hclustfunc(distfunc(t(xx)))
    
    # extract cluster assignments; i.e. k=8 (rows) k=5 (columns)
    gr.row <- cutree(cl.row, clrn)
    gr.col <- cutree(cl.col, clcn)
  }else{
    # set default values
    gr.row <- rep(1, nrow(xx))
    gr.col <- rep(1, ncol(xx))
  }
  if(!is.na(filter_col)){
    xx <- xx[,which(colMaxs(as.matrix(abs(xx)))>filter_col)]
  }
  #pdf(filename,width=25,height=25)
  #pdf(filename)
  #$par(oma=c(10,4,4,10) + 0.1)
  
  breakSeq <- seq(quantile(unlist(xx), na.rm = TRUE, probs = (1-legend.limit)), quantile(unlist(xx), na.rm = TRUE, probs = legend.limit), break_step)
  if(length(breakSeq)==1){
    if(sign(breakSeq)==1){
      breakSeq <- seq(0,breakSeq,break_step)
    }else{
      breakSeq <- seq(breakSeq,0,break_step)
    }
  }
  color <- vector()
  if(-1%in%sign(unlist(xx))){
    color_down <- colorRampPalette(c("blue", "white"))(length(seq(quantile(unlist(xx), na.rm = TRUE, probs = (1-legend.limit)), 0, break_step))) #blue(n=length(breakSeq)-1) #, low="blue", mid = "white", high="red")
    color <- c(color,color_down)
  }
  if(1%in%sign(unlist(xx))){
    color_up <- colorRampPalette(c("white", "red"))(length(seq(0, quantile(unlist(xx), na.rm = TRUE, probs = legend.limit), break_step)))
    color <- c(color,color_up)
  }
  if((nrow(xx)>1 | ncol(xx)>1) & plot.fig == T){
    heatmap.plot <- pheatmap(xx, cluster_cols=colClust, cluster_rows=rowClust, clustering_distance_rows = distMethod, clustering_distance_cols = distMethod,
                             clustering_method = clusterMethod, annotation_col=annCol, annotation_row = annRow, 
                             breaks = breakSeq, color = color, annotation_colors = annotation_colors,
                             fontsize_row = fontsize_row, fontsize_col = fontsize_col, cutree_rows = clrn, filename = filename,
                             height = height, width = width, border_color = border_col, display_numbers = display_numbers)
    system(paste("inkscape -l ", filename, ".svg ", filename, sep = ""))
    return(list(row_cluster=gr.row, col_cluster=gr.col, plot=heatmap.plot))
  }else{
    return(list(row_cluster=gr.row, col_cluster=gr.col))
  }
  #dev.off()
}
