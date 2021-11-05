library(pheatmap)
library(ComplexHeatmap)

# Heatmap 
plotHeatmap <- function(x, filename = "no_name_set.pdf", row_subset = NA, distMethod = "euclidean", clusterMethod = "complete", clrn = 1, clcn = 1,
                        rowClust = T, colClust = T, rowNames = T, colNames = T, fontsize_row = 0.8, fontsize_col = 10, annCol = NA, annRow = NA, border_col = "grey60", plot.fig = T,
                        filter_col = NA, annotation_colors = NA, height = NA, width = NA, break_step = 0.1, display_numbers = F, cellwidth = NULL, cellheight = NULL, 
                        legend.cut = 1, legend.limit.up = NA, legend.limit.down = NA, ...){
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
  
  # Calculate legend limits, if not defined
  if(is.na(legend.limit.up)){
    legend.limit.up <- quantile(unlist(xx), na.rm = TRUE, probs = legend.cut)
  }
  if(is.na(legend.limit.down)){
    legend.limit.down <- quantile(unlist(xx), na.rm = TRUE, probs = (1-legend.cut))  
  }
  
  color <- vector()
  legend.limit <- max(legend.limit.up, abs(legend.limit.down))
  if(-1%in%sign(unlist(xx))){
    color_down <- colorRampPalette(c("blue", "white"))(length(seq(-legend.limit, 0, break_step))) #blue(n=length(breakSeq)-1) #, low="blue", mid = "white", high="red")
    color <- c(color,color_down)
  }
  if(1%in%sign(unlist(xx))){
    color_up <- colorRampPalette(c("white", "red"))(length(seq(0, legend.limit, break_step)))
    color <- c(color,color_up[-1])
  }
  col_fun = circlize::colorRamp2(c(-legend.limit, 0, legend.limit), c("blue", "white", "red"))
  
  if(!is.null(cellwidth) & !is.null(cellheight)){
    cellwidth <- ncol(xx)*unit(cellwidth, "pt")
    cellheight <- nrow(xx)*unit(cellheight, "pt")
  }
  
  col_split <- sub("_.*","",colnames(xx))
  ca <- columnAnnotation(foo = anno_empty(border = FALSE, width = max_text_width(unlist(col_split)) + unit(2, "mm"), height = unit(1, "mm")))
  class(as.matrix(xx)[1,1])
  if((nrow(xx)>1 | ncol(xx)>1) & plot.fig == T){
    pdf(filename)
    draw(Heatmap(as.matrix(xx), show_row_names = rowNames, column_names_side = "top", cluster_columns = colClust, column_split = factor(col_split, levels=unique(col_split)), 
                 top_annotation = ca, show_column_names = colNames, col = col_fun, width = cellwidth, height = cellheight, column_names_rot = 0,
                 cluster_column_slices = F, column_title_gp = gpar(fontsize = fontsize_col), column_labels = sub("h","",sub(".*_","",colnames(xx))),
                 heatmap_legend_param = list(title = NULL, labels_gp = gpar(fontsize = 13)), column_names_gp = gpar(fontsize = fontsize_col), row_names_gp = gpar(fontsize = fontsize_row), ...))
    for(i in 1:length(unique(col_split))){
      decorate_annotation("foo", slice = i, {
             grid.rect(x = 0, width = 1, height = unit(0.1, "mm"), gp = gpar(fill = 1, col = NA), just = "left")
         })
    }
    dev.off()
    # heatmap.plot <- pheatmap(xx, cluster_cols=colClust, cluster_rows=rowClust, clustering_distance_rows = distMethod, clustering_distance_cols = distMethod,
    #                          clustering_method = clusterMethod, annotation_col=annCol, annotation_row = annRow, 
    #                          breaks = breakSeq, color = color, annotation_colors = annotation_colors,
    #                          fontsize_row = fontsize_row, fontsize_col = fontsize_col, cutree_rows = clrn, filename = filename,
    #                          height = height, width = width, border_color = border_col, display_numbers = display_numbers, ...)
    # system(paste("inkscape -l ", filename, ".svg ", filename, sep = ""))
    return(list(row_cluster=gr.row, col_cluster=gr.col))
  }else{
    return(list(row_cluster=gr.row, col_cluster=gr.col))
  }
  #dev.off()
}
