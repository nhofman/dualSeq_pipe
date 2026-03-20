library(pheatmap)
library(ComplexHeatmap)

# Plot complex heatmap for a data.frame of expression values (rows = genes, columns = virus and time point)
plotHeatmap <- function(x, filename = "no_name_set.pdf", row_subset = NA, distMethod = "euclidean", clusterMethod = "complete", clrn = 1, clcn = 1, split.dend = NULL, row.dend = T,
                        rowClust = T, colClust = T, rowNames = T, colNames = T, fontsize_row = 0.8, fontsize_col = 10, family = "sans", annCol = NA, annRow = NA, border_col = "grey60", plot.fig = T,
                        filter_col = NA, annotation_colors = NA, height = NA, width = NA, break_step = 0.1, display_numbers = F, cellwidth = NULL, cellheight = NULL, 
                        legend.cut = 1, legend.limit.up = NA, legend.limit.down = NA, ink = "-l", ...){
  if(!dir.exists(dirname(filename))){
    dir.create(dirname(filename), recursive = T)
  }

  if(is.na(row_subset[1])){
    xx <- x
  }else{
    xx <- x[row_subset,, drop = F]
  }
  xx <- na.omit(xx)
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
    color_down <- colorRampPalette(c("blue", "white"))(length(seq(-legend.limit, 0, break_step))) 
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

  if((nrow(xx)>1 | ncol(xx)>1)){
    heat <- draw(Heatmap(as.matrix(xx), col = col_fun, width = cellwidth, height = cellheight, #top_annotation = ca, 
                 show_row_names = rowNames, cluster_rows = rowClust, show_row_dend = row.dend, split = split.dend, cluster_row_slices = F, clustering_distance_rows = distMethod, clustering_method_rows = clusterMethod,
                 show_column_names = colNames, column_names_side = "top", cluster_columns = colClust, column_split = factor(col_split, levels=unique(col_split)), column_names_rot = 0,
                 clustering_distance_columns = distMethod, clustering_method_columns = clusterMethod, column_names_centered = T,
                 cluster_column_slices = F, column_title_gp = gpar(fontsize = fontsize_col*1.5, fontface = "bold"), column_labels = sub("h","",sub(".*_","",colnames(xx))),
                 column_names_gp = gpar(fontsize = fontsize_col, family = family), row_names_gp = gpar(fontsize = fontsize_row, family = family), 
                 heatmap_legend_param = list(title = NULL, labels_gp = gpar(fontsize = fontsize_col)), ...))

    if(plot.fig == T){
      pdf(filename, family = family)
      draw(heat)
      dev.off()
      system(paste0("inkscape ", ink, filename, ".svg ", filename))
    }
    return(heat)
  }else{
    return(list(row_cluster=gr.row, col_cluster=gr.col))
  }
}
