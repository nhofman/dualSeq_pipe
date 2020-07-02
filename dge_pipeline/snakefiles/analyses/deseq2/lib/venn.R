library("VennDiagram")
library("gplots")
library("ggplot2")
library("openxlsx")

draw_venn <- function(venn.list, out_name, h = 3000, w = 3000, cex = 2, cat.cex = 2, imagetype = "tiff", margin = 0.05, 
                      main = "", fill = c("blue", "green", "red", "orange","purple")){
  if(!dir.exists(dirname(out_name))){
    dir.create(dirname(out_name))
  }
  list.intersect <- list()
  venn.plot <- venn.diagram(venn.list, fill = fill, filename = NULL, main = main,
                            imagetype = imagetype, height = h, width = w, cex = cex, cat.cex = cat.cex, margin = margin, ext.text = F)
  ggsave(paste(basename(out_name), ".", imagetype, sep = ""), venn.plot, "svg", dirname(out_name))
  #system(paste("inkscape -A ", out_name, ".pdf ", out_name, ".", imagetype, sep = ""))
  list.intersect <- attr(venn(venn.list, show.plot = F), "intersections")
  intersections <- sapply(list.intersect, "[", seq(max(sapply(list.intersect, length))))
  write.csv(unique(intersections), file = paste0(out_name,".csv"), row.names = F)
  write.xlsx(unique(intersections), file = paste0(out_name,".xlsx"))
  return(list.intersect)
}
