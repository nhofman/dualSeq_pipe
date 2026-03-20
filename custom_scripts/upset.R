library(UpSetR)

upset_json <- function(file, out.dir, name, header = 0, sep = ",", start, end){
  json.list <- list("file"=paste0(out.dir,"/",file,".csv"), 
                    "name"=name, 
                    "header"=0, 
                    "separator"=",", 
                    "skip"=0, 
                    "meta"=data.frame("type"="id", "index"=0, "name"="SYMBOL"), 
                    "sets"=data.frame("format"="binary", "start"=start, "end"=end))
  dge.json <- jsonlite::toJSON(json.list, pretty = T, auto_unbox = T)
  write(dge.json, paste0(out.dir,"/",file,".json"))
}

list2binary <- function(data.list, filename){
  data.binary <- Reduce(function(x,y)merge(x,y,by="SYMBOL",all=T),sapply(names(data.list),function(n){
    if(length(data.list[[n]])>0){
      data.df <- data.frame(data.list[[n]],1)
    }else{
      data.df <- data.frame(matrix(ncol=2,nrow=0)) 
    }
    colnames(data.df) <- c("SYMBOL", n)
    return(data.df)
  }, USE.NAMES = T, simplify = F))
  data.binary[is.na(data.binary)] <- 0
  if(!missing(filename)){
    write.csv(data.binary, filename, row.names = F)
  }
  return(data.binary)
}
